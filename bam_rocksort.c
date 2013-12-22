#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "bam.h"
#include <rocksdb/c.h>

static int strnum_cmp(const char *_a, const char *_b)
{
	const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
	const unsigned char *pa = a, *pb = b;
	while (*pa && *pb) {
		if (isdigit(*pa) && isdigit(*pb)) {
			while (*pa == '0') ++pa;
			while (*pb == '0') ++pb;
			while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
			if (isdigit(*pa) && isdigit(*pb)) {
				int i = 0;
				while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
				return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
			} else if (isdigit(*pa)) return 1;
			else if (isdigit(*pb)) return -1;
			else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
		} else {
			if (*pa != *pb) return (int)*pa - (int)*pb;
			++pa; ++pb;
		}
	}
	return *pa? 1 : *pb? -1 : 0;
}

static int change_SO(bam_header_t *h, const char *so)
{
	char *p, *q, *beg = 0, *end = 0, *newtext;
	if (h->l_text > 3) {
		if (strncmp(h->text, "@HD", 3) == 0) {
			if ((p = strchr(h->text, '\n')) == 0) return -1;
			*p = '\0';
			if ((q = strstr(h->text, "\tSO:")) != 0) {
				*p = '\n'; // change back
				if (strncmp(q + 4, so, p - q - 4) != 0) {
					beg = q;
					for (q += 4; *q != '\n' && *q != '\t'; ++q);
					end = q;
				} else return 0; // no need to change
			} else beg = end = p, *p = '\n';
		}
	}
	if (beg == 0) { // no @HD
		h->l_text += strlen(so) + 15;
		newtext = malloc(h->l_text + 1);
		sprintf(newtext, "@HD\tVN:1.3\tSO:%s\n", so);
		strcat(newtext, h->text);
	} else { // has @HD but different or no SO
		h->l_text = (beg - h->text) + (4 + strlen(so)) + (h->text + h->l_text - end);
		newtext = malloc(h->l_text + 1);
		strncpy(newtext, h->text, beg - h->text);
		sprintf(newtext + (beg - h->text), "\tSO:%s", so);
		strcat(newtext, end);
	}
	free(h->text);
	h->text = newtext;
	return 0;
}

/* compare function for variable-length, little-endian sequence numbers
   we append to RocksDB key names to ensure uniqueness and sort stability */
static int compare_sequence_numbers(const char *seqnoa, size_t alen, const char *seqnob, size_t blen) {
	uint64_t ai = 0, bi = 0;
	for(; alen; alen--) {
		ai <<= 8;
		ai |= (uint8_t) *(seqnoa++);
	}
	for(; blen; blen--) {
		bi <<= 8;
		bi |= (uint8_t) *(seqnob++);
	}
	if (ai < bi) {
		return -1;
	} else if (ai > bi) {
		return 1;
	} else {
		/* shouldn't happen - the sequence numbers are unique... */
		return 0;
	}
}

/* RocksDB comparator for genome positions encoded as uint64_t's (as in
   bam1_lt in bam_sort.c) */
static int pos_rocksdb_comparator(void *ctx, const char *a, size_t alen, const char *b, size_t blen) {
	register uint64_t ai = *((uint64_t*) a);
	register uint64_t bi = *((uint64_t*) b);
	if (ai < bi) {
		return -1;
	} else if (ai > bi) {
		return 1;
	}
	/* rarely - break tie with sequence numbers */
	return compare_sequence_numbers(a+sizeof(uint64_t),alen-sizeof(uint64_t),
		                            b+sizeof(uint64_t),blen-sizeof(uint64_t));
}
static const char* pos_rocksdb_comparator_name(void *ctx) {
	return "samtools_rocksort_pos";
}

/* RocksDB comparator for qnames. The key encoding consists of the null-
   terminated qname string followed immediately by flag:uint32_t. */
static int qname_rocksdb_comparator(void *ctx, const char *a, size_t alen, const char *b, size_t blen) {
	register int t;
	uint32_t aflag, bflag;
	size_t qnamealen, qnameblen;
	if ((t = strnum_cmp(a,b)) != 0) {
		return t;
	}
	qnamealen = strlen(a);
	qnameblen = strlen(b);
	aflag = *(uint32_t*)(a + qnamealen + 1);
	bflag = *(uint32_t*)(b + qnameblen + 1);
	if (aflag < bflag) {
		return -1;
	} else if (aflag > bflag) {
		return 1;
	}
	/* rarely - break tie with sequence numbers */
	return compare_sequence_numbers(a+qnamealen+1+sizeof(uint32_t), alen-qnamealen-1-sizeof(uint32_t),
		                            b+qnameblen+1+sizeof(uint32_t), blen-qnameblen-1-sizeof(uint32_t));
}
static const char* qname_rocksdb_comparator_name(void *ctx) {
	return "samtools_rocksort_qname";
}

/* No-op RocksDB comparator destructor */
static void nop_rocksdb_comparator_destructor(void *c) {
}

/* Load contents of BAM fp into RocksDB rdb, keyed by leftmost genome position
   (encoded as a uint64_t) or by qname (if is_by_qname) */
static int bam_to_rocksdb(bamFile fp, rocksdb_t *rdb, const int is_by_qname, unsigned long long *count) {
	int ret = 0;
	bam1_t *b = 0;
	size_t buflen = 1024;
	rocksdb_writeoptions_t *rdbwropts = 0;
	char *rdberr = 0;
	char *key = 0;
	size_t keybufsz = 0;
	size_t keylen = 0;
	uint64_t seqno;

	key = calloc(2, sizeof(uint64_t));
	if (!key) {
		ret = -4;
		goto cleanup;
	}
	keybufsz = 2*sizeof(uint64_t);

	*count=0;

	rdbwropts = rocksdb_writeoptions_create();
	b = (bam1_t*) malloc(sizeof(bam1_t)+buflen);
	if (!rdbwropts || !b) {
		ret = -4;
		goto cleanup;
	}
	/* disable RocksDB write-ahead log since durability (to power loss, etc.)
	   is unnecessary in this context */
	rocksdb_writeoptions_disable_WAL(rdbwropts, 1);
	memset(b, 0, sizeof(bam1_t));
	
	/* for each BAM record */
	while ((ret = bam_read1(fp,b)) >= 0) {
		/* formulate key for insertion into RocksDB */
		if (is_by_qname) {
			/* null-terminated qname followed by flags byte */
			keylen = b->core.l_qname + sizeof(uint32_t);
			if (keylen+sizeof(uint64_t) > keybufsz) {
				keybufsz = keylen+sizeof(uint64_t);
				kroundup32(keybufsz);
				key = realloc(key, keybufsz);
				if (!key) {
					ret = -4;
					goto cleanup;
				}
			}
			memcpy(key, bam1_qname(b), b->core.l_qname);
			*(uint32_t*)(key+b->core.l_qname) = b->core.flag;
		} else {
			/* uint64_t encoded genome position (tid,pos,strand) */
			*(uint64_t*)key = ((uint64_t)b->core.tid<<32|(b->core.pos+1)<<1|bam1_strand(b));
			keylen = sizeof(uint64_t);
		}
		/* append the variable-length, little-endian sequence number to the key,
		   ensuring uniqueness & sort stability */
		seqno = ++(*count);
		while (seqno > 0) {
			*(uint8_t*)(key+keylen) = (uint8_t) (seqno & 0xFF);
			seqno >>= 8;
			keylen++;
		}

		/* Copy b->data into the buffer directly following the bam1_t at b,
		   so that we can then stick the whole thing into RocksDB. It would
		   be nice to avoid this copying, but bam_read1 may call realloc on
		   b->data so unfortunately we can't just put it anywhere we want. */
		if (buflen < b->l_data) {
			// expand the buffer
			buflen = b->l_data;
			kroundup32(buflen);
			b = (bam1_t*) realloc(b, sizeof(bam1_t)+buflen);
			if (!b) {
				ret = -4;
				goto cleanup;
			}
		}
		memcpy(b+1, b->data, b->l_data);

		/* insert this key & value into RocksDB */
		rocksdb_put(rdb, rdbwropts, key, keylen, (char*) b, sizeof(bam1_t) + b->l_data, &rdberr);
		if (rdberr) {
			fprintf(stderr, "[bam_rocksort_core] %s\n", rdberr);
			ret = -1;
			goto cleanup;
		}
	}
	if (ret != -1)
		fprintf(stderr, "[bam_rocksort_core] truncated file. Continue anyway.\n");
	ret = 0;

cleanup:
	if (rdbwropts) rocksdb_writeoptions_destroy(rdbwropts);
	if (b) {
		if (b->data) free(b->data);
		free(b);
	}
	if (keybufsz) {
		free(key);
	}
	if (rdberr) {
		free(rdberr);
	}
	return ret;
}

/* Traverse the RocksDB and output a BAM file */
static int rocksdb_to_bam(rocksdb_t *rdb, const bam_header_t *header, const char *fnout, const int n_threads, const int level, unsigned long long *count) {
	int ret = 0;
	rocksdb_readoptions_t *rdbrdopts = 0;
	rocksdb_iterator_t *rdbiter = 0;
	char *rdberr = 0;
	bamFile fp = 0;
	bam1_t *b;
	size_t vsz = 0;
	char mode[8];

	*count = 0;

	if (!(rdbrdopts = rocksdb_readoptions_create())) {
		ret = -4;
		goto cleanup;
	}
	rocksdb_readoptions_set_verify_checksums(rdbrdopts, 0);
	if (!(rdbiter = rocksdb_create_iterator(rdb, rdbrdopts))) {
		ret = -4;
		goto cleanup;
	}

	strcpy(mode, "w");
	if (level >= 0) sprintf(mode + 1, "%d", level < 9? level : 9);
	if (!(fp = strcmp(fnout, "-") ? bam_open(fnout, mode) : bam_dopen(fileno(stdout), mode))) {
		fprintf(stderr, "[bam_rocksort_core] fail to create the output file %s\n",fnout);
		ret = -1;
		goto cleanup;
	}
	bam_header_write(fp, header);
	if (n_threads > 1) bgzf_mt(fp, n_threads, 256);

	/* For each record in RocksDB */
	for (rocksdb_iter_seek_to_first(rdbiter); rocksdb_iter_valid(rdbiter); rocksdb_iter_next(rdbiter)) {
        /* Extract bam1_t, knowing that the RocksDB value consists of a bam1_t
		   with a garbage data pointer, immediately followed by the data */
		b = (bam1_t*) rocksdb_iter_value(rdbiter, &vsz);
		b->data = (uint8_t*)(b+1);
		b->m_data = vsz - sizeof(bam1_t);
		/* Write to output BAM */
		bam_write1(fp, b);
		(*count)++;
	}

	rocksdb_iter_get_error(rdbiter, &rdberr);
	if (rdberr) {
		fprintf(stderr, "[bam_rocksort_core] %s\n", rdberr);
		ret = -1;
		goto cleanup;
	}

cleanup:
	if (fp) bam_close(fp);
	if (rdbiter) rocksdb_iter_destroy(rdbiter);
	if (rdberr) free(rdberr);
	if (rdbrdopts) rocksdb_readoptions_destroy(rdbrdopts);

	return ret;
}

/* If TMPDIR environment variable is defined, create a temp directory there.
   Otherwise create a directory name based on prefix.*/
char* choose_rocksdb_path(const char *prefix) {
	char *tmpdir = getenv("TMPDIR");
	char *ans = 0;
	const char *template_const = "/samtools_rocksort_XXXXXX";
	const char *suffix = ".rocksort";

	if (tmpdir && *tmpdir) {
		ans = malloc(strlen(tmpdir)+strlen(template_const)+1);
		strcpy(ans, tmpdir);
		strcat(ans, template_const);
		if (!mkdtemp(ans)) {
			free(ans);
			return NULL;
		}
	} else {
		ans = malloc(strlen(prefix)+strlen(suffix)+1);
		strcpy(ans, prefix);
		strcat(ans, suffix);
	}

	return ans;
}

/*!
  @abstract Sort an unsorted BAM file based on the chromosome order
  and the leftmost position of an alignment

  @param  is_by_qname whether to sort by query name
  @param  fn       name of the file to be sorted
  @param  prefix   prefix of the temporary files (prefix.NNNN.bam are written)
  @param  fnout    name of the final output file to be written
  @param  max_mem  approxiate maximum memory (very inaccurate)
  @return 0 for successful sorting, negative on errors

  @discussion It may create multiple temporary subalignment files
  and then merge them by calling bam_merge_core(). This function is
  NOT thread safe.
 */
int bam_rocksort_core_ext(int is_by_qname, const char *fn, const char *prefix, const char *fnout, const size_t _max_mem, const int n_threads, const int level, const int keep_db) {
	int ret = 0;
	bamFile fp = 0;
	bam_header_t *header = 0;
	rocksdb_env_t *rdbenv = 0;
	rocksdb_t *rdb = 0;
	rocksdb_options_t *rdbopts = 0;
	rocksdb_comparator_t *rdbcomp = 0;
	rocksdb_universal_compaction_options_t *rdbucopts = 0;
	rocksdb_cache_t *rdbcache = 0;
	char *rdberr = 0;
	char *rdbpath = 0;
	unsigned long long count1 = 0, count2 = 0;

	if (!(fp = strcmp(fn, "-")? bam_open(fn, "r") : bam_dopen(fileno(stdin), "r"))) {
		fprintf(stderr, "[bam_rocksort_core] fail to open file %s\n", fn);
		ret = -1;
		goto cleanup;
	}
	if (!(header = bam_header_read(fp))) {
		fprintf(stderr, "[bam_rocksort_core] fail to parse file %s\n", fn);
		ret = -1;
		goto cleanup;
	}
	if (is_by_qname) change_SO(header, "queryname");
	else change_SO(header, "coordinate");

	/* Configure & open RocksDB */
	if (is_by_qname) {
		rdbcomp = rocksdb_comparator_create(0, nop_rocksdb_comparator_destructor, qname_rocksdb_comparator, qname_rocksdb_comparator_name);
	} else {
		rdbcomp = rocksdb_comparator_create(0, nop_rocksdb_comparator_destructor, pos_rocksdb_comparator, pos_rocksdb_comparator_name);
	}
	rdbenv = rocksdb_create_default_env();
	rdbucopts = rocksdb_universal_compaction_options_create();
	rdbopts = rocksdb_options_create();
	rdbcache = rocksdb_cache_create_lru(_max_mem);
	if (!rdbcomp || !rdbopts || !rdbenv || !rdbucopts || !rdbcache) {
		ret = -4;
		goto cleanup;
	}

	/* ROCKSDB CONFIGURATION: ESSENTIALS */
	rocksdb_options_set_create_if_missing(rdbopts, 1);
	rocksdb_options_set_error_if_exists(rdbopts, 1);
	rocksdb_options_set_comparator(rdbopts, rdbcomp);
	rocksdb_options_set_num_levels(rdbopts, 1);
	rocksdb_options_set_compression(rdbopts, rocksdb_snappy_compression);
	rocksdb_options_set_cache(rdbopts, rdbcache);
	rocksdb_options_set_max_open_files(rdbopts, 256);
	/* concurrency */
	rocksdb_options_set_max_background_compactions(rdbopts, n_threads);
	rocksdb_env_set_background_threads(rdbenv, n_threads);
	rocksdb_options_set_env(rdbopts, rdbenv);
	/* turn off some durability and throttling features (unnecessary here) */
	rocksdb_options_set_disable_data_sync(rdbopts, 1);
	rocksdb_options_set_paranoid_checks(rdbopts, 0);
	rocksdb_options_set_disable_seek_compaction(rdbopts, 1);
	rocksdb_options_set_level0_slowdown_writes_trigger(rdbopts, 1<<30);
	rocksdb_options_set_level0_stop_writes_trigger(rdbopts, 1<<30);

	/* TUNE IN-MEMORY COMPACTION */

	/* make the 'open' write buffer at most only 16MB, since insertions
	   into this buffer are on the unparallelized critical path */
	rocksdb_options_set_write_buffer_size(rdbopts, 16<<20);
	/* keep as many of these buffers in memory as possible */
	int max_write_buffer_number = _max_mem * n_threads / (16<<20);
	if (max_write_buffer_number < 16) max_write_buffer_number = 16;
	rocksdb_options_set_max_write_buffer_number(rdbopts, max_write_buffer_number);
	/* background-flush the buffers to disk in batches of at least 1/3 the
	   total number in memory, with a 1MB compression block size */
	rocksdb_options_set_min_write_buffer_number_to_merge(rdbopts, max_write_buffer_number / 3);
	rocksdb_options_set_block_size(rdbopts, 1<<20);

	/* TUNE ON-DISK COMPACTION  */

	/* use 'universal' compaction which has lower write amplification and
	   higher read amplification than level compaction - this is the right
	   tradeoff given that we're ultimately just going to do one sequential
	   scan over the DB and then delete it. */
	rocksdb_options_set_compaction_style(rdbopts, rocksdb_universal_compaction);
	/* background merge each group of 16 files. this uses background cpu and
	   disk bandwidth to reduce the merging work that has to be done in the
	   unparallelized critical path of writing out the final BAM file */
	rocksdb_options_set_level0_file_num_compaction_trigger(rdbopts, 16);
	rocksdb_universal_compaction_options_set_min_merge_width(rdbucopts, 16);
	rocksdb_universal_compaction_options_set_max_merge_width(rdbucopts, 16);
	rocksdb_universal_compaction_options_set_max_size_amplification_percent(rdbucopts, 1<<30);
	rocksdb_options_set_universal_compaction_options(rdbopts, rdbucopts);

	/* useless:
	rocksdb_options_set_memtable_vector_rep(rdbopts);
	rocksdb_options_set_target_file_size_base(rdbopts, 4*16777216);
	rocksdb_options_set_source_compaction_factor(rdbopts, 10);
	rocksdb_universal_compaction_options_set_stop_style(rdbucopts, rocksdb_similar_size_compaction_stop_style);
	*/

	rdbpath = choose_rocksdb_path(prefix);
	if (!rdbpath) {
		fprintf(stderr, "[bam_rocksort] Failed creating temporary directory; check given prefix and TMPDIR environment variable.\n");
		ret = -1;
		goto cleanup;
	}
	if (!(rdb = rocksdb_open(rdbopts, rdbpath, &rdberr))) {
		fprintf(stderr, "[bam_rocksort] %s\n", rdberr ? rdberr : "failed creating RocksDB");
		ret = -4;
		goto cleanup;
	}

	/* Load input BAM into RocksDB */ 
	fprintf(stderr, "[bam_rocksort_core] Loading in %s (you can change this with the TMPDIR environment variable)...\n", rdbpath);
	if ((ret = bam_to_rocksdb(fp, rdb, is_by_qname, &count1)) != 0) {
		goto cleanup;
	}

	/* TODO: halt or wait for any ongoing background compactions */

	/* Export sorted BAM from RocksDB */
	fprintf(stderr, "[bam_rocksort_core] Writing %llu records to %s...\n", count1, fnout);
	if ((ret = rocksdb_to_bam(rdb, header, fnout, n_threads, level, &count2)) != 0) {
		goto cleanup;
	}

	if (count1 != count2) {
		fprintf(stderr, "[bam_rocksort_core] BUG: read %llu, wrote %llu records!\n", count1, count2);
		ret = -1;
		goto cleanup;
	}

	fprintf(stderr, "[bam_rocksort_core] OK\n");

cleanup:
	if(fp) bam_close(fp);
	if(header) bam_header_destroy(header);
	if(rdb) {
		rocksdb_close(rdb);
		if (!keep_db) rocksdb_destroy_db(rdbopts, rdbpath, &rdberr);
	}
	if(rdbcomp) rocksdb_comparator_destroy(rdbcomp);
	if(rdbopts) rocksdb_options_destroy(rdbopts);
	if(rdbenv) rocksdb_env_destroy(rdbenv);
	if(rdbucopts) rocksdb_universal_compaction_options_destroy(rdbucopts);
	if(rdbcache) rocksdb_cache_destroy(rdbcache);
	if(rdberr) free(rdberr);
	if(rdbpath) free(rdbpath);

	if (ret != 0) {
		fprintf(stderr, "[bam_rocksort_core] Exit status %d\n", ret);
	}
	return ret;
}

int bam_rocksort_core(int is_by_qname, const char *fn, const char *prefix, size_t max_mem)
{
	int ret;
	char *fnout = calloc(strlen(prefix) + 4 + 1, 1);
	sprintf(fnout, "%s.bam", prefix);
	ret = bam_rocksort_core_ext(is_by_qname, fn, prefix, fnout, max_mem, 0, -1, 0);
	free(fnout);
	return ret;
}

int bam_rocksort(int argc, char *argv[])
{
	size_t max_mem = 768<<20; // 512MB
	int c, is_by_qname = 0, is_stdout = 0, ret = 0,
	    n_threads = 0, level = -1, full_path = 0, keep_db = 0;
	char *fnout;
	while ((c = getopt(argc, argv, "fnokm:@:l:")) >= 0) {
		switch (c) {
		case 'f': full_path = 1; break;
		case 'o': is_stdout = 1; break;
		case 'n': is_by_qname = 1; break;
		case 'k': keep_db=1; break;
		case 'm': {
				char *q;
				max_mem = strtol(optarg, &q, 0);
				if (*q == 'k' || *q == 'K') max_mem <<= 10;
				else if (*q == 'm' || *q == 'M') max_mem <<= 20;
				else if (*q == 'g' || *q == 'G') max_mem <<= 30;
				break;
			}
		case '@': n_threads = atoi(optarg); break;
		case 'l': level = atoi(optarg); break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   samtools rocksort [options] <in.bam> <out.prefix>\n\n");
		fprintf(stderr, "Options: -n        sort by read name\n");
		fprintf(stderr, "         -f        use <out.prefix> as full file name instead of prefix\n");
		fprintf(stderr, "         -o        final output to stdout\n");
		fprintf(stderr, "         -k        keep RocksDB instead of deleting it when done\n");
		fprintf(stderr, "         -l INT    compression level, from 0 to 9 [-1]\n");
		fprintf(stderr, "         -@ INT    number of sorting and compression threads [1]\n");
		fprintf(stderr, "         -m INT    max memory per thread; suffix K/M/G recognized [768M]\n");
		fprintf(stderr, "\n");
		return 1;
	}

	if (is_stdout) fnout = strdup("-");
	else {
		fnout = calloc(strlen(argv[optind+1]) + 4 + 1, 1);
		sprintf(fnout, "%s%s", argv[optind+1], full_path? "" : ".bam");
	}

	if (bam_rocksort_core_ext(is_by_qname, argv[optind], argv[optind+1], fnout, max_mem, n_threads, level, keep_db) < 0) ret = 1;
	free(fnout);
	return ret;
}
