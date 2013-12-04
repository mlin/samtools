#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "bam.h"
#include <leveldb/c.h>

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
   we append to LevelDB key names to ensure uniqueness and sort stability */
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

/* LevelDB comparator for genome positions encoded as uint64_t's (as in
   bam1_lt in bam_sort.c) */
static int pos_leveldb_comparator(void *ctx, const char *a, size_t alen, const char *b, size_t blen) {
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
static const char* pos_leveldb_comparator_name(void *ctx) {
	return "samtools_sort_lsm_pos";
}

/* LevelDB comparator for qnames. The key encoding consists of the null-
   terminated qname string followed immediately by flag:uint32_t. */
static int qname_leveldb_comparator(void *ctx, const char *a, size_t alen, const char *b, size_t blen) {
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
static const char* qname_leveldb_comparator_name(void *ctx) {
	return "samtools_sort_lsm_qname";
}

/* No-op LevelDB comparator destructor */
static void nop_leveldb_comparator_destructor(void *c) {
}

/* Load contents of BAM fp into LevelDB ldb, keyed by leftmost genome position
   (encoded as a uint64_t) or by qname (if is_by_qname) */
static int bam_to_leveldb(bamFile fp, leveldb_t *ldb, const int is_by_qname, unsigned long long *count) {
	int ret = 0;
	bam1_t *b = 0;
	size_t buflen = 1024;
	leveldb_writeoptions_t *ldbwropts = 0;
	char *ldberr = 0;
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

	ldbwropts = leveldb_writeoptions_create();
	b = (bam1_t*) malloc(sizeof(bam1_t)+buflen);
	if (!ldbwropts || !b) {
		ret = -4;
		goto cleanup;
	}
	memset(b, 0, sizeof(bam1_t));
	
	/* for each BAM record */
	while ((ret = bam_read1(fp,b)) >= 0) {
		/* formulate key for insertion into LevelDB */
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
		   so that we can then stick the whole thing into LevelDB. It would
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

		/* insert this key & value into LevelDB */
		leveldb_put(ldb, ldbwropts, key, keylen, (char*) b, sizeof(bam1_t) + b->l_data, &ldberr);
		if (ldberr) {
			fprintf(stderr, "[bam_sort_lsm_core] %s\n", ldberr);
			ret = -1;
			goto cleanup;
		}
	}
	if (ret != -1)
		fprintf(stderr, "[bam_sort_lsm_core] truncated file. Continue anyway.\n");
	ret = 0;

cleanup:
	if (ldbwropts) leveldb_writeoptions_destroy(ldbwropts);
	if (b) {
		if (b->data) free(b->data);
		free(b);
	}
	if (keybufsz) {
		free(key);
	}
	if (ldberr) {
		free(ldberr);
	}
	return ret;
}

/* Traverse the LevelDB and output a BAM file */
static int leveldb_to_bam(leveldb_t *ldb, const bam_header_t *header, const char *fnout, const int n_threads, const int level, unsigned long long *count) {
	int ret = 0;
	leveldb_readoptions_t *ldbrdopts = 0;
	leveldb_iterator_t *ldbiter = 0;
	char *ldberr = 0;
	bamFile fp = 0;
	bam1_t *b;
	size_t vsz = 0;
	char mode[8];

	*count = 0;

	if (!(ldbrdopts = leveldb_readoptions_create())) {
		ret = -4;
		goto cleanup;
	}
	leveldb_readoptions_set_verify_checksums(ldbrdopts, 0);
	if (!(ldbiter = leveldb_create_iterator(ldb, ldbrdopts))) {
		ret = -4;
		goto cleanup;
	}

	strcpy(mode, "w");
	if (level >= 0) sprintf(mode + 1, "%d", level < 9? level : 9);
	if (!(fp = strcmp(fnout, "-") ? bam_open(fnout, mode) : bam_dopen(fileno(stdout), mode))) {
		fprintf(stderr, "[bam_sort_lsm_core] fail to create the output file %s\n",fnout);
		ret = -1;
		goto cleanup;
	}
	bam_header_write(fp, header);
	if (n_threads > 1) bgzf_mt(fp, n_threads, 256);

	/* For each record in LevelDB */
	for (leveldb_iter_seek_to_first(ldbiter); leveldb_iter_valid(ldbiter); leveldb_iter_next(ldbiter)) {
        /* Extract bam1_t, knowing that the LevelDB value consists of a bam1_t
		   with a garbage data pointer, immediately followed by the data */
		b = (bam1_t*) leveldb_iter_value(ldbiter, &vsz);
		b->data = (uint8_t*)(b+1);
		b->m_data = vsz - sizeof(bam1_t);
		/* Write to output BAM */
		bam_write1(fp, b);
		(*count)++;
	}

	leveldb_iter_get_error(ldbiter, &ldberr);
	if (ldberr) {
		fprintf(stderr, "[bam_sort_lsm_core] %s\n", ldberr);
		ret = -1;
		goto cleanup;
	}

cleanup:
	if (fp) bam_close(fp);
	if (ldbiter) leveldb_iter_destroy(ldbiter);
	if (ldberr) free(ldberr);
	if (ldbrdopts) leveldb_readoptions_destroy(ldbrdopts);

	return ret;
}

/* If TMPDIR environment variable is defined, create a temp directory there.
   Otherwise create a directory name based on prefix.*/
char* choose_leveldb_path(const char *prefix) {
	char *tmpdir = getenv("TMPDIR");
	char *ans = 0;
	const char *template_const = "/samtools_sort_lsm_XXXXXX";
	const char *suffix = ".sort_lsm";

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
int bam_sort_lsm_core_ext(int is_by_qname, const char *fn, const char *prefix, const char *fnout, const size_t _max_mem, const int n_threads, const int level) {
	int ret = 0;
	bamFile fp = 0;
	bam_header_t *header = 0;
	leveldb_t *ldb = 0;
	leveldb_options_t *ldbopts = 0;
	leveldb_comparator_t *ldbcomp = 0;
	char *ldberr = 0;
	char *ldbpath = 0;
	unsigned long long count1 = 0, count2 = 0;

	if (!(fp = strcmp(fn, "-")? bam_open(fn, "r") : bam_dopen(fileno(stdin), "r"))) {
		fprintf(stderr, "[bam_sort_lsm_core] fail to open file %s\n", fn);
		ret = -1;
		goto cleanup;
	}
	if (!(header = bam_header_read(fp))) {
		fprintf(stderr, "[bam_sort_lsm_core] fail to parse file %s\n", fn);
		ret = -1;
		goto cleanup;
	}
	if (is_by_qname) change_SO(header, "queryname");
	else change_SO(header, "coordinate");

	/* Configure & open LevelDB */
	if (is_by_qname) {
		ldbcomp = leveldb_comparator_create(0, nop_leveldb_comparator_destructor, qname_leveldb_comparator, qname_leveldb_comparator_name);
	} else {
		ldbcomp = leveldb_comparator_create(0, nop_leveldb_comparator_destructor, pos_leveldb_comparator, pos_leveldb_comparator_name);
	}
	ldbopts = leveldb_options_create();
	if (!ldbcomp || !ldbopts) {
		ret = -4;
		goto cleanup;
	}
	leveldb_options_set_create_if_missing(ldbopts, 1);
	leveldb_options_set_error_if_exists(ldbopts, 1);
	leveldb_options_set_paranoid_checks(ldbopts, 0);
	leveldb_options_set_verify_compactions(ldbopts, 0);
	leveldb_options_set_write_buffer_size(ldbopts, _max_mem/2);
	leveldb_options_set_block_size(ldbopts, 16777216);
	leveldb_options_set_compression(ldbopts, leveldb_snappy_compression);
	leveldb_options_set_comparator(ldbopts, ldbcomp);

	ldbpath = choose_leveldb_path(prefix);
	if (!ldbpath) {
		fprintf(stderr, "[bam_sort_lsm] Failed creating temporary directory; check given prefix and TMPDIR environment variable.\n");
		ret = -1;
		goto cleanup;
	}
	if (!(ldb = leveldb_open(ldbopts, ldbpath, &ldberr))) {
		fprintf(stderr, "[bam_sort_lsm] %s\n", ldberr ? ldberr : "failed creating LevelDB");
		ret = -4;
		goto cleanup;
	}

	/* Load input BAM into LevelDB */ 
	fprintf(stderr, "[bam_sort_lsm_core] Sorting in %s (you can change this with the TMPDIR environment variable)...\n", ldbpath);
	if ((ret = bam_to_leveldb(fp, ldb, is_by_qname, &count1)) != 0) {
		goto cleanup;
	}

	/* Export sorted BAM from LevelDB */
	fprintf(stderr, "[bam_sort_lsm_core] Sorted %llu records. Writing out %s...\n", count1, fnout);
	if ((ret = leveldb_to_bam(ldb, header, fnout, n_threads, level, &count2)) != 0) {
		goto cleanup;
	}

	if (count1 != count2) {
		fprintf(stderr, "[bam_sort_lsm_core] BUG: read %llu, wrote %llu records!\n", count1, count2);
		ret = -1;
		goto cleanup;
	}

	fprintf(stderr, "[bam_sort_lsm_core] OK\n");

cleanup:
	if(fp) bam_close(fp);
	if(header) bam_header_destroy(header);
	if(ldb) {
		leveldb_close(ldb);
		leveldb_destroy_db(ldbopts, ldbpath, &ldberr);
	}
	if(ldbcomp) leveldb_comparator_destroy(ldbcomp);
	if(ldbopts) leveldb_options_destroy(ldbopts);
	if(ldberr) free(ldberr);
	if(ldbpath) free(ldbpath);

	if (ret != 0) {
		fprintf(stderr, "[bam_sort_lsm_core] Exit status %d\n", ret);
	}
	return ret;
}

int bam_sort_lsm_core(int is_by_qname, const char *fn, const char *prefix, size_t max_mem)
{
	int ret;
	char *fnout = calloc(strlen(prefix) + 4 + 1, 1);
	sprintf(fnout, "%s.bam", prefix);
	ret = bam_sort_lsm_core_ext(is_by_qname, fn, prefix, fnout, max_mem, 0, -1);
	free(fnout);
	return ret;
}

int bam_sort_lsm(int argc, char *argv[])
{
	size_t max_mem = 768<<20; // 512MB
	int c, is_by_qname = 0, is_stdout = 0, ret = 0, n_threads = 0, level = -1, full_path = 0;
	char *fnout;
	while ((c = getopt(argc, argv, "fnom:@:l:")) >= 0) {
		switch (c) {
		case 'f': full_path = 1; break;
		case 'o': is_stdout = 1; break;
		case 'n': is_by_qname = 1; break;
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
		fprintf(stderr, "Usage:   samtools sort [options] <in.bam> <out.prefix>\n\n");
		fprintf(stderr, "Options: -n        sort by read name\n");
		fprintf(stderr, "         -f        use <out.prefix> as full file name instead of prefix\n");
		fprintf(stderr, "         -o        final output to stdout\n");
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

	if (bam_sort_lsm_core_ext(is_by_qname, argv[optind], argv[optind+1], fnout, max_mem, n_threads, level) < 0) ret = 1;
	free(fnout);
	return ret;
}
