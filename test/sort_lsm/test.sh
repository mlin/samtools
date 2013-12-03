#!/bin/bash -e

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ] ; do SOURCE="$(readlink "$SOURCE")"; done
srcdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

htslibdir="${srcdir}/../../../htslib"
samtoolsdir="${srcdir}/../.."
samtools="${samtoolsdir}/samtools"

if [ -z "$TMPDIR" ]; then
	export TMPDIR="/tmp"
fi

# Fetch a test BAM file
if [ -z "$TEST_BAM" ]; then
	TEST_BAM="${TMPDIR}/samtools_sort_lsm_test.bam"
	if ! [ -f "$TEST_BAM" ]; then
		wget -O "$TEST_BAM" "https://test.galaxyproject.org/library_common/download_dataset_from_folder?library_id=04c38ea196069a5d&cntrller=library&use_panels=False&id=4a29c5f588caf7fb"
	fi
fi

# Shuffle it
testdir=$(mktemp -d --tmpdir samtools_sort_lsm_test.XXXXXX)
echo "Temporary test directory: $testdir"
cd "$testdir"
shuffled_bam="${testdir}/samtools_sort_lsm_test_shuffled.bam"
echo "Shuffling $TEST_BAM..."
$samtools bamshuf $TEST_BAM "${testdir}/samtools_sort_lsm_test_shuffled"

# Compile and test bamsorted.c (FIXME ugly)
cc -o bamsorted -O2 ${srcdir}/bamsorted.c -I"$htslibdir" -I"$samtoolsdir" "$htslibdir/libhts.a" "$samtoolsdir/libbam.a" "$htslibdir/faidx.o" -lz -lpthread
echo -n "Positive control: "
./bamsorted "$TEST_BAM"
echo -n "Negative control: "
code=1
./bamsorted "$shuffled_bam" || code=$?
if [ "$code" -eq "0" ]; then
	echo "should have failed: bamsorted $shuffled_bam"
	exit 1
fi

echo "Running samtools sort_lsm..."
time $samtools sort_lsm -l 1 -@ 2 -m 16M "$shuffled_bam" "${testdir}/samtools_sort_lsm_test_sorted"
sorted_bam="${testdir}/samtools_sort_lsm_test_sorted.bam"

echo -n "Verifying product: "
./bamsorted "$sorted_bam"

# Clean up
rm -rf $testdir
echo "Test passed. Cleaned up $testdir"
