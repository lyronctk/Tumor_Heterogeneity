#!/bin/bash

## This script will pull out the number of different species in each window across the selector and build a table for use in R, Matlab, or other statistical programming tool
## Writen by Joanne Soo and Lyron Co Ting Keh (4/27/17)

if [ -z "$7" ]
then
  echo "    ****"
  echo "    Usage: $0 <bam_file> <germline_file> <freq_file> <cosmic_file> <selector> <window_length>"
  echo "    bam_file: sorted bam file for sample to be analyzed"
  echo "    germline_file: germline freq file for patient of sample to be analyzed"
  echo "    freq_file: sample freq file of sample to be analyzed"
  echo "    cosmic_file: sample cosmic.split file of sample to be analyzed"
  echo "    selector: selector for the sample used"
  echo "    window_length: desired length of window"
  echo "    window_slide: desired shift in sliding window"
  echo "****"
  exit 1
fi

file=$1
germline_file=$2
freq_file=$3
cosmic_file=$4
selector=$5
window_length=$6
window_slide=$7
ref_genome=/data/indexes/hg19.fa

echo '----Processing bam file...' $file
file_name=$(echo $file | cut -d '_' -f 2- | cut -d '.' -f 1)

echo '----Limit to selector space...'
samtools view -b -L $selector $file > $file_name.test.bam

echo '----Remove reference bases...'
samtools view -b $file_name.test.bam | samtools fillmd -e -b - $ref_genome > $file_name.refremoved.bam
samtools view $file_name.refremoved.bam | cut -d "  " -f3-4,6,10 > $file_name.refremoved.txt
rm $file_name.test.bam

echo '----Subset mutated rows...'
length_raw=$(less $file_name.refremoved.txt | head -1 | cut -f 3)
length=${length_raw%?}
dummy=$(head -c $length < /dev/zero | tr '\0' '\141' | sed -e 's,a,=,g')
awk -v dummy="$dummy" '$4!=dummy' $file_name.refremoved.txt > $file_name.mutatedrows.txt 

echo '----Summarizing cosmic file and depth file...'
cut -f13-14 $cosmic_file | tail -n +2 > $file_name.summarized.cosmic.split.txt
cut -f1-3 $file_name.refremoved.txt > $file_name.depth.txt

# clang++ -std=c++11 -stdlib=libc++ generateWindows.cpp
echo "$file_name.mutatedrows.txt $selector $germline_file $freq_file $file_name.depth.txt $cosmic_file $window_length $window_slide $file_name.window_counts.txt" | ./a.out

# clean up intermediate files
rm $file_name.refremoved.bam $file_name.refremoved.txt $file_name.mutatedrows.txt $file_name.summarized.cosmic.split.txt $file_name.depth.txt