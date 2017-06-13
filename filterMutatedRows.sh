#!/bin/bash

##This script will pull out the number of different species in each window across the selector and build a table for use in R, Matlab, or other statistical programming tool
##Writen by Joanne Soo 4/27/17

# if [ -z "$6" ]
# then
#   echo "  "
#   echo "Usage: $0 <bam_file> <germline_file> <freq_file> <cosmic_file> <selector> <window_length>"
#   echo "bam_file: sorted bam file for sample to be analyzed"
#   echo "germline_file: germline freq file for patient of sample to be analyzed"
#   echo "freq_file: sample freq file of sample to be analyzed"
#   echo "cosmic_file: sample cosmic.split file of sample to be analyzed"
#   echo "selector: selector for the sample used"
#   echo "window_length: desired length of window"
#   echo "  "
#   exit 1
# fi

# file=$1
# germline_file=$2
# freq_file=$3
# cosmic_file=$4
# selector=$5
# window_length=$6
# ref_genome=/data/indexes/hg19.fa

# echo 'processing bam file' $file

# echo 'limit to selector space'
# samtools view -b -L $selector $file > $file_name.test.bam

# echo 'remove reference bases'
# samtools view -b $file_name.test.bam | samtools fillmd -e -b - $ref_genome > $file_name.refremoved.bam
# samtools view $file_name.refremoved.bam | cut -d "  " -f3-4,6,10 > $file_name.refremoved.txt
# rm $file_name.test.bam

# echo 'subset mutated rows'
# length_raw=$(less $file_name.refremoved.txt | head -1 | cut -f 3)
# length=${length_raw%?}
# dummy=$(head -c $length < /dev/zero | tr '\0' '\141' | sed -e 's,a,=,g')
# awk -v dummy="$dummy" '$4!=dummy' $file_name.refremoved.txt > $file_name.mutatedrows.txt 

clang++ -std=c++11 -stdlib=libc++ generateWindows.cpp
echo "DLBCL021-Tumor.mutatedrows.txt selector.bed Sample_DLBCL021_Normal.singleindex-deduped.sorted.freq.paired.Q30.txt Sample_DLBCL021_Tumor.singleindex-deduped.sorted.freq.paired.Q30.txt 50 windows.txt" | ./a.out