if [ -z "$3" ]
then
        echo '    Usage: sh clonal_heterogeneity_wrapper.sh <selector> <window_size> <window_slide>'
        exit 1
fi
selector=$1
window_size=$2
window_slide=$3

#1. list files in folder
ls *.sorted.bam > bamfiles.txt
ls *Normal* > normfiles.txt
ls *cosmic* > cosmicfiles.txt
ls *Tumor*freq.paired.Q30.txt > freqfiles.txt 

#2. create windows for selector
echo '----Creating windows for selector...'
while read -r a b c
do
        window_end=$((c-window_size+1))
        for (( n=$b; n <= $window_end; n += $window_size));
        do
                echo $a $n  >> windows.txt
        done
done < $selector

parallel --jobs 8 --xapply sh filterMutatedRows.sh bamfiles.txt :::: normfiles.txt :::: freqfiles.txt :::: cosmicfiles.txt ::: $selector ::: $window_size ::: $window_slide

rm bamfiles.txt normfiles.txt cosmicfiles.txt freqfiles.txt windows.txt