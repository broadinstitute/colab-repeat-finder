set -ex

# run a slightly modified version of PERF @ https://github.com/bw2/perf
# detect all exact repeats on chr22 with motif sizes 1-6bp

# -a output html
# -l min length
# -u  min number of repeats
# -i input sequence

gunzip -c ../chr22.fa.gz > chr22.fa
time PERF -a  -l 9  -u 3  -i chr22.fa
rm chr22.fa

mv chr22_perf.bed chr22.perf.bed
rm chr22_perf.tsv chr22_perf.html

set +x
echo PERF output: 
wc -l chr22.perf.bed
