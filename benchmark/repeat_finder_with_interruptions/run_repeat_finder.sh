set -ex

# run a slightly modified version of PERF @ https://github.com/bw2/perf
# detect all exact repeats on chr22 with motif sizes 1-6bp

# -a output html
# -l min length
# -u  min number of repeats
# -i input sequence

gunzip -c ../chr22.fa.gz > chr22.fa
#time python3 ../../python/repeat_finder.py  --min-motif-size 1 --max-motif-size 6 --min-repeats 3 --min-span 9  --max-interruptions 1  -o chr22_repeats chr22.fa --show-progress-bar
time python3 ../../python/repeat_finder.py  --min-motif-size 1 --max-motif-size 6 --min-repeats 3 --min-span 9  -i chr22:15226420-15229480  --plot repeats.svg ../chr22.fa.gz --show-progress-bar  --max-interruptions 1 -o chr22_repeats
rm chr22.fa

set +x
echo output: 
wc -l chr22_repeats.bed

exit 0


