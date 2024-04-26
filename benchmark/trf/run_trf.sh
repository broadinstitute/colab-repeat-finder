set -ex

# run TandemRepeatFinder to detect only perfect repeats (installed from https://github.com/Benson-Genomics-Lab/TRF)
gunzip -c ../chr22.fa.gz > chr22.fa
time trf409.macosx chr22.fa 2 10000000 10000000 80 10 2 2000 -h -l 6 -ngs > chr22.trf.dat
rm chr22.fa

# use str-analysis script to convert .dat to .bed   (installed via  python3 -m pip install str-analysis) 
python3 -m str_analysis.convert_dat_to_bed chr22.trf.dat

# filter to 1-6bp motifs
cat chr22.trf.bed | awk -F $'\t' '{ if( length($4) <= 6 ) { print $0 } }' > chr22.trf.filtered.bed
mv chr22.trf.filtered.bed chr22.trf.bed

set +x
echo TRF output:
wc -l chr22.trf.bed
