This repo contains:

- [**perfect_repeat_finder.py**](python/perfect_repeat_finder.py) -  a tool that takes a nucleotide sequence or FASTA file as input, finds all perfect tandem repeats (ie. those without interruptions) that pass user-defined criteria, and outputs their exact genomic coordinates and repeat motifs to a BED file.
- *other tools are under development*

---
Example command-line:

```
python3 python/perfect_repeat_finder.py \
  --min-span 9 \
  --min-repeats 3 \
  --min-motif-size 2 \
  --max-motif-size 6 \
  --interval chr1:1-10000000 \
  --show-progress-bar \
  /path/to/hg38.fa
```

It takes 55 seconds and detects all 63,738 perfect repeats in the first 10Mb of chr1 that pass the following criteria:  
- span at least 9bp from start to end
- include at least 3 perfect repeats of some motif
- have 2bp ≤ motif size ≤ 6bp

**NOTE:** run time is proportional to the length of the input sequence and the range of motif sizes included in the output. 

---
Example BED output file:

```
...
chr1	10397	10442	CCCTAA
chr1	10440	10468	CCCTAA
chr1	10485	10498	GCCC
chr1	10629	10635	GC
chr1	10652	10658	AG
chr1	10658	10664	GC
...
```
---
All command-line options:

```
usage: perfect_repeat_finder.py [-h] [-min MIN_MOTIF_SIZE] [-max MAX_MOTIF_SIZE] [--min-repeats MIN_REPEATS] [--min-span MIN_SPAN] [-i INTERVAL] [-p PLOT] [-o OUTPUT_PREFIX] [--verbose] [--debug] [--show-progress-bar] input_sequence

positional arguments:
  input_sequence        The nucleotide sequence, or a FASTA file path

optional arguments:
  -h, --help            show this help message and exit
  -i INTERVAL, --interval INTERVAL
                        Only consider sequence from this interval (chrom:start_0based-end). (default: None)
  -p PLOT, --plot PLOT  Write out a plot with this filename. (default: None)
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        The output filename prefix for the output TSV file. If the input is a FASTA file, a BED file will also be generated. (default: None)
  --verbose             Print verbose output. (default: False)
  --debug               Print debugging output. (default: False)
  --show-progress-bar   Show progress bar. (default: False)

Repeat Filters:
  -min MIN_MOTIF_SIZE, --min-motif-size MIN_MOTIF_SIZE
                        Minimum motif size in base pairs. (default: 1)
  -max MAX_MOTIF_SIZE, --max-motif-size MAX_MOTIF_SIZE
                        Maximum motif size in base pairs. (default: 50)
  --min-repeats MIN_REPEATS
                        The minimum number of repeats to look for. (default: 3)
  --min-span MIN_SPAN   The repeats should span at least this many consecutive bases in the input sequence. (default: 9)
```
