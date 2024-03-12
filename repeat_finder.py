# copied from https://colab.research.google.com/drive/1wa_96-zPbsJpQpEyVnQMJ2bwMtYZ5Mq8#scrollTo=36bfbda3
import argparse
import pyfastx
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import os
import re

plt.rcParams['figure.figsize'] = [16.5, 5]
plt.rcParams['font.size'] = 12


class RepeatTracker:
    """This class tracks repeats of a single motif size in a given input sequence. It outputs all perfect repeats of
    this that pass filter criteria while scanning the input sequence from left to right"""

    def __init__(self, motif_size, min_repeats, min_span, input_sequence, output_intervals):
        self.motif_size = motif_size
        self.min_repeats = min_repeats
        self.min_span = min_span
        self.input_sequence = input_sequence
        self.output_intervals = output_intervals

        self.current_position = 0
        self.run_length = 1

    def advance(self):
        seq = self.input_sequence
        i = self.current_position
        period = self.motif_size

        if i >= len(seq) - period:
            return

        if seq[i] == seq[i + period]:
            self.current_position += 1
            self.run_length += 1
            return

        #print(f"Mismatch for {period} bp motif at {i} {seq[i:i+period]}")
        self.output_interval_if_it_passes_filters()
        self.current_position += 1
        self.run_length = 1

    def done(self):
        self.output_interval_if_it_passes_filters()

    def output_interval_if_it_passes_filters(self):
        seq = self.input_sequence
        i = self.current_position
        period = self.motif_size

        # found mismatch, check if previous run is worth outputing
        start_0based = i - self.run_length + 1
        motif = seq[start_0based:start_0based + period]
        if motif == "N":
            return

        #if len(motif) == 2:
        #    print(f"start_0based={start_0based}, i={i}, run_length={self.run_length}, motif={motif}")
        if self.run_length + period - 1 >= self.min_span and self.run_length + period - 1 >= self.min_repeats * period:
            while i < len(seq) - 1 and seq[i+1] == seq[i+1 - period]:
                #print("Increment by 1")
                self.run_length += 1
                i += 1

        if self.run_length >= self.min_span and self.run_length >= self.min_repeats * period:
            end = i + 1

            previous_motif = self.output_intervals.get((start_0based, end))
            if previous_motif is None or len(motif) < len(previous_motif):
                #print(f"Adding: {(end - start_0based)/len(motif):0.1f} x {motif} [{start_0based}-{end}]")
                self.output_intervals[(start_0based, end)] = motif
            #else:
            #    print(f"Skipping: {(end - start_0based)/len(motif):0.1f} x {motif} [{start_0based}-{end}]")


def shift_string_by(string, shift):
    """Shift a string to the right by a given number of characters (eg. "AGTTT" shifted by 2 becomes "TTAGT")."""

    return string[-shift:] + string[:-shift]


def get_period_matrix(min_motif_size, max_motif_size, query):
    min_motif_size = max(min_motif_size, 1)
    max_motif_size = min(max_motif_size, len(query) // 2)
    matrix = [[0 for _ in range(len(query))] for _ in range(max_motif_size)]
    
    for row in range(min_motif_size - 1, max_motif_size):
        period = row + 1
        for column in range(len(query) - period):
            if query[column] == query[column + period]:
                # set the matrix value to a hash that stays the same as long as the motif is the same
                value = abs(hash(shift_string_by(query[column:column + period], column % period)))
                matrix[row][column] = value

    return matrix


def detect_repeats(input_sequence, filter_settings, verbose=False):
    """Detect repeats in a given input sequence.

    Args:
        input_sequence (str): The input sequence.
        filter_settings (Namespace): Filter settings from command-line args.

    Return:
        list: A list of (start_0based, end, motif) tuples reprsenting all detected repeats in the input sequence.
    """
    if not getattr(filter_settings, "min_motif_size") or filter_settings.min_motif_size < 1:
        raise ValueError(f"min_motif_size is set to {filter_settings.min_motif_size}. It must be at least 1.")
    if not getattr(filter_settings, "max_motif_size") or filter_settings.max_motif_size < filter_settings.min_motif_size:
        raise ValueError(f"max_motif_size is set to {filter_settings.max_motif_size}. It must be at least min_motif_size.")
    if not getattr(filter_settings, "min_repeats") or filter_settings.min_repeats < 1:
        raise ValueError(f"min_repeats is set to {filter_settings.min_repeats}. It must be at least 1.")
    if not getattr(filter_settings, "min_span") or filter_settings.min_span < 1:
        raise ValueError(f"min_span is set to {filter_settings.min_span}. It must be at least 1.")

    # dictionary to store detected repeats. The key is (start_0based, end) and the value is the detected motif.
    output_intervals = {}

    # generate all intervals
    run_trackers = []
    for motif_size in range(filter_settings.min_motif_size, filter_settings.max_motif_size + 1):
        run_tracker = RepeatTracker(
            motif_size=motif_size,
            min_repeats=filter_settings.min_repeats,
            min_span=filter_settings.min_span,
            input_sequence=input_sequence,
            output_intervals=output_intervals,
        )
        run_trackers.append(run_tracker)

    position_range = range(len(input_sequence))
    if verbose:
        import tqdm
        position_range = tqdm.tqdm(position_range, unit=" bp", unit_scale=True, total=len(input_sequence))

    for i in position_range:
        for run_tracker in run_trackers:
            run_tracker.advance()

    for run_tracker in run_trackers:
        run_tracker.done()

    return [(start_0based, end, motif) for (start_0based, end), motif in sorted(output_intervals.items())]


def plot_results(input_sequence, output_intervals, max_motif_size, output_filename):

    matrix = [[0 for _ in range(len(input_sequence))] for _ in range(max_motif_size)]
    for start_0based, end, motif in output_intervals:
        row = len(motif) - 1
        for i in range(start_0based, end + 1):
            matrix[row][i] = abs(hash(motif)) % 10 + 1

    #%%
    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10, 8))
    cmap = ListedColormap(['#F3F3F3'] + list(plt.get_cmap("Pastel2", 12).colors))
    ax1.set_xticks(range(0, len(input_sequence)), minor=True)
    ax1.matshow(matrix, aspect="auto", cmap=cmap)
    ax1.set_xlabel("Query position")
    ax1.set_ylabel("Period")

    scores = [sum([1 for _ in row if _ > 0]) for row in matrix]
    periods = list(range(len(scores)))
    scores = [scores[i]/(len(input_sequence) - periods[i] + 1) for i in range(len(scores))]
    ax2.bar(periods, scores)
    ax2.set_ylabel("Fraction of matches")
    ax2.set_xlabel("Period");

    #plt.show()
    plt.savefig(output_filename)
    print(f"Wrote {output_filename}")


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_argument_group("Repeat Filters")
    group.add_argument("-min", "--min-motif-size", default=1, type=int, help="Minimum motif size in base pairs.")
    group.add_argument("-max", "--max-motif-size", default=50, type=int, help="Maximum motif size in base pairs.")
    group.add_argument("--min-repeats", default=3, type=int, help="The minimum number of repeats to look for.")
    group.add_argument("--min-span", default=9, type=int, help="The repeats should span at least this many consecutive "
                                                               "bases in the input sequence.")
    parser.add_argument("-i", "--interval", help="Only consider sequence from this interval (chrom:start_0based-end).")
    parser.add_argument("-p", "--plot", help="Write out a plot with this filename.")
    parser.add_argument("-o", "--output-prefix", help="The output filename prefix for the output TSV file. If the input "
                                                      "is a FASTA file, a BED file will also be generated.")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output.")
    parser.add_argument("input_sequence", help="The nucleotide sequence, or a FASTA file path")

    args = parser.parse_args()

    interval_sequence = None
    if os.path.isfile(args.input_sequence):
        if not args.output_prefix:
            args.output_prefix = re.sub(".fa(sta)?(.gz)?", "", args.input_sequence)

        output_bed_path = f"{args.output_prefix}.bed"
        fasta_entries = pyfastx.Fasta(args.input_sequence)
        if args.interval:
            interval = re.split("[:-]", args.interval)
            if len(interval) != 3:
                parser.error(f"Invalid --plot-interval format. Must be chrom:start_0based-end")
            interval_chrom, interval_start_0based, interval_end = interval
            interval_start_0based = int(interval_start_0based)
            interval_end = int(interval_end)

            # iterate over chromosomes in the FASTA file
            if interval_chrom not in fasta_entries:
                parser.error(f"Chromosome {interval_chrom} not found in the input FASTA file")

            interval_sequence = fasta_entries[interval_chrom][interval_start_0based:interval_end].seq
            fasta_entries = [
                argparse.Namespace(name=args.interval, seq=interval_sequence)
            ]

        with open(output_bed_path, "wt") as bed_file:
            for fasta_entry in fasta_entries:
                seq = fasta_entry.seq
                chrom = fasta_entry.name
                print(f"Processing {chrom} ({len(seq):,d} bp)")
                output_intervals = detect_repeats(seq, args, verbose=args.verbose)
                print(f"Found {len(output_intervals):,d} repeats")
                for start_0based, end, motif in output_intervals:
                    bed_file.write("\t".join([chrom, str(start_0based), str(end), motif]) + "\n")

    elif set(args.input_sequence.upper()) <= set("ACGTN"):
        interval_sequence = args.input_sequence
        if not args.output_prefix:
            args.output_prefix = "repeats"
        output_tsv_path = f"{args.output_prefix}.tsv"

        output_intervals = detect_repeats(args.input_sequence, args)
        print(f"Found {len(output_intervals):,d} repeats")

        with open(output_tsv_path, "wt") as tsv_file:
            tsv_file.write("\t".join(["start_0based", "end", "motif"]) + "\n")
            for start_0based, end, motif in output_intervals:
                row = "\t".join([str(start_0based), str(end), motif]) + "\n"
                tsv_file.write(row)
    else:
        parser.error(f"Invalid input: {args.input_sequence}. This should be a FASTA file path or a string of nucleotides.")

    # detect repeats
    if args.plot and interval_sequence:
        if len(interval_sequence) > 5_000:
            print(f"Warning: The input sequence is too long ({len(sequence_to_plot):,d} bp). Skipping plot...")
        else:
            plot_results(interval_sequence, output_intervals, args.max_motif_size, args.plot)


if __name__ == "__main__":
    main()