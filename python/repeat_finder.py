import argparse
import pyfastx
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import os
import re
import tqdm

from utils.pure_repeat_tracker import PureRepeatTracker
from utils.repeat_tracker import RepeatTracker

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


def detect_repeats(input_sequence, filter_settings, verbose=False, show_progress_bar=False, debug=False):
    """Detect repeats in a given input sequence.

    Args:
        input_sequence (str): The input sequence.
        filter_settings (Namespace): Filter settings from command-line args.
        verbose (bool): Print verbose output for debugging.
        show_progress_bar (bool): Show a progress bar while traversing the input sequence.
        debug (bool): Print debug output.

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

    input_sequence = input_sequence.upper()

    # generate all intervals
    output_intervals = {}
    repeat_trackers = {}
    for motif_size in range(filter_settings.min_motif_size, filter_settings.max_motif_size + 1):

        max_interruptions = filter_settings.max_interruptions_by_motif_size.get(motif_size, 0)
        if max_interruptions > 0:
            repeat_tracker = RepeatTracker(
                motif_size=motif_size,
                min_repeats=filter_settings.min_repeats,
                min_span=filter_settings.min_span,
                max_interruptions=max_interruptions,
                allow_multibase_homopolymer_motifs=filter_settings.allow_multibase_homopolymer_motifs,
                allow_ending_with_different_motif=filter_settings.allow_ending_with_different_motif,
                input_sequence=input_sequence,
                output_intervals=output_intervals,
                verbose=filter_settings.verbose,
                debug=filter_settings.debug,
            )
        else:
            repeat_tracker = PureRepeatTracker(
                motif_size=motif_size,
                min_repeats=filter_settings.min_repeats,
                min_span=filter_settings.min_span,
                input_sequence=input_sequence,
                output_intervals=output_intervals,
            )
        repeat_trackers[motif_size] = repeat_tracker

    if all(max_interruptions == 0 for max_interruptions in filter_settings.max_interruptions_by_motif_size.values()):
        if show_progress_bar:
            input_sequence = tqdm.tqdm(input_sequence, unit=" bp", unit_scale=True, total=len(input_sequence))

        repeat_trackers = list(repeat_trackers.values())
        for _ in input_sequence:
            for repeat_tracker in repeat_trackers:
                repeat_tracker.advance()
    else:
        progress_bar = None
        if show_progress_bar:
            progress_bar = iter(tqdm.tqdm(range(0, len(input_sequence)), unit=" bp", unit_scale=True, total=len(input_sequence)))
            progress_bar_position = 0

        while repeat_trackers:
            repeat_trackers_that_are_done = []
            for i, (motif_size, repeat_tracker) in enumerate(repeat_trackers.items()):
                advanced_to_next_position = repeat_tracker.advance()
                if not advanced_to_next_position:
                    repeat_tracker.done()
                    repeat_trackers_that_are_done.append(motif_size)
                    if verbose: print(f"Done with motif size {motif_size}bp")

                if progress_bar:
                    while progress_bar_position < repeat_tracker.current_position:
                        next(progress_bar)
                        progress_bar_position += 1

            for motif_size in repeat_trackers_that_are_done:
                del repeat_trackers[motif_size]

    return [(start_0based, end, motif) for (start_0based, end), motif in sorted(output_intervals.items())]


def plot_results(input_sequence, output_intervals, max_motif_size, output_path):
    """Plot the repeats detected in the given input sequence. This code was copied from
    https://colab.research.google.com/drive/1wa_96-zPbsJpQpEyVnQMJ2bwMtYZ5Mq8#scrollTo=36bfbda3

    Args:
        input_sequence (str): The input sequence.
        output_intervals (list): A list of (start_0based, end, motif) tuples representing all repeats detected in the input sequence.
        max_motif_size (int): The maximum motif size in base pairs.
        output_path (str): The output image filename for the plot.

    """

    plt.rcParams['figure.figsize'] = [16.5, 5]
    plt.rcParams['font.size'] = 12

    matrix = [[0 for _ in range(len(input_sequence))] for _ in range(max_motif_size)]
    for start_0based, end, motif in output_intervals:
        row = len(motif) - 1
        for i in range(start_0based, end + 1):
            matrix[row][i] = abs(hash(motif)) % 10 + 1

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

    plt.savefig(output_path)
    print(f"Wrote {output_path}")


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_argument_group("Repeat Filters")
    group.add_argument("-min", "--min-motif-size", default=1, type=int, help="Minimum motif size in base pairs.")
    group.add_argument("-max", "--max-motif-size", default=50, type=int, help="Maximum motif size in base pairs.")
    group.add_argument("--min-repeats", default=3, type=int, help="The minimum number of repeats to look for.")
    group.add_argument("--min-span", default=9, type=int, help="The repeats should span at least this many consecutive "
                                                               "bases in the input sequence.")
    group.add_argument("--max-interruptions", default="0", help="How many bases within a motif are allowed to vary "
        "across repeats. For example, setting this to 1 would be equivalent to matching repeats with patterns like "
        "(GCN)* or (AANAG)* or (ANAAG)*, where a single base can vary. If the value is an integer, it is applied to "
        "all motif sizes greater than 2. If a config file path is provided instead, it should be a TSV table with "
        "2 columns and the following header: 'MotifSize\tMaxInterruptions'. The values in the MotifSize columns should "
        "be integers greater than 2, and the values in the MaxInterruptions column should be non-negative integers.")

    group.add_argument("--allow-multibase-homopolymer-motifs", action="store_true", help="When --min-motif-size is "
        "greater than 1, allow motifs like 'AA' or 'TTTT'. By default, they are discarded. Similarly, when interruptions "
        "are enabled, allow repeats like (ANAAA)* where all bases are either a homopolymer or can vary across repeats.")
    group.add_argument("--allow-ending-with-different-motif", action="store_true", help="By default, interrupted "
        "repeats must still end with at least 1 copy of the exact same motif that they started with. This option "
        "disables that rule.")

    parser.add_argument("-i", "--interval", help="Only consider sequence from this interval (chrom:start_0based-end).")
    parser.add_argument("-p", "--plot", help="Write out a plot with this filename.")
    parser.add_argument("-o", "--output-prefix", help="The output filename prefix for the output TSV file. If the input "
                                                      "is a FASTA file, a BED file will also be generated.")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output.")
    parser.add_argument("--debug", action="store_true", help="Print debugging output.")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show progress bar.")
    parser.add_argument("input_sequence", help="The nucleotide sequence, or a FASTA file path")

    args = parser.parse_args()

    if args.min_motif_size < 1:
        parser.error(f"--min-motif-size is set to {args.min_motif_size}. It must be at least 1.")
    if args.max_motif_size < args.min_motif_size:
        parser.error(f"--max-motif-size is set to {args.max_motif_size}. It must be at least --min-motif-size.")
    if args.min_repeats < 1:
        parser.error(f"--min-repeats is set to {args.min_repeats}. It must be at least 1.")
    if args.min_span < 1:
        parser.error(f"--min-span is set to {args.min_span}. It must be at least 1.")

    args.max_interruptions_by_motif_size = {}
    try:
        max_interruptions = int(args.max_interruptions)
        for motif_size in range(args.min_motif_size, args.max_motif_size + 1):
            if motif_size < 3 or motif_size <= max_interruptions:
                args.max_interruptions_by_motif_size[motif_size] = 0
            else:
                args.max_interruptions_by_motif_size[motif_size] = max_interruptions
    except ValueError:
        pass

    if not args.max_interruptions_by_motif_size and os.path.isfile(args.max_interruptions):
        for motif_size in range(args.min_motif_size, args.max_motif_size + 1):
            args.max_interruptions_by_motif_size[motif_size] = 0

        with (open(args.max_interruptions, "rt") as tsv_file):
            header = tsv_file.readline().strip().split("\t")
            if header != ["MotifSize", "MaxInterruptions"]:
                parser.error(f"Invalid header in {args.max_interruptions}. Expected 2 columns named: 'MotifSize\tMaxInterruptions'")

            for i, line in enumerate(tsv_file):
                fields = line.strip().split("\t")
                if len(fields) != 2:
                    parser.error(f"Invalid line {i+1} in {args.max_interruptions}. Expected 2 tab-delimited columns.")

                motif_size, max_interruptions = fields
                try:
                    motif_size = int(motif_size)
                    max_interruptions = int(max_interruptions)
                except ValueError:
                    parser.error(f"Invalid value in {args.max_interruptions} line {i+1}: {motif_size}\t{max_interruptions}. Both values must be integers.")

                if motif_size < 1:
                    parser.error(f"Invalid motif size in {args.max_interruptions} line {i+1}: {motif_size}. It must be greater than 0.")
                if motif_size < 3 and max_interruptions > 0:
                    parser.error(f"Invalid setting in {args.max_interruptions} line {i+1}: Interruptions are not supported in motif sizes less than 3bp.")
                if max_interruptions >= motif_size:
                    parser.error(f"Invalid setting in {args.max_interruptions} line {i+1}: MaxInterruptions must be less than MotifSize.")

                if motif_size >= args.min_motif_size and motif_size <= args.max_motif_size:
                    args.max_interruptions_by_motif_size[motif_size] = max_interruptions

    if not args.max_interruptions_by_motif_size:
        parser.error(f"Invalid --max-interruptions value: {args.max_interruptions}. It must be an integer or a TSV file path.")

    if args.verbose:
        for motif_size, max_interruptions in args.max_interruptions_by_motif_size.items():
            print(f"interruptions limit = {max_interruptions} per repeat for {motif_size}bp motifs")

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
                argparse.Namespace(name=interval_chrom, seq=interval_sequence)
            ]

        with open(output_bed_path, "wt") as bed_file:
            for fasta_entry in fasta_entries:
                seq = fasta_entry.seq
                chrom = fasta_entry.name
                print(f"Processing {chrom} ({len(seq):,d} bp)")
                output_intervals = detect_repeats(seq, args, verbose=args.verbose, show_progress_bar=args.show_progress_bar)
                print(f"Found {len(output_intervals):,d} repeats")
                for start_0based, end, motif in output_intervals:
                    if args.interval:
                        start_0based += interval_start_0based
                        end += interval_start_0based
                    bed_file.write("\t".join([chrom, str(start_0based), str(end), motif]) + "\n")
        print(f"Wrote results to {output_bed_path}")

    elif set(args.input_sequence.upper()) <= set("ACGTN"):
        if args.interval:
            parser.error("The --interval option is only supported for FASTA files.")

        interval_sequence = args.input_sequence
        if not args.output_prefix:
            args.output_prefix = "repeats"
        output_tsv_path = f"{args.output_prefix}.tsv"

        output_intervals = detect_repeats(args.input_sequence, args, verbose=args.verbose, show_progress_bar=args.show_progress_bar)
        print(f"Found {len(output_intervals):,d} repeats")

        with open(output_tsv_path, "wt") as tsv_file:
            tsv_file.write("\t".join(["start_0based", "end", "motif"]) + "\n")
            for start_0based, end, motif in output_intervals:
                row = "\t".join([str(start_0based), str(end), motif]) + "\n"
                tsv_file.write(row)
        print(f"Wrote results to {output_tsv_path}")
    else:
        parser.error(f"Invalid input: {args.input_sequence}. This should be a FASTA file path or a string of nucleotides.")

    if args.plot and interval_sequence:
        if len(interval_sequence) > 5_000:
            print(f"Warning: The input sequence is too long ({len(sequence_to_plot):,d} bp). Skipping plot...")
        else:
            plot_results(interval_sequence, output_intervals, args.max_motif_size, args.plot)


if __name__ == "__main__":
    main()