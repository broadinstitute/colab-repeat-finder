import argparse
import pyfastx
import os
import re
import tqdm

from utils.perfect_repeat_tracker import PerfectRepeatTracker
from utils.repeat_tracker import RepeatTracker
from utils.plot_utils import plot_results

# campture Ctrl-C and print a newline
#repeat_trackers = None

#import signal
#signal.signal(signal.SIGINT, lambda x, y: all(r.log(f"Ctrl-C: motif size = {r.motif_size}bp {r.current_position:,d}", force=True) for r in repeat_trackers.values()))


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
    if not getattr(filter_settings, "max_interruptions_by_motif_size"):
        raise ValueError(f"max_interruptions_by_motif_size is not set.")

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
                input_sequence=input_sequence,
                output_intervals=output_intervals,
                verbose=filter_settings.verbose,
                debug=filter_settings.debug,
            )
        else:
            repeat_tracker = PerfectRepeatTracker(
                motif_size=motif_size,
                min_repeats=filter_settings.min_repeats,
                min_span=filter_settings.min_span,
                input_sequence=input_sequence,
                output_intervals=output_intervals,
            )
        repeat_trackers[motif_size] = repeat_tracker

    if all(max_interruptions == 0 for max_interruptions in filter_settings.max_interruptions_by_motif_size.values()):
        # if all repeats must be perfect, the repeat trackers can just advance linearly through the sequence
        if verbose: print("Running perfect repeat trackers for motifs", list(repeat_trackers.keys()))
        if show_progress_bar:
            input_sequence = tqdm.tqdm(input_sequence, unit=" bp", unit_scale=True, total=len(input_sequence))

        repeat_trackers = list(repeat_trackers.values())
        for _ in input_sequence:
            for repeat_tracker in repeat_trackers:
                repeat_tracker.advance()

        for repeat_tracker in repeat_trackers:
            assert not repeat_tracker.advance(), f"{repeat_tracker.motif_size}bp motif RepeatTracker did not reach end of the sequence"

            repeat_tracker.done()
    else:
        # if interruptions are allowed for at least some motifs, some repeat trackers may need to jump backwards
        # before proceeding fowards again through the sequence, so the iteration must deal with this
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
                current_position = min(repeat_tracker.current_position for repeat_tracker in repeat_trackers.values())
                while progress_bar_position < current_position:
                    next(progress_bar)
                    progress_bar_position += 1

            for motif_size in repeat_trackers_that_are_done:
                del repeat_trackers[motif_size]

    return [(start_0based, end, motif) for (start_0based, end), motif in sorted(output_intervals.items())]


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

        with open(args.max_interruptions, "rt") as tsv_file:
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
        # process FASTA input file
        if not args.output_prefix:
            args.output_prefix = re.sub(".fa(sta)?(.gz)?", "", args.input_sequence)

        output_bed_path = f"{os.path.basename(args.output_prefix)}.bed"
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
                output_intervals = detect_repeats(seq, args, verbose=args.verbose, show_progress_bar=args.show_progress_bar, debug=args.debug)
                print(f"Found {len(output_intervals):,d} repeats")
                for start_0based, end, motif in output_intervals:
                    if args.interval:
                        start_0based += interval_start_0based
                        end += interval_start_0based
                    bed_file.write("\t".join([chrom, str(start_0based), str(end), motif]) + "\n")
        print(f"Wrote results to {output_bed_path}")

    elif set(args.input_sequence.upper()) <= set("ACGTN"):
        # process nucleotide sequence specified on the command line
        if args.interval:
            parser.error("The --interval option is only supported for FASTA files.")

        interval_sequence = args.input_sequence

        output_intervals = detect_repeats(args.input_sequence, args, verbose=args.verbose, show_progress_bar=args.show_progress_bar)
        print(f"Found {len(output_intervals):,d} repeat(s):")
        for start_0based, end, motif in output_intervals:
            print(f"{start_0based}-{end} = {motif} x {(end - start_0based)/len(motif):0.1f}")

        if args.output_prefix:
            output_tsv_path = f"{args.output_prefix}.tsv"
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
            print(f"Generating plot {args.plot}")
            plot_results(interval_sequence, output_intervals, args.max_motif_size, args.plot)


if __name__ == "__main__":
    main()