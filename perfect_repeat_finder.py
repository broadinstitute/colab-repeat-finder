import argparse
import pyfastx
import os
import re
import tqdm

from utils.perfect_repeat_tracker import PerfectRepeatTracker
from utils.plot_utils import plot_results

def detect_repeats(input_sequence, filter_settings, verbose=False, show_progress_bar=False, debug=False):
    """Detect repeats in a given input sequence.

    Args:
        input_sequence (str): The input sequence.
        filter_settings (Namespace): Filter settings from command-line args.
        verbose (bool): Print verbose output for debugging.
        show_progress_bar (bool): Show a progress bar while traversing the input sequence.
        debug (bool): Currently identical to verbose; print verbose output for debugging.

    Returns:
        list: A list of (start_0based, end, motif) tuples representing all detected repeats in the input sequence.
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

    interval_start_0based = getattr(filter_settings, "interval_start_0based", 0)
    interval_end = getattr(filter_settings, "interval_end", len(input_sequence))
    input_sequence_length = len(input_sequence)

    # skip any Ns at the beginning and end of the interval
    while interval_start_0based < interval_end and input_sequence[interval_start_0based] == "N":
        interval_start_0based += 1
    while interval_end > interval_start_0based and input_sequence[interval_end - 1] == "N":
        interval_end -= 1
        input_sequence_length -= 1

    input_sequence = input_sequence[interval_start_0based:input_sequence_length]

    # generate all intervals
    output_intervals = {}
    repeat_trackers = {}
    for motif_size in range(filter_settings.min_motif_size, filter_settings.max_motif_size + 1):
        repeat_tracker = PerfectRepeatTracker(
            motif_size=motif_size,
            min_repeats=filter_settings.min_repeats,
            min_span=filter_settings.min_span,
            input_sequence=input_sequence,
            output_intervals=output_intervals,
        )
        repeat_trackers[motif_size] = repeat_tracker

    end_position = interval_end - interval_start_0based
    position_iter = range(len(input_sequence))
    if show_progress_bar:
        position_iter = tqdm.tqdm(position_iter, unit=" bp", unit_scale=True, total=end_position)

    for position in position_iter:
        for repeat_tracker in repeat_trackers.values():
            repeat_tracker.advance()

        # keep processing past the end of the interval until every repeat that started inside the interval has ended
        if position + 1 >= end_position and all(
                repeat_tracker.current_run_start_0based() >= end_position for repeat_tracker in repeat_trackers.values()):
            break

    for repeat_tracker in repeat_trackers.values():
        if interval_end == input_sequence_length:
            assert not repeat_tracker.advance(), f"{repeat_tracker.motif_size}bp motif RepeatTracker did not reach end of the sequence"
        repeat_tracker.done()

    # drop repeats that start after the interval (they're only seen because processing continued past its end)
    repeats = [
        (start_0based + interval_start_0based, end + interval_start_0based, motif)
        for (start_0based, end), motif in sorted(output_intervals.items()) if start_0based < end_position
    ]
    if verbose or debug:
        for start_0based, end, motif in repeats:
            print(f"Found {(end - start_0based) / len(motif):0.1f} x {motif} at {start_0based}-{end}")

    return repeats

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group = parser.add_argument_group("Repeat Filters")
    group.add_argument("-min", "--min-motif-size", default=1, type=int, help="Minimum motif size in base pairs.")
    group.add_argument("-max", "--max-motif-size", default=50, type=int, help="Maximum motif size in base pairs.")
    group.add_argument("--min-repeats", default=3, type=int, help="The minimum number of repeats to look for.")
    group.add_argument("--min-span", default=9, type=int, help="The repeats should span at least this many consecutive "
                                                               "bases in the input sequence.")
    parser.add_argument("-i", "--interval", help="Only consider sequence from this interval "
                                                "(chrom:start_0based-end). A repeat that starts inside the interval is "
                                                "reported in full even if it extends past the interval's end, but a "
                                                "repeat that begins before the interval is truncated to the portion "
                                                "from the interval's start onward, and may be reported with a "
                                                "different (rotated) motif as a result.")
    parser.add_argument("-p", "--plot", help="Write out a plot with this filename.")
    parser.add_argument("-o", "--output-prefix", help="The output filename prefix. Results are written to <prefix>.bed "
                                                      "if the input is a FASTA file, or to <prefix>.tsv if it's a "
                                                      "nucleotide sequence.")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output.")
    parser.add_argument("--debug", action="store_true", help="Currently identical to --verbose.")
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

    interval_sequence = None
    if os.path.isfile(args.input_sequence):
        if args.plot and not args.interval:
            parser.error("--plot requires --interval when the input is a FASTA file.")
        if not args.output_prefix:
            args.output_prefix = os.path.basename(re.sub(r"\.(fa|fasta|fna)(\.b?gz)?$", "", args.input_sequence))

        output_bed_path = f"{args.output_prefix}.bed"
        if args.interval:
            interval = re.split("[:-]", args.interval)
            if len(interval) != 3:
                parser.error(f"Invalid --interval format. Must be chrom:start_0based-end")
            args.interval_chrom, args.interval_start_0based, args.interval_end = interval
            args.interval_start_0based = int(args.interval_start_0based)
            args.interval_end = int(args.interval_end)

            # build an index so the requested chromosome can be looked up directly, without scanning preceding contigs
            fasta_index = pyfastx.Fasta(args.input_sequence)
            if args.interval_chrom not in fasta_index:
                parser.error(f"Chromosome {args.interval_chrom} not found in the input FASTA file")
            chrom_sequence = fasta_index[args.interval_chrom].seq

            fasta_entries = [(args.interval_chrom, chrom_sequence)]
        else:
            # read (name, sequence) tuples without building an index, since contigs are processed in file order
            fasta_entries = pyfastx.Fasta(args.input_sequence, build_index=False)

        with open(output_bed_path, "wt") as bed_file:
            for chrom, seq in fasta_entries:
                seq_len = len(seq)
                if args.interval:
                    args.interval_end = min(args.interval_end, seq_len)
                    seq_len = args.interval_end - args.interval_start_0based
                print(f"Processing {chrom} ({seq_len:,d} bp)")
                output_intervals = detect_repeats(
                    seq, args, verbose=args.verbose, show_progress_bar=args.show_progress_bar, debug=args.debug)
                print(f"Found {len(output_intervals):,d} repeats")
                for start_0based, end, motif in output_intervals:
                    bed_file.write("\t".join([chrom, str(start_0based), str(end), motif]) + "\n")

        print(f"Wrote results to {output_bed_path}")

        if args.interval:
            # for plotting, make repeat coordinates relative to the interval and clip them to it
            interval_sequence = chrom_sequence[args.interval_start_0based:args.interval_end]
            output_intervals = [
                (start_0based - args.interval_start_0based, min(end, args.interval_end) - args.interval_start_0based, motif)
                for start_0based, end, motif in output_intervals
            ]

    elif set(args.input_sequence.upper()) <= set("ACGTN"):
        # process nucleotide sequence specified on the command line
        if args.interval:
            parser.error("The --interval option is only supported for FASTA files.")

        interval_sequence = args.input_sequence
        if not args.output_prefix:
            args.output_prefix = "repeats"
        output_tsv_path = f"{args.output_prefix}.tsv"

        output_intervals = detect_repeats(
            args.input_sequence, args, verbose=args.verbose, show_progress_bar=args.show_progress_bar, debug=args.debug)
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
            print(f"Warning: The input sequence is too long ({len(interval_sequence):,d} bp). Skipping plot...")
        else:
            plot_results(interval_sequence, output_intervals, args.max_motif_size, args.plot)


if __name__ == "__main__":
    main()