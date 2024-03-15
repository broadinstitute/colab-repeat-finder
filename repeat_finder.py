import argparse
import pyfastx
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import os
import re


class RepeatTracker:
    """This class tracks repeats of a single motif size in a given input sequence. It outputs all perfect repeats of
    this motif size that pass filter criteria while scanning the input sequence from left to right.
    If the max_interruptions parameter is > 0 then it will also look for interrupted repeats where a fixed number of
    positions within the motif can vary across repeats. In this mode, the algorithm traverses left to right, but
    may occasionally jump backward to revisit some previously-evaluated positions."""

    def __init__(self, motif_size, min_repeats, min_span, max_interruptions, input_sequence, output_intervals, verbose=False, debug=False):
        """Initialize a RepeatTracker object.

        Args:
            motif_size (int): The motif size (in base pairs) that will be tracked by this RepeatTracker.
            min_repeats (int): Only add repeats to the output when there's at least this many repeats in a row.
            min_span (int):  Only add repeats to the output when they span at least this many base pairs.
            max_interruptions (int): How many bases within a motif are allowed to vary across repeats.
            input_sequence (str): The input sequence.
            output_intervals (dict): A dictionary to store detected repeats. The key is (start_0based, end) and the
                value is the detected motif.
            verbose (bool): Print verbose output for debugging.
        """

        self.motif_size = motif_size
        self.min_repeats = min_repeats
        self.min_span = min_span
        self.max_interruptions = max_interruptions
        self.input_sequence = input_sequence
        self.output_intervals = output_intervals
        self.verbose = verbose

        self.current_position = 0   # 0-based position in the input sequence
        self.run_length = 0

        self.current_interrupted_positions_in_motif = set()
        self.position_of_first_interruption = None

        self.previous_output_interval = None  # 3-tuple (start_0based, end, motif) of the most recent interval added to output_intervals

        # a list of strings for debugging. Each entry represents traversal paths through the input sequence, starting with a
        # position integer, followed by 1 or more nucleotides, ending at the end of a repeat interval or at the end of the input sequence
        self.debug_strings = ["0 - "] if debug else None

    def log(self, message):
        """Print a message if verbose output is enabled."""
        if not self.verbose:
            return

        start_0based = max(0, self.current_position - self.run_length)
        motif = self.input_sequence[start_0based : start_0based + self.motif_size]
        print(f"{message:100s}  || "
              f"{len(self.input_sequence):,d}bp  [{start_0based}:{self.current_position}], "
              f"run={self.run_length}, "
              f"i0={self.current_position}: "
              f"{(self.current_position - start_0based)/len(motif):0.2f} x {motif} "
              f"==> {self.input_sequence[start_0based : self.current_position]}")

    def advance(self):
        """Increment current position within the input sequence while updating internal state and
        recording any detected repeats.

        Return False if the end of the input sequence has been reached.
        
        Pure repeats: 
        
        For pure repeats, the algorithm moves forward 1 base at a time, incrementing run_length until it reaches  
        a position in the sequence where the current base differs from the base at curret position + motif size.
        Then, it checks if the base(s) between the few extra bases between the curent position and positions 
        ahead of the interruption at current position + motif size are still consistent with the current repeat run. 
        If yes, it adds them. Then, it checks if the current repeat run satisfies filter critera for minimum span and
        minimum number of repeats. If yes, it adds the repeat run to the output_intervals dictionary. 
        
        Interrupted repeats:
        
        For interrupted repeats, the algorithm is the same as above. However, when it encounters an interruption, it
        saves the currrent position within the motif.
               |
        AACAACAGCAACAACAACAAC
           
        """

        seq = self.input_sequence
        if self.current_position >= len(seq) - self.motif_size:
            self.log(f"advance() reached end of sequence")
            return False

        if self.debug_strings:
            self.debug_strings[-1] += seq[self.current_position]

        if seq[self.current_position] == seq[self.current_position + self.motif_size]:
            self.log(f"Adding position {self.current_position} to repeat run")
            self.run_length += 1
            self.current_position += 1
            return True

        if self.verbose:
            print(f"Position {self.current_position} ({seq[self.current_position]}) doesn't match position {self.current_position + self.motif_size} ({seq[self.current_position + self.motif_size]})")

        if self.run_length > 0:
            if self.position_of_first_interruption is None:
                # save the current position
                self.log(f"Saving position {self.current_position} of the first interruption")
                self.position_of_first_interruption = self.current_position

            current_position_in_motif = self.run_length % self.motif_size
            if len(self.current_interrupted_positions_in_motif) < self.max_interruptions:
                self.log(f"Allowing base #{current_position_in_motif + 1} within the motif to vary across repeats.")
                self.current_interrupted_positions_in_motif.add(current_position_in_motif)

            if current_position_in_motif in self.current_interrupted_positions_in_motif:
                if self.verbose:
                    print(f"Continuing to position {self.current_position + 1} despite interruption at position {current_position_in_motif}")
                self.run_length += 1
                self.current_position += 1
                return True

        self.output_interval_if_it_passes_filters()

        self.run_length = 0
        self.current_position += 1
        return True

    def done(self):
        """Output the last interval if it passes filters"""
        self.output_interval_if_it_passes_filters()

        if self.debug_strings:
            self.print_debug_strings()

    def output_interval_if_it_passes_filters(self):
        """Check internal state to see if enough repeats have accumulated to output an interval. If yes, add it to
        the output. self.current_position will be the 0-based position in the input sequence where the interval ends.
        """
        if self.run_length + self.motif_size < self.min_span or self.run_length + self.motif_size < self.min_repeats * self.motif_size:
            return

        self.log(f"Checking whether to output current repeat run")
        seq = self.input_sequence

        start_0based = self.current_position - self.run_length
        motif = seq[start_0based : start_0based + self.motif_size]
        if "N" in motif:
            return

        # extend the interval to the right, up to self.motif_size - 1 bases
        while self.current_position < len(seq) and (seq[self.current_position] == seq[self.current_position - self.motif_size] or (self.run_length % self.motif_size) in self.current_interrupted_positions_in_motif):
            self.run_length += 1
            self.current_position += 1
            if self.debug_strings and self.current_position < len(seq):
                self.debug_strings[-1] += seq[self.current_position]
            self.log(f"Extend {motif} repeat to {start_0based}-{self.current_position}: {seq[start_0based:self.current_position]}")

        if self.run_length >= self.min_span and self.run_length >= self.min_repeats * self.motif_size:
            end = self.current_position

            motif_previously_detected_at_this_interval = self.output_intervals.get((start_0based, end))
            if (
                # if another motif size has already been record for this exact interval by another RepeatTracker,
                # keep the short motif (eg. replace AAGAAG with AAG)
                motif_previously_detected_at_this_interval is None or len(motif) < len(motif_previously_detected_at_this_interval)
            ) and (
                # only output overlapping intervals if they extend beyond the previously added interval by at least one
                # repeat unit. This becomes important when allowing interruptions.
                self.previous_output_interval is None or end - self.previous_output_interval[1] >= self.motif_size
            ):
                final_motif = motif
                for k in self.current_interrupted_positions_in_motif:
                    final_motif = final_motif[:k] + "N" + final_motif[k+1:]

                if self.verbose: print(f"==> Yes! Adding repeat run to output.")
                self.output_intervals[(start_0based, end)] = final_motif
                self.previous_output_interval = (start_0based, end, final_motif)
            else:
                if self.verbose: print(f"==> No, this interval signicantly overlaps another, previously detected repeat sequence.")
        else:
            if self.verbose: print(f"==> No, this run doesn't span enough bases or have enough repeats of a {self.motif_size} bp motif.")

        # start searching again from the position where the 1st interruption was detected

        if self.position_of_first_interruption is not None:
            if self.verbose:
                print(f"Jumping back from position {self.current_position} to position of first interruption "
                  f"({self.position_of_first_interruption - 1})")

            if self.debug_strings:
                self.debug_strings.append(f"{self.position_of_first_interruption} - ")

            self.current_position = self.position_of_first_interruption
            self.position_of_first_interruption = None

        self.run_length = 0

        self.current_interrupted_positions_in_motif.clear()

    def print_debug_strings(self):
        """Print debug strings for debugging."""
        if self.debug_strings:
            print("=======")
            print(f"Debug strings for {self.motif_size}bp motif:")
            for i, s in enumerate(self.debug_strings):
                print("  "*i, s)


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
    run_trackers = {}
    for motif_size in range(filter_settings.min_motif_size, filter_settings.max_motif_size + 1):
        run_tracker = RepeatTracker(
            motif_size=motif_size,
            min_repeats=filter_settings.min_repeats,
            min_span=filter_settings.min_span,
            max_interruptions=filter_settings.max_interruptions_by_motif_size.get(motif_size, 0),
            input_sequence=input_sequence,
            output_intervals=output_intervals,
            verbose=filter_settings.verbose,
            debug=debug,
        )
        run_trackers[motif_size] = run_tracker

    progress_bar = None
    if show_progress_bar:
        import tqdm
        progress_bar = tqdm.tqdm(unit=" bp", unit_scale=True, total=len(input_sequence))
        progress_bar.n = len(input_sequence)

    while run_trackers:
        run_trackers_that_are_done = []
        for motif_size, run_tracker in run_trackers.items():
            advanced = run_tracker.advance()
            if not advanced:
                run_tracker.done()
                run_trackers_that_are_done.append(motif_size)
                if verbose: print(f"Done with motif size {motif_size}bp")

            if progress_bar:
                progress_bar.n = run_tracker.current_position
                progress_bar.refresh()

        for motif_size in run_trackers_that_are_done:
            del run_trackers[motif_size]

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