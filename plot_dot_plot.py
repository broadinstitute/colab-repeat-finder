"""
This script takes a DNA sequence of length N and generates an NxN square image - aka. dot plot - where
pixel (i,j) is dark only if positions i and j have the same base

Additional details:
- to reduce noise, pixels are only shown if they are part of a diagonal stretch of dark pixels longer than some threshold
- the image can optionally be downsampled to smaller width x height
"""

import argparse
import gzip
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os
from numpy import kaiser
import pyfaidx

def parse_interval(interval_string):
    """Parses interval string like "chr1:12345-54321" and returns 3-tuple (chrom, start, end)"""

    try:
        tokens = interval_string.split(":")
        chrom = ":".join(tokens[:-1])  # some super-contig names have : in them
        start, end = map(int, tokens[-1].split("-"))
    except Exception as e:
        raise ValueError(f"Unable to parse interval: '{interval_string}': {e}")

    return chrom, start, end



def generate_matrix(sequence):
    """Return a square 2D list where entry (i,j) == 1 if s[i]==s[j], else 0.

    Args:
        sequence (list): A DNA sequence of length N

    Returns:
        list: 2D square matrix of length NxN where each entry is either 0 or 1
    """

    return [[int(base_i == base_j) for base_i in sequence] for base_j in sequence]


def is_noise(matrix, i, j, min_diagonal_run=3):
    """Looks along both diagonals starting at (i, j) and checks whether there is a continuous stretch of cells == 1
    along at least one of the diagonals starting at (i, j). If there is a run of at least "min_diagonal_run" cells,
    returns False, otherwise returns True .

    Args:
          matrix (list): 2D square matrix of length NxN where each entry is either 0 or 1
          i (int): position i in the matrix
          j (int): position j in the matrix
          min_diagonal_run (int): threshold for min continuous cells == 1 along at least one of the diagonals starting at (i, j)

    Returns:
        bool: whether (i, j) represents noise in the dot plot matrix and so can be filtered out
    """

    matrix_size = len(matrix)
    for diagonal_direction_j in 1, -1:  # determines whether i and j move in the same direction
        run_length = 0
        for diagonal_direction in 1, -1:
            current_i = i
            current_j = j
            while current_i >= 0 and current_j >= 0 and current_i < matrix_size and current_j < matrix_size and matrix[current_i][current_j] > 0:
                current_i += diagonal_direction
                current_j += diagonal_direction * diagonal_direction_j
                run_length += 1
                if run_length >= min_diagonal_run:
                    return False

    return True


def filter_out_noise(matrix, min_diagonal_run=3, set_noise_to=0):
    """Filters the matrix in place to reduce noise

    Args:
        matrix (list): 2D square matrix of length NxN
        min_diagonal_run (int): threshold for min continuous cells == 1 along at least one of the diagonals
        set_noise_to (int): if a pixel is determined to be noise, set it to this value
    """
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if matrix[i][j] > 0 and is_noise(matrix, i, j, min_diagonal_run=min_diagonal_run):
                matrix[i][j] = set_noise_to


def plot_dot_plot(matrix, save_path=None, show=False, figure_size=None):
    """Plot the matrix as a 2D square image.

    Args:
        matrix (list): 2D square matrix of length NxN
        save_path (string): path to save the image
        show (bool): whether to show the image
        figure_size (int): size of the figure. If not specified, figure size will be set so that there's one screen pixel per matrix cell
    """
    if figure_size is None:
        dpi = 150
        figure_size = 5 * len(matrix) / dpi

    vmax = max(max(matrix[r]) for r in range(len(matrix)))
    fig, ax = plt.subplots(figsize=(figure_size, figure_size))
    ax.imshow(matrix, cmap=ListedColormap(["white", "black", "red"][:vmax + 1]), interpolation="nearest", vmin=0, vmax=vmax)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        #spine.set_visible(False)
        #spine.set_linestyle("--")
        spine.set_color("#CCCCCC")
        spine.set_linewidth(0.5)

    if show:
        plt.show()

    if save_path:
        fig.savefig(save_path, dpi=fig.get_dpi(), bbox_inches="tight", pad_inches=0)

    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot dot plot for the given input sequence.")
    parser.add_argument("--image-size", "-w", type=float, help="Image width and height in inches")
    parser.add_argument("--output-dir", "-d", default=".", help="optional output directory")
    parser.add_argument("--output-path", "-o", default="dot_plot.png", help="Output image path (PNG). Only used if a single input sequence is provided.")
    parser.add_argument("--filter-threshold", "-t", type=int, help="Minimum motif size", default=3)
    parser.add_argument("--show-filtered-pixels", action="store_true", help="Instead of hiding filtered-out pixels, show them in a different color")
    parser.add_argument("--show-plot", action="store_true", help="Show the image window before saving it to a file")
    parser.add_argument("-R", "--reference-fasta", help="Path to reference genome FASTA file")
    parser.add_argument("--verbose", "-v", action="store_true", help="Print verbose output")
    parser.add_argument("input_sequence", nargs="+", help="Nucleotide sequence(s) to plot, or BED file path(s), or interval(s) like chrom:start0based-end.")
    args = parser.parse_args()

    input_sequences = []
    output_filenames = []
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    for i, seq in enumerate(args.input_sequence):
        invalid_chars = set(seq) - set("ACGT")
        if not invalid_chars:
            # this is a literal nucleotide sequence

            input_sequences.append(seq)
            if len(args.input_sequence) == 1:
                output_filenames.append(os.path.join(args.output_dir, args.output_path))
            else:
                output_filenames.append(os.path.join(args.output_dir, f"dot_plot_{i+1:03d}_of_{len(args.input_sequence)}.{len(seq)}bp_sequence.png"))

        else:
            if not "bed" in seq and not ":" in seq and not "-" in seq:
                parser.error(f"Error: {seq} is not a valid nucleotide sequence, BED file path, or interval")

            if not args.reference_fasta:
                parser.error("Error: --reference-fasta is required when the input is a BED file or interval")

            intervals = []
            if "bed" in seq:
                bed_path = seq
                if not os.path.isfile(bed_path):
                    parser.error(f"Error: {bed_path} file not found")
                            
                fopen = gzip.open if bed_path.endswith("gz") else open
                with fopen(bed_path, "rt") as bed_file:
                    for line_i, line in enumerate(bed_file):
                        fields = line.strip().split("\t")
                        if len(fields) < 3:
                            parser.error(f"Error: {bed_path} line #{line_i+1} is invalid: '{line.strip()}'")
                        
                        chrom, start_0based, end = fields[:3]
                        intervals.append((chrom, int(start_0based), int(end)))
            elif ":" in seq and "-" in seq:
                chrom, start_0based, end = parse_interval(seq)
                intervals.append((chrom, int(start_0based), int(end)))
            else:
                parser.error(f"Error: {seq} is not a valid interval or BED file")

            fasta_entries = pyfaidx.Fasta(args.reference_fasta, as_raw=True, one_based_attributes=False, sequence_always_upper=True)
            for i, (chrom, start_0based, end) in enumerate(intervals):
                chrom = chrom.replace("chr", "")
                seq = fasta_entries[f"chr{chrom}"][start_0based:end]                
                input_sequences.append(seq)
                output_filenames.append(os.path.join(args.output_dir, f"dot_plot_{i+1:03d}_of_{len(intervals)}.chr{chrom}_{start_0based}-{end}.{len(seq)}bp_sequence.png"))

    print(f"Loaded {len(intervals):,d} interval(s) from {args.reference_fasta}")

    for i, (seq, output_filename) in enumerate(zip(input_sequences, output_filenames)):
        if args.verbose:
            if len(input_sequences) > 1:
                print(f"Generating dot plot for {len(seq)}bp sequence {i+1} of {len(input_sequences)}")
            else:
                print(f"Generating dot plot for {len(seq)}bp sequence")

        matrix = generate_matrix(seq)

        filter_out_noise(matrix, set_noise_to=2 if args.show_filtered_pixels else 0, min_diagonal_run=args.filter_threshold)

        plot_dot_plot(matrix, save_path=output_filename, show=args.show_plot, figure_size=args.image_size)

        if args.verbose or len(input_sequences) == 1:
            print(f"Saved image to {output_filename} (sequence length {len(seq)}bp)")

    print(f"Done generating {len(input_sequences):,d} dot plot(s)", f"in {args.output_dir} directory" if args.output_dir != "." else "")


