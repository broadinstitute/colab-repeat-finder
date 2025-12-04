"""
This script takes a DNA sequence of length N and generates an NxN square image - aka. dot plot - where
pixel (i,j) is dark only if positions i and j have the same base

Additional details:
- to reduce noise, pixels are only shown if they are part of a diagonal stretch of dark pixels longer than some threshold
- the image can optionally be downsampled to smaller width x height
"""

import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

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
        print(f"Saved image to {save_path}")

    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot dot plot for the given input sequence.")
    parser.add_argument("--image-size", "-w", type=float, help="Image width and height in inches")
    parser.add_argument("--output-path", "-o", default="dot_plot.png", help="Output image path (PNG)")
    parser.add_argument("--filter-threshold", "-t", type=int, help="Minimum motif size", default=3)
    parser.add_argument("--show-filtered", action="store_true", help="Instead of hiding filtered-out pixels, show them in a different color")
    parser.add_argument("--show-plot", action="store_true", help="Show the image window")
    parser.add_argument("input_sequence", help="Sequence to plot")
    args = parser.parse_args()

    seq = args.input_sequence
    invalid_chars = set(seq) - set("ACGT")
    if invalid_chars:
        parser.error(f"Error: input sequence contains invalid characters {', '.join(sorted(invalid_chars))}; only A,C,G,T allowed.")

    print(f"Generating dot plot for {len(seq)}bp sequence")
    matrix = generate_matrix(seq)

    filter_out_noise(matrix, set_noise_to=2 if args.show_filtered else 0, min_diagonal_run=args.filter_threshold)

    plot_dot_plot(matrix, save_path=args.output_path, show=args.show_plot, figure_size=args.image_size)

