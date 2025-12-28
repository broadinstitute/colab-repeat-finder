"""Plot the periodicity matrix for a given input sequence.

This script generates a visualization that shows the periodicity patterns in a DNA sequence.
The periodicity matrix displays which positions in the sequence match at different motif sizes,
helping to identify tandem repeat structures.
"""

import argparse
from utils.plot_utils import get_period_matrix, plot_periodicity_matrix

parser = argparse.ArgumentParser()
parser.add_argument("input_sequence", help="The input sequence.")
parser.add_argument("--min-motif-size", default=1, type=int, help="The minimum motif size in base pairs.")
parser.add_argument("--max-motif-size", default=50, type=int, help="The maximum motif size in base pairs.")
parser.add_argument("-o", "--output-path", default="periodicity_matrix.png", help="The output path for the periodicity matrix.")
args = parser.parse_args()

periodicity_matrix = get_period_matrix(args.min_motif_size, args.max_motif_size, args.input_sequence)
plot_periodicity_matrix(periodicity_matrix, args.output_path)

