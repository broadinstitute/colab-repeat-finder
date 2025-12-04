"""Plot the periodicity matrix for a given input sequence."""

import argparse
from utils.plot_utils import get_period_matrix, plot_periodicity_matrix

p = argparse.ArgumentParser()
p.add_argument("input_sequence", help="The input sequence.")
p.add_argument("--min-motif-size", default=1, type=int, help="The minimum motif size in base pairs.")
p.add_argument("--max-motif-size", default=50, type=int, help="The maximum motif size in base pairs.")
p.add_argument("-o", "--output-path", default="periodicity_matrix.png", help="The output path for the periodicity matrix.")
args = p.parse_args()

periodicity_matrix = get_period_matrix(args.min_motif_size, args.max_motif_size, args.input_sequence)
plot_periodicity_matrix(periodicity_matrix, args.output_path)

