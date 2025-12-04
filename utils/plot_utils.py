import zlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


def shift_string_by(string, shift):
	"""Shift a string to the right by a given number of characters (eg. "AGTTT" shifted by 2 becomes "TTAGT")."""

	return string[-shift:] + string[:-shift]


def get_period_matrix(min_motif_size, max_motif_size, input_sequence):
	min_motif_size = max(min_motif_size, 1)
	max_motif_size = min(max_motif_size, len(input_sequence) // 2)
	matrix = [[0 for _ in range(len(input_sequence))] for _ in range(max_motif_size)]

	for row in range(min_motif_size - 1, max_motif_size):
		period = row + 1
		for column in range(len(input_sequence) - period):
			if input_sequence[column] == input_sequence[column + period]:
				# set the matrix value to a hash that stays the same as long as the motif is the same
				value = abs(hash(bytes(shift_string_by(input_sequence[column:column + period], column % period), encoding="utf-8")))
				matrix[row][column] = value

	return matrix



def plot_periodicity_matrix(periodicity_matrix, output_path):
	
	matrix_width = len(periodicity_matrix[0])
	fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10, 8))
	cmap = ListedColormap(['#F3F3F3'] + list(plt.get_cmap("Pastel2", 12).colors))
	ax1.matshow(periodicity_matrix, aspect="auto", cmap=cmap)
	#ax1.set_xticks(range(0, matrix_width + 1, 10))
	#ax1.set_xticks(range(5, matrix_width + 1, 5), minor=True)
	#ax1.set_xticklabels([str(i + 1) for i in range(10, matrix_width, 10)])
	#ax1.set_xticks(range(matrix_width), minor=True)
	#ax1.set_yticks(range(len(periodicity_matrix)))
	#ax1.set_yticklabels([str(i + 1) for i in range(len(periodicity_matrix))])
	ax1.set_xlabel("Input sequence position")
	ax1.set_ylabel("Period")

	scores = [sum([1 for _ in row if _ > 0]) for row in periodicity_matrix]
	periods = list(range(1, len(scores) + 1))
	scores = [scores[i]/(matrix_width - periods[i] + 1) for i in range(len(scores))]
	ax2.bar(periods, scores)
	ax2.set_xticks(periods)
	ax2.set_ylabel("Fraction of matches")
	ax2.set_xlabel("Period")

	plt.savefig(output_path)
	print(f"Wrote {output_path}")


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
			#matrix[row][i] = zlib.adler32(bytes(motif, encoding="utf-8")) % 10 + 1
			matrix[row][i] = abs(hash(bytes(motif, encoding="utf-8")) % 10 + 1)

	plot_periodicity_matrix(matrix, output_path)
