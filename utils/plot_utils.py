import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


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

