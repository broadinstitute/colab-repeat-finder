"""Original code for the RepeatTracker class before it was refactored to allow detection of interrupted repeats"""

class PureRepeatTracker:
	"""This class tracks repeats of a single motif size in a given input sequence. It outputs all perfect repeats of
	this that pass filter criteria while scanning the input sequence from left to right"""

	def __init__(self, motif_size, min_repeats, min_span, input_sequence, output_intervals):
		"""Initialize a RepeatTracker object.

		Args:
			motif_size (int): The motif size (in base pairs) that will be tracked by this RepeatTracker.
			min_repeats (int): Only add repeats to the output when there's at least this many repeats in a row.
			min_span (int):  Only add repeats to the output when they span at least this many base pairs.
			input_sequence (str): The input sequence.
			output_intervals (dict): A dictionary to store detected repeats. The key is (start_0based, end) and the
				value is the detected motif.
		"""

		self.motif_size = motif_size
		self.min_repeats = min_repeats
		self.min_span = min_span
		self.input_sequence = input_sequence
		self.output_intervals = output_intervals

		self.current_position = 0
		self.run_length = 1

	def advance(self):
		"""Increment current position within the input sequence while updating internal state and recording any detected repeats"""

		seq = self.input_sequence
		i = self.current_position
		period = self.motif_size

		if i >= len(seq) - period:
			return False

		if seq[i] == seq[i + period]:
			self.current_position += 1
			self.run_length += 1
			return True

		self.output_interval_if_it_passes_filters()
		self.current_position += 1
		self.run_length = 1
		return True

	def done(self):
		"""Output the last interval if it passes filters"""
		self.output_interval_if_it_passes_filters()

	def output_interval_if_it_passes_filters(self):
		"""Check internal state to see if enough repeats have accumulated to output an interval. If yes, add it to
		the output.
		"""

		seq = self.input_sequence
		i = self.current_position
		period = self.motif_size

		# found mismatch, check if previous run is worth outputing
		start_0based = i - self.run_length + 1
		motif = seq[start_0based:start_0based + period]
		if "N" in motif:
			return

		if self.run_length + period - 1 >= self.min_span and self.run_length + period - 1 >= self.min_repeats * period:
			while i < len(seq) - 1 and seq[i+1] == seq[i+1 - period]:
				self.run_length += 1
				i += 1

		if self.run_length >= self.min_span and self.run_length >= self.min_repeats * period:
			end = i + 1

			previous_motif = self.output_intervals.get((start_0based, end))
			if previous_motif is None or len(motif) < len(previous_motif):
				self.output_intervals[(start_0based, end)] = motif

