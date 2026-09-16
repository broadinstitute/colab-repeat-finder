"""Original code for the RepeatTracker class before it was refactored to allow detection of interrupted repeats"""

class PerfectRepeatTracker:
	"""This class tracks repeats of a single motif size in a given input sequence. It outputs all perfect repeats
	that pass filter criteria while scanning the input sequence from left to right"""

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

		self._current_position = 0   # 0-based position in the input sequence
		self._run_length = 1  # number of bases added so far to the current repeat interval

	def advance(self):
		"""Increment current position within the input sequence while updating internal state and recording any detected repeats"""

		seq = self.input_sequence
		i = self._current_position
		period = self.motif_size

		if i >= len(seq) - period:
			return False

		if seq[i] == seq[i + period] and seq[i] != "N":
			self._current_position += 1
			self._run_length += 1
			return True

		self.output_interval_if_it_passes_filters()
		self._current_position += 1
		self._run_length = 1
		return True

	def current_run_start_0based(self):
		"""Return the 0-based start position of the run that the tracker is currently extending. Any repeat this
		tracker outputs in the future will start at or after this position."""
		return self._current_position - self._run_length + 1

	def done(self):
		"""Output the last interval if it passes filters"""
		self.output_interval_if_it_passes_filters()

	def output_interval_if_it_passes_filters(self):
		"""Check internal state to see if enough repeats have accumulated to output an interval. If yes, add it to
		the output.
		"""

		seq = self.input_sequence
		i = self._current_position
		period = self.motif_size

		# found mismatch, check if previous run is worth outputing
		start_0based = i - self._run_length + 1
		motif = seq[start_0based:start_0based + period]
		if "N" in motif:
			return

		if self._run_length + period - 1 >= self.min_span and self._run_length + period - 1 >= self.min_repeats * period:
			while i < len(seq) - 1 and i + 1 - period >= 0 and seq[i+1] == seq[i+1 - period]:
				self._run_length += 1
				i += 1

		if self._run_length >= self.min_span and self._run_length >= self.min_repeats * period:
			end = i + 1

			previous_motif = self.output_intervals.get((start_0based, end))
			if previous_motif is not None and len(motif) > len(previous_motif):
				return

			if consists_of_perfect_repeats(motif):
				return

			self.output_intervals[(start_0based, end)] = motif

	@property
	def current_position(self):
		return self._current_position


def consists_of_perfect_repeats(sequence):
	"""Check whether the given nucleotide sequence consists entirely of repeats of some smaller repeat unit
    (eg. CAGCAGCAGCAGCAG = 5xCAG). If yes, return the shorter repeat unit. Otherwise, return None.

    Args:
        sequence (str): nucleotide sequence

    Return:
         the repeat unit or None if the input sequence does not consist of perfect repeats of a smaller repeat unit
    """

	# quickly check if the sequence consists of exact repeats of a 1bp, 2bp, or 3bp motif
	if len(sequence) >= 2 and sequence == sequence[0]*len(sequence):
		return sequence[0]
	elif len(sequence) >= 4 and len(sequence) % 2 == 0 and sequence == sequence[:2]*(len(sequence)//2):
		return sequence[:2]
	elif len(sequence) >= 6 and len(sequence) % 3 == 0 and sequence == sequence[:3]*(len(sequence)//3):
		return sequence[:3]

	# find the smallest repeat unit that covers the entire sequence
	repeat_unit_length = 4
	while repeat_unit_length <= len(sequence)/2:
		if len(sequence) % repeat_unit_length != 0:
			repeat_unit_length += 1
			continue

		repeat_unit = sequence[:repeat_unit_length]
		num_repeats = sequence.count(repeat_unit)
		if num_repeats * repeat_unit_length == len(sequence):
			return repeat_unit

		repeat_unit_length += 1

	# no repeat unit found in sequence
	return None
