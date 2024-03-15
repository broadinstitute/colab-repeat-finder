"""RepeatTracker class for detecting both pure and interrupted repeats of a given motif size in a sequence."""

class RepeatTracker:
	"""This class tracks repeats of a single motif size in a given input sequence. It outputs all perfect repeats of
	this motif size that pass filter criteria while scanning the input sequence from left to right.
	If the max_interruptions parameter is > 0 then it will also look for interrupted repeats where a fixed number of
	positions within the motif can vary across repeats. In this mode, the algorithm traverses left to right, but
	may occasionally jump backward to revisit some previously-evaluated positions."""

	def __init__(self,
				 motif_size,
				 min_repeats,
				 min_span,
				 max_interruptions,
				 input_sequence,
				 output_intervals,
				 allow_multibase_homopolymer_motifs=False,
				 allow_ending_with_different_motif=False,
				 verbose=False,
				 debug=False):
		"""Initialize a RepeatTracker object.

		Args:
			motif_size (int): The motif size (in base pairs) that will be tracked by this RepeatTracker.
			min_repeats (int): Only add repeats to the output when there's at least this many repeats in a row.
			min_span (int):  Only add repeats to the output when they span at least this many base pairs.
			max_interruptions (int): How many bases within a motif are allowed to vary across repeats.
			allow_multibase_homopolymer_motifs (bool): When the motif size is greater than 1, allow motifs like 'AA' or 'TTTT'.
			allow_ending_with_different_motif (bool): By default, interrupted repeats must still end with at least 1 copy of the exact same motif that they started with.
			input_sequence (str): The input sequence.
			output_intervals (dict): A dictionary to store detected repeats. The key is (start_0based, end) and the
				value is the detected motif.
			verbose (bool): Print verbose output for debugging.
		"""

		self.motif_size = motif_size
		self.min_repeats = min_repeats
		self.min_span = min_span
		self.max_interruptions = max_interruptions
		self.allow_multibase_homopolymer_motifs = allow_multibase_homopolymer_motifs
		self.allow_ending_with_different_motif = allow_ending_with_different_motif
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
			  f"{len(self.input_sequence):,d}bp  [{start_0based}:{self.current_position+1}], "
			  f"run={self.run_length}, "
			  f"i0={self.current_position}: "
			  f"{(self.current_position - start_0based)/len(motif):0.2f} x {motif} "
			  f"==> {self.input_sequence[start_0based : self.current_position+1]}")

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
		while self.current_position < len(seq) and (
			seq[self.current_position] == seq[self.current_position - self.motif_size]
			or (self.run_length % self.motif_size) in self.current_interrupted_positions_in_motif
		):
			if self.debug_strings:
				self.debug_strings[-1] += seq[self.current_position]
			self.run_length += 1
			self.current_position += 1
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
				final_motif_bases = {b for b in final_motif if b != "N"}
				if self.motif_size > 1 and len(final_motif_bases) == 1 and not self.allow_multibase_homopolymer_motifs:
					if self.verbose: print(f"==> No, the motif is equivalent to a homopolymer: {final_motif}.")
					return

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

