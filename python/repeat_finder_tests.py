import argparse
import unittest
from repeat_finder import detect_repeats
from utils.plot_utils import shift_string_by


class RepeatFinderTests(unittest.TestCase):

	def test_shift_string_by(self):
		self.assertEqual(shift_string_by("A", 1), "A")
		self.assertEqual(shift_string_by("TTTCG", 0), "TTTCG")
		self.assertEqual(shift_string_by("TTTCG", 1), "GTTTC")
		self.assertEqual(shift_string_by("TTTCG", 2), "CGTTT")
		self.assertEqual(shift_string_by("TTTCG", 3), "TCGTT")
		self.assertEqual(shift_string_by("TTTCG", 4), "TTCGT")
		self.assertEqual(shift_string_by("TTTCG", 5), "TTTCG")


	def test_detect_repeats(self):
		filter_settings = argparse.Namespace(
			min_motif_size=1,
			max_motif_size=7,
			min_repeats=3,
			max_interruptions_by_motif_size={i: 0 for i in range(1, 8)},
			min_span=6,
			verbose=False,
			debug=False)

		for motif in "A", "CA", "CAG", "CAGA", "CAGAT", "CAGATT", "CAGATTA", "CAGATTAG":
			seq = 6*motif
			repeats = detect_repeats(seq, filter_settings)
			if len(motif) <= 7:
				self.assertEqual(repeats, [(0, len(seq), motif)], f"Error on sequence: {seq}")
			else:
				self.assertEqual(repeats, [], f"Error on sequence:= {seq}. Got {repeats}")

		for motif in "A", "CA", "CAG", "CAGA", "CAGAT", "CAGATT", "CAGATAT", "CAGATTAG":
			seq = 6*motif + 10*"TA"
			repeats = detect_repeats(seq, filter_settings)
			expected_TA_repeat_result = (6*len(motif) - 1, len(seq), "AT") if motif.endswith("A") else (6*len(motif), len(seq), "TA")
			if len(motif) <= 7:
				self.assertEqual(repeats, [
					(0, 6*len(motif), motif),
					expected_TA_repeat_result,
				], f"Error on sequence: {seq}. Got {repeats}")
			else:
				self.assertEqual(repeats, [
					expected_TA_repeat_result,
				], f"Error on sequence: {seq}. Got {repeats}")

		seq = "A"*9 + "C"*11 + "G"*10 + "T"*9
		repeats = detect_repeats(seq, filter_settings)
		self.assertEqual(repeats, [(0, 9, "A"), (9, 20, "C"), (20, 30, "G"), (30, 39, "T")], f"Error on sequence: {seq}. Got {repeats}")

		seq = "CA"*9 + "GT"*11 + "CG"*10 + "TA"*9
		repeats = detect_repeats(seq, filter_settings)
		self.assertEqual(repeats, [(0, 18, "CA"), (18, 40, "GT"), (40, 60, "CG"), (60, 78, "TA")], f"Error on sequence: {seq}. Got {repeats}")

		motif = "A"
		filter_settings.min_motif_size = 2
		repeats = detect_repeats(7*motif, filter_settings)
		self.assertEqual(repeats, [(0, 7, "AA")])

		motif1 = "A"
		motif2 = "AGAC"
		motif3 = "AAC"
		seq = 7*motif1 + 7*motif2 + 10*motif3
		self.assertEqual(len(seq), 65)

		filter_settings.min_motif_size = 2
		filter_settings.max_motif_size = 10
		repeats = detect_repeats(seq, filter_settings)
		self.assertEqual(repeats, [(0, 8, "AA"), (7, 36, "AGAC"), (33, 65, "ACA")], f"Error on sequence: {seq}. Got {repeats}")

		filter_settings.min_motif_size = 1
		repeats = detect_repeats(seq, filter_settings)
		self.assertEqual(repeats, [(0, 8, "A"), (7, 36, "AGAC"), (33, 65, "ACA")], f"Error on sequence: {seq}. Got {repeats}")

