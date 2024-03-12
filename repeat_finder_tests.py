import argparse
import unittest
from repeat_finder import shift_string_by, detect_repeats


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
			min_span=6)

		for motif in "A", "CA", "CAG", "CAGA", "CAGAT", "CAGATT", "CAGATTA", "CAGATTAG":
			repeats = detect_repeats(6*motif, filter_settings)
			if len(motif) <= 7:
				self.assertEqual(repeats, [(0, 6*len(motif), motif)])
			else:
				self.assertEqual(repeats, [])

		motif = "A"
		filter_settings.min_motif_size = 2
		repeats = detect_repeats(7*motif, filter_settings)
		self.assertEqual(repeats, [(0, 7, "AA")])

		motif1 = "A"
		motif2 = "AAC"
		filter_settings.min_motif_size = 2
		repeats = detect_repeats(7*motif1 + 10*motif2, filter_settings)
		self.assertEqual(repeats, [(0, 9, "AA"), (7, 37, "AAC")])
