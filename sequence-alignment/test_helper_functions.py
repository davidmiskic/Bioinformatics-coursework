import unittest

from helper_functions import local_alignment


class TestLocalAlignment(unittest.TestCase):
    def test_local_alignment(self):
        seq1 = "ATATTAGGTTTTTACCTACCCAGGAAAAGCCAACCAACCTCGATCTCTTGTAGATCTGTTC"
        seq2 = "TATCTCACTTCGTTCTCT"

        def scoring_function(x, y):
            if x == "-" or y == "-":
                return -2
            elif x != y:
                return -1
            else:
                return 2

        valid_alignments = [
            ("CAACCAACCTCGATCTCT", "TATCTCACTTCGTTCTCT"),
            ("CAAC-CAACCTCGATCTCT", "TATCTC-ACTTCGTTCTCT"),
            ("CAAC-CAACCTCGATCTCT", "TATCTCA-CTTCGTTCTCT"),
        ]

        align1, align2, score = local_alignment(seq1, seq2, scoring_function)
        self.assertEqual(score, 19)
        self.assertIn((align1, align2), valid_alignments)


if __name__ == "__main__":
    unittest.main()
