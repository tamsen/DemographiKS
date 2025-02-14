import unittest
import modules.summary_stats as summary_stats


class TestSummaryStates(unittest.TestCase):

    def test_sequence_identity(self):

        sequence_a='AAACCCTTT'
        sequence_b='AAACCCTTT'
        result = summary_stats.compute_sequence_identity(sequence_a, sequence_b)
        self.assertEqual(result, 1)  # add assertion here

        sequence_b='CCCTTTGGG'
        result = summary_stats.compute_sequence_identity(sequence_a, sequence_b)
        self.assertEqual(result, 0)  # ad

        sequence_b = 'AAACCCGGG'
        result = summary_stats.compute_sequence_identity(sequence_a, sequence_b)
        self.assertEqual(result, (2.0/3.0))  # ad

    def test_compute_num_sequence_diffs(self):

        sequence_a='AAACCCTTT'
        sequence_b='AAACCCTTT'
        result = summary_stats.compute_num_sequence_diffs(sequence_a, sequence_b)
        self.assertEqual(result, 0)  # add assertion here

        sequence_b='CCCTTTGGG'
        result = summary_stats.compute_num_sequence_diffs(sequence_a, sequence_b)
        self.assertEqual(result, 9)  # ad

        sequence_b = 'AAACCCGGG'
        result = summary_stats.compute_num_sequence_diffs(sequence_a, sequence_b)
        self.assertEqual(result, 3)  # ad


    def test_nucleotide_diverisity(self):

        sequence_a ='AAACCCTTT'
        sequence_b ='CCCTTTGGG'
        sequence_c = 'AAACCCGGG'
        seq_list1 =[sequence_a,sequence_b]
        pi = summary_stats.pi_nuc_diversity(seq_list1)
        self.assertEqual(pi, 1)

        seq_list2 =[sequence_c, sequence_c]
        pi = summary_stats.pi_nuc_diversity(seq_list2)
        self.assertEqual(pi, 0)  # add assertion here

        seq_list3 =["AACC", "AAGG"]
        pi = summary_stats.pi_nuc_diversity(seq_list3)
        self.assertEqual(pi, 0.5)

        seq_list3 =["AACC", "AAGG", "AAGG"]
        pi = summary_stats.pi_nuc_diversity(seq_list3)
        #2,2,0, so a total of 4 difference out of 4*3=12 comparisons = 4 / 12 = 1/3
        self.assertEqual(pi, (1/3))

if __name__ == '__main__':
    unittest.main()
