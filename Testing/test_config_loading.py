import unittest
import os
import demographiKS
class TestConfigLoading(unittest.TestCase):
    def test_config_loading(self):

        config_folder="../sample_configs"
        in_file_name = os.path.join(config_folder, "test-case.xml")
        mock_arguments=['cmd',in_file_name]
        conf = demographiKS.setup(mock_arguments)

        self.assertEqual(conf.bottleneck_Ne, 10000)  # add assertion here
        self.assertEqual(conf.ancestral_Ne, 100000)

        self.assertEqual(conf.DIV_time_Ge, 50)  # add assertion here
        self.assertEqual(conf.WGD_time_Ge, 10)

        self.assertAlmostEqual(conf.mean_WGD_life_span_in_GE, 5770.78,2)  # add assertion here
        self.assertEqual(conf.recombination_rate, 1.26e-06)
        self.assertEqual(conf.mutation_rate, 0.001)
        
        self.assertEqual(conf.burnin_time, 478)

        self.assertEqual(conf.mig_start, 5)
        self.assertEqual(conf.mig_stop, 10)
        self.assertEqual(conf.mig_rate, 0.5)

        self.assertEqual(conf.homoeologous_exchange_rate, 0.123)

        self.assertEqual(conf.sample_ancestral_genomes_for_Tc[0],[1,10])
        self.assertEqual(conf.sample_ancestral_genomes_for_Tc[1],[5,7])

        self.assertEqual(conf.DemographiKS_random_seed, 7)
        self.assertEqual(conf.Msprime_random_seed, 8)
        self.assertEqual(conf.SLiM_rep, 9)

        print("mig rate=" + str(conf.mig_rate))
        print("dij rate=" + str(conf.homoeologous_exchange_rate))

if __name__ == '__main__':
    unittest.main()
