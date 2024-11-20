import unittest
import os
import demographiKS
class TestConfigLoading(unittest.TestCase):
    def test_config_loading(self):

        config_folder="../sample_configs"
        in_file_name = os.path.join(config_folder, "short-run.xml")
        mock_arguments=['cmd',in_file_name]
        conf = demographiKS.setup(mock_arguments)

        self.assertEqual(conf.bottleneck_Ne, 100)  # add assertion here
        self.assertEqual(conf.ancestral_Ne, 1000)

        self.assertEqual(conf.DIV_time_Ge, 75)  # add assertion here
        self.assertEqual(conf.WGD_time_Ge, 25)

if __name__ == '__main__':
    unittest.main()
