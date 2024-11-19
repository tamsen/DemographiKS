import os.path
import unittest
import process_wrapper


class MyDataFetcher(unittest.TestCase):

    def test_fetch_data(self):


        #///usr/scratch2/userdata2/tdunn/DemographiKS_Output
        #//home/tamsen/Data/DemographiKS/demographiKS_output/allotetraploid_bottleneck_hist.png
        #

        output_root_folder="/usr/scratch2/userdata2/tdunn/DemographiKS_Output"
        run_name="DGKS_1000_100_m08d29y2024_h15m27s59"
        #DGKS_1000_1000_m08d29y2024_h16m08s09
        #DGKS_100_100_m08d29y2024_h16m09s38
        #DGKS_1000_100_m08d29y2024_h15m27s59
        #DGKS_10_10_m08d29y2024_h16m13s27
        #DGKS_1000_10_m08d29y2024_h16m12s08

        #run_name_splat=run_name.split("_")
        run_specific_folder=run_name + "/demographiKS_output"
        #nickname="_".join(run_name_splat[0],run_name_splat[1],run_name_splat[2])
        ks_data_file="allotetraploid_bottleneck.csv"
        config_file="*.used.xml"

        me_at_remote_URL =  'tdunn@mesx.sdsu.edu'
        local_output_folder = os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
                                           run_name)

        if not os.path.exists(local_output_folder):
            os.makedirs(local_output_folder)

        remote_csv_file = os.path.join(output_root_folder,run_specific_folder,ks_data_file)
        remote_config_file = os.path.join(output_root_folder,run_name,config_file)

        cmd1 = ['scp', '-r', me_at_remote_URL + ':' + remote_config_file, local_output_folder]
        print(" ".join(cmd1))
        out_string, error_string = process_wrapper.run_and_wait_with_retry(
            cmd1, local_output_folder, "Connection reset by peer", 2, 5)

        cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_csv_file, local_output_folder]
        print(" ".join(cmd2))
        out_string, error_string = process_wrapper.run_and_wait_with_retry(
            cmd2, local_output_folder, "Connection reset by peer", 2, 5)

        expected_output_file=os.path.join(local_output_folder,ks_data_file)
        self.assertEqual(os.path.exists(expected_output_file), True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
