import os.path
import unittest
import process_wrapper


class MyDGXDataFetcher(unittest.TestCase):

    def test_fetch_demographics_data(self):

        run_list=['DGKS_5000_5000_m5_BI_40K_RC7_m01d18y2025_h08m07s39',
                        'DGKS_5000_5000_m5_BI_40K_RC7_m01d18y2025_h08m07s39_test2']
        run_collection_name="GE"

        #me_at_remote_URL =  'tdunn@mesx.sdsu.edu'
        me_at_remote_URL =  'mesx_cluster'

        if 'DGKS' in run_list[0]:
            output_root_folder = os.path.join("/usr/scratch2/userdata2/tdunn/DemographiKS_Output")
        else:
            output_root_folder=os.path.join("/usr/scratch2/userdata2/tdunn/DemographiKS_output",
                                        run_collection_name)



        for i in range(0,len(run_list)):
            run_name = run_list[i]
            local_output_folder = os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
                                               run_name)
            self.pull_down_run_data(local_output_folder, me_at_remote_URL, output_root_folder, run_name)

        self.assertEqual(True, True)  # add assertion here

    def pull_down_run_data(self, local_output_folder, me_at_remote_URL, output_root_folder, run_name):

        demographics_folder=run_name + "/demographiKS_output"

        if not os.path.exists(local_output_folder):
            os.makedirs(local_output_folder)

        file_needed_from_demographics_folder=["*.csv","*.png", "*.used.xml"]
        file_needed_from_run_folder=["*log*", "*.used.xml","*.png","*.slim"]

        for file in file_needed_from_demographics_folder:
            remote_file = os.path.join(output_root_folder, demographics_folder, file)
            cmd1 = ['scp', '-r', me_at_remote_URL + ':' + remote_file, local_output_folder]
            print(" ".join(cmd1))
            out_string, error_string = process_wrapper.run_and_wait_with_retry(
                cmd1, local_output_folder, "Connection reset by peer", 2, 5)

        for file in file_needed_from_run_folder:
            remote_file = os.path.join(output_root_folder, run_name, file )
            cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_file, local_output_folder]
            print(" ".join(cmd2))
            out_string, error_string = process_wrapper.run_and_wait_with_retry(
                cmd2, local_output_folder, "Connection reset by peer", 2, 5)


if __name__ == '__main__':
    unittest.main()
