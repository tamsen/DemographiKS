import os.path
import unittest
#import process_wrapper
import config
import process_wrapper


class MyDataFetcher(unittest.TestCase):

    def test_fetch_data(self):

        RC_run_list=['RC10_m12d18y2024_h14m35s15',
                     'RC10_m12d18y2024_h14m35s17',
                     'RC10_m12d18y2024_h13m36s12','RC10_m12d18y2024_h14m35s19',
                     'RC10_m12d18y2024_h14m30s54','RC10_m12d18y2024_h14m39s19',
                     'RC10_m12d18y2024_h14m30s58']

        me_at_remote_URL =  'tdunn@mesx.sdsu.edu'
        output_root_folder="/usr/scratch2/userdata2/tdunn/DemographiKS_output/RC/"

        #run_name="RC10_m12d18y2024_h11m43s13"
        #TODO, add fetching the log, so I can get the time stamps
        #for i in range(3,4):
        run_name = "RC10_m12d18y2024_h14m35s19"#RC_run_list[2]
        local_output_folder = os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
                                               run_name)
        self.pull_down_run_data(local_output_folder, me_at_remote_URL, output_root_folder, run_name)

        self.assertEqual(True, True)  # add assertion here

    def pull_down_run_data(self, local_output_folder, me_at_remote_URL, output_root_folder, run_name):

        run_specific_folder=run_name + "/demographiKS_output"

        if not os.path.exists(local_output_folder):
            os.makedirs(local_output_folder)

        file_needed=["*.csv","*.png", "*.used.xml"]
        for file in file_needed:
            remote_file = os.path.join(output_root_folder, run_specific_folder, file)
            cmd1 = ['scp', '-r', me_at_remote_URL + ':' + remote_file, local_output_folder]
            print(" ".join(cmd1))
            out_string, error_string = process_wrapper.run_and_wait_with_retry(
                cmd1, local_output_folder, "Connection reset by peer", 2, 5)

        remote_log_file = os.path.join(output_root_folder, run_name, '*log*')
        cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_log_file, local_output_folder]
        print(" ".join(cmd2))
        out_string, error_string = process_wrapper.run_and_wait_with_retry(
            cmd2, local_output_folder, "Connection reset by peer", 2, 5)


if __name__ == '__main__':
    unittest.main()
