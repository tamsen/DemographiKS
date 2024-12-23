import os.path
import unittest
#import process_wrapper
import config
import process_wrapper


class MyDataFetcher(unittest.TestCase):

    def test_fetch_data(self):

        TE_run_list=['TE01_m12d20y2024_h14m16s21','TE03_m12d20y2024_h14m26s56',
                     'TE02_m12d20y2024_h14m22s23']#,'TE04_m12d20y2024_h14m31s31


        RC_run_list=['RC10_m12d18y2024_h14m35s15','RC10_m12d18y2024_h14m35s19']
        #             'RC10_m12d18y2024_h14m35s17',
        #             'RC10_m12d18y2024_h13m36s12','RC10_m12d18y2024_h14m35s19',
        #             'RC10_m12d18y2024_h14m30s54','RC10_m12d18y2024_h14m39s19',
        #             'RC10_m12d18y2024_h14m30s58']

        BI_run_list=["BI2_m12d19y2024_h10m57s18","BI4_m12d19y2024_h11m02s24",
                    "BI7_m12d19y2024_h11m02s31","BI3_m12d19y2024_h11m02s23",
                     "BI5_m12d19y2024_h11m02s27"]


        Ne_run_list=["Ne3_m12d18y2024_h16m08s40"]
        #"Ne1_m12d18y2024_h16m08s35","Ne5_m12d18y2024_h16m10s54",
        #             "Ne2_m12d18y2024_h16m08s38",
        #             "Ne3_m12d18y2024_h16m08s40",
        #             "Ne4_m12d18y2024_h16m10s52"]
        #"Ne6_m12d18y2024_h16m10s57",,"Ne6_m12d19y2024_h11m14s34",
        # "Ne7_m12d19y2024_h11m14s41","Ne7_m12d18y2024_h16m11s01"

        GE_run_list=["GE4_m12d19y2024_h11m47s58","GE6_m12d19y2024_h11m48s02",
                     "GE8_m12d19y2024_h13m34s13","GE5_m12d19y2024_h11m47s58",
                     "GE7_m12d19y2024_h13m30s32"]

        run_list=TE_run_list
        run_collection_name="TE"
        me_at_remote_URL =  'tdunn@mesx.sdsu.edu'
        output_root_folder=os.path.join("/usr/scratch2/userdata2/tdunn/DemographiKS_output",
                                        run_collection_name)

        #run_name="RC10_m12d18y2024_h11m43s13"
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
