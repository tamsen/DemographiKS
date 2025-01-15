import os.path
import unittest
#import process_wrapper
import config
import process_wrapper


class MySPXDataFetcher(unittest.TestCase):

    def test_fetch_specks_data(self):

        #/ usr / scratch2 / userdata2 / tdunn / DemographiKS_output / SP / SP1MY_m12d24y2024_h09m00s34
        TE5_run_list=['specks_TE05_m12d30y2024_h11m50s03',
                      'specks_TE08_m12d30y2024_h12m10s13',
                      'specks_TE07_m12d30y2024_h12m10s15','specks_TE09_m12d30y2024_h12m10s11']
        TE5_run_list = ['specks_TE05_m12d31y2024_h09m10s39','specks_TE07_m12d31y2024_h09m10s28',
                      'specks_TE08_m12d31y2024_h09m10s32',
                      'specks_TE09_m12d31y2024_h09m10s34']
        TE5_run_list = ['specks_TE10_m12d31y2024_h09m30s26']#specks_TE11_m12d31y2024_h09m30s22']#specks_TE10_m12d31y2024_h09m30s26']
        #TE5_run_list = ['specks_TE1000_m01d13y2025_h13m17s53',
        #                'specks_TE100_m01d13y2025_h13m17s56','specks_TE10_m01d13y2025_h13m18s28',
        #                'specks_TE5000_m01d13y2025_h13m18s40']
        TE5_run_list = ["specks_TE08_m01d14y2025_h09m50s16"]
            #"specks_TE05_m01d14y2025_h09m50s16" ,
        #"specks_TE07_m01d14y2025_h09m50s16",
        #"specks_TE08_m01d14y2025_h09m50s16" ,
        #"specks_TE09_m01d14y2025_h09m50s16" ]
        run_list=TE5_run_list
        run_collection_name="TE"
        #me_at_remote_URL =  'tdunn@mesx.sdsu.edu'
        me_at_remote_URL =  'mesx_cluster'
        output_root_folder=os.path.join("/usr/scratch2/userdata2/tdunn/SpecKS_Output/",
                                        run_collection_name)

        for i in range(0,len(run_list)):
            run_name = run_list[i]
            local_output_folder = os.path.join("/home/tamsen/Data/Specks_output_from_mesx",
                                               run_name)
            self.pull_down_run_data(local_output_folder, me_at_remote_URL, output_root_folder, run_name)

        self.assertEqual(True, True)  # add assertion here

    def pull_down_run_data(self, local_output_folder, me_at_remote_URL, output_root_folder, run_name):

        demographics_folder=run_name + "/*TE*"

        if not os.path.exists(local_output_folder):
            os.makedirs(local_output_folder)

        file_needed_from_allo_folder=["8_final_results/*.csv",
                                      "2_randomized_gene_trees/*.png", "2_randomized_gene_trees/*.txt"]
        file_needed_from_run_folder=["*log*", "*.used.xml"]

        for file in file_needed_from_allo_folder:
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
