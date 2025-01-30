import os.path
import unittest
import process_wrapper


class MyDGXDataFetcher(unittest.TestCase):

    def test_fetch_demographics_data(self):

        #/usr/scratch2/userdata2/tdunn/DemographiKS_Output
        run_list = ['Auto_100Na_1KNb_m01d29y2025_h18m22s20']
        #    'Auto100Na_100Nb_m01d29y2025_h18m22s20','Auto_100Na_500Nb_m01d29y2025_h18m22s20',
        #            'Auto_100Na_1KNb_m01d29y2025_h18m22s20','Auto_100Na_5KNb_m01d29y2025_h18m22s20']

        #/usr/scratch2/userdata2/tdunn/DemographiKS_Output
        #run_list = ['Auto_Ne1000_m01d29y2025_h17m52s03']
        #'Auto_Ne10_m01d29y2025_h17m47s29','Auto_Ne100_m01d29y2025_h17m47s29']
        #            'Auto_Ne1000_m01d29y2025_h17m52s03','Auto_Ne5000_m01d29y2025_h17m47s29']

        #/usr/scratch2/userdata2/tdunn/DemographiKS_Output
        #run_list = [ 'Auto_Twgd500000_m01d29y2025_h17m22s23','Auto_Twgd1000000_m01d29y2025_h17m22s24']
        #'Auto_Twgd10000_m01d29y2025_h17m22s23', 'Auto_Twgd100000_m01d29y2025_h17m25s36',

        run_collection_name="TE"
        me_at_remote_URL =  'mesx_cluster'
        output_root_folder = os.path.join("/usr/scratch2/userdata2/tdunn/DemographiKS_Output")
        #if False:
        #    output_root_folder = os.path.join("/usr/scratch2/userdata2/tdunn/DemographiKS_Output")
        #else:
        #    output_root_folder=os.path.join("/usr/scratch2/userdata2/tdunn/DemographiKS_output",
        #                                run_collection_name)



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
