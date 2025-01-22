import unittest

from data_aggregation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestKsByNe(unittest.TestCase):

    def test_Ks_for_varying_RC(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'
        demographics_run_list = [False, 'KSvsRC6_m01d21y2025_h20m06s34','KSvsRC7_m01d21y2025_h12m16s26',
                                 'KSvsRC8_m01d21y2025_h08m21s04','KSvsRC9_m01d21y2025_h14m24s03',
                                         'KSvsRC11_m01d21y2025_h20m31s14']

        specks_TE5_run_list = [False,False,False,False,False,False,False]

        bin_sizes_Tc = [200, 200, 200, 200, 200, 200, 200]
        bin_sizes_Ks = [0.002, 0.002, 0.002,0.002, 0.002, 0.002, 0.002]
        xmax_Ks = [False for f in demographics_run_list] #[0.025, 0.025, 0.025, 0.025, 0.025]  # 0.001  # max(demographiKS_ks_results)
        xmax_Tc = [15000 for f in demographics_run_list]

        run_list_name = "Ks_for_varying_varying_RC"
        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate
        Ks_per_YR = 0.833*10**-5
        ymax_KS = [False for f in demographics_run_list]
        show_KS_predictions = [False, False, False]
        total_num_genes = 1000
        suptitle = "SLiM and SpecKS Ks histograms\n"
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, suptitle,
                                     show_KS_predictions)

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
