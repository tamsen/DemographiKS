import unittest

from data_aggregation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestKsByMig(unittest.TestCase):

    # Would like to see the impace of migration,
    # Mig rate =0, 10% directly after TDIV,   10% directly after 1/2 TDIV,  50% directly after 1/2 TDIV,
    def test_Ks_for_varying_Migration(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Mg'

        demographics_run_list = [False,
                                 'Mig9_m02d06y2025_h16m59s59','Mig11_m02d06y2025_h17m09s33']
        #'Mig7_m02d06y2025_h16m03s11', 'Mig6_m02d06y2025_h16m03s14',
        #'Mig4_m02d06y2025_h13m45s58', 'Mig5_m02d06y2025_h13m45s56',
        #'Mig2_m02d06y2025_h12m42s08', 'Mig3_m02d06y2025_h12m51s06',
        specks_TE5_run_list = [False, False, False, False, False,False,False]


        #xmax_Ks = [0.01,0.05,0.05,0.1,0.5,1]
        xmax_Ks = [0.02 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i/50 for xmax_KS_i in xmax_Ks]

        #xmax_Tc = [5000,5000,5000,10000,50000,100000]
        xmax_Tc = [5000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/50 for xmax_Tc_i in xmax_Tc]


        #ymax_KS = [False for f in demographics_run_list]
        #ymax_Tc = [False for f in demographics_run_list]
        ymax_KS = [100 for f in demographics_run_list]
        ymax_Tc = [200 for f in demographics_run_list]

        run_list_name = "Ks_for_varying_varying_Mg_constantRC"
        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate

        show_KS_predictions = [False, False, False]
        suptitle = "Ks histograms\n"
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, demographiKS_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions)

        self.assertEqual(True, True)  # add assertion here

    def test_Ks_for_varying_RC_1KNe(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_RC/save_Ne_1K'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'


        #full, w/Ne1000
        #KSvsRC9_m01d21y2025_h14m24s03

        #demographics_run_list = [False, 'KSvsRC6_m01d21y2025_h20m06s34', 'KSvsRC7_m01d21y2025_h12m16s26',
        #                           'KSvsRC8_m01d23y2025_h10m39s34', 'KSvsRC9_m01d21y2025_h14m24s03',
        #                           'KSvsRC10_m01d24y2025_h09m16s53', 'KSvsRC11_m01d23y2025_h11m04s33']

        #ks_hist_by_Ks_for_varying_varying_RC_test_Ne1000_long_burnin_save_for_paper.png
        demographics_run_list = [False, 'KSvsRC6_m01d21y2025_h20m06s34','KSvsRC7_m01d21y2025_h12m16s26',
                                 'KSvsRC8_m01d23y2025_h10m39s34','KSvsRC9_m01d24y2025_h08m51s36',
                                 'KSvsRC10_m01d24y2025_h09m16s53','KSvsRC11_m01d23y2025_h11m04s33']

        #problems only
        #demographics_run_list = [False,  'KSvsRC8_m01d23y2025_h10m39s34',
        #                             'KSvsRC10_m01d24y2025_h10m40s24','KSvsRC11_m01d23y2025_h11m04s33']
        #demographics_run_list_old= [False,  'KSvsRC8_m01d21y2025_h08m21s04','KSvsRC11_m01d22y2025_h11m14s42']

        specks_TE5_run_list = [False,False,False,False,False,False,False]

        bin_sizes_Tc = [200, 200, 200, 200, 200, 200, 200]
        bin_sizes_Ks = [0.002, 0.002, 0.002,0.002, 0.002, 0.002, 0.002]
        xmax_Ks = [0.2 for f in demographics_run_list] #[0.025, 0.025, 0.025, 0.025, 0.025]  # 0.001  # max(demographiKS_ks_results)
        xmax_Tc = [15000 for f in demographics_run_list]
        ymax_KS = [140 for f in demographics_run_list]
        ymax_Tc = [100 for f in demographics_run_list]

        run_list_name = "Ks_for_varying_varying_RC"
        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate

        show_KS_predictions = [False, False, False]
        suptitle = "SLiM and SpecKS Ks histograms\n"
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions)

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
