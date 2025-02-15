import unittest

from data_aggregation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestKsByMig(unittest.TestCase):

    # Would like to see the impace of migration,
    # Mig rate =0, 10% directly after TDIV,   10% directly after 1/2 TDIV,  50% directly after 1/2 TDIV,
    def test_Ks_for_varying_Migration(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Mg'

        #When Migration is halfway between DIV and NOW, for 100 years.
        run_list_name = "Ks_for_100yr_of_Mig_with_varying_rates"
        demographics_run_list = [False, 'Mig14_m02d11y2025_h10m52s03', 'Mig13_m02d11y2025_h10m52s04'
                                 ]
        #                         'Mig14v2_m02d08y2025_h12m59s40',
        #                         'Mig13v2_m02d09y2025_h09m19s48',

        #                         'Mig14v3_m02d10y2025_h17m03s26',
        #                         'Mig13v3_m02d10y2025_h17m04s04',
        #                         'Mig12v3_m02d10y2025_h17m04s44']

        #                         'Mig14v2_m02d08y2025_h12m59s40',
        #                         'Mig13v2_m02d09y2025_h09m19s48','Mig12v2_m02d09y2025_h12m17s34',
        #                         'Mig11v2_m02d09y2025_h12m21s19']

        # run_list_name = "Ks_for_5yr_of_Mig_with_varying_rates"
        #demographics_run_list = [False,
        #                         'Mig14v2_m02d08y2025_h12m59s40',
        #                         'Mig15v2_m02d09y2025_h19m13s36',
        #                         'Mig16v2_m02d09y2025_h19m13s40',
        #                         'Mig17v2_m02d09y2025_h19m13s44']

        run_list_name = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration"
        demographics_run_list = [False,
                                 'Mig25v2_m02d14y2025_h09m48s26',
                                  'Mig26v2_m02d14y2025_h09m48s27',
                                  'Mig27v2_m02d14y2025_h09m48s30',
                                  'Mig28v2_m02d14y2025_h09m48s32',
                                  'Mig29v2_m02d14y2025_h09m48s35']

        specks_TE5_run_list = [False, False, False, False, False,False,False]


        #xmax_Ks = [0.01,0.05,0.05,0.1,0.5,1]
        xmax_Ks = [0.02 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i/25 for xmax_KS_i in xmax_Ks]

        #xmax_Tc = [5000,5000,5000,10000,50000,100000]
        xmax_Tc = [10000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/25 for xmax_Tc_i in xmax_Tc]


        #ymax_KS = [False for f in demographics_run_list]
        #ymax_Tc = [False for f in demographics_run_list]
        ymax_KS = [500 for f in demographics_run_list]
        ymax_Tc = [800 for f in demographics_run_list]

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
