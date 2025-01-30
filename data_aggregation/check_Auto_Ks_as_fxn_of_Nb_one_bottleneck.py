import unittest

from data_aggregation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestKsForAuto(unittest.TestCase):

    def test_Ks_for_varying_Nb_one_bottleneck_Auto(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'

        demographics_run_list = [False,'Auto1KNe_m01d29y2025_h12m27s12',
                                 'Auto1KNe_m01d29y2025_h11m39s31',
                                 'KSvs1KNe_m01d26y2025_h18m18s41']
        specks_TE5_run_list = [False,False,False,False,False,False,False]


        xmax_Ks = [0.10 for f in demographics_run_list ]
        bin_sizes_Ks = [xmax_KS_i/50 for xmax_KS_i in xmax_Ks]

        xmax_Tc = [10000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/50 for xmax_Tc_i in xmax_Tc]


        ymax_KS = [100 for f in demographics_run_list]
        ymax_Tc = [100 for f in demographics_run_list]

        run_list_name = "Ks_for_Auto"
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


    def test_Ks_for_varying_Nb_one_bottleneck_Na_fixed_at_10K(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Nb/Na10K'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'

        #full, w/Ne 10K
        demographics_run_list = [False,
                                 'KSvs10KNa_100Nb_m01d27y2025_h16m40s00','KSvs10KNa_1KNb_m01d27y2025_h16m40s00',
                                 'KSvs10KNa_500Nb_m01d27y2025_h16m40s00','KSvs10KNa_5KNb_m01d27y2025_h16m41s35']
        specks_TE5_run_list = [False,False,False,False,False,False,False]


        xmax_Ks = [0.8 for f in demographics_run_list ]
        bin_sizes_Ks = [xmax_KS_i/50 for xmax_KS_i in xmax_Ks]

        xmax_Tc = [80000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/50 for xmax_Tc_i in xmax_Tc]


        ymax_KS = [False for f in demographics_run_list]
        ymax_Tc = [False for f in demographics_run_list]

        run_list_name = "Ks_for_varying_Nb_constantRC8andNa10K"
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
