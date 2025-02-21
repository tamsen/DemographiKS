import unittest

from data_aggregation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots
from data_aggregation.ks_plot_aggregations_auto_vs_allo import make_Tc_Ks_Allo_vs_Auto_fig_with_subplots


class TestKsForAuto(unittest.TestCase):

    def test_Ks_for_varying_Nb_one_bottleneck_Auto_Na_fixed_at_100(self):

        demographiKS_allo_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Nb/Na100'
        demographics_auto_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto'

        demographics_allo_run_list = [False,
         'KSvs100Na_100Nb_m01d27y2025_h16m51s00', 'KSvs100Na_500Nb_m01d27y2025_h16m51s10',
         'KSvs100Na_1KNb_m01d27y2025_h16m52s55', 'KSvs100Na_5KNb_m01d27y2025_h16m52s44']

        demographics_auto_run_list   = [False,
                                      'Auto100Na_100Nb_m01d29y2025_h18m22s20',
                                      'Auto_100Na_500Nb_m01d29y2025_h18m22s20',
                                      'Auto_100Na_1KNb_m01d29y2025_h18m22s20',
                                      'Auto_100Na_5KNb_m01d29y2025_h18m22s20'
                                      ]


        xmax_Ks = [0.10 for f in demographics_auto_run_list ]
        bin_sizes_Ks = [xmax_KS_i/50 for xmax_KS_i in xmax_Ks]

        xmax_Tc = [5000 for f in demographics_auto_run_list]
        bin_sizes_Tc = [xmax_Tc_i / 50 for xmax_Tc_i in xmax_Tc]

        ymax_KS = [False for f in demographics_auto_run_list]
        ymax_Tc = [100 for f in demographics_auto_run_list]
        run_list_name = "Ks_for_Allo_and_Auto_varying_varying_Nb_100Na"


        show_KS_predictions = [False, False, False]
        suptitle = "Auto and Allo Ks histograms\n"
        make_Tc_Ks_Allo_vs_Auto_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_allo_out_path, demographics_allo_run_list, run_list_name,
                                     demographics_auto_run_list, demographics_auto_out_path ,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions)

        self.assertEqual(True, True)  # add assertion here


    def test_Ks_for_varying_Nb_one_bottleneck_Na_fixed_at_10K(self):

        #full, w/Ne 10K
        demographiKS_allo_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto/Auto_vs_Nb/Na10K'
        demographics_auto_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto/Auto_vs_Nb/Na10K'

        demographics_allo_run_list = [
            False,
               'KSvs10KNa_100Nb_m01d27y2025_h16m40s00', 'KSvs10KNa_1KNb_m01d27y2025_h16m40s00',
               'KSvs10KNa_500Nb_m01d27y2025_h16m40s00', 'KSvs10KNa_5KNb_m01d27y2025_h16m41s35']

        demographics_auto_run_list = [
            False,
            'Auto10KNa_100Nb_m02d20y2025_h17m32s04',
            'Auto_10KNa_500Nb_m02d20y2025_h17m32s04',
            'Auto_10KNa_1KNb_m02d20y2025_h17m32s04',
            'Auto_10KNa_5KNb_m02d20y2025_h17m32s04']

        #xmax_Ks = [0.8  for f in demographics_auto_run_list ]
        xmax_Ks = [[0.05, 0.05,0.05,0.6,0.8 ],[0.4, 0.4,0.5,0.6,0.8 ]]
        bin_sizes_Ks = [[xmax_KS_i/50 for xmax_KS_i in xmax_Ks[0]],
                        [xmax_KS_i / 50 for xmax_KS_i in xmax_Ks[1]]]

        xmax_Tc = [80000 for f in demographics_auto_run_list]
        bin_sizes_Tc = [xmax_Tc_i / 50 for xmax_Tc_i in xmax_Tc]

        ymax_KS = [[False for f in demographics_auto_run_list],
                   [100 for f in demographics_auto_run_list]]

        ymax_Tc = [False for f in demographics_auto_run_list]
        run_list_name = "Ks_for_Allo_and_Auto_varying_varying_Nb_10K_Na"

        show_KS_predictions = [False, False, False]
        suptitle = "Auto and Allo Ks histograms\n"
        make_Tc_Ks_Allo_vs_Auto_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
        demographiKS_allo_out_path, demographics_allo_run_list, run_list_name,
        demographics_auto_run_list, demographics_auto_out_path,
        xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
        suptitle, show_KS_predictions)

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
