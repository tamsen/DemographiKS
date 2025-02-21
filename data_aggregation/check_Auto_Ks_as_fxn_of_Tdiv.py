import unittest

from data_aggregation.ks_plot_aggregations_auto_vs_allo import make_Tc_Ks_Allo_vs_Auto_fig_with_subplots


class TestAlloVsAuto_Tdiv(unittest.TestCase):

    def test_allo_vs_auto_Ks_for_varying_varying_Tdiv_times(self):


        demographiKS_auto_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto'
        demographiKS_allo_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Tdiv'

        #specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Tdiv'
        demographics_allo_run_list=[False,
          'TE05fix__m01d06y2025_h11m17s15','TE07_fix__m01d08y2025_h15m09s22',
           'TE08_fix_m01d14y2025_h09m17s20','TE09_fix_m01d13y2025_h14m13s12']

        demographics_auto_run_list=[False,
                                    'Auto_Twgd10000_m01d29y2025_h17m22s23',
                                    'Auto_Twgd100000_m01d29y2025_h17m25s36',
                                    'Auto_Twgd500000_m01d29y2025_h17m22s23',
                                    'Auto_Twgd100000_m01d29y2025_h17m25s36']

        bin_sizes_Tc = [200,200, 200, 200,200]

        xmax_Ks = [0.025,0.025,0.025,0.025,0.025] #0.001  # max(demographiKS_ks_results)
        xmax_Tc = [False,False,False,False,False]
        #xmax_Ks = [0.0005 for f in demographics_allo_run_list]
        bin_sizes_Ks = [xmax_Ks_i/100.0 for xmax_Ks_i in xmax_Ks]

        run_list_name="Ks_for_Allo_and_Auto_varying_varying_Tdiv"
        ymax_KS = [800,800,800,800,800]
        ymax_Tc = [False for f in demographics_allo_run_list]
        #show_KS_predictions=[True,True,True]
        show_KS_predictions = [False, False, False]
        suptitle = "Allo and Auto Ks histograms\n" + \
                                  "Recombination rate = 8e-9, Ne and BI constant"

        make_Tc_Ks_Allo_vs_Auto_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                          demographiKS_allo_out_path, demographics_allo_run_list, run_list_name,
                                          demographics_auto_run_list, demographiKS_auto_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS,  ymax_Tc,suptitle,
                                     show_KS_predictions)

        self.assertEqual(True, True)  # add assertion here

if __name__ == '__main__':
    unittest.main()

