import unittest

from data_aggregation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots
from data_aggregation.ks_plot_aggregations_auto_vs_allo import make_Tc_Ks_Allo_vs_Auto_fig_with_subplots


class TestAutoKsByNe(unittest.TestCase):

    def test_Auto_Ks_for_varying_Ne_Tdiv_1000(self):

        demographiKS_auto_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto'
        demographiKS_allo_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Ne'


        demographiKS_allo_run_list = [False, "DGKS_10_10_m5_RC7_m01d14y2025_h09m14s47",
                                    "DGKS_100_100_m5_RC7_m01d14y2025_h09m14s18",
                                    "DGKS_1000_1000_m5_BI40_RC7_m01d13y2025_h15m36s22",
                                    "DGKS_5000_5000_m5_BI_40K_RC7_m01d14y2025_h09m13s52"]

        demographiKS_auto_run_list = [False,
                                      'Auto_Ne10_m01d29y2025_h17m47s29', 'Auto_Ne100_m01d29y2025_h17m47s29',
                                      'Auto_Ne1000_m01d29y2025_h17m52s03','Auto_Ne5000_m01d29y2025_h17m47s29'
                                      ]

        xmax_Ks = [0.4,0.4,0.4,0.4,0.4] #[0.01,0.01,0.01,0.1,0.2]#False#0.08  # for mut rate e-5
        bin_sizes_Ks = [xmax_Ks_i / 50.0 for xmax_Ks_i in xmax_Ks]
        xmax_Tc = [2000,2000,2000,20000,40000]
        bin_sizes_Tc = [xmax_Tc_i / 50.0 for xmax_Tc_i in xmax_Tc ]

        run_list_name="Ks_for_Allo_and_Auto_varying_varying_Ne"
        ymax_Ks = [50 for f in demographiKS_allo_run_list ]
        ymax_Tc = [False for f in demographiKS_allo_run_list ]

        suptitle = "Auto vs Allo, Tcoal and Ks\n" + \
                   "Recombination rate = 1.26e-7, mut rate 1.0e-5"
        show_KS_predictions=[False,False,False]


        make_Tc_Ks_Allo_vs_Auto_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_allo_out_path, demographiKS_allo_run_list, run_list_name,
                                     demographiKS_auto_run_list, demographiKS_auto_out_path ,
                                     xmax_Ks, xmax_Tc, ymax_Ks, ymax_Tc,
        suptitle, show_KS_predictions)

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
