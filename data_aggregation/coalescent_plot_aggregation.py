import datetime
import math
import os
import unittest
import process_wrapper
import numpy as np
from datetime import datetime
from matplotlib import pyplot as plt


class TestTCoalAggregation(unittest.TestCase):
    
    TE_run_list = ['TE01_m12d20y2024_h14m16s21', 'TE03_m12d20y2024_h14m26s56',
                   'TE02_m12d20y2024_h14m22s23']  # ,'TE04_m12d20y2024_h14m31s31
    
    #"SpecKS like simulations, gradually increasing the population size
    def test_TCoal_with_varying_TE(self):

        #make 1 histogram.
        aggragate_output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx"

        TE_run_list = ['TE01_m12d20y2024_h14m16s21', 'TE02_m12d20y2024_h14m22s23',
                       'TE03_m12d20y2024_h14m26s56']  # ,'TE04_m12d20y2024_h14m31s31

        burnin_times_in_generations=[5e6, 5e7, 5e7]
        Ne_in_generations=[10, 100, 1000]
       
        num_runs=len(TE_run_list)
        bin_sizes = [2,20, 200]
        xmax = False#1000
        ymax=False
        png_out=os.path.join(aggragate_output_folder,"mrcsa_by_TE_tests.png")
        fig, ax = plt.subplots(1, num_runs,figsize=(20,5))
        fig.suptitle("SLiM Tcoal by gene in ancestral species at Tdiv\n" + \
            "Recombination rate = 8e-9, Ne varies, BI varies")

        for i in range(0,num_runs):
            run_name=TE_run_list[i]
            local_output_folder = os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
                                               run_name)

            run_duration_in_m = get_run_time_in_minutes(local_output_folder)
            slim_csv_file=os.path.join(local_output_folder,"simulated_ancestral_gene_mrcas.csv")
            loci, slim_mrcas_by_gene=read_data_csv(slim_csv_file)
            
            theory_output_file = os.path.join(local_output_folder, "theoretical_ancestral_gene_mrcas.csv")
            loci, theory_mrcas_by_gene=read_data_csv(theory_output_file)
            plot_title="burnin time=" + str(burnin_times_in_generations[i]) + " gen\n"\
                        + "Ne=" + str(Ne_in_generations[i])
            plot_mrca(ax[i],slim_mrcas_by_gene, [],theory_mrcas_by_gene,
                run_duration_in_m, plot_title,bin_sizes[i],xmax, ymax)
        
        ax[0].set(ylabel="# genes in bin")
        plt.tight_layout()
        plt.savefig(png_out, dpi=550)
        plt.clf()
        plt.close()
        
        self.assertEqual(True, True)  # add assertion here

    
    def test_TCoal_with_varying_BI(self):

        #make 1 histogram.
        aggragate_output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx"

        BI_run_list=["BI2_m12d19y2024_h10m57s18","BI3_m12d19y2024_h11m02s23",
                     "BI4_m12d19y2024_h11m02s24","BI5_m12d19y2024_h11m02s27"]
                    #"BI7_m12d19y2024_h11m02s31"]

        burnin_times_in_generations=[5e2, 5e3, 5e4, 5e5, 5e6] #, 5e7]#8e-8, <- longer BI times still running..
        theory_output_file=os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
                                        BI_run_list[0],"theoretical_ancestral_gene_mrcas.csv")
        num_runs=len(BI_run_list)
        bin_size = 10000
        xmax = 160_000
        png_out=os.path.join(aggragate_output_folder,"mrcsa_by_burnin_time.png")
        fig, ax = plt.subplots(1, num_runs,figsize=(20,5))
        fig.suptitle("SLiM Tcoal by gene in ancestral species at Tdiv\n" + \
            "Recombination rate = 8e-10, Ne=1e4")

        for i in range(0,num_runs):
            run_name=BI_run_list[i]
            local_output_folder = os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
                                               run_name)

            run_duration_in_m = get_run_time_in_minutes(local_output_folder)
            slim_csv_file=os.path.join(local_output_folder,"simulated_ancestral_gene_mrcas.csv")
            loci, slim_mrcas_by_gene=read_data_csv(slim_csv_file)
            loci, theory_mrcas_by_gene=read_data_csv(theory_output_file)
            plot_title="burnin time=" + str(burnin_times_in_generations[i]) + " generations"
            plot_mrca(ax[i],slim_mrcas_by_gene, [],theory_mrcas_by_gene,
                run_duration_in_m, plot_title,bin_size,xmax)
        
        ax[0].set(ylabel="# genes in bin")
        plt.tight_layout()
        plt.savefig(png_out, dpi=550)
        plt.clf()
        plt.close()
        
        self.assertEqual(True, True)  # add assertion here

    def test_TCoal_with_varying_RC(self):

        # make 1 histogram.
        aggragate_output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx"

        RC_run_list = ["RC08_m12d18y2024_h14m30s58",
            "RC09_m12d18y2024_h14m30s54", "RC10_m12d18y2024_h13m36s12",
                       "RC11_m12d18y2024_h14m35s15", "RC12_m12d18y2024_h14m35s17",
                       "RC13_m12d18y2024_h14m35s19", "RC14_m12d18y2024_h14m39s19"]
        recombination_rates = [8e-8,8e-9, 8e-10, 8e-11, 8e-12, 8e-13, 8e-14]  # 8e-8,
        theory_output_file = os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
                                          RC_run_list[0], "theoretical_ancestral_gene_mrcas.csv")
        num_runs = len(RC_run_list)
        bin_size = 10000
        xmax = 160_000
        png_out = os.path.join(aggragate_output_folder, "mrcsa_by_recombination_rate.png")
        fig, ax = plt.subplots(1, num_runs, figsize=(20, 5))
        fig.suptitle("SLiM Tcoal by gene in ancestral species at Tdiv\n" + \
                     "burnin_time = 5000000 generations")

        for i in range(0, num_runs):
            run_name = RC_run_list[i]
            local_output_folder = os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
                                               run_name)
            run_duration_in_m = get_run_time_in_minutes(local_output_folder)
            slim_csv_file = os.path.join(local_output_folder, "simulated_ancestral_gene_mrcas.csv")
            loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file)
            loci, theory_mrcas_by_gene = read_data_csv(theory_output_file)
            plot_title="recombination rate=" + str(recombination_rates[i])
            plot_mrca(ax[i], slim_mrcas_by_gene, [], theory_mrcas_by_gene,
                      run_duration_in_m,plot_title, bin_size, xmax)

        ax[0].set(ylabel="# genes in bin")
        plt.tight_layout()
        plt.savefig(png_out, dpi=550)
        plt.clf()
        plt.close()

        self.assertEqual(True, True)  # add assertion here


    def test_TCoal_with_varying_Ne(self):

        # make 1 histogram.
        aggragate_output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx"

        Ne_run_list = ["Ne1_m12d18y2024_h16m08s35",
                     "Ne2_m12d18y2024_h16m08s38",
                     "Ne3_m12d18y2024_h16m08s40",
                     "Ne4_m12d18y2024_h16m10s52","Ne5_m12d18y2024_h16m10s54",]
        effective_population_size = [1e1,1e2, 1e3, 1e4,1e5]

        num_runs = len(Ne_run_list )

        xmax = False#160_000
        ymax = 2_000
        png_out = os.path.join(aggragate_output_folder, "mrcsa_by_ancestral_pop_size_Ne.png")
        fig, ax = plt.subplots(1, num_runs, figsize=(20, 5))
        fig.suptitle("SLiM Tcoal by gene in ancestral species at Tdiv, by Ne\n" + \
                     "burnin_time = 5000000 generations, recombination rate = 8e-9")

        for i in range(0, num_runs):
            bin_size = effective_population_size[i]/10.0
            run_name = Ne_run_list[i]
            local_output_folder = os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
                                               run_name)
            run_duration_in_m = get_run_time_in_minutes(local_output_folder)
            slim_csv_file = os.path.join(local_output_folder, "simulated_ancestral_gene_mrcas.csv")
            theory_output_file = os.path.join(local_output_folder, "theoretical_ancestral_gene_mrcas.csv")
            loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file)
            loci, theory_mrcas_by_gene = read_data_csv(theory_output_file)
            plot_title="pop size (Ne)=" + str(effective_population_size[i])
            plot_mrca(ax[i], slim_mrcas_by_gene, [], theory_mrcas_by_gene,
                      run_duration_in_m,plot_title, bin_size, xmax, ymax)

        ax[0].set(ylabel="# genes in bin")
        plt.tight_layout()
        plt.savefig(png_out, dpi=550)
        plt.clf()
        plt.close()

        self.assertEqual(True, True)  # add assertion here

    def test_TCoal_with_varying_Ge(self):

        # make 1 histogram.
        aggragate_output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx"

        GE_run_list =["GE4_m12d19y2024_h11m47s58","GE5_m12d19y2024_h11m47s58",
                      "GE6_m12d19y2024_h11m48s02","GE7_m12d19y2024_h13m30s32",
                      "GE8_m12d19y2024_h13m34s13"]
        genome_size = [1e4,1e5, 1e6, 1e7,1e8]
        #bin_sizes= [10,10,100,1000,10000]
        bin_sizes = [1000, 1000, 1000, 1000, 1000]
        num_runs = len(GE_run_list )
        #xmaxes = [20,160,160_0,160_00,160_000]# False#160_000
        xmaxes = [75_000, 75_000, 75_000, 75_000, 75_000]
        ymax = False#2_000
        png_out = os.path.join(aggragate_output_folder, "mrcsa_by_genome_size_GE.png")
        fig, ax = plt.subplots(1, num_runs, figsize=(20, 5))
        fig.suptitle("SLiM Tcoal by gene in ancestral species at Tdiv, by genome size\n" + \
                     "burnin_time = 5000000 generations, recombination rate = 8e-9")

        for i in range(0, num_runs):

            bin_size = bin_sizes[i]
            #ymax = 10**i
            run_name = GE_run_list[i]
            local_output_folder = os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
                                               run_name)
            run_duration_in_m = get_run_time_in_minutes(local_output_folder)
            slim_csv_file = os.path.join(local_output_folder, "simulated_ancestral_gene_mrcas.csv")
            theory_output_file = os.path.join(local_output_folder, "theoretical_ancestral_gene_mrcas.csv")
            loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file)
            loci, theory_mrcas_by_gene = read_data_csv(theory_output_file)
            plot_title="genome size = " + str(genome_size[i])
            plot_mrca(ax[i], slim_mrcas_by_gene, [], theory_mrcas_by_gene,
                      run_duration_in_m,plot_title, bin_size, xmaxes[i], ymax)

        ax[0].set(ylabel="# genes in bin")
        plt.tight_layout()
        plt.savefig(png_out, dpi=550)
        plt.clf()
        plt.close()

        self.assertEqual(True, True)  # add assertion here



def read_data_csv(csv_file):

    loci=[]
    mrcas=[]

    with open(csv_file, "r") as f:

        while True:
            line = f.readline()
            if "Start" in line:
                continue
            if len(line)==0:
                break

            data = line.strip().split(",")
            loci.append(int(data[0]))
            mrcas.append(float(data[1]))

    return loci,mrcas


def plot_mrca(this_ax,slim_mrcas_by_gene, slim_mrcas_by_tree, theoretical_mrcas_by_gene,
              run_duration_in_m, title, Ne, Ks_per_YR, bin_size,xmax, ymax):

    #fig = plt.figure(figsize=(10, 10), dpi=350
    #Co.T=(1/2N)*e^-((t-1)/2N))

    #max_mrca = max(slim_mrcas_by_gene)
    num_slim_genes=len(slim_mrcas_by_gene)
    num_theory_genes=len(theoretical_mrcas_by_gene)
    avg_slim_Tc=sum(slim_mrcas_by_gene)/num_slim_genes
    avg_theory_Tc=sum(theoretical_mrcas_by_gene)/num_theory_genes
    num_segments=len(slim_mrcas_by_tree)

    if not xmax:
        xmax = max(slim_mrcas_by_gene)
    bins = np.arange(0, xmax , bin_size)

    two_Ne=2.0*Ne
    print("Tc plot bin_size_in_time: " + str(bin_size))
    kingman = [min(num_theory_genes,
                   (bin_size*num_theory_genes/two_Ne) * math.e ** ((-1 * i) / two_Ne))
               for i in bins]


    if len(slim_mrcas_by_gene) > 0:
        this_ax.hist(slim_mrcas_by_gene, bins=bins, facecolor='b', alpha=0.25,
                                label='SLiM Tcoal by gene\n'
                                + "(" +str(num_slim_genes) + " genes in genome,\n"
                                 +"avg Tc " +str(int(avg_slim_Tc)) + " generations.",
                                density=False)
    #label = 'SLiM Tcoal by gene (total: ' + str(num_genes) + ')',

    if len(slim_mrcas_by_tree) > 0:
        this_ax.hist(slim_mrcas_by_tree, bins=bins, facecolor='r', alpha=0.25,
                 label='SLiM Tcoal by segment (total: '+str(num_segments)+')',
                 density=False)
    
    if len(theoretical_mrcas_by_gene) > 0:
        this_ax.hist(theoretical_mrcas_by_gene, bins=bins, facecolor='c', alpha=0.25,
                 label='Theoretical Tcoal by gene\n'
                                + "(" +str(num_theory_genes) + " genes in genome,\n"
                                 +"avg Tc " +str(int(avg_theory_Tc)) + " generations.",
                 density=False)

    this_ax.plot(bins,kingman,c='red', label='Expectations under Kingman')
    #Tc_to_Ks = avg_theory_Tc * Ks_per_YR
    Tc_to_Ks = avg_slim_Tc * Ks_per_YR
    Tc_info='   mean Tc by Kingman = ' + \
        str(2.0*Ne)

    ks_info='   simulated mean Ks at Tdiv = ' + \
        "{:.2E}".format(Tc_to_Ks )

    ks_info_2='   2*Ne*Ks_per_YR = ' + \
        "{:.2E}".format(2.0*Ne*Ks_per_YR)

    mut_info = '   simulated num mutations per gene = ' + \
        "{:.2E}".format(Tc_to_Ks*3.0*1000.0)
    this_ax.text(0,0, "\n".join([Tc_info,ks_info,ks_info_2,mut_info])+"\n", fontsize = 12)
    x_axis_label="MRCA time\n" + "run time: " +str(round(run_duration_in_m,2)) + " min"

    if ymax:
        this_ax.set(ylim=[0, ymax])
    this_ax.set(xlim=[0, xmax])
    this_ax.set(xlabel=x_axis_label)
    this_ax.set(title=title)
    #this_ax.set(title="Recom. rate: " + str(recombination_rate))
    this_ax.legend()
    #plt.ylabel("# genes in bin")
    #plt.savefig(png_out)

def get_run_time_in_minutes(local_output_folder):
        out_string, error_string = process_wrapper.run_and_wait_with_retry(
            ['ls'], local_output_folder, "Connection reset by peer", 2, 5)

        print("folder: \n" + local_output_folder)
        print("ls: \n" + out_string)
        log_file = [s for s in out_string.split() if 'log.txt' in s][0]
        print("log file: " + log_file)
        head_string, error_string = process_wrapper.run_and_wait_with_retry(
            ['head', log_file], local_output_folder, "Connection reset by peer", 2, 5)
        tail_string, error_string = process_wrapper.run_and_wait_with_retry(
            ['tail', log_file], local_output_folder, "Connection reset by peer", 2, 5)
        first_line = head_string.split("\n")[0].split("\t")[0]
        last_line = tail_string.split("\n")[-4].split("\t")[0]
        datetime_start = datetime.strptime(first_line, '%d/%m/%Y,%H:%M:%S:')
        datetime_end = datetime.strptime(last_line, '%d/%m/%Y,%H:%M:%S:')
        difference = datetime_end - datetime_start
        duration_in_m = difference.total_seconds() / 60.0
        return duration_in_m

if __name__ == '__main__':
    unittest.main()
