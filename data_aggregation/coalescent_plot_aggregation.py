import datetime
import glob
import math
import os
import unittest
import process_wrapper
import numpy as np
from datetime import datetime
from matplotlib import pyplot as plt
import config

def make_Tc_fig_with_subplots(bin_sizes_Tc,
                                 demographiKS_out_path, demographics_TE9_run_list, run_list_name,
                                 Ks_per_YR,
                                 xmax_Tc, suptitle, total_num_genes):
    ymax_Tc = False
    num_runs = len(demographics_TE9_run_list)
    png_out = os.path.join(demographiKS_out_path, "ks_hist_for_{0}.png".format(run_list_name))
    fig, ax = plt.subplots(1, 4, figsize=(20, 5))
    fig.suptitle(suptitle)
    for i in range(0, num_runs):
        dgx_run_name = demographics_TE9_run_list[i]

        if dgx_run_name:

            dgx_run_path = os.path.join(demographiKS_out_path, dgx_run_name)
            print("dgx_run_path: " +dgx_run_path )
            glob_results=glob.glob(dgx_run_path + '/*.used.xml')
            input_xml_file = glob_results[0]
            config_used = config.DemographiKS_config(input_xml_file)
            dgx_run_duration_in_m = get_run_time_in_minutes(dgx_run_path)

        else:
            config_used = False
            dgx_run_duration_in_m = 0



        slim_csv_file = os.path.join(dgx_run_path, "simulated_ancestral_gene_mrcas.csv")
        loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file)

        theory_output_file = os.path.join(dgx_run_path, "theoretical_ancestral_gene_mrcas.csv")
        loci, theory_mrcas_by_gene = read_data_csv(theory_output_file)
        burnin_times_in_generations=config_used.burnin_time
        plot_title = "Tcoal at Tdiv\nburnin time=" + str(burnin_times_in_generations) + " gen, " \
                     + "Ne=" + str(config_used.ancestral_Ne)

        plot_mrca(ax[i], slim_mrcas_by_gene, False, theory_mrcas_by_gene,
                  dgx_run_duration_in_m, plot_title, config_used.ancestral_Ne,
                  Ks_per_YR, bin_sizes_Tc[i], xmax_Tc[i], ymax_Tc, total_num_genes[i])

    ax[0].set(ylabel="# genes in bin")
    plt.tight_layout()
    plt.savefig(png_out, dpi=550)
    plt.clf()
    plt.close()



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
            if len(data) > 1:
                loci.append(int(data[0]))
                mrcas.append(float(data[1]))
            else:
                mrcas.append(float(data[0]))

    return loci,mrcas


def plot_mrca(this_ax, slim_mrcas_by_gene, specks_mrcas_by_gene, theoretical_mrcas_by_gene,
              dgks_run_duration_in_m, title, Ne, Ks_per_YR, bin_size, xmax, ymax, total_num_genes):

    #fig = plt.figure(figsize=(10, 10), dpi=350
    #Co.T=(1/2N)*e^-((t-1)/2N))

    #max_mrca = max(slim_mrcas_by_gene)
    num_slim_genes=len(slim_mrcas_by_gene)
    if num_slim_genes==0:
        avg_slim_Tc = 0
    else:
        avg_slim_Tc=sum(slim_mrcas_by_gene)/num_slim_genes
    #num_segments=len(slim_mrcas_by_tree)

    if not xmax:
        xmax = max(slim_mrcas_by_gene)
    bins = np.arange(0, xmax , bin_size)

    two_Ne=2.0*Ne
    print("Tc plot bin_size_in_time: " + str(bin_size))
    kingman = [min(total_num_genes,
                   (bin_size*total_num_genes/two_Ne) * math.e ** ((-1 * i) / two_Ne))
               for i in bins]


    if slim_mrcas_by_gene:
        this_ax.hist(slim_mrcas_by_gene, bins=bins, facecolor='b', alpha=0.25,
                                label='SLiM Tcoal by gene\n'
                                + "(" +str(num_slim_genes) + " genes in genome,\n"
                                 +"avg Tc " +str(int(avg_slim_Tc)) + " generations)",
                                density=False)
    #label = 'SLiM Tcoal by gene (total: ' + str(num_genes) + ')',

    if specks_mrcas_by_gene:
        num_specks_genes = len(specks_mrcas_by_gene)
        specks_mrcas_by_gene_in_YRs=[m*10.0**6 for m in specks_mrcas_by_gene]
        this_ax.hist(specks_mrcas_by_gene_in_YRs, bins=bins, facecolor='c', alpha=0.25,
                                label='SpecKS Tcoal by gene\n'
                                + "(" +str(num_specks_genes) + " genes in genome)\n",
                 density=False)
    
    if theoretical_mrcas_by_gene:
        num_theory_genes = len(theoretical_mrcas_by_gene)
        avg_theory_Tc = sum(theoretical_mrcas_by_gene) / num_theory_genes
        this_ax.hist(theoretical_mrcas_by_gene, bins=bins, facecolor='r', alpha=0.25,
                 label='Theoretical Tcoal by gene\n'
                                + "(" +str(num_theory_genes) + " genes in genome,\n"
                                 +"avg Tc " +str(int(avg_theory_Tc)) + " generations)",
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

    annotation_txt= "\n".join([Tc_info,ks_info,ks_info_2,mut_info])+"\n"
    this_ax.annotate(annotation_txt, (0, 0), (0, -60), xycoords='axes fraction', textcoords='offset points', va='top')

    x_axis_label="MRCA time\n" + "demographiKS run time: " + str(round(dgks_run_duration_in_m, 2)) + " min"

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
