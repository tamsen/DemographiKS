import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import shutil
from modules import ks_calculator


def get_Ks_from_file(paml_out_file):

    ks_results =[]
    ordered_ortholog_list = []
    with open(paml_out_file, "r") as f:
        lines = f.readlines()

    row_index=0
    for l in lines[1:len(lines)]:
            data = l.split()
            if len(data)==0:
                continue

            #print("line: " + l)
            ordered_ortholog_list.append(data[0])
            for col_index in range(1,len(data)):
                d=data[col_index]
                new_ks_result=ks_data(float(d),row_index,col_index-1)
                ks_results.append(new_ks_result)
            row_index = row_index+1

    for ks_result in ks_results:
        ks_result.set_orthologs(ordered_ortholog_list)

    return ks_results



def extract_K_values(csv_file_out, res_files):
    KS_values = []
    with open(csv_file_out, 'w') as f:

        #f.writelines("SpecKS " + version_string + "\n")  # version_info.to_string())
        #f.writelines("GeneTree,leaf names,Ks,path to original output file\n")
        for paml_out_file in res_files:
            print(paml_out_file)
            base_name = os.path.basename(paml_out_file)
            Ks_for_og = get_Ks_from_file(paml_out_file)
            for Ks_data in Ks_for_og:
                Ks_value = Ks_data.ks_between_ortholog_and_LCA
                ortholog_names_str = str(Ks_data.ortholog_pair).replace(",", " ")
                f.writelines(base_name + "," + ortholog_names_str + "," +
                             str(Ks_value) + "," + paml_out_file + "\n")
                KS_values.append(Ks_value)
    return KS_values

def plot_Ks_histogram(PAML_hist_out_file, species_name, Ks_results, WGD_as_Ks, SPEC_as_Ks,
                      max_Ks, max_y, alg_name, color, bin_size):

    fig = plt.figure(figsize=(10, 10), dpi=100)
    x = Ks_results

    if max_Ks:
        bins = np.arange(0, max_Ks + 0.1, bin_size)
        #log.write_to_log(PAML_hist_out_file)
        n, bins, patches = plt.hist(x, bins=bins, facecolor=color, alpha=0.25, label='histogram data')
        plt.xlim([0, max_Ks * (1.1)])
    else:
        nBins=50
        n, bins, patches = plt.hist(x, bins=nBins, facecolor=color, alpha=0.25, label='histogram data')
    if max_y:
        plt.ylim([0, max_y])

    if WGD_as_Ks:
        plt.axvline(x=WGD_as_Ks, color='b', linestyle='-', label="WGD time as Ks")
    if SPEC_as_Ks:
        plt.axvline(x=SPEC_as_Ks, color='r', linestyle='--', label="SPEC time as Ks")

    plt.legend()
    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")
    plt.title("Ks histogram for " + species_name+ "\n" +
              "algorithm: PAML " + alg_name + " num paralogs: " + str(len(Ks_results)))
    plt.savefig(PAML_hist_out_file)
    plt.clf()
    plt.close()

class ks_data():

    round_trip_ks_between_orthologs=0
    ks_between_ortholog_and_LCA=0
    ortholog_pair=[]
    row_index=-1
    col_index=-1
    def __init__(self, round_trip_ks_data_result, row_index, col_index):
        self.round_trip_ks_between_orthologs = round_trip_ks_data_result
        self.ks_between_ortholog_and_LCA = self.round_trip_ks_between_orthologs * 0.5
        self.ortholog_pair = []
        self.row_index = row_index
        self.col_index = col_index

    def set_orthologs(self, ordered_ortholog_list):

        ortholog1=ordered_ortholog_list[self.row_index]
        ortholog2=ordered_ortholog_list[self.col_index]
        self.ortholog_pair = [ortholog1,ortholog2]