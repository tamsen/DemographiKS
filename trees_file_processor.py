import math
import os
from scipy.stats import expon
import numpy as np
import tskit
from matplotlib import pyplot as plt


def plot_coalescent(trees_file, genome_index_1,genome_index_2, config, output_folder):

        ts = tskit.load(trees_file)
        reduced_ts = ts.simplify([genome_index_1,genome_index_2])  # simplify to the first 10 samples
        genome_length=config.total_num_bases
        gene_length=3*config.num_codons_in_a_gene
        num_genes=int(genome_length/gene_length)

        #get tree for every gene
        gene_starts=[i*gene_length for i in range(0,num_genes)]
        gene_centers=[g+int(gene_length/2)for g in gene_starts]
        mrcas=[]
        for gene_center in gene_centers:
            tree_for_gene=reduced_ts.at(gene_center)
            mrca = reduced_ts.node(tree_for_gene.root)
            mrcas.append(mrca.time)
            #print("gene {0}, mrca time: {1}".format(gene_center,mrca.time))

        csv_file_out1 = os.path.join(output_folder, "simulated_ancestral_mrcas.csv")
        save_mrca_values(csv_file_out1, mrcas,gene_starts)
        theoretical_mrcas=theoretical_coalescent(num_genes,config.ancestral_Ne)
        csv_file_out2 = os.path.join(output_folder, "theoretical_ancestral_mrcas.csv")
        save_mrca_values(csv_file_out2,theoretical_mrcas,gene_starts)
        png_out1 = os.path.join(output_folder, "mrca_hist1.png")
        plot_mrca(mrcas,theoretical_mrcas, png_out1)
        png_out2 = os.path.join(output_folder, "mrca_hist2.png")
        plot_mrca(mrcas,[], png_out2)

def theoretical_coalescent(num_genes,N):
    #Co.T=(1/2N)*e^-((t-1)/2N))
    loc=0
    xscale=2.0*N
    random_draws_from_distribution = expon.rvs(loc=loc, size=num_genes, scale=xscale)
    return random_draws_from_distribution

def plot_mrca(slim_mrcas, theoretical_mrcas,png_out):
    fig = plt.figure(figsize=(10, 10), dpi=350)
    
    #Co.T=(1/2N)*e^-((t-1)/2N))
    x = slim_mrcas
    label = "Coalescent Times For Genes From Sampled Two Ancestral Genomes"
    print(label)
    max_mrca = max(slim_mrcas)
    num_genes=len(slim_mrcas)
    bin_size = 100
    bins = np.arange(0, max_mrca, bin_size)

    plt.hist(x, bins=bins, facecolor='b', alpha=0.25,
                                label='SLiM genes (total: '+str(num_genes)+')',
                                density=False)
    plt.hist(theoretical_mrcas, bins=bins, facecolor='c', alpha=0.25,
                                label='Theoretical (total: '+str(num_genes)+')',
                                density=False)
    plt.xlim([0, max_mrca * (1.1)])
    plt.title(label)
    # plt.xlim([0, max_mrca])
    # plt.ylim([0, 300])
    plt.legend()
    plt.xlabel("MRCA time")
    plt.ylabel("# genes in bin")
    plt.savefig(png_out)
    plt.clf()
    plt.close()
    
def save_mrca_values(csv_file_out, mrcas_values,gene_starts):

    with open(csv_file_out, 'w') as f:

        f.writelines("GeneStart,TimeToMRCA\n")
        for i in range(0,len(mrcas_values)):
            mrcas_value = mrcas_values[i]
            gene_start= gene_starts[i]
            f.writelines(str(gene_start) + "," + str(mrcas_value))
