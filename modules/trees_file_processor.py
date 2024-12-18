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
        mrcas_by_gene=[]
        mrcas_by_tree=[]
        for gene_center in gene_centers:
            tree_for_gene=reduced_ts.at(gene_center)
            if not tree_for_gene.has_multiple_roots:
                mrca = reduced_ts.node(tree_for_gene.root)
                mrcas_by_gene.append(mrca.time)

        tree_interval_starts=[]
        for tree in reduced_ts.trees():
            if not tree.has_multiple_roots:
                tree_interval_start=tree.interval.left
                mrca = reduced_ts.node(tree.root)
                mrcas_by_tree.append(mrca.time)
                tree_interval_starts.append(tree_interval_start)

        csv_file_out1 = os.path.join(output_folder, "simulated_ancestral_gene_mrcas.csv")
        save_mrca_values(csv_file_out1, mrcas_by_gene,gene_starts)
        theoretical_mrcas=theoretical_coalescent(num_genes,config.ancestral_Ne)
        csv_file_out2 = os.path.join(output_folder, "theoretical_ancestral_gene_mrcas.csv")
        save_mrca_values(csv_file_out2,theoretical_mrcas,gene_starts)
        csv_file_out3 = os.path.join(output_folder, "simulated_ancestral_tree_mrcas.csv")
        save_mrca_values(csv_file_out3,mrcas_by_tree,tree_interval_starts)
        png_out1 = os.path.join(output_folder, "mrca_hist1.png")
        plot_mrca(mrcas_by_gene,mrcas_by_tree,theoretical_mrcas, png_out1)
        png_out2 = os.path.join(output_folder, "mrca_hist2.png")
        plot_mrca(mrcas_by_gene,[],[], png_out2)

def theoretical_coalescent(num_genes,N):
    #Co.T=(1/2N)*e^-((t-1)/2N))
    loc=0
    xscale=2.0*N
    random_draws_from_distribution = expon.rvs(loc=loc, size=num_genes, scale=xscale)
    return random_draws_from_distribution

def plot_mrca(slim_mrcas_by_gene, slim_mrcas_by_tree, theoretical_mrcas_by_gene,
              png_out):
    fig = plt.figure(figsize=(10, 10), dpi=350)
    
    #Co.T=(1/2N)*e^-((t-1)/2N))
    x = slim_mrcas_by_gene
    label = "Coalescent Times For Genes From Sampled Two Ancestral Genomes"
    max_mrca = max(slim_mrcas_by_gene)
    num_genes=len(slim_mrcas_by_gene)
    num_segments=len(slim_mrcas_by_tree)
    bin_size = 100
    bins = np.arange(0, max_mrca, bin_size)

    plt.hist(x, bins=bins, facecolor='b', alpha=0.25,
                                label='SLiM Tcoal by gene (total: '+str(num_genes)+')',
                                density=False)

    if len(slim_mrcas_by_tree) > 0:
        plt.hist(slim_mrcas_by_tree, bins=bins, facecolor='r', alpha=0.25,
                 label='SLiM Tcoal by segment (total: '+str(num_segments)+')',
                 density=False)
    
    if len(theoretical_mrcas_by_gene) > 0:
        plt.hist(theoretical_mrcas_by_gene, bins=bins, facecolor='c', alpha=0.25,
                 label='Theoretical Tcoal by gene (total: '+str(num_genes)+')',
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

        f.writelines("IntervalStart,TimeToMRCA\n")
        for i in range(0,len(mrcas_values)):
            mrcas_value = mrcas_values[i]
            gene_start= gene_starts[i]
            f.writelines(str(gene_start) + "," + str(mrcas_value)+"\n")
