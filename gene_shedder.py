import random
from scipy.stats import expon


def shed_genes(avg_WGD_gene_lifespan, sequences_by_paralog_name_dict, time_since_WGD):

    fraction_WGD_genes_remaining_at_time_since_WGD = avg_WGD_gene_lifespan * expon.pdf(
        time_since_WGD, loc=0, scale=avg_WGD_gene_lifespan)
    all_sequences = list(sequences_by_paralog_name_dict.keys())
    num_original_sequences = float(len(all_sequences))
    num_genes_to_shed = int(num_original_sequences * (1.0 - fraction_WGD_genes_remaining_at_time_since_WGD))
    # random.seed(polyploid.general_sim_config.specks_random_seed)
    random.seed(42)
    gene_trees_to_loose_a_duplicate_gene = random.sample(all_sequences, num_genes_to_shed)

    print("original num genes: " + str(num_original_sequences))
    print("num to shed genes: " + str(num_genes_to_shed))
    print("genes shed: " + str(list(gene_trees_to_loose_a_duplicate_gene)))

    for gene in gene_trees_to_loose_a_duplicate_gene:
        del sequences_by_paralog_name_dict[gene]

    print("remaining genes: " + str(len(sequences_by_paralog_name_dict)))
    return sequences_by_paralog_name_dict
