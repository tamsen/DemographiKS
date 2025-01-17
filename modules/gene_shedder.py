import random
from scipy.stats import expon
import log

def decide_genes_to_shed(paralog_names, config):

    seed=config.DemographiKS_random_seed
    avg_WGD_gene_lifespan=config.mean_WGD_life_span_in_GE
    time_since_WGD=config.WGD_time_Ge
    log.write_to_log("time since WGD (generations): " + str(time_since_WGD))
    log.write_to_log("avg_WGD_gene_lifespan (generations): " + str(avg_WGD_gene_lifespan))
    fraction_WGD_genes_remaining_at_time_since_WGD = avg_WGD_gene_lifespan * expon.pdf(
        time_since_WGD, loc=0, scale=avg_WGD_gene_lifespan)

    all_sequences = paralog_names #list(sequences_by_paralog_name_dict.keys())
    num_original_sequences = float(len(all_sequences))
    num_genes_to_shed = int(num_original_sequences * (1.0 - fraction_WGD_genes_remaining_at_time_since_WGD))
    random.seed(seed + 1)
    gene_trees_to_loose_a_duplicate_gene = random.sample(all_sequences, num_genes_to_shed)

    log.write_to_log("original num genes: " + str(num_original_sequences))
    log.write_to_log("num to shed genes: " + str(num_genes_to_shed))
    log.write_to_log("genes shed: " + str(list(gene_trees_to_loose_a_duplicate_gene)))

    return gene_trees_to_loose_a_duplicate_gene


def shed_genes(sequences_by_paralog_name_dict, config):

    seed=config.DemographiKS_random_seed
    avg_WGD_gene_lifespan=config.mean_WGD_life_span_in_GE
    time_since_WGD=config.WGD_time_Ge
    log.write_to_log("time since WGD (generations): " + str(time_since_WGD))
    log.write_to_log("avg_WGD_gene_lifespan (generations): " + str(avg_WGD_gene_lifespan))
    fraction_WGD_genes_remaining_at_time_since_WGD = avg_WGD_gene_lifespan * expon.pdf(
        time_since_WGD, loc=0, scale=avg_WGD_gene_lifespan)

    all_sequences = list(sequences_by_paralog_name_dict.keys())
    num_original_sequences = float(len(all_sequences))
    num_genes_to_shed = int(num_original_sequences * (1.0 - fraction_WGD_genes_remaining_at_time_since_WGD))
    random.seed(seed + 1)
    gene_trees_to_loose_a_duplicate_gene = random.sample(all_sequences, num_genes_to_shed)

    log.write_to_log("original num genes: " + str(num_original_sequences))
    log.write_to_log("num to shed genes: " + str(num_genes_to_shed))
    log.write_to_log("genes shed: " + str(list(gene_trees_to_loose_a_duplicate_gene)))

    for gene in gene_trees_to_loose_a_duplicate_gene:
        del sequences_by_paralog_name_dict[gene]

    print("remaining genes: " + str(len(sequences_by_paralog_name_dict)))
    return sequences_by_paralog_name_dict
