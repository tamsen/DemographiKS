import sys
import random
from utilities import birth_and_death_with_escape

def run():
    args = sys.argv
    path_to_input_csv=args[1]
    path_to_output_csv = path_to_input_csv.replace(".csv", "_with_SSD_model.csv")

    rseed=10
    leaf_to_leaf_Ks=2.0

    fraction_genes_maintained = 0.6
    if len(sys.argv)>2:
        fraction_genes_maintained=float(args[2])

    # Note:
    #num SSD genes to pull from SSD model = int(fraction_genes_maintained * len(demographiKS_ks_results))

    run_results = read_Ks_csv(path_to_input_csv, False)
    run_results = [k * leaf_to_leaf_Ks for k in run_results]
    demographiKS_ks_results =  run_results
    demographiKS_plus_selection_results = add_selection(demographiKS_ks_results, fraction_genes_maintained, rseed)
    with open(path_to_output_csv, "w", encoding="utf-8") as f:
        for result in demographiKS_plus_selection_results:
            f.write(str(result)+ "\n")
    print("SSDs added")


def add_selection(demographiKS_ks_results, fraction_genes_maintained, rseed):
    step_size = 0.001
    decay_constant = 33  # =(1/mean life expectancy of SSD gene), in Ks space
    rate_escape = 0.3  # % of SSD are maintained

    # where decay constant comes from:
    #
    #    #(Ferris and Whitt, 1977; Nadeau and Sankoff, 1997; Lynch and Conery, 2000;
    #    # #Blanc and Wolfe, 2004b; Soltis et al.,2016; Cheng et al., 2018).
    #    # if a SSD lasts a few million years (say 3x10^6)
    #    # In Ks space, that depends mutation rate.
    #    # if Mutation rate is 10^-8, then
    #    # SSD lasts a few million years (say 3x10^-2 ) in KS
    #    # so, decay constant is 1/(3x10^-2)

    # where rate escape comes from
    #  Proportion of retained duplicates (what portion of the overall genome is retained duplicates):
    #  between 0.3, 0.4  human and mouse https://pmc.ncbi.nlm.nih.gov/articles/PMC1413713/
    #  Not the same as plants, but at least we are sure those numbers not confounded w/ WGD

    # Retention of duplicated genes in evolution (20-40% of the whole genome)
    # so, a genome of size 100 has 20-40 retained duplcatate SSDs
    # if they were born randomly over 574 MY (~KS=4), that's 0.035-0.070 gene retained per MY
    # from a starting genome of size 100 genes.
    # So, thats (0.035 - 0.07)/100, so 0.00035 to 0.0007 dup retained per gene million year.
    # Gene birth rate is ~0.00162 per g per million years (Tiley)
    # If 0.00162 are born, and 0.00035 to 0.0007 survive, then
    # ~0.21-0.43 survive. Wow. more than I expected.
    # refs ()

    random.seed(rseed)
    include_debugging_plots = True
    my_pdf, xs = birth_and_death_with_escape.gene_birth_death_with_escape_pdf(
        step_size, decay_constant, rate_escape, rseed, include_debugging_plots)

    num_genes_needed = int(fraction_genes_maintained * len(demographiKS_ks_results))
    SSD_ks = birth_and_death_with_escape.draw_SSDs_from_pdf(my_pdf, xs,
                                                            num_genes_needed,
                                                            rseed, include_debugging_plots)

    print("Num SSD genes simulated: " + str(len(SSD_ks)))

    demographiKS_plus_selection_results = demographiKS_ks_results + SSD_ks

    return demographiKS_plus_selection_results


def read_Ks_csv(csv_file, expect_header):

    ks_results = []
    with open(csv_file, "r") as f:

        reading_header=expect_header
        while True:
            line = f.readline()
            if "ersion" in line:
                continue
            if "Git" in line:
                continue
            if "leaf names" in line:
                continue
            if not line:
                break
            if len(line)==0:
                break
            if reading_header:
                reading_header=False
                continue
            data = line.split(",")
            if len(data[2])>0:
                ks_value=float(data[2])
                #if ks_value > 0:
                ks_results.append(ks_value)

    return ks_results

if __name__ == '__main__':
    run()
