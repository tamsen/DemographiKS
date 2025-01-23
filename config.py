# eventually this should be read in from an xml config or similar
import math
import xml.etree.ElementTree as ET


class DemographiKS_config:

    output_folder_root = "/home/DemographiKS_output"
    sim_name = "allotetraploid_bottleneck"
    pre_existing_trees_file =  False
    test_burnin = False
    stop_codons=["TAA","TGA","TAG"]
    num_codons_in_a_gene=1000
    len_codon=3
    recombination_rate=1.25*10**(-6)
    mutation_rate=1*10**(-5)
    gene_length_in_bases= num_codons_in_a_gene * len_codon #nucleotides
    stop_codons=["TAA","TGA","TAG"]
    max_num_paralogs_to_process=20 #per genome
    log_file_name = "log.txt"
    total_num_bases = 2000
    ancestral_Ne = 100
    bottleneck_Ne = 20
    DIV_time_Ge=4000
    WGD_time_Ge=2000

    # typically 2*10*Ne, but you should check how long it takes for SLiM
    # ancestral population to be at a steady state
    burnin_time=2*10*100 #

    #In specs WGD_half_life_MY was 31 my.
    #A half-life of 31 MY = mean life span of 44.723 MY
    #So, thats 44.723*10**6 GE if one year = 1 generation
    mean_WGD_life_span_in_GE= 44.723*10**6

    #randomness
    SLiM_rep =1
    Msprime_random_seed = 42
    DemographiKS_random_seed = 17

    #debugging
    stop_at_step = 999
    keep_intermediary_files=False

    def __init__(self, config_file):

        mytree = ET.parse(config_file)
        myroot = mytree.getroot()

        for top_layer in myroot:

                incoming_tag = top_layer.tag.strip()
                incoming_txt = top_layer.text.strip()

                if (incoming_tag == "StopAtStep"):
                    self.stop_at_step = int(incoming_txt)

                if (incoming_tag == "LogFileName"):
                    if incoming_txt.upper() == "FALSE":
                        self.log_file_name = False
                    else:
                        self.log_file_name = incoming_txt
                if (incoming_tag == "KeepIntermediaryFiles"):
                    if incoming_txt.upper() == "FALSE":
                        self.keep_intermediary_files = False
                    else:
                        self.keep_intermediary_files = incoming_txt

                if (incoming_tag == "Paths"):
                    for inner_layer in top_layer:
                        incoming_txt = inner_layer.text.strip()
                        incoming_tag = inner_layer.tag.strip()
                        if (incoming_tag == "output_folder_root"):
                            self.output_folder_root = incoming_txt
                        if (incoming_tag == "pre_existing_trees_file"):
                            if incoming_txt.upper() == "FALSE":
                                self.pre_existing_trees_file=False
                            else:
                                self.pre_existing_trees_file = incoming_txt

                if (incoming_tag == "Chromosome"):
                    for inner_layer in top_layer:
                        incoming_txt = inner_layer.text.strip()
                        incoming_tag = inner_layer.tag.strip()

                        if (incoming_tag == "recombination_rate"):
                            self.recombination_rate = float(incoming_txt)
                        if (incoming_tag == "total_num_bases"):
                            self.total_num_bases = int(incoming_txt)
                        if (incoming_tag == "num_codons_in_a_gene"):
                            self.num_codons_in_a_gene = int(incoming_txt)
                            self.gene_length_in_bases = self.num_codons_in_a_gene * self.len_codon
                        if (incoming_tag == "max_num_paralogs_to_process"):
                            self.max_num_paralogs_to_process = parse_int_or_false(incoming_txt)

                if (incoming_tag == "SequenceEvolution"):
                    for inner_layer in top_layer:
                        incoming_txt = inner_layer.text.strip()
                        incoming_tag = inner_layer.tag.strip()
                        if (incoming_tag == "mutation_rate"):
                            self.mutation_rate = float(incoming_txt)
                        if (incoming_tag == "WGD_gene_half_life_in_GE"):
                            self.mean_WGD_life_span_in_GE = half_life_to_mean_life(incoming_txt)
                        if (incoming_tag == "num_codons"):
                            self.num_codons = int(incoming_txt)
                        if (incoming_tag == "Ks_per_Myr"):
                            self.Ks_per_Myr = float(incoming_txt)
                        if (incoming_tag == "per_site_evolutionary_distance"):
                            self.per_site_evolutionary_distance = float(incoming_txt)
                        if (incoming_tag == "evolver_random_seed"):
                            self.evolver_random_seed = int(incoming_txt)

                if (incoming_tag == "Population"):
                    for inner_layer in top_layer:
                        incoming_txt = inner_layer.text.strip()
                        incoming_tag = inner_layer.tag.strip()
                        if (incoming_tag == "ancestral_ne"):
                            self.ancestral_Ne= int(incoming_txt)
                        if (incoming_tag == "bottleneck_ne"):
                            self.bottleneck_Ne = int(incoming_txt)
                        if (incoming_tag == "burnin_time"):
                            self.burnin_time = int(incoming_txt)

                if (incoming_tag == "Speciation"):
                    for inner_layer in top_layer:
                        incoming_txt = inner_layer.text.strip()
                        incoming_tag = inner_layer.tag.strip()
                        if (incoming_tag == "DIV_time_Ge"):
                            self.DIV_time_Ge = float(incoming_txt)
                        if (incoming_tag == "WGD_time_Ge"):
                            self.WGD_time_Ge = float(incoming_txt)

                if (incoming_tag == "Randomization"):
                    for inner_layer in top_layer:
                        incoming_txt = inner_layer.text.strip()
                        incoming_tag = inner_layer.tag.strip()
                        if (incoming_tag == "SLiM_rep"):
                            self.SLiM_rep = int(incoming_txt)
                        if (incoming_tag == "Msprime_random_seed"):
                            self.Msprime_random_seed = int(incoming_txt)
                        if (incoming_tag == "DemographiKS_random_seed"):
                            self.DemographiKS_random_seed = int(incoming_txt)

        #we are going with Ks_per_YR = mut rate * (1/1.2)
        # we multiply by 1/1.2 since thats syn / total mut rate
        self.Ks_per_YR = 0.833 * self.mutation_rate
        self.mean_Ks_from_Tc = 2.0 * self.ancestral_Ne * self.Ks_per_YR
        self.num_genes = int(self.total_num_bases / self.gene_length_in_bases)
        self.t_div_as_ks = self.DIV_time_Ge * self.Ks_per_YR

def parse_tuple_string(tuple_string):
    if tuple_string.upper() == "FALSE":
        return False
    else:
        splat = tuple_string.replace("(", "").replace(")", "").split(",")
        data = [float(s) for s in splat]
        return data

def parse_float_or_false(input_string):
    if input_string.upper() == "FALSE":
        return False
    else:
        return float(input_string)

def parse_int_or_false(input_string):
    if input_string.upper() == "FALSE":
        return False
    else:
        return int(input_string)

def half_life_to_mean_life(input_string):

    half_life=parse_float_or_false(input_string)
    if not half_life:
        return False
    else:
        ln2=math.log(2)
        return half_life/ln2

def parse_comma_separated_values(input_string):
    if input_string.upper() == "FALSE":
        return False
    else:
        cleaned_string=input_string.replace("(", "").replace(")", "")
        splat = cleaned_string.split(",")
        return splat
