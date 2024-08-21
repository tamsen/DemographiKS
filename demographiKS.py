import os
import subprocess, pyslim
import msprime
import numpy as np
import tskit
from cairosvg import svg2png
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from ks_calculator import sequences_to_codeml_in, run_codeml
from ks_histogramer import get_Ks_from_file, extract_K_values, plot_Ks_histogram


def run_sim():

    #conf = setup(sys.argv)
    #if not conf:
    #    return

    num_codons_in_a_gene=1000
    len_codon=3
    gene_length=num_codons_in_a_gene*len_codon #nucleotides
    stop_codons=["TAA","TGA","TAG"]
    max_num_paralogs_to_process=20 #per genome
    focal_genomes=["n11","n245"] #pick two randomly

    slim_out_folder="SLiM_output"
    demographics_out_folder="demographiKS_output"
    #sim_name="diploid_snm"
    #trees_file = os.path.join(slim_out_folder,"diploid_trees.txt")
    sim_name = "allotetraploid_bottleneck"
    trees_file = os.path.join(slim_out_folder,"allotetraploid_trees.txt")
    my_SLiM_script= os.path.join("SLiM_scripts", "allotetraploid_bottleneck_trees.slim")
    out_fasta=os.path.join(demographics_out_folder,sim_name + ".fa")
    out_csv=os.path.join(demographics_out_folder,sim_name + ".csv")
    out_png=os.path.join(demographics_out_folder,sim_name + "_hist.png")

    #TODO:
    # mk slim_out_folder, demographics_out_folder if they dont exist already

    # Run the SLiM model
    print("Running SLiM:\t" + str(my_SLiM_script))
    subprocess.check_output(["slim", "-m", "-s", "0", my_SLiM_script])

    print("Loading:\t" + str(trees_file))
    ts = tskit.load(trees_file)
    metadata=ts.metadata["SLiM"]
    n_populations=ts.nodes_population[1:10]
    n_individuals=ts.nodes_individual[1:10]
    print("SLiM metadata dict:\t" + str(metadata))
    print("size SLiM population:\t" + str((ts.individuals_population.size)))
    print("size SLiM samples:\t" + str((ts.num_samples)))
    #print("size SLiM individuals:\t" + str((ts.individuals_population)))
    #print("Tree Seq Max Root Time:\t" + str(ts.max_root_time))
    individuals_to_population_map= dict(zip([i for i in range(0,len(ts.individuals_population))], ts.individuals_population))
    print("individuals_to_population_map:\t" + str(individuals_to_population_map))
    #overlays neutral mutations
    mts = msprime.sim_mutations(ts, rate=1e-5, random_seed=42, keep=True)
    v_list = [v for v in mts.variants()]
    print(str(len(v_list)) + " mutations added.")

    #giant string...
    result =mts.as_fasta(reference_sequence=tskit.random_nucleotides(mts.sequence_length))
    with open(out_fasta, "w") as f:
        f.write(result)

    print("Sequences written to FASTA file: " + out_fasta + ".")

    sequences_by_paralog_name_dict = write_per_genome_per_paralog_fastas(demographics_out_folder, focal_genomes,
                                                                         gene_length, max_num_paralogs_to_process,
                                                                         out_fasta, sim_name)

    print("Removing STOP codons. PAML needs sequences that code for AA only")

    problem_codon_indexes_by_paralog_name_dict = get_index_of_any_STOP_codons(num_codons_in_a_gene,
                                                                              sequences_by_paralog_name_dict,
                                                                              stop_codons)

    cleaned_sequences_by_paralog_name_dict = set_STOP_codons_to_NNN(len_codon,
                                                                    problem_codon_indexes_by_paralog_name_dict,
                                                                    sequences_by_paralog_name_dict)

    print("Runnning CODEML on paralogs.")
    paml_out_files = run_CODEML_on_paralogs(cleaned_sequences_by_paralog_name_dict, demographics_out_folder)

    print("Extracting Ks values from PAML.")
    results = extract_K_values(out_csv, paml_out_files)
    plot_Ks_histogram(out_png, sim_name,results ,
                      None,None,None,None,"ML","b", 0.001)
    print(results)

    print("Done.")
    return


def write_per_genome_per_paralog_fastas(demographics_out_folder, focal_genomes, gene_length,
                                        max_num_paralogs_to_process, out_fasta, sim_name):
    sequences_by_paralog_name_dict = {}
    with open(out_fasta) as f:

        for seq_record in SeqIO.parse(f, 'fasta'):
            seq_record.id = seq_record.description = seq_record.id.replace('.seq', '')
            if seq_record.id in focal_genomes:
                start_index_in_sequence = 0
                num_paralogs_processed = 0
                genome_name = sim_name + "_" + seq_record.id
                seq = seq_record.seq
                full_seq_length = len(seq)

                while True:

                    paralog_name = genome_name + "_paralog_" + str(start_index_in_sequence)
                    end_index_in_sequence = start_index_in_sequence + gene_length

                    if end_index_in_sequence >= full_seq_length:
                        break
                    if num_paralogs_processed >= max_num_paralogs_to_process:
                        break

                    if start_index_in_sequence not in sequences_by_paralog_name_dict:
                        sequences_by_paralog_name_dict[start_index_in_sequence] = {}

                    subsequence = seq[start_index_in_sequence:end_index_in_sequence]
                    sequences_by_paralog_name_dict[start_index_in_sequence][paralog_name] = subsequence

                    # print("Subsequence : " + subsequence)
                    start_index_in_sequence = start_index_in_sequence + gene_length
                    out_per_genome_per_paralog_fasta = os.path.join(demographics_out_folder, paralog_name + ".fa")
                    # print("Writing data for : " + out_per_genome_fasta + ".")
                    record = SeqRecord(subsequence,
                                       id=seq_record.id, name=paralog_name,
                                       description="simulated paralogous gene")
                    SeqIO.write(record, out_per_genome_per_paralog_fasta, "fasta")
                    num_paralogs_processed = num_paralogs_processed + 1
    return sequences_by_paralog_name_dict


def set_STOP_codons_to_NNN(len_codon, problem_codon_indexes_by_paralog_name_dict, sequences_by_paralog_name_dict):
    cleaned_sequences_by_paralog_name_dict = {}
    for paralog_key, sequences_dict in sequences_by_paralog_name_dict.items():
        cleaned_sequences_by_paralog_name_dict[paralog_key] = {}
    for paralog_key, cleaned_sequences_dict in cleaned_sequences_by_paralog_name_dict.items():
        problem_codon_indexes = problem_codon_indexes_by_paralog_name_dict[paralog_key]
        for seq_key, sequence in sequences_by_paralog_name_dict[paralog_key].items():
            # print("problem codons " + str(problem_codon_indexes))
            # print("original sequence:\t" + str(sequence))
            revised_seq = str(sequence)
            for problem_codon_index in problem_codon_indexes:
                revised_seq = revised_seq[:problem_codon_index * len_codon] + "NNN" + revised_seq[(
                                                                                                              problem_codon_index + 1) * len_codon:]
            # print("revised sequence :\t" + revised_seq)
            cleaned_sequences_by_paralog_name_dict[paralog_key][seq_key] = revised_seq
    return cleaned_sequences_by_paralog_name_dict


def get_index_of_any_STOP_codons(num_codons_in_a_gene, sequences_by_paralog_name_dict, stop_codons):
    problem_codon_indexes_by_paralog_name_dict = {}
    for paralog_key, sequences_dict in sequences_by_paralog_name_dict.items():
        print("checking paralog " + str(paralog_key))
        problem_codon_indexes = []
        for seq_key, sequence in sequences_dict.items():
            for i in range(0, num_codons_in_a_gene):
                codon = str(sequence[3 * i:3 * (i + 1)])
                if codon in stop_codons:
                    problem_codon_indexes.append(i)
        problem_codon_indexes_by_paralog_name_dict[paralog_key] = problem_codon_indexes
    return problem_codon_indexes_by_paralog_name_dict


def run_CODEML_on_paralogs(cleaned_sequences_by_paralog_name_dict, demographics_out_folder):
    paml_out_files = []
    for paralog in cleaned_sequences_by_paralog_name_dict:

        paralog_folder = os.path.join(demographics_out_folder, "paralog_" + str(paralog))
        if not os.path.exists(paralog_folder):
            os.makedirs(paralog_folder)

        paralog_fa_file = os.path.join(paralog_folder, "paralog_" + str(paralog) + ".fa")
        codeml_input_fa_file = sequences_to_codeml_in(cleaned_sequences_by_paralog_name_dict[paralog], paralog_fa_file)
        print("codeml_input_fa_file written to " + codeml_input_fa_file)

        # need "conda install -c bioconda paml"
        result = run_codeml(codeml_input_fa_file, paralog_folder)
        paml_out_files.append(result.ML_dS_file)
    return paml_out_files


if __name__ == '__main__':
    run_sim()
