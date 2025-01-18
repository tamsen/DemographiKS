import datetime
import os
import shutil
import msprime
import sys
import tskit
import random
from io import StringIO
import config
import modules.FASTA_extracta
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import log
from datetime import datetime

from modules.FASTA_extracta import get_sequences_by_paralog_name
from modules.gene_shedder import shed_genes,decide_genes_to_shed
from modules import SLiM_runner, ks_calculator, FASTA_extracta, ks_histogramer
from modules.trees_file_processor import plot_coalescent


def run(conf):

    # start the log
    log.write_start_to_log(conf.output_folder, conf.log_file_name, conf.version_info)
    log.write_to_log('Command Arguments Given: %s' % sys.argv)

    slim_out_folder = os.path.join(conf.output_folder, "SLiM_output")
    demographics_out_folder = os.path.join(conf.output_folder, "demographiKS_output")
    final_trees_file = os.path.join(slim_out_folder, conf.sim_name + "_trees.txt")
    trees_file_at_div = os.path.join(slim_out_folder, conf.sim_name + "_trees_at_div.txt")
    my_SLiM_script = os.path.join("SLiM_scripts", "allotetraploid_bottleneck_trees.slim")

    out_fasta = os.path.join(demographics_out_folder, conf.sim_name + ".fa")
    out_csv = os.path.join(demographics_out_folder, conf.sim_name + ".csv")
    out_png = os.path.join(conf.output_folder, conf.sim_name + "_hist.png")

    folders_needed = [conf.output_folder, demographics_out_folder, slim_out_folder]
    for f in folders_needed:
        if not os.path.exists(f):
            os.mkdir(f)

    # Run the SLiM model
    log.write_to_log("Step 1: running SLiM")
    if conf.pre_existing_trees_file:
        final_trees_file = conf.pre_existing_trees_file
        log.write_to_log("Using pre-exisiting trees file:\t" + str(conf.pre_existing_trees_file))
    else:
        path_to_current_py_script = os.path.abspath(__file__)
        full_slim_script = os.path.join(os.path.dirname(path_to_current_py_script), my_SLiM_script)
        log.write_to_log("Running SLiM script:\t" + str(full_slim_script))
        SLiM_runner.run_slim(conf, final_trees_file, trees_file_at_div, full_slim_script)

    # random ancestral genomes:
    genome_index_1 = 1
    genome_index_2 = 5

    plot_coalescent(trees_file_at_div, genome_index_1, genome_index_2,
                    conf, demographics_out_folder)

    if conf.stop_at_step < 2:
        return

    log.write_to_log("Step 2: Generating paralogs from trees file")
    log.write_to_log("Loading:\t" + str(final_trees_file))
    ts = tskit.load(final_trees_file)
    metadata = ts.metadata["SLiM"]
    # log.write_to_log("SLiM metadata dict:\t" + str(metadata))
    log.write_to_log("size SLiM population:\t" + str((ts.individuals_population.size)))
    log.write_to_log("size SLiM samples:\t" + str((ts.num_samples)))

    # pick a random polyploid individual (ie, two random subgenomes from the two populations of parental subgenomes)
    num_genomes = conf.bottleneck_Ne * 2  # because diploid individuals
    random.seed(conf.DemographiKS_random_seed)
    focal_genomes = ["n" + str(random.randint(1, num_genomes)),
                     "n" + str(random.randint(1 + num_genomes, 2 * num_genomes))]
    log.write_to_log("random focal polyploid individual:\t" + str(focal_genomes))

    # overlays neutral mutations
    mts = msprime.sim_mutations(ts, rate=conf.mutation_rate, random_seed=conf.Msprime_random_seed, keep=True)
    v_list = [v for v in mts.variants()]
    log.write_to_log(str(len(v_list)) + " mutations added.")

    paralog_names = get_sequences_by_paralog_name(conf.gene_length, conf.total_num_bases,
                                                  conf.max_num_paralogs_to_process)
    genes_to_loose_a_duplicate = decide_genes_to_shed(paralog_names, conf)

    out_csv_2 = os.path.join(demographics_out_folder, conf.sim_name + "_2.csv")
    random_nuceotides_seed=42

    log.write_to_log("Getting paralog sequences from TS data.")
    cleaned_sequences_by_paralog_name_dict = FASTA_extracta.extract_paralog_sequences(demographics_out_folder,
                                                                                      focal_genomes,
                                                                                      conf, mts, out_fasta)



    #write out the fasta for out focal genomes
    fasta_string = mts.as_fasta(reference_sequence=tskit.random_nucleotides(mts.sequence_length,
                                                                            seed=random_nuceotides_seed))
    fasta_io = StringIO(fasta_string)
    SeqDict = SeqIO.to_dict(SeqIO.parse(fasta_io , "fasta"))
    Ks_values=[]
    for paralog_ID in paralog_names:

        raw_sequences_for_paralog={}
        indexes_of_concern=[]
        for subgenome in focal_genomes:
            genome_name = conf.sim_name + "_" + subgenome
            subsequence=SeqDict[subgenome][paralog_ID:paralog_ID+conf.gene_length].seq
            paralog_name = genome_name + "_paralog_" + str(paralog_ID)
            out_per_genome_per_paralog_fasta = os.path.join(demographics_out_folder, paralog_name + "_foo.fa")
            print("Writing data for : " + paralog_name + ".")
            idx_of_stop_codons=FASTA_extracta.get_nucleotide_index_of_any_STOP_codons_in_seq(
                conf.num_codons_in_a_gene,str(subsequence), conf.stop_codons)

            print("stop codons: " + ",".join([str(i) for i in idx_of_stop_codons]))
            #print("fixed_subsequence: " + fixed_subsequence)
            record = SeqRecord(subsequence,
                       id=subgenome, name=paralog_name,
                       description="simulated paralogous gene")
            SeqIO.write(record, out_per_genome_per_paralog_fasta, "fasta")
            raw_sequences_for_paralog[subgenome]=str(subsequence)
            indexes_of_concern= indexes_of_concern+idx_of_stop_codons
    
        final_sequences_for_paralog={}
        for subgenome in focal_genomes:
            raw_seq= raw_sequences_for_paralog[subgenome]
            fixed_subsequence = FASTA_extracta.replace_str_indexes(raw_seq,indexes_of_concern, "NNN")
            final_sequences_for_paralog[subgenome] = fixed_subsequence

        if paralog_ID in genes_to_loose_a_duplicate:
            log.write_to_log("Step 4:\tshedding paralog" + str(paralog_ID))
            continue

        log.write_to_log("Step 4:\tRunnning CODEML on paralog" + str(paralog_ID))
        codeml_ML_dS_file = ks_calculator.run_CODEML_by_paralog(paralog_ID,final_sequences_for_paralog,
                                                                    demographics_out_folder)

        Ks_values_for_paralog, file_lines_for_paralog = ks_histogramer.extract_Ks_values_by_file(codeml_ML_dS_file)
        with open(out_csv_2, 'a') as f:
            f.writelines(file_lines_for_paralog)

        Ks_values=Ks_values+Ks_values_for_paralog
        # else: the paralog was shed


    ks_histogramer.plot_Ks_histogram(out_png, conf, Ks_values,
                                     None, None, None, None, "ML", "b", 0.001)

    log.write_end_to_log()

    return
