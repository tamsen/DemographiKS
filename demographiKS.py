import datetime
import math
import os
import shutil
import msprime
import sys
import tskit
import random
from datetime import datetime
import config
import version
import log
from modules import SLiM_runner, ks_calculator, FASTA_extracta, ks_histogramer, gene_shedder, trees_file_processor, \
    summary_stats


def run():

    conf = setup(sys.argv)
    if not conf:
        return

    # start the log
    log.write_start_to_log(conf.output_folder, conf.log_file_name, conf.version_info)
    log.write_to_log('Command Arguments Given: %s' % sys.argv)

    slim_out_folder = os.path.join(conf.output_folder, "SLiM_output")
    demographics_out_folder = os.path.join(conf.output_folder, "demographiKS_output")
    final_trees_file = os.path.join(slim_out_folder, conf.sim_name + "_trees.txt")
    trees_file_at_div = os.path.join(slim_out_folder, conf.sim_name + "_trees_at_div.txt")
    my_SLiM_allo_script = os.path.join("SLiM_scripts", "allotetraploid_bottleneck_trees.slim")
    my_SLiM_auto_script = os.path.join("SLiM_scripts", "autotetraploid_bottleneck_trees.slim")
    my_SLiM_allo_with_migration_script = os.path.join("SLiM_scripts", "allotetraploid_migration_trees.slim")
    my_SLiM_allo_with_assortative_mating_script = os.path.join("SLiM_scripts", "allotetraploid_bottleneck_inbreeding_trees.slim")

    out_fasta = os.path.join(demographics_out_folder, conf.sim_name + ".fa")
    ks_csv = os.path.join(demographics_out_folder, conf.sim_name + ".csv")
    ss_csv = os.path.join(demographics_out_folder, "summary_stats.csv")
    out_png = os.path.join(conf.output_folder, conf.sim_name + "_hist.png")

    folders_needed = [conf.output_folder, demographics_out_folder, slim_out_folder]
    for f in folders_needed:
        if not os.path.exists(f):
            os.mkdir(f)

    # Run the SLiM model
    log.write_to_log("Step 1: running SLiM")
    if conf.pre_existing_trees_file:
        final_trees_file = conf.pre_existing_trees_file
        log.write_to_log("Using pre-existing trees file:\t" + str(conf.pre_existing_trees_file))
    else:
        path_to_current_py_script = os.path.abspath(__file__)
        if conf.DIV_time_Ge:
            full_slim_script = os.path.join(os.path.dirname(path_to_current_py_script), my_SLiM_allo_script)
            if str(conf.mig_rate) != str(False):
                full_slim_script = os.path.join(os.path.dirname(path_to_current_py_script), my_SLiM_allo_with_migration_script)
            if str(conf.assortative_mating_coefficient) != str(False):
                full_slim_script = os.path.join(os.path.dirname(path_to_current_py_script),
                                                my_SLiM_allo_with_assortative_mating_script)

        else: #if there was no parental divergence, then this must be an autopolyploid
            full_slim_script = os.path.join(os.path.dirname(path_to_current_py_script), my_SLiM_auto_script)
        log.write_to_log("Running SLiM script:\t" + str(full_slim_script))
        SLiM_runner.run_slim(conf, final_trees_file, trees_file_at_div, full_slim_script)

    # select random ancestral genomes to calculate the T coalescent:
    genome_index_1 = 1
    genome_index_2 = 5
    trees_file_processor.plot_coalescent(trees_file_at_div, genome_index_1, genome_index_2,
                    conf, demographics_out_folder)

    if conf.stop_at_step < 2:
        return

    log.write_to_log("Step 2: Generating paralogs from trees file")
    log.write_to_log("Loading:\t" + str(final_trees_file))
    ts = tskit.load(final_trees_file)

    # log.write_to_log("SLiM metadata dict:\t" + str(metadata))
    num_individuals_if_diploid = ts.individuals_population.size
    num_polyploids = num_individuals_if_diploid / 2.0
    num_genomes = ts.num_samples
    num_genome_per_individual= num_genomes / num_polyploids
    log.write_to_log("size SLiM population:\t" + str(num_individuals_if_diploid))
    log.write_to_log("size SLiM samples:\t" + str(num_genomes))
    log.write_to_log("num polyploids in population:\t" + str(num_polyploids))
    log.write_to_log("num subgenomes being simulated:\t" + str(num_genomes))
    log.write_to_log("ploidy per polyploid:\t" + str(int(num_genome_per_individual)))

    # pick a random polyploid individual (ie, two random subgenomes from the two populations of parental subgenomes)
    num_genomes = conf.bottleneck_Ne * 2  # because diploid individuals, and thats SLiM default
    random.seed(conf.DemographiKS_random_seed)
    focal_genomes_as_int = [(random.randint(1, num_genomes)),
                            (random.randint(1 + num_genomes, 2 * num_genomes))]
    focal_genomes_as_str = ["n" + str(i) for i in focal_genomes_as_int]
    log.write_to_log("random focal polyploid individual:\t" + str(focal_genomes_as_str))

    # overlays neutral mutations
    mts = msprime.sim_mutations(ts, rate=conf.mutation_rate, random_seed=conf.Msprime_random_seed, keep=True)
    v_list = [v for v in mts.variants()]
    log.write_to_log(str(len(v_list)) + " mutations added.")

    log.write_to_log("Getting paralog names.")
    log.write_to_log("Gene length:\t" + str(conf.gene_length_in_bases))
    log.write_to_log("Total num bases:\t" + str(conf.total_num_bases))
    log.write_to_log("Max num paralogs:\t" + str(conf.max_num_paralogs_to_process))
    paralog_names = FASTA_extracta.get_paralog_names(conf.gene_length_in_bases, conf.total_num_bases,
                                                     conf.max_num_paralogs_to_process)
    log.write_to_log("Num paralogs before shedding:\t" + str(len(paralog_names)))
    genes_to_loose_a_duplicate = gene_shedder.decide_genes_to_shed(paralog_names, conf)


    #write out the fasta for our focal genomes
    log.write_to_log("Getting paralog sequences from TS data.")
    #https://tskit.dev/tskit/docs/stable/python-api.html
    wrap_width=0
    fasta_string = mts.as_fasta(reference_sequence=tskit.random_nucleotides(mts.sequence_length,
                                                                            seed=conf.Msprime_random_seed+1),
                                wrap_width=wrap_width)


    if conf.keep_intermediary_files:
        with open(out_fasta, "w") as f:
            f.write(fasta_string )
        log.write_to_log("Sequences written to FASTA file: " + out_fasta + ".")
    log.write_to_log("Final genomes complete")

    #memory intensive, so do it a different way..
    #fasta_io = StringIO(fasta_string)
    #SeqDict = SeqIO.to_dict(SeqIO.parse(fasta_io , "fasta"))
    Ks_values=[]

    #the different way...index straight into the fasta
    seq_dict={}
    for k in range(0,len(focal_genomes_as_int)):
        i=focal_genomes_as_int[k]
        num_places=int(math.log10(i))
        terms=[9*(m+1)*10**m for m in range(0,num_places)]
        index_update=1+sum(terms)+(i-10**num_places)*(num_places+1)+(i)*4
        next_to_add= 5 + num_places
        str_inx=i*conf.total_num_bases+index_update+next_to_add
        seq_dict[focal_genomes_as_str[k]]=fasta_string[str_inx:str_inx+conf.total_num_bases]

    if conf.stop_at_step < 3:
            return

    log.write_to_log("Step 3-5: Entering paralog-processing loop")
    with open(ss_csv, 'a') as f:
        f.write("paralog_ID\tpi\n")

    for paralog_ID in paralog_names:

        log.write_to_log("\tStep 3: Gene shedding for " + str(paralog_ID))

        #These would be shed the end of the sim anyway, so by deciding early on
        # what to shed, we dont loose any time processing them and then killing them.
        if paralog_ID in genes_to_loose_a_duplicate:
            log.write_to_log("\tGene shed: " + str(paralog_ID))
            log.write_to_log("\tSteps 4 & 5 skipped for shed gene " + str(paralog_ID))
            continue

        log.write_to_log("\tGene retained: " + str(paralog_ID))
        if conf.stop_at_step < 4:
            continue
        log.write_to_log("\tStep 4: Writing paralog sequences for " + str(paralog_ID))
        raw_sequences_for_paralog={}
        indexes_of_concern=[]
        for subgenome in focal_genomes_as_str:
            genome_name = conf.sim_name + "_" + subgenome
            subsequence = seq_dict[subgenome][paralog_ID:paralog_ID + conf.gene_length_in_bases]
            paralog_name = genome_name + "_paralog_" + str(paralog_ID)

            log.write_to_log("Writing data for: " + paralog_name + ".")
            idx_of_stop_codons=FASTA_extracta.get_nucleotide_index_of_any_STOP_codons_in_seq(
                conf.num_codons_in_a_gene,str(subsequence), conf.stop_codons)

            log.write_to_log("stop codons: " + ",".join([str(i) for i in idx_of_stop_codons]))
            if conf.keep_intermediary_files:
                out_per_genome_per_paralog_fasta = os.path.join(demographics_out_folder, paralog_name + ".fa")
                with open(out_per_genome_per_paralog_fasta, "w") as f:
                    f.writelines([">" + subgenome + " simulated paralogous gene\n",subsequence])
                    print("written:" + subsequence)

            raw_sequences_for_paralog[paralog_name]=str(subsequence)
            indexes_of_concern= indexes_of_concern+idx_of_stop_codons
    
        final_sequences_for_paralog={}
        unique_paralog_names=list(raw_sequences_for_paralog.keys())
        for paralog_name in unique_paralog_names:
            raw_seq= raw_sequences_for_paralog[paralog_name]
            fixed_subsequence = FASTA_extracta.replace_str_indexes(raw_seq,indexes_of_concern, "NNN")
            final_sequences_for_paralog[paralog_name] = fixed_subsequence

        #we have the sequences in memory at this point, we can calculate summary stats
        pi = summary_stats.pi_nuc_diversity(list(final_sequences_for_paralog.values()))
        with open(ss_csv, 'a') as f:
            summary_stats_data=str(paralog_ID) + "\t" + str(pi) + "\n"
            f.write(str(summary_stats_data))

        paralog_folder = os.path.join(demographics_out_folder, "paralog_" + str(paralog_ID))
        if not os.path.exists(paralog_folder):
            os.makedirs(paralog_folder)

        if conf.stop_at_step < 5:
            continue
        log.write_to_log("Step 5:\tRunnning CODEML on paralog " + str(paralog_ID))
        codeml_ML_dS_file = ks_calculator.run_CODEML_by_paralog(paralog_ID,final_sequences_for_paralog,
                                                                    paralog_folder)

        Ks_values_for_paralog, file_lines_for_paralog = ks_histogramer.extract_Ks_values_by_file(codeml_ML_dS_file)

        if not conf.keep_intermediary_files:
            shutil.rmtree(paralog_folder)

        with open(ks_csv, 'a') as f:
            f.writelines(file_lines_for_paralog)

        Ks_values=Ks_values+Ks_values_for_paralog

    log.write_to_log("Step 3-5: Paralog-processing loop complete")
    if conf.stop_at_step < 6:
        return

    log.write_to_log("Step 6:\tPlotting Ks histogram.")
    ks_histogramer.plot_Ks_histogram(out_png, conf, Ks_values, "ML", "b")

    log.write_end_to_log()
    return


def setup(arguments):

    print('Command Arguments Given: %s' % arguments)
    if len(arguments) < 2:
        print('Please give an input file path.')
        return False

    config_file=arguments[1]
    now = datetime.now()
    date_time = now.strftime("m%md%dy%Y_h%Hm%Ms%S")
    conf = config.DemographiKS_config(config_file)
    conf.output_folder = conf.output_folder_root + "_" + date_time
    conf.log_file_name = date_time + "_" + conf.log_file_name
    conf.version_info = version.version_info()
    cwd=os.getcwd()

    print('Config file: %s' % config_file)
    print("Current environment: %s" + str(os.environ))
    print("Current Working Directory:\t" + cwd)
    if conf.output_folder[0:2]== "./":
        conf.output_folder = os.path.join(os.getcwd(),conf.output_folder.replace("./",""))

    if (conf.mig_rate and conf.assortative_mating_coefficient):
        print("For simplicty, DemographiKS does not currently simulate migration and inbreeding.")
        print("Please choose one or the other. Thank you.")
        return False

    config_file_used=os.path.basename(config_file).replace(".xml",".used.xml")
    print("Output folder:\t" + conf.output_folder)
    if not os.path.exists(conf.output_folder):
        os.makedirs(conf.output_folder)

    #move a copy of the config file into the output folder so we remember what was run
    dst = os.path.join(conf.output_folder,config_file_used)
    shutil.copyfile(config_file, dst)

    return conf

if __name__ == '__main__':
    run()
