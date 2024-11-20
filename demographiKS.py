import datetime
import os
import shutil
import msprime
import sys
import tskit
import random
import config
import version
import log
from datetime import datetime
from scipy.stats import expon

from gene_shedder import shed_genes
from modules import SLiM_runner,ks_calculator, FASTA_extracta,ks_histogramer

def run_sim():

    conf = setup(sys.argv)
    if not conf:
        return

    #start the log
    log.write_start_to_log(conf.output_folder,conf.log_file_name, conf.version_info)
    log.write_to_log('Command Arguments Given: %s' % sys.argv)

    slim_out_folder=os.path.join(conf.output_folder,"SLiM_output")
    demographics_out_folder=os.path.join(conf.output_folder,"demographiKS_output")
    trees_file = os.path.join(slim_out_folder,conf.sim_name + "_trees.txt")
    my_SLiM_script= os.path.join("SLiM_scripts", "allotetraploid_bottleneck_trees.slim")

    out_fasta=os.path.join(demographics_out_folder,conf.sim_name + ".fa")
    out_csv=os.path.join(demographics_out_folder,conf.sim_name + ".csv")
    out_png=os.path.join(conf.output_folder,conf.sim_name + "_hist.png")

    folders_needed=[conf.output_folder,demographics_out_folder,slim_out_folder]
    for f in folders_needed:
        if not os.path.exists(f):
            os.mkdir(f)

    # Run the SLiM model

    if conf.pre_existing_trees_file:
        trees_file=conf.pre_existing_trees_file
        log.write_to_log("Using pre-exisitng trees file:\t" + str(conf.pre_existing_trees_file))
    else:
        path_to_current_py_script = os.path.abspath(__file__)
        full_slim_script = os.path.join( os.path.dirname(path_to_current_py_script), my_SLiM_script)
        log.write_to_log("Running SLiM:\t" + str(full_slim_script))
        SLiM_runner.run_slim(conf,trees_file, full_slim_script)

    log.write_to_log("Loading:\t" + str(trees_file))
    ts = tskit.load(trees_file)
    metadata=ts.metadata["SLiM"]
    #log.write_to_log("SLiM metadata dict:\t" + str(metadata))
    log.write_to_log("size SLiM population:\t" + str((ts.individuals_population.size)))
    log.write_to_log("size SLiM samples:\t" + str((ts.num_samples)))

    # pick a random polyploid indiviual (ie, two random subgenomes from the two populations of parental subgenomes)
    num_genomes=conf.bottleneck_Ne*2 #because diploid individuals
    random.seed(conf.DemographiKS_random_seed)
    focal_genomes = ["n" +str(random.randint(1,num_genomes)),
                     "n" +str(random.randint(1+num_genomes, 2*num_genomes))]
    log.write_to_log("random focal polyploid individual:\t" + str(focal_genomes))

    #overlays neutral mutations
    mts = msprime.sim_mutations(ts, rate=conf.mutation_rate, random_seed=conf.Msprime_random_seed, keep=True)
    v_list = [v for v in mts.variants()]
    log.write_to_log(str(len(v_list)) + " mutations added.")

    log.write_to_log("Getting paralog sequences from TS data.")
    cleaned_sequences_by_paralog_name_dict = FASTA_extracta.extract_paralog_sequences(demographics_out_folder,
                                                                                      focal_genomes,
                                                                       conf, mts, out_fasta)

    log.write_to_log("Shedding paralogs lost to dpiloidization & fractionation")
    remaining_sequences_by_paralog_name_dict  = shed_genes(cleaned_sequences_by_paralog_name_dict,conf)

    log.write_to_log("Runnning CODEML on remaining paralogs.")
    paml_out_files = ks_calculator.run_CODEML_on_paralogs(
        remaining_sequences_by_paralog_name_dict, demographics_out_folder)

    log.write_to_log("Extracting Ks values from PAML.")
    results = ks_histogramer.extract_K_values(out_csv, paml_out_files)
    ks_histogramer.plot_Ks_histogram(out_png, conf,results ,
                      None,None,None,None,"ML","b", 0.001)

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

    config_file_used=os.path.basename(config_file).replace(".xml",".used.xml")
    print("Output folder:\t" + conf.output_folder)
    if not os.path.exists(conf.output_folder):
        os.makedirs(conf.output_folder)

    #move a copy of the config file into the output folder so we remember what was run
    dst = os.path.join(conf.output_folder,config_file_used)
    shutil.copyfile(config_file, dst)

    return conf

if __name__ == '__main__':
    run_sim()
