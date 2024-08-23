import datetime
import os
import shutil
import msprime
import sys
import tskit
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import config
import version
from datetime import datetime
from modules import SLiM_runner,ks_calculator, FASTA_extracta,ks_histogramer

def run_sim():

    conf = setup(sys.argv)
    if not conf:
        return

    #num_codons_in_a_gene=1000
    #len_codon=3
    #gene_length=num_codons_in_a_gene*len_codon #nucleotides
    #stop_codons=["TAA","TGA","TAG"]
    #max_num_paralogs_to_process=20 #per genome
    focal_genomes=["n11","n245"] #pick two randomly from each parent
    slim_out_folder=os.path.join(conf.output_folder,"SLiM_output")
    demographics_out_folder=os.path.join(conf.output_folder,"demographiKS_output")
    sim_name = "allotetraploid_bottleneck"
    trees_file = os.path.join(slim_out_folder,"allotetraploid_trees.txt")
    my_SLiM_script= os.path.join("SLiM_scripts", "allotetraploid_bottleneck_trees.slim")
    out_fasta=os.path.join(demographics_out_folder,sim_name + ".fa")
    out_csv=os.path.join(demographics_out_folder,sim_name + ".csv")
    out_png=os.path.join(demographics_out_folder,sim_name + "_hist.png")

    folders_needed=[conf.output_folder,demographics_out_folder,slim_out_folder]
    for f in folders_needed:
        if not os.path.exists(f):
            os.mkdir(f)

    # Run the SLiM model
    print("Running SLiM:\t" + str(my_SLiM_script))
    SLiM_runner.run_slim(conf,trees_file, my_SLiM_script)

    print("Loading:\t" + str(trees_file))
    ts = tskit.load(trees_file)
    metadata=ts.metadata["SLiM"]
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

    print("Getting paralog sequences from TS data.")
    cleaned_sequences_by_paralog_name_dict = FASTA_extracta.extract_paralog_sequences(demographics_out_folder,
                                                                                      focal_genomes,
                                                                       conf, mts, out_fasta, sim_name)

    print("Runnning CODEML on paralogs.")
    paml_out_files = ks_calculator.run_CODEML_on_paralogs(cleaned_sequences_by_paralog_name_dict, demographics_out_folder)

    print("Extracting Ks values from PAML.")
    results = ks_histogramer.extract_K_values(out_csv, paml_out_files)
    ks_histogramer.plot_Ks_histogram(out_png, sim_name,results ,
                      None,None,None,None,"ML","b", 0.001)
    print(results)

    print("Done.")
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
