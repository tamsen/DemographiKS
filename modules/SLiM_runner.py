import os
import shutil
import log
import process_wrapper


def print_slim_version(config):
    cmd = ["slim","-v", ]

    log.write_to_log("\t cmd: " + " ".join(cmd))
    log.write_to_log("\t cwd: " + config.output_folder)
    out_string, error_string = process_wrapper.run_and_wait_on_process(cmd, config.output_folder)
    log.write_to_log("\t slim -v out_string: " + out_string)
    log.write_to_log("\t slim -v error_string: " + error_string)
    return
def run_slim(config,trees_file_name, my_SLiM_script):

    #TODO, expose parameter: config.recombination_rate

    print_slim_version(config)

    log.write_to_log("copy slim script:\t" + my_SLiM_script)
    shutil.copy(my_SLiM_script,config.output_folder)
    full_path_to_slim_script_destination = os.path.join(config.output_folder, os.path.basename(my_SLiM_script))
    log.write_to_log("full_path_to_slim_script_destination:\t" + full_path_to_slim_script_destination)
	#Parameters:
	#	nuBot: Proportion of the ancestral population size remaining after bottleneck.
	#	T1: The amount of time in dadi units (# of 2N generations) that the parents
	#		are isolated before forming the allotetraploid.
	#	T2: The amount of time the allotetraploid lineage has existed before we sample
	#		it.
	#	rep: Simulation replicate number (for running things in a for loop or
	#		 an array job on an HPC).

    #slim -d "nuBot=0.1" -d "T1=0.5" -d "T2=0.25" -d "rep=1" allotetraploid_bottleneck.slim

    delta_t=  ( float(config.DIV_time_Ge) - float(config.WGD_time_Ge) )
    burnin_time = 2 * 10 * config.ancestral_Ne
    cmd = ["slim",
           "-d", "trees_file_name='"+str(trees_file_name)+"'",
           "-d", "L=" + str(config.total_num_bases),
           "-d", "Na=" + str(config.ancestral_Ne),
           "-d", "Nb=" + str(config.bottleneck_Ne),
           "-d", "delta_t=" + str(delta_t),
           "-d", "Tdiv_gen=" + str(config.DIV_time_Ge),
           "-d", "BurninTime=" + str(burnin_time),
           "-d", "recombination_rate=" + str(config.recombination_rate),
           "-d", "rep=" + str(config.SLiM_rep),
           "-m", "-s", "0", full_path_to_slim_script_destination]


    log.write_to_log("\t cmd: " + " ".join(cmd))
    log.write_to_log("\t cwd: " + config.output_folder)
    out_string,error_string = process_wrapper.run_and_wait_on_process(cmd, config.output_folder)
    return