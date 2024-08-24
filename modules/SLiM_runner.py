import os
import shutil

import process_wrapper


def run_slim(config,trees_file_name, my_SLiM_script):

    #TODO - path issiue here
    shutil.copy(my_SLiM_script,config.output_folder,)
    #full_path_to_slim_script=os.path.join(os.getcwd(),my_SLiM_script)
    full_path_to_slim_script = os.path.join(config.output_folder, my_SLiM_script)

	#Parameters:
	#	nuBot: Proportion of the ancestral population size remaining after bottleneck.
	#	T1: The amount of time in dadi units (# of 2N generations) that the parents
	#		are isolated before forming the allotetraploid.
	#	T2: The amount of time the allotetraploid lineage has existed before we sample
	#		it.
	#	rep: Simulation replicate number (for running things in a for loop or
	#		 an array job on an HPC).

    cmd = ["slim",
           "-d", "trees_file_name='"+str(trees_file_name)+"'",
           "-d", "L=" + str(config.total_num_bases),
           "-d", "nuBot=" + str(0.1),
           "-d", "T2=" + str(0.5),
           "-d", "T1=" + str(0.25),
           "-d", "rep=" + str(1),
           "-m", "-s", "0", full_path_to_slim_script]

    print("\t cmd: " + " ".join(cmd))
    print("\t cwd: " + config.output_folder)
    out_string,error_string = process_wrapper.run_and_wait_on_process(cmd, config.output_folder)
    return