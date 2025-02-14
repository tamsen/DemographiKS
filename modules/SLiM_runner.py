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
def run_slim(config,trees_file_name, trees_file_name_at_div, my_SLiM_script):

    print_slim_version(config)
    log.write_to_log("copy slim script:\t" + my_SLiM_script)
    shutil.copy(my_SLiM_script,config.output_folder)
    full_path_to_slim_script_destination = os.path.join(config.output_folder, os.path.basename(my_SLiM_script))
    log.write_to_log("full_path_to_slim_script_destination:\t" + full_path_to_slim_script_destination)
    delta_t=  ( float(config.DIV_time_Ge) - float(config.WGD_time_Ge) )

    #burnin_time = 2 * 10 * config.ancestral_Ne <- reccomended by SLiM manual, but not really always enough.
    if not config.DIV_time_Ge:
        div_time_string="-1"
    else:
        div_time_string=str(config.DIV_time_Ge)

    burnin_time = config.burnin_time
    cmd = ["slim",
           "-d", "trees_file_name='"+str(trees_file_name)+"'",
           "-d", "trees_file_name_at_div='" + str(trees_file_name_at_div) + "'",
           "-d", "L=" + str(config.total_num_bases),
           "-d", "Na=" + str(config.ancestral_Ne),
           "-d", "Nb=" + str(config.bottleneck_Ne),
           "-d", "delta_t=" + str(delta_t),
           "-d", "Tdiv_gen=" + div_time_string,
           "-d", "Twgd_gen=" + str(config.WGD_time_Ge),
           "-d", "BurninTime=" + str(burnin_time),
           "-d", "recombination_rate=" + str(config.recombination_rate),
           "-d", "dij=" + str(config.homoeologous_exchange_rate),
           "-d", "rep=" + str(config.SLiM_rep)]

    if str(config.mig_rate) != str(False):
        extra_commands_for_migration=[
            "-d", "MigStart="+str(config.mig_start),
            "-d", "MigStop="+str(config.mig_stop),
            "-d", "MigRate="+str(config.mig_rate)]
        cmd = cmd  + extra_commands_for_migration

    if str(config.assortative_mating_coefficient) != str(False):
        extra_commands_for_inbreeding=[
            "-d", "AMCoef="+str(config.assortative_mating_coefficient)]
        cmd = cmd  + extra_commands_for_inbreeding

    cmd= cmd + ["-m", "-s", "0", full_path_to_slim_script_destination]
    log.write_to_log("\t cmd: " + " ".join(cmd))
    log.write_to_log("\t cwd: " + config.output_folder)
    out_string,error_string = process_wrapper.run_and_wait_on_process(cmd, config.output_folder)
    return