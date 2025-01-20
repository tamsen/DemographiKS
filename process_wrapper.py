import os
import subprocess
import time
import log

def run_and_wait_with_retry(cmd, folder, excuse, num_retries_allowed, sleepy_time):

    num_tries=0
    while True:

        out_string,error_string = run_and_wait_on_process(cmd, folder)
        num_tries = num_tries+1
        if num_tries > num_retries_allowed:
            break

        if excuse in error_string:
            print("Got "+ excuse +". Retrying " + str(num_tries) + " time.")
            print("Wait for " + str(sleepy_time) + " secs.")
            time.sleep(sleepy_time)
        else:
            break


    return out_string,error_string

def run_and_wait_on_process(cmd, cwd_folder):

    program=cmd[0]
    log.write_to_log(" ".join(cmd) )
    process_completed_result = subprocess.run(cmd, capture_output=True, cwd=cwd_folder)
    stderr_string=process_completed_result.stderr.decode()
    out_string=process_completed_result.stdout.decode()

    #https://stackoverflow.com/questions/287871/how-do-i-print-colored-text-to-the-terminal
    colored_error_string = '\033[93m' + "WARNING:  " + stderr_string + '\x1b[0m'
    if len(stderr_string)> 0:
        log.write_to_log(colored_error_string )

    with open(os.path.join(cwd_folder, program + "_stderr.txt"), 'w') as f:
        f.writelines(stderr_string)
    with open(os.path.join(cwd_folder, program + "_stdout.txt"), 'w') as f:
        f.writelines(out_string)

    return out_string,stderr_string