import os
import shutil
from pathlib import Path
import process_wrapper



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


def sequences_to_codeml_in(sequences, fa_out_file):

    if len(sequences)== 0:
        return []

    sequence_names = sequences.keys()
    with open(fa_out_file, 'w') as f:

        for seq_name in sequence_names:
            f.writelines(">" + seq_name + "\n")
            f.writelines(sequences[seq_name] + "\n")

    return fa_out_file

def write_codeml_control_file(template_ctl_file, sequence_file):
    lines_to_write = []
    base_name = os.path.basename(sequence_file).replace(".codonalign", "").replace(".fa", "")
    new_ctl_file = sequence_file.replace(".fa", ".ctl")

    with open(template_ctl_file, 'r') as f:

        while (True):

            line = f.readline()
            new_line = line

            if "DUMMY.codonalign.fa" in line:
                #new_line = line.replace("DUMMY.codonalign.fa", sequence_file)
                new_line = line.replace("DUMMY.codonalign.fa",os.path.basename(sequence_file))
                print(new_line)

            if "mlcTree_DUMMY.out" in line:
                new_line = line.replace("DUMMY", base_name)

            lines_to_write.append(new_line)

            if not line:
                break

    with open(new_ctl_file, 'w') as f:

        for line in lines_to_write:
            f.writelines(line)

    return new_ctl_file

def run_codeml(fa_out_file,out_folder):


    template_codeml_ctl_file = get_codeml_ctl_template()

    control_file = write_codeml_control_file(template_codeml_ctl_file, fa_out_file)

    cmd = ["codeml",os.path.basename(control_file)]
    print("\t cmd: " + " ".join(cmd))
    print("\t cwd: " + out_folder)
    print("\t calculating Ks.. ")
    process_wrapper.run_and_wait_on_process(cmd, out_folder)
    print("\t Ks determined...")
    result= codeml_result(out_folder)

    return result


def get_codeml_ctl_template():
    par_dir = Path(__file__).parent.parent
    template_codeml_ctl_file = os.path.join(par_dir, "paml_input_templates",
                                            "codeml_input_example.ctl")
    return template_codeml_ctl_file


class codeml_result():
    ML_dS_file= ""
    ML_dN_file= ""
    #ML_file=""
    def __init__(self, gene_tree_subfolder):
        self.ML_dS_file = os.path.join(gene_tree_subfolder, "2ML.dS")
        self.ML_dN_file = os.path.join(gene_tree_subfolder, "2ML.dN")
        #self.NG_file = os.path.join(gene_tree_subfolder,"2NG.dS")