import os
import subprocess, pyslim
import msprime
import numpy as np
import tskit
from cairosvg import svg2png
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from ks_calculator import sequences_to_codeml_in, run_codeml
from ks_histogramer import get_Ks_from_file, extract_K_values


def run_sim():

    #conf = setup(sys.argv)
    #if not conf:
    #    return

    num_codons=5
    gene_length=num_codons*3 #nucleotides
    slim_out_folder="SLiM_output"
    demographics_out_folder="demographics_output"
    sim_name="diploid_snm"
    trees_file = os.path.join(slim_out_folder,"diploid_trees.txt")
    out_fasta=os.path.join(demographics_out_folder,sim_name + ".fa")

    # Run the SLiM model
    #subprocess.check_output(["slim", "-m", "-s", "0", my_SLiM_script])


    ts = tskit.load(trees_file)

    #overlays neutral mutations
    mts = msprime.sim_mutations(ts, rate=1e-5, random_seed=42, keep=True)
    v_list = [v for v in mts.variants()]
    print(str(len(v_list)) + " mutations added.")

    #giant string...
    result =mts.as_fasta(reference_sequence=tskit.random_nucleotides(mts.sequence_length))
    with open(out_fasta, "w") as f:
        f.write(result)

    print("Sequences written to FASTA file: " + out_fasta + ".")

    sequences_we_care_about=["n11","n45"] #pick two randomly

    max_num_paralogs_to_process=5
    sequences_by_paralog_name_dict={}
    with open(out_fasta) as f:

        for seq_record in SeqIO.parse(f, 'fasta'):
            seq_record.id = seq_record.description = seq_record.id.replace('.seq','')
            if seq_record.id in sequences_we_care_about:
                start_index_in_sequence = 0
                num_paralogs_porcessed=0
                genome_name = sim_name + "_" + seq_record.id
                seq=seq_record.seq
                full_seq_length=len(seq)

                while True:

                    paralog_name = genome_name + "_paralog_" + str(start_index_in_sequence)
                    end_index_in_sequence=start_index_in_sequence+gene_length

                    if end_index_in_sequence >= full_seq_length:
                        break
                    if num_paralogs_porcessed >= max_num_paralogs_to_process:
                        break

                    if start_index_in_sequence not in sequences_by_paralog_name_dict:
                        sequences_by_paralog_name_dict[start_index_in_sequence] = {}

                    subsequence=seq[start_index_in_sequence:end_index_in_sequence]
                    sequences_by_paralog_name_dict[start_index_in_sequence][paralog_name]=subsequence

                    print("Subsequence : " + subsequence)
                    start_index_in_sequence = start_index_in_sequence+gene_length
                    out_per_genome_fasta = os.path.join(demographics_out_folder, paralog_name + ".fa")
                    print("Writing data for : " + out_per_genome_fasta + ".")
                    record = SeqRecord(subsequence,
                                       id=seq_record.id, name=paralog_name,
                                       description="simulated paralogous gene")
                    SeqIO.write(record, out_per_genome_fasta , "fasta")
                    num_paralogs_porcessed = num_paralogs_porcessed + 1

    print("Sorting paralogs.")
    for paralog in sequences_by_paralog_name_dict:
        paralog_fa_file = os.path.join(demographics_out_folder, "paralog_" + str(paralog) +".fa")

        codeml_input_fa_file=sequences_to_codeml_in(sequences_by_paralog_name_dict[paralog], paralog_fa_file)
        print("codeml_input_fa_file written to " + codeml_input_fa_file)
        # need "conda install -c bioconda paml"
        codeml_out_folder=""
        result=run_codeml(codeml_input_fa_file, demographics_out_folder)
        #paml_out_file=result.ML_dS_file
    #results = extract_K_values(out_csv, [paml_out_file])
    #print(results)

    print("Done.")
    return

if __name__ == '__main__':
    run_sim()
