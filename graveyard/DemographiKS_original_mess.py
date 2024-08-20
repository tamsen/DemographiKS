import os
import subprocess, pyslim
import msprime
import numpy as np
import tskit
from cairosvg import svg2png

from ks_calculator import sequences_to_codeml_in, run_codeml
from ks_histogramer import get_Ks_from_file, extract_K_values


def run_sim_old():

    #conf = setup(sys.argv)
    #if not conf:
    #    return

    print("hello!!")
    # Run the SLiM model

    slim_folder="/home/tamsen/Data/SLiM_scripts/"
    my_SLiM_script= os.path.join("SLiM_scripts", "TreeSeqRecording.slim")
    # trees_file=os.path.join(slim_folder,"Tex1_TS.trees")"SLiM_output/diploid_trees.txt"
    trees_file = os.path.join("SLiM_output/diploid_trees.txt")
    # sequence_file=os.path.join(slim_folder,"ex1_TS_overlaid.trees")
    sequence_file = os.path.join("SLiM_output/TS_overlaid_trees.txt")

    out_fasta=os.path.join("SLiM_output/ex1_out.fa")
    out_png=os.path.join("SLiM_output/ex1_out.png")
    out_csv=os.path.join(slim_folder,"ex1_out.csv")
    subprocess.check_output(["slim", "-m", "-s", "0", my_SLiM_script])
    # Overlay neutral mutations
    ts = tskit.load(trees_file)

    #https://tskit.dev/msprime/docs/stable/mutations.html

    #overlays neutral mutations
    mts = msprime.sim_mutations(ts, rate=1e-5, random_seed=42, keep=True)
    print("Mutations added")

    #writes them out here:
    mts.dump(sequence_file)
    #print(str(num_variants) + " mutations written over the trees to " + sequence_file)
    print("Mutations written over the trees to " + sequence_file)

    v_list = [v for v in mts.variants()]
    print(str(len(v_list)) + " mutations written (dumped)")

    #giant string...
    result =mts.as_fasta(reference_sequence=tskit.random_nucleotides(mts.sequence_length))
    with open(out_fasta, "w") as f:
        f.write(result)

    print("Sequences written (to FASTA " + out_fasta + " )")
    print("done")
    return

    #visualization that isnt very good:
    #result = mts.draw_svg()
    #svg2png(bytestring=result,write_to=out_png)

    i=0
    ancestral_seq = ["AAA"] * 100
    mutated_seq = ancestral_seq.copy()
    samp_ids = mts.samples()
    print("  ID of diploid individual: \t\t\t", " ".join([f"{ts.node(s).individual:3}" for s in samp_ids]))
    print("       ID of (sample) node: \t\t\t", " ".join([f"{s:3}" for s in samp_ids]))
    for var in mts.variants():
        #print(str(i))
        site=var.site
        alleles = np.array(var.alleles)
        print(f"Site {site.position} (ancestral state '{site.ancestral_state}')\t", alleles[var.genotypes])
        original_codon=ancestral_seq[int(site.position)][0:2]+site.ancestral_state
        print("original_codon="+original_codon)
        new_codon=ancestral_seq[int(site.position)][0:2]+ alleles[var.genotypes[0]]
        print("new_codon="+new_codon)
        ancestral_seq[int(site.position)] = original_codon
        mutated_seq[int(site.position)] = new_codon
        i=i+1

    print(ancestral_seq)
    print(mutated_seq)
    sequences={"ancestral":"".join(ancestral_seq),
                "derived":"".join(mutated_seq)}

    codeml_input_fa_file=sequences_to_codeml_in(sequences, out_fasta)

    # need "conda install -c bioconda paml"
    result=run_codeml(codeml_input_fa_file, slim_folder)
    paml_out_file=result.ML_dS_file
    results = extract_K_values(out_csv, [paml_out_file])
    print(results)

if __name__ == '__main__':
    run_sim()
