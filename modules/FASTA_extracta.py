import os
import tskit
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import log
def extract_paralog_sequences(demographics_out_folder, focal_genomes, config, mts, out_fasta):
    # giant string...

    result = mts.as_fasta(reference_sequence=tskit.random_nucleotides(mts.sequence_length))
    with open(out_fasta, "w") as f:
        f.write(result)
    log.write_to_log("Sequences written to FASTA file: " + out_fasta + ".")
    sequences_by_paralog_name_dict = write_per_genome_per_paralog_fastas(demographics_out_folder, focal_genomes,
                                                                         config.gene_length,
                                                                         config.max_num_paralogs_to_process,
                                                                         out_fasta, config.sim_name)
    log.write_to_log("Removing STOP codons. PAML needs sequences that code for AA only")
    problem_codon_indexes_by_paralog_name_dict = get_index_of_any_STOP_codons(config.num_codons_in_a_gene,
                                                                              sequences_by_paralog_name_dict,
                                                                              config.stop_codons)
    cleaned_sequences_by_paralog_name_dict = set_STOP_codons_to_NNN(config.len_codon,
                                                                    problem_codon_indexes_by_paralog_name_dict,
                                                                    sequences_by_paralog_name_dict)
    return cleaned_sequences_by_paralog_name_dict


def get_sequences_by_paralog_name(gene_length,genome_size,
                                        max_num_paralogs_to_process):
    paralog_names = []
    expected_num_paralogs=int(genome_size/gene_length)
    max_num_paralogs=min(expected_num_paralogs,max_num_paralogs_to_process)
    for i in range(0,max_num_paralogs):
        start_index=i*gene_length
        paralog_names.append(start_index)
    return paralog_names


def write_per_genome_per_paralog_fastas(demographics_out_folder, focal_genomes, gene_length,
                                        max_num_paralogs_to_process, out_fasta, sim_name):
    sequences_by_paralog_name_dict = {}
    with open(out_fasta) as f:

        for seq_record in SeqIO.parse(f, 'fasta'):
            seq_record.id = seq_record.description = seq_record.id.replace('.seq', '')
            if seq_record.id in focal_genomes:
                start_index_in_sequence = 0
                num_paralogs_processed = 0
                genome_name = sim_name + "_" + seq_record.id
                seq = seq_record.seq
                full_seq_length = len(seq)

                while True:

                    paralog_name = genome_name + "_paralog_" + str(start_index_in_sequence)
                    end_index_in_sequence = start_index_in_sequence + gene_length

                    if end_index_in_sequence >= full_seq_length:
                        break
                    if max_num_paralogs_to_process:
                        if num_paralogs_processed >= max_num_paralogs_to_process:
                            break

                    if start_index_in_sequence not in sequences_by_paralog_name_dict:
                        sequences_by_paralog_name_dict[start_index_in_sequence] = {}

                    subsequence = seq[start_index_in_sequence:end_index_in_sequence]
                    sequences_by_paralog_name_dict[start_index_in_sequence][paralog_name] = subsequence

                    # print("Subsequence : " + subsequence)
                    start_index_in_sequence = start_index_in_sequence + gene_length
                    out_per_genome_per_paralog_fasta = os.path.join(demographics_out_folder, paralog_name + ".fa")
                    # print("Writing data for : " + out_per_genome_fasta + ".")
                    record = SeqRecord(subsequence,
                                       id=seq_record.id, name=paralog_name,
                                       description="simulated paralogous gene")
                    SeqIO.write(record, out_per_genome_per_paralog_fasta, "fasta")
                    num_paralogs_processed = num_paralogs_processed + 1
    return sequences_by_paralog_name_dict


def set_STOP_codons_to_NNN(len_codon, problem_codon_indexes_by_paralog_name_dict, sequences_by_paralog_name_dict):
    cleaned_sequences_by_paralog_name_dict = {}
    for paralog_key, sequences_dict in sequences_by_paralog_name_dict.items():
        cleaned_sequences_by_paralog_name_dict[paralog_key] = {}
    for paralog_key, cleaned_sequences_dict in cleaned_sequences_by_paralog_name_dict.items():
        problem_codon_indexes = problem_codon_indexes_by_paralog_name_dict[paralog_key]
        for seq_key, sequence in sequences_by_paralog_name_dict[paralog_key].items():
            # print("problem codons " + str(problem_codon_indexes))
            # print("original sequence:\t" + str(sequence))
            revised_seq = str(sequence)
            for problem_codon_index in problem_codon_indexes:
                revised_seq = revised_seq[:problem_codon_index * len_codon] + "NNN" + revised_seq[(
                                                                                                              problem_codon_index + 1) * len_codon:]
            # print("revised sequence :\t" + revised_seq)
            cleaned_sequences_by_paralog_name_dict[paralog_key][seq_key] = revised_seq
    return cleaned_sequences_by_paralog_name_dict


def get_index_of_any_STOP_codons(num_codons_in_a_gene, sequences_by_paralog_name_dict, stop_codons):
    problem_codon_indexes_by_paralog_name_dict = {}
    for paralog_key, sequences_dict in sequences_by_paralog_name_dict.items():
        log.write_to_log("checking paralog " + str(paralog_key))
        problem_codon_indexes = []
        for seq_key, sequence in sequences_dict.items():
            for i in range(0, num_codons_in_a_gene):
                codon = str(sequence[3 * i:3 * (i + 1)])
                if codon in stop_codons:
                    problem_codon_indexes.append(i)
        problem_codon_indexes_by_paralog_name_dict[paralog_key] = problem_codon_indexes
    return problem_codon_indexes_by_paralog_name_dict

