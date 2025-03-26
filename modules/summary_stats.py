import itertools

#https://stackoverflow.com/questions/72531894/calculate-sequence-identity-of-two-sequences-of-equal-length
def compute_sequence_identity(sequence_a, sequence_b):
    if len(sequence_a) == len(sequence_b):
        return sum([sequence_a[i] == sequence_b[i] for i in range(len(sequence_a))]) / len(sequence_a)
    else:
        print("Sequences are not of equal length.")
        return None

def compute_num_sequence_diffs(sequence_a, sequence_b):
    if len(sequence_a) == len(sequence_b):
        return sum([sequence_a[i] != sequence_b[i] for i in range(len(sequence_a))])
    else:
        print("Sequences are not of equal length.")
        #print("SequenceA:\t" + sequence_a)
        #print("Len SequenceA:\t" + str(len(sequence_a)))
        #print("SequenceB:\t" + sequence_b)
        #print("Len SequenceB:\t" + str(len(sequence_b)))
        return None

#sequences = ['AUUGCAUG', 'CGUGGCUA']
#sequence_identity = compute_sequence_identity(sequences[0], sequences[1])
#print("Sequence identity: " + str(sequence_identity) + "%")

#Inbreeding Coefficient measures the excess heterozygosity
#Fis = (Het_exp - Het_obs) / Het_exp

# https://pmc.ncbi.nlm.nih.gov/articles/PMC11192967/#:~:text=2.,(2)
#Het_exp = 4 * Ne * m / ( 4 * Ne * m + 1)

def expected_heterozygosity(Ne,mutation_rate):
    numerator =4.0*Ne*mutation_rate
    denominator = numerator + 1.0
    return  numerator/ denominator

def inbreeding_coefficient(het_exp,het_obs):
    numerator =het_exp - het_obs
    return  numerator/ het_exp

# https://www.nature.com/articles/hdy201643
# heterozygosity excess
# Ne = π / (4μ);
# and π is nucleotide diversity
# (average number of nuc differences in a population or sample)
# https://evolutionarygenetics.github.io/Chapter7.html

# seems like a nice package
# https://scikit-allel.readthedocs.io/en/stable/stats/diversity.html
def observed_Ne(pi,mu):
    return  pi/ (4.0 * mu)

def pi_nuc_diversity(seq_list):
    all_pairs_of_seqs = list(itertools.combinations(seq_list, 2))
    print("pairs: " + str(all_pairs_of_seqs))
    diffs = [compute_num_sequence_diffs(p[0],p[1]) for p in all_pairs_of_seqs  ]
    print("diffs: " + str(diffs))
    num_comparisons = (len(all_pairs_of_seqs)*len(seq_list[0]))
    diversity = sum(diffs) / num_comparisons
    return diversity
