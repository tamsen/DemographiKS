

#https://stackoverflow.com/questions/72531894/calculate-sequence-identity-of-two-sequences-of-equal-length
def compute_sequence_identity(sequence_a, sequence_b):
    if len(sequence_a) == len(sequence_b):
        return sum([sequence_a[i] == sequence_b[i] for i in range(len(sequence_a))]) / len(sequence_a)
    else:
        print("Sequences are not of equal length.")
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

def inbreeding_coefficient(Het_exp,Het_obs):
    numerator =4.0*Ne*mutation_rate
    denominator = numerator + 1.0
    return  numerator/ denominator
