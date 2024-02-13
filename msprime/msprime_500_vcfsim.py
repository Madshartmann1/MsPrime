#python script for creating fasta files of individuals via bash execution 
import os
import msprime
import sys
import demes

#Checking the arguments given to the program
if len(sys.argv) != 2:
    print("Usage: python3 THIS_PROGRAM.py <demes_file_path> > OUTPUT.fasta")
    sys.exit(1)

demes_file_path = sys.argv[1]
if not os.path.exists(demes_file_path):
    print(f"File {demes_file_path} not found.")
    sys.exit(1)


#parameters
mutation_rate = 1e-8 #static rate for mutation across the sequence
sq_len = 248387328    #248387328 length of sequence (chrom 1 here)
n_individuals = 500 #number of individuals
rate_recomb = 1e-8 #average recombination rate in humans
static_seed = 5555 #eliminate randomness

#load deme and simulation populations
graph = demes.load(demes_file_path)
demography = msprime.Demography.from_demes(graph)
ts = msprime.sim_ancestry({"CEU": n_individuals,
                           "YRI": n_individuals,
                           "CHB": n_individuals,
                           "KAR": n_individuals},
                          sequence_length=sq_len,
                          demography=demography,
                          random_seed=static_seed,
                          recombination_rate=rate_recomb)
mts = msprime.sim_mutations(ts,
                            rate=mutation_rate,
                            random_seed=static_seed)




#make a vcf
mts.write_vcf(sys.stdout)
