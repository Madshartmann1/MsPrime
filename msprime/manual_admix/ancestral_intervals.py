#!/usr/bin/env python

import numpy as np
import random




class Individual:
    def __init__(self, mom_id=None, dad_id=None, chrom_pair=[[], []]):
        self.mom_id = mom_id
        self.dad_id = dad_id
        self.chrom_pair = chrom_pair

class Sequence:
    def __init__(self, start_pos=None, end_pos=None, id=None):
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.id = id




# Function to generate expanded ancestral haplotype IDs with _M and _P suffixes
def generate_ancestral_haplotype_ids_with_suffixes():
    #these numbers match  the rations from the DEMES file deciding the admixed population ancestery
    max_numbers = {"KAR": 200, "YRI": 150, "CEU": 50}
    expanded_ancestral_haplotype_id = []
    for key, max_num in max_numbers.items():
        for i in range(0, max_num):
            expanded_ancestral_haplotype_id.append(f"{key}#{i}_M")
            expanded_ancestral_haplotype_id.append(f"{key}#{i}_P")
    return expanded_ancestral_haplotype_id


def get_segments(start_pos, end_pos, chrom):
    i = 0
    segments = []
    while True:
        if start_pos > chrom[i].end_pos:
            i += 1
        elif end_pos > chrom[i].end_pos:
            segments.append(Sequence(start_pos, chrom[i].end_pos, chrom[i].id))
            start_pos = chrom[i].end_pos+1
            i += 1
        else:
            segments.append(Sequence(start_pos, end_pos, chrom[i].id))
            break
    return segments


def get_chrom(chrom_pair, recomb_rate):
    num_recomb = np.random.binomial(n=seq_length, p=recomb_rate)
    
    break_pos = [0] * (num_recomb+1)
    for i in range(num_recomb):
        break_pos[i] = random.randint(a=0, b=seq_length-1)
    break_pos[-1] = seq_length
    break_pos.sort()

    chrom_index = random.getrandbits(1)

    seq = get_segments(0, break_pos[0]-1, chrom_pair[chrom_index])
    chrom_index = not chrom_index

    for i in range(num_recomb):
        segment = get_segments(break_pos[i], break_pos[i+1]-1, chrom_pair[chrom_index])
        if segment[0].id == seq[-1].id:
            segment[0].start_pos = seq.pop().start_pos
        seq += segment
        chrom_index = not chrom_index
        
    return seq


def get_population_size(init_population_size, num_generations, growth_rate):
    population_sizes = [0] * num_generations
    population_sizes[0] = init_population_size
    for i in range(1, num_generations):
        population_sizes[i] = population_sizes[i-1] * (growth_rate+1)
    return list(map(int, population_sizes))


def run_wright_fisher(init_population_size, num_generations, growth_rate, seq_length, recomb_rate, ancestor_pool):
    population_sizes = get_population_size(init_population_size, num_generations, growth_rate)
    population = [[0 for _ in range(population_sizes[i])] for i in range(num_generations)]

    #setting up the initial generation
    for i in range(init_population_size):

        # Pick a random ancestor ID from the list
        ancestor_id = random.choice(ancestor_pool)
        # Create two Sequence objects with this ID
        chrom = [Sequence(0, seq_length-1, ancestor_id)]
        # Assign the individual with chromosomes to the population
        population[0][i] = Individual(chrom_pair=[chrom, chrom])


    #simulation of following generations
    for i in range(1, num_generations):
        for j in range(population_sizes[i]):
            mom_id = random.randint(0, population_sizes[i-1] - 1)
            dad_id = random.randint(0, population_sizes[i-1] - 1)

            mom_seq = get_chrom(population[i-1][mom_id].chrom_pair, recomb_rate)
            dad_seq = get_chrom(population[i-1][dad_id].chrom_pair, recomb_rate)

            population[i][j] = Individual(mom_id, dad_id, [mom_seq, dad_seq])
    return population



#start of simulation
seed = 5555
random.seed(seed)
np.random.seed(seed)

#paramters
init_population_size = 300
num_generations = 12
seq_length = 248387328
recomb_rate = 1e-8
growth_rate = 0

# Generation of list of haplotypes IDs
ancestral_haplotype_id = generate_ancestral_haplotype_ids_with_suffixes()

# Simulate
population = run_wright_fisher(init_population_size,
                               num_generations,
                               growth_rate,
                               seq_length,
                               recomb_rate,
                               ancestral_haplotype_id)


# Assuming we're focusing on the last generation
last_generation_index = num_generations - 1  # Index of the last generation

# Print table header
#print("Individual\tM\t\t\tP\n")

for j in range(len(population[last_generation_index])):
    # Lists to hold segment information strings for maternal (M) and paternal (P) chromosomes
    segments_info_m = []  # Maternal chromosome segments
    segments_info_p = []  # Paternal chromosome segments
    
    # Assuming the first chromosome in chrom_pair is maternal (M) and the second is paternal (P)
    maternal_chrom = population[last_generation_index][j].chrom_pair[0]
    paternal_chrom = population[last_generation_index][j].chrom_pair[1]

    # Process maternal chromosome segments
    for segment in maternal_chrom:
        segment_info_m = f"{segment.start_pos}, {segment.end_pos}, {segment.id}"
        segments_info_m.append(segment_info_m)

    # Process paternal chromosome segments
    for segment in paternal_chrom:
        segment_info_p = f"{segment.start_pos}, {segment.end_pos}, {segment.id}"
        segments_info_p.append(segment_info_p)

    # Print individual's row with segments info
    # Assuming you want to print lists of segments info directly under "M" and "P"
    print(f"Individual {j}")
    print("\tM\t" + "\n\t\t".join(segments_info_m))
    print("\tP\t" + "\n\t\t".join(segments_info_p) + "\n")






