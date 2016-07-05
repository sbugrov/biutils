# biutils

Pretty basic functions that are useful in bioinformatics.

all_kmers(k): return a list of all possible k-mers

hamming_distance(a, b): Compute the Hamming distance between two strings.
 
kmers(k, genome): Returns a list of the k-mers ordered in lexicographic order.

pattern_count(genome, pattern): Finds the number of Patterns in a given genome or region of the Genome.
  
most_frequent_kmer(genome, k): Finds the most frequent pattern(s) of predetermined length K in

reverse_complement(strand): Returns reverse complement of a given sequence

starting_positions(pattern, genome): Finds all starting positions where Pattern appears as a substring of a sequence.
  
clump_finding(genome, k, L, t): Finds all patterns of length K forming (L, T)-clumps in a sequence.

skew(genome): Skew function

positions_of_min_skew(genome): Finds the position(s) where skew function minimum(s) is attained.
  
debruijn(edges): Construct the de Bruijn graph from a set of k-mers.
