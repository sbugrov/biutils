
def pattern_count(genome, pattern):
    ''' Finds the number of Patterns in a given genome or region of the Genome.

    # Arguments

        genome: genome, a string
        pattern: some hidden message, a string

    # Returns

        A number of Patterns in a given genome or region of the Genome.

    # Note

        Patterns may overlap.
        For various biological processes, certain nucleotide strings often
        appear surprisingly often in small regions of the genome.
    '''
    genome = genome.upper()
    pattern = pattern.upper()

    return sum(pattern == genome[i:i+len(pattern)] for i in xrange(len(genome)-(len(pattern)-1)))

def most_frequent_kmer(genome, k):
    ''' Finds the most frequent pattern(s) of predetermined length K in
        a given Genome or region of the genome.

    # Arguments

        genome: genome, a string
        k: length of some hidden message (k-mer), an integer number

    # Returns

        A list of the most frequent pattern(s) of predetermined length K.

    # Note

        Patterns may overlap.
        For various biological processes, certain nucleotide strings often
        appear surprisingly often in small regions of the genome.
    '''
    genome = genome.upper()
    unique_kmers = set([genome[i:i + k] for i in xrange(len(genome) - (k - 1))])
    kmer_counts = [(b, pattern_count(genome, b)) for b in unique_kmers]
    max_count = max([count for (kmer, count) in kmer_counts])

    return [kmer for (kmer, count) in kmer_counts if count == max_count]

def reverse_complement(strand):
    ''' Returns reverse complement of a given genome Strand

    # Arguments

        strand: a template strand, a string
    '''

    mapping = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    strand = strand.upper()

    return ''.join(map(lambda c: mapping[c], strand[::-1]))

def starting_positions(pattern, genome):
    ''' Finds all starting positions where Pattern appears as a substring of Genome.

    #Returns

        A list of integers specifying all starting positions where Pattern appears as a substring of Genome.

    # Note

        Patterns may overlap.
    '''
    genome = genome.upper()
    pattern = pattern.upper()

    return [i for i in xrange(len(genome)-(len(pattern)-1)) if (pattern == genome[i:i+len(pattern)])]

def clump_finding(genome, k, L, t):
    ''' Finds all patterns of length K forming (L, T)-clumps in a Genome.

     # Arguments

        genome: a string Genome
        k: length of a pattern
        L: length of a window
        t: minimum number of times a pattern appears in a given Genome

     # Returns
        All distinct K-mers forming (L, T)-clumps in Genome.
    '''
    genome = genome.upper()
    windows = [genome[i:i + L] for i in xrange(len(genome) - (L - 1))]
    kmers_w = [[window[i:i + k] for i in xrange(len(window) - (k - 1))] for window in windows]
    kmer_counts_w = [[(b, pattern_count(windows[i], b)) for b in set(kmers_w[i])] for i in xrange(len(kmers_w))]
    kmer_counts = [[(v, k) for v, k in kmers if k >= t] for kmers in kmer_counts_w]

    return list(set([item[0] for sublist in kmer_counts for item in sublist]))

def skew(genome):
    ''' Skew function
    
    # Arguments

        genome: a string Genome
    '''
    def cumulative_sum(values):
        start=0
        for v in values:
            start += v

            yield start
      
    mapping = {'A': 0, 'T': 0, 'G': 1, 'C': -1}
    genome = genome.upper()
    mapped_genome = map(lambda c: mapping[c], genome)
    
    return [0] + list(cumulative_sum(mapped_genome))
    
