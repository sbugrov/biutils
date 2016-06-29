from collections import Counter

def pattern_count(genome, pattern):
    ''' Finds the number of patterns in a given genome or region of the genome.

    # Arguments

        genome: genome, a string
        pattern: some hidden message, a string

    # Returns

        A number of patterns in a given genome or region of the genome.

    # Note

        Patterns may overlap.
        For various biological processes, certain nucleotide strings often
        appear surprisingly often in small regions of the genome.
    '''

    return sum(pattern == genome[i:i+len(pattern)] for i in xrange(len(genome)-(len(pattern)-1)))

def most_frequent_kmer(genome, k):
    ''' Finds the most frequent pattern(s) of predetermined length k in 
        a given genome or region of the genome.

    # Arguments

        genome: genome, a string
        k: length of some hidden message (k-mer), an integer number

    # Returns

        A list of the most frequent pattern(s) of predetermined length k.

    # Note

        Patterns may overlap.
        For various biological processes, certain nucleotide strings often
        appear surprisingly often in small regions of the genome.
    '''
    
    kmers = dict(Counter([genome[i:i+k] for i in range(len(genome)-k-1)]))
    max_count = max(kmers.values())

    return [kmer for kmer, count in kmers.iteritems() if count == max_count]
