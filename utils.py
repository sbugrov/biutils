from collections import Counter

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

    kmers = dict(Counter([genome[i:i+k] for i in range(len(genome)-k-1)]))
    max_count = max(kmers.values())

    return [kmer for kmer, count in kmers.iteritems() if count == max_count]

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
