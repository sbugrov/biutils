def kmers(k, genome):
    ''' Solve the String Composition Problem.
    # Arguments

        genome: genome, a string
        k: length of needed k-mer

    # Returns
        List of the k-mers ordered in lexicographic order.
    '''
    
    result = [genome[i:i+k] for i in xrange(len(genome) - k + 1)]
    result.sort()
  
    return result
    
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
    
def positions_of_min_skew(genome):
    ''' Finds the position(s) where skew function minimum(s) is attained.
    
    # Arguments

        genome: a string Genome
    
     # Returns
        List of the position(s) of the minimum(s).
    '''
    
    skew_list = skew(genome)
    minimum = min(skew_list)
    return [i for i in range(len(skew_list)) if (skew_list[i]==minimum)]
    
def debruijn(edges):
    '''
    Construct the de Bruijn graph from a set of k-mers.
    Input: A collection of k-mers Patterns.
    Output: De Bruijn graph. Format of the graph:
  
    de_bruijn_graph: [nodes = [], edges = []}
    node: [node, edges = [], degree = |edges|, degreeUnvisited = |edges.visited == false|]
    edge: [edge, nodes = [], visited = False]
    '''
  
    ins = [edge[:-1] for edge in edges]
    unque_ins = set(ins)
    nodes = [[unque_in, [edge[1:] for edge in edges if edge[:-1] == unque_in], len(['dummy' for edge in edges if (edge[:-1] == unque_in or edge[1:] == unque_in)])] for unque_in in unque_ins]
    nodes = [[node[0], node[1], node[2], node[2]] for node in nodes]
    edges = [[edge, [edge[:-1], edge[1:]], False] for edge in edges]
  
    de_bruijn = [nodes, edges]
    
    return de_bruijn
