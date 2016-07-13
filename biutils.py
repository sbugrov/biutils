from itertools import product
import numpy as np

def rna_to_aa(rna):
  ''' Translate an RNA string into an amino acid string.
      Input: An RNA string Pattern.
      Output: The translation of Pattern into an amino acid string Peptide.
  '''
  rna_codon_table = { "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L", "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
                      "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*", "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
                      "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L", "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
                      "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
                      "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M", "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
                      "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
                      "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V", "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
                      "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  
  rna = rna.upper()
  rna = [rna[i:3+i] for i in range(0, len(rna), 3)]
    
  return ''.join([rna_codon_table[codon] for codon in rna])

def linear_spectrum(peptide):
  
  '''Generate the theoretical spectrum of a linear peptide.
     Input: An amino acid string Peptide.
     Output: Cyclospectrum(Peptide).
  '''
  
  amino_acid_mass = { 'G' : 57,  'A' : 71,  'S' : 87,  'P' : 97, 'V' : 99,  'T' : 101, 'C' : 103, 'I' : 113, 'L' : 113, 'N' : 114,
                      'D' : 115, 'K' : 128, 'Q' : 128, 'E' : 129,'M' : 131, 'H' : 137, 'F' : 147, 'R' : 156, 'Y' : 163, 'W' : 186}
  
  prefix_mass = [0]
  
  for i in xrange(1, len(peptide) + 1):
    prefix_mass.append(prefix_mass[i-1] + amino_acid_mass[peptide[i-1]])
  
  for i in xrange(1, len(peptide)+1):
    for j in xrange(i + 1, len(peptide)+1):
      prefix_mass.append(prefix_mass[j] - prefix_mass[i])
      
  prefix_mass.sort()
  
  return prefix_mass
  
def cyclic_spectrum(peptide):
  
  '''Generate the theoretical spectrum of a cyclic peptide.
     Input: An amino acid string Peptide.
     Output: Cyclospectrum(Peptide).
  '''
  
  amino_acid_mass = { 'G' : 57,  'A' : 71,  'S' : 87,  'P' : 97, 'V' : 99,  'T' : 101, 'C' : 103, 'I' : 113, 'L' : 113, 'N' : 114,
                      'D' : 115, 'K' : 128, 'Q' : 128, 'E' : 129,'M' : 131, 'H' : 137, 'F' : 147, 'R' : 156, 'Y' : 163, 'W' : 186}
  
  prefix_mass = [0]
  
  peptide_length = len(peptide)
  
  for i in xrange(1, peptide_length + 1):
    prefix_mass.append(prefix_mass[i-1] + amino_acid_mass[peptide[i-1]])
  
  total_peptide_mass = prefix_mass[peptide_length]
  
  for i in xrange(1, peptide_length+1):
    for j in xrange(i + 1, peptide_length+1):
      prefix_mass.append(prefix_mass[j] - prefix_mass[i])
      if i > 0 and j < peptide_length:
        prefix_mass.append(total_peptide_mass - (prefix_mass[j] - prefix_mass[i]))
        
  prefix_mass.sort()
  
  return prefix_mass

def score_spectrum(peptide, spectrum):
  ''' Returns the number of matches between the theoretical
      spectrum of a given cyclic Peptide and the given experimental Spectrum.
  '''
  
  score = 0
  peptide = cyclic_spectrum(peptide)
  for i in xrange(len(peptide)):
    for j in xrange(len(spectrum)):
      if peptide[i] == spectrum[j]:
        score += 1
        del spectrum[j]
        break
  return score

def all_kmers(k):
  '''
  return a list of all possible k-mers
  '''
  return [''.join(a) for a in product(['A','C','T','G'], repeat=k)]
  
def hamming_distance(a, b):
  '''Compute the Hamming distance between two strings.
  # Arguments
      a: a string
      b: a string
  # Returns
      The Hamming distance between strings A and B.
  '''
  length = len(a)
  distance = length
  for i in range(length):
      distance -= a[i] == b[i]
  return distance

def pattern_matching_with_difference (pattern, sequence, difference):
  ''' Returns all starting positions where Pattern appears 
  as an approximate substring of Sequence with a at most 'd' character difference
  '''
  p_len = len(pattern)
  s_len = len(sequence)
  match = []
  for i in range(s_len - p_len + 1):
    if hamming_distance(sequence[i:i+p_len],pattern) <= d:
      match.append(i)
            
  return match

def kmers(k, genome):
  ''' Solve the String Composition Problem.
  # Arguments
      genome: genome, a string
      k: length of needed k-mer
  # Returns
      List of the k-mers.
  '''

  result = [genome[i:i+k] for i in xrange(len(genome) - k + 1)]

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

def build_de_bruijn_graph(edges):
  '''
  Construct the de Bruijn graph from a set of k-mers.
  Input: A collection of k-mers Patterns.
  Output: De Bruijn graph. Format of the graph:
  
  de_bruijn_graph: {node: [outs]}
  '''
  
  ins = [edge[:-1] for edge in edges]
  unque_ins = set(ins)
  de_bruijn_graph = {unque_in: [edge[1:] for edge in edges if edge[:-1] == unque_in] for unque_in in unque_ins}
 
  return de_bruijn_graph
  
def find_eulerian_path(graph):
  '''
  '''
  
  eulerian_path = graph.values()
  eulerian_path = [item for sublist in eulerian_path for item in sublist]
  unque_nodes = set(eulerian_path)
  visits_left = {node: eulerian_path.count(node) for node in unque_nodes}
  cycle = []
  graph['AGA'] = graph['AGA'] + ['AAG'] # figure out the first and the last node instead
  print graph
  
  def find_new_starting_point():
    if len(cycle) == 0:
      cycle.append('AAG') # append the first node
      return 'AAG' # return the first node
    else:
      for node in cycle:
        if len(graph[node]) > 0:
          return node
      return 0
    
  def append_node(cycle, position, new_node): # needs improvement
    cycle.append(new_node)
    return cycle
    
    
  while(True):
    new_startring_point = find_new_starting_point()
    
    if (new_startring_point == 0):
      break
          
    while(True):
      next_step = graph[new_startring_point][0]
      
      if len(graph[new_startring_point]) > 1:
        graph[new_startring_point] = graph[new_startring_point][1:]
      else:
        graph[new_startring_point] = []
      
      if len(graph[next_step]) == 0:
        cycle = append_node(cycle, 'here must be something, a position', next_step)
        break
        
      cycle = append_node(cycle, 'here must be something a position', next_step)
      new_startring_point = next_step
        
  eulerian_path = cycle
  
  return eulerian_path

def best_match(motif, sequence):
  ''' Returns a substring of the Sequence that has the lowest
      number of differences with the Motif
  '''

  def hamming_distance(a, b):
    length = len(a)
    distance = length
    for i in range(length):
      distance -= a[i] == b[i]
    return distance

  p_len = len(motif)
  s_len = len(sequence)
  best_match = sequence[0:p_len]
  best_distance = p_len

  for i in range(s_len - p_len + 1):
    h_d = hamming_distance(sequence[i:i+p_len], motif)
    if h_d < best_distance:
      match = sequence[i:i+p_len]
      best_distance = h_d

  return match

def profile(motifs):
  '''
  '''
  
  mapping = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1]}
  n_motifs = len(motifs)
  motif_len = len(motifs[0])
  profile = np.ones([4, motif_len])
  
  for i in xrange(n_motifs):
    for j in xrange(motif_len):
      profile[:, j] += mapping[motifs[i][j]]
      
  return profile

print profile(motifs)
