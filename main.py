from graph import *
from error_correction import *
import sys


def reads_from_file(file):

    reads = []
    for line in file:
        if not line.startswith('>'):
            reads.append(line.strip())
    file.close()
    return reads


# Set length of k-mer
k = 18
# Name of input reads
name = 'reads5'
# Change recursion limit
sys.setrecursionlimit(2500)

if len(sys.argv) == 3:
    input = sys.argv[1]
    output = sys.argv[2]
elif len(sys.argv) == 2:
    input = sys.argv[1]
else:
    input = './reads/%s.fasta' % name

if 'output' not in globals():
    output = './%s_contigs.fasta' % input

reads = reads_from_file(open(input, 'r'))
khist = kmerHist(reads, k)
new_reads = remove_errors(reads, k, 1)
wrong_kmers = remove_rare(new_reads, k)

graph = DeBruijnGraph(new_reads, k, wrong_kmers, mean(khist.values()) - stdev(khist.values()), name, output)

# if you want to save graph into dot file
# graph.to_dot()
