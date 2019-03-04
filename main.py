from graph import *
from error_correction import *
import sys


def reads_from_file(file):
    """
    Saving reads from the fasta file to a list.
    :param file: object of the open fasta file
    :return: list of the reads
    """

    reads = []
    for line in file:
        if not line.startswith('>'):
            reads.append(line.strip())
    file.close()
    return reads


def build_up(reads, k):
    """
    Removing errors from the given reads, establishing list of tentative kmers.
    :param reads: list of input reads
    :param k: size of k-mers
    :return:
    new_reads - list of corrected reads,
    wrong_kmers - list of tentative kmers,
    thresh - threshold of kmers coverage
    """

    khist, thresh = kmerHist(reads, k)
    new_reads = remove_errors(reads, k, khist)
    wrong_kmers = rare_kmers(khist, thresh)

    return new_reads, wrong_kmers, thresh


# Change recursion limit
sys.setrecursionlimit(8000)

if len(sys.argv) == 3:
    input = sys.argv[1]
    output = sys.argv[2]
elif len(sys.argv) == 2:
    input = sys.argv[1]
else:
    read = 5
    input = './reads/reads%d.fasta' % read

name = input.split('/')[-1].split('.fasta')[0]
if 'output' not in globals():
    output = './%s_contigs.fasta' % name

n = 0
bestk = 0
reads = reads_from_file(open(input, 'r'))
ref = 0
for r in reads:
    ref += len(r)
ref /= 5

# looking for the optimal k value
for k in range(15, 24):
    new_reads, wrong_kmers, thresh = build_up(reads, k)
    try:
        graph = DeBruijnGraph(new_reads, k, wrong_kmers, thresh, name, output)
    except RecursionError:
        print('Recursion error for k = %d' % k)
        continue
    mark = sum([len(c.seq) for c in graph.contigs])/(ref*math.log(4+len(graph.contigs), 5))
    print('Found %d contigs for k = %d, mark = %.3f' % (len(graph.contigs), k, mark))
    if mark >= n:
        n = mark
        bestk = k
        bestgraph = graph

print('Best k = %d' % bestk)
bestgraph.contigs_to_file()
print('Contigs saved!')


# if you want to save the best graph into dot file
# prev_graph.to_dot()
