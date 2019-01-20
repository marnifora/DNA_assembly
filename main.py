from DeBruijn import *
import DeBruijn_copy as copy
from error_correction import *
import sys


def reads_from_file(file, l):

    reads = []
    for i, line in enumerate(file):
        if i >= 2*l:
            break
        if not line.startswith('>'):
            reads.append(line.strip())
    file.close()
    return reads


def build_graph(reads, k, name, wrong_kmers, thresh):

    graph = DeBruijnGraph(reads, k, wrong_kmers, thresh, name)
    # graph.toDot(weights=True)

    # path = list(graph.eulerianPath())
    # assembly = path[0] + ''.join(map(lambda x: x[-1], path[1:]))

    # print(path, assembly)


def remove_errors(reads, k, thresh):
    khist = kmerHist(reads, k)
    corrected_reads = []
    for i, read in enumerate(reads):
        corrected = correct1mm(read, k, khist, 'ACGT', thresh)
        if corrected != read:
            print('%s => %s' % (read, corrected))
        corrected_reads.append(corrected)
    return corrected_reads


k = 25
name = 'reads5'

reads = reads_from_file(open('./reads/%s.fasta' % name, 'r'), 1001)
khist = kmerHist(reads, k)
new_reads = remove_errors(reads, k, 1)
wrong_kmers = remove_rare(new_reads, k)

build_graph(new_reads, k, name, wrong_kmers, mean(khist.values()))

# DeBruijnGraph.simplification()
