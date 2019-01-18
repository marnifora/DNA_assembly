from DeBruijn import *
from error_correction import *


def reads_from_file(file, l):

    reads = []
    for i, line in enumerate(file):
        if i >= 2*l:
            break
        if not line.startswith('>'):
            reads.append(line.strip())
    return reads


def test_graph(reads, k, output):

    graph = DeBruijnGraph(reads, k)
    graph.toDot(output)

    path = list(graph.eulerianPath())
    assembly = path[0] + ''.join(map(lambda x: x[-1], path[1:]))

    print(path, assembly)


def test_error(reads, k, thresh):
    khist = kmerHist(reads, k)
    corrected_reads = []
    for i, read in enumerate(reads):
        corrected = correct1mm(read, k, khist, 'ACGT', thresh)
        if corrected != read:
            print('%s => %s' % (read, corrected))
        corrected_reads.append(corrected)
    return corrected_reads


k = 3
'''
reads = reads_from_file(open('./reads/reads0.fasta', 'r'), 1001)
khist = kmerHist(reads, k)
print(len(khist))
print(len([v for v in khist.values() if v == 1]))
new_reads = test_error(reads, k, 1)
new_khist = kmerHist(new_reads, k)
print(len(new_khist))
print(len([v for v in new_khist.values() if v == 1]))

# reads = ['AAABBBA']
'''
reads = ['to_every_thing_turn_turn_turn_there_is_a_season']
output = open('test.txt', 'w')
test_graph(reads, k, output)
