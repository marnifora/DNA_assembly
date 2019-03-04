from statistics import stdev, mean
import math


def neighbors1mm(kmer, alpha):
    """ Generate all neighbors at Hamming distance 1 from kmer """
    neighbors = []
    for j in range(len(kmer)-1, -1, -1):
        oldc = kmer[j]
        for c in alpha:
            if c == oldc:
                continue
            neighbors.append(kmer[:j] + c + kmer[j+1:])
    return neighbors


def kmerHist(reads, k):
    """ Return k-mer histogram and average # k-mer occurrences """
    kmerhist = {}
    for read in reads:
        for kmer in [read[i:i+k] for i in range(0, len(read)-(k-1))]:
            kmerhist[kmer] = kmerhist.get(kmer, 0) + 1
    return kmerhist, max(math.floor(mean(kmerhist.values()) - stdev(kmerhist.values())), 0)


def correct1mm(read, k, kmerhist, alpha, thresh):
    """ Return an error-corrected version of read.  k = k-mer length.
        kmerhist is kmer count map.  alpha is alphabet.  thresh is
        count threshold above which k-mer is considered correct. """
    # Iterate over k-mers in read
    for i in range(0, len(read)-(k-1)):
        kmer = read[i:i+k]
        # If k-mer is infrequent...
        if kmerhist.get(kmer, 0) <= thresh:
            # Look for a frequent neighbor
            for newkmer in neighbors1mm(kmer, alpha):
                if kmerhist.get(newkmer, 0) > thresh:
                    # Found a frequent neighbor; replace old kmer
                    # with neighbor
                    read = read[:i] + newkmer + read[i+k:]
                    break
    # Return possibly-corrected read
    return read


def rare_kmers(khist, thresh):
    """ Finding rare kmers, which should be omitted. """
    wrong_kmers = []
    for k in khist.keys():
        if khist[k] <= thresh:
            wrong_kmers.append(k)
    return wrong_kmers


def remove_errors(reads, k, khist):
    """ Removing errors from all given reads. """
    corrected_reads = []
    for read in reads:
        corrected = correct1mm(read, k, khist, 'ACGT', 1)
        corrected_reads.append(corrected)
    return corrected_reads
