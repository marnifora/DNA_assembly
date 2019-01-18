def neighbors1mm(kmer, alpha):
    """ Generate all neighbors at Hamming distance 1 from kmer """
    neighbors = []
    for j in range(len(kmer)-1, -1, -1):
        oldc = kmer[j]
        for c in alpha:
            if c == oldc: continue
            neighbors.append(kmer[:j] + c + kmer[j+1:])
    return neighbors


def kmerHist(reads, k):
    """ Return k-mer histogram and average # k-mer occurrences """
    kmerhist = {}
    for read in reads:
        for kmer in [read[i:i+k] for i in range(0, len(read)-(k-1))]:
            kmerhist[kmer] = kmerhist.get(kmer, 0) + 1
    return kmerhist


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
