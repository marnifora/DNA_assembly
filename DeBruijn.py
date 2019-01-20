from collections import deque
from statistics import mean


class DeBruijnGraph:
    """ A de Bruijn multigraph built from a collection of strings.
        User supplies strings and k-mer length k.  Nodes of the de
        Bruijn graph are k-1-mers and edges correspond to the k-mer
        that joins a left k-1-mer to a right k-1-mer. """

    @staticmethod
    def chop(st, k):
        """ Chop a string up into k mers of given length """
        for i in range(0, len(st)-(k-1)):
            yield (st[i:i+k], st[i:i+k-1], st[i+1:i+k])

    class Node:
        """ Node in a de Bruijn graph, representing a k-1 mer.  We keep
            track of # of incoming/outgoing edges so it's easy to check
            for balanced, semi-balanced. """
        
        def __init__(self, km1mer):
            self.km1mer = km1mer
            self.children = {}
            self.parents = {}
            self.weights = []
        
        def __str__(self):
            return self.km1mer
    
    def __init__(self, strIter, k, wrong_kmers, thresh, name):
        """ Build de Bruijn multigraph given string iterator and k-mer
            length k """

        self.name = name
        self.k = k
        self.G = {}     # multimap from nodes to neighbors
        self.nodes = {} # maps k-1-mers to Node objects
        for st in strIter:
            for kmer, km1L, km1R in self.chop(st, k):
                if km1L in self.nodes:
                    nodeL = self.nodes[km1L]
                else:
                    nodeL = self.nodes[km1L] = self.Node(km1L)
                if km1R in self.nodes:
                    nodeR = self.nodes[km1R]
                else:
                    nodeR = self.nodes[km1R] = self.Node(km1R)

                nodeL.children[nodeR] = nodeL.children.setdefault(nodeR, 0) + 1
                nodeR.parents[nodeL] = nodeR.parents.setdefault(nodeL, 0) + 1

        # removing wrong kmers, only if they aren't head or tail
        removed = self.remove_rare_kmers(wrong_kmers)
        print('Number of removed rare kmers = %d' % removed)

        print('Number of nodes: %d' % len(self.nodes))

        notconnected = 0
        for n in list(self.nodes.values()).copy():
            if len(n.children) == 0 and len(n.parents) == 0:
                del (self.nodes[n.km1mer])
                notconnected += 1

        print('Number of not connected nodes: %d' % notconnected)

        self.head = [n for n in self.nodes.values() if len(n.parents) == 0 and len(n.children) > 0]
        print('Real heads = %d' % len(self.head))
        self.tail = [n for n in self.nodes.values() if len(n.children) == 0 and len(n.parents) > 0]
        print('Real tails = %d' % len(self.tail))

        self.done = []
        for h in self.head:
            self.merge_down(h)

        for h in self.tail:
            self.merge_up(h)

        i = 0
        for n in self.nodes.values():
            if len(n.km1mer) >= 300:
                i += 1
            print('%s\t%d\t%d\t%s' % (n.km1mer, len(n.parents), len(n.children), n in self.done))
        print('Number of nodes after merging: %d' % len(self.nodes))
        print(len(self.done))
        print('Number of nodes longer than 300: %d' % i)

        self.contigs = []
        self.best_contigs()
        print(self.contigs)
        print(len(self.contigs))
        print([len(l) for l in self.contigs])
        self.contigs_to_file()

        # removing nodes which are tips
        # print(self.remove_tips(thresh))

    def contigs_to_file(self):
        o = open('%s_contigs.fasta' % self.name, 'w')
        for i, contig in enumerate(self.contigs):
            o.write('>contig%d\n%s\n' % (i, contig))
        o.close()

    def to_csv(self):

        file = open('graph.csv', 'w')
        file.write('Source\tTarget\tweight\n')
        edges = {}
        for g in self.G:
            for gg in self.G[g]:
                if g not in edges:
                    edges[g] = [gg, 1]
                else:
                    edges[g][1] += 1

        for k in edges.keys():
            file.write('%s\t%s\t%d\n' % (k, edges[k][0], edges[k][1]))
        file.close()

    def remove_rare_kmers(self, wrong_kmers):

        def side_kmers(kmers, side):
            nodes = set()
            for k in kmers:
                if side == 'right':
                    nodes.add(k[1:])
                elif side == 'left':
                    nodes.add(k[:-1])
            return nodes

        removed = 0
        for lista, side in [[side_kmers(wrong_kmers, 'left'), 'left'], [side_kmers(wrong_kmers, 'right'), 'right']]:
            for s in lista:
                try:
                    node = self.nodes[s]
                except KeyError:  # if the same km1mer was in the other group (right/left)
                    continue
                if (side == 'left' and len(node.parents) > 0) or (side == 'right' and len(node.children) > 0):
                    for n in node.children.keys():
                        del (n.parents[node])
                    for n in node.parents.keys():
                        del (n.children[node])
                    del (self.nodes[s])
                    removed += 1
        return removed

    def remove_tips(self, thresh):
        """ Remove tips from graph. A node is considered a tip if it is disconnected
            on one of its ends, the length of the information stored in the node is
            shorter than 2k (important for simplified graph) and the weight of arc
            leading to this node is < thresh. """
        removed = []
        for node in self.nodes.values():
            if node.nin == 0 or node.nout == 0:
                if node not in self.G or len(self.G[node]) < thresh:
                    removed.append(node)
        for node in removed:
            del(self.nodes[str(node)])
            if node in self.G:
                del(self.G[node])
        return 'Tips removed = %d' % len(removed)

    def merge_down(self, node):

        if node not in self.done:
            self.done.append(node)
            if len(node.parents) > 1:
                for n in node.children.keys():
                    self.merge_down(n)
                    return 0
            if len(node.children) == 1:
                ch = next(iter(node.children.keys()))
                if len(ch.parents) == 1:
                    del (self.nodes[ch.km1mer])
                    del (self.nodes[node.km1mer])
                    node.weights.append(node.children[ch])
                    ch.km1mer = node.km1mer + ch.km1mer[self.k-2:]
                    ch.weights = node.weights
                    ch.parents = node.parents
                    self.nodes[ch.km1mer] = ch
                    self.merge_down(ch)
                else:
                    self.merge_down(ch)
                    return 0
            else:
                for ch in node.children.keys():
                    self.merge_down(ch)
                    return 0
        else:
            return 0

    def merge_up(self, node):

        if node not in self.done:
            self.done.append(node)
            if len(node.children) > 1:
                for n in node.parents.keys():
                    self.merge_up(n)
                    return 0
            if len(node.parents) == 1:
                ch = next(iter(node.parents.keys()))
                if len(ch.children) == 1:
                    del (self.nodes[ch.km1mer])
                    del (self.nodes[node.km1mer])
                    node.weights.append(node.parents[ch])
                    ch.km1mer = ch.km1mer[:-self.k+2] + node.km1mer
                    ch.weights = node.weights
                    ch.children = node.children
                    self.nodes[ch.km1mer] = ch
                    self.merge_up(ch)
                else:
                    self.merge_up(ch)
                    return 0
            else:
                for ch in node.parents.keys():
                    self.merge_up(ch)
                    return 0
        else:
            return 0

    def find_contigs(self, node, contig):

        if node not in self.done:
            self.done.append(node)
            contig.append(node)
            if len(node.children) == 0:
                self.all.append(contig)
                return 0
            elif len(node.children) == 1:
                self.find_contigs(next(iter(node.children.keys())), contig)
            else:
                for ch in node.children.keys():
                    self.find_contigs(ch, contig)
        elif contig:
            self.all.append(contig)

    def best_contigs(self):

        for h in self.head:
            self.all = []
            self.done = []
            self.find_contigs(h, [])
            approved = []
            for contig in self.all:
                l = contig[0].km1mer
                for el in contig[1:]:
                    l += el.km1mer[self.k-2:]
                if len(l) > 300:
                    approved.append(contig)
            if len(approved) > 1:
                maks = 0
                best = None
                for contig in approved:
                    cov = []
                    for el in contig:
                        cov += el.weights
                    cov = mean(cov)
                    if cov > maks:
                        maks = cov
                        best = contig
            elif len(approved) == 1:
                best = approved[0]
            else:
                continue
            l = best[0].km1mer
            for el in best[1:]:
                l += el.km1mer[self.k - 2:]
            self.contigs.append(l)

        
    def toDot(self, weights=False):
        """ Write dot representation to given filehandle.  If 'weights'
            is true, label edges corresponding to distinct k-1-mers
            with weights, instead of writing a separate edge for each
            copy of a k-1-mer. """
        dotFh = open('%s_graph.dot' % self.name, 'w')
        dotFh.write("digraph \"Graph\" {\n")
        dotFh.write("  bgcolor=\"transparent\";\n")
        for node in self.nodes():
            lab = node.km1mer
            dotFh.write("  %s [label=\"%s\"] ;\n" % (lab, lab))
        for src in self.nodes:
            dsts = src.children.keys()
            srclab = src.km1mer
            if weights:
                weightmap = {}
                if weights:
                    for dst in dsts:
                        weightmap[dst] = weightmap.get(dst, 0) + 1
                for dst, v in weightmap.items():
                    dstlab = dst.km1mer
                    dotFh.write("  %s -> %s [label=\"%d\"] ;\n" % (srclab, dstlab, v))
            else:
                for dst in dsts:
                    srclab = src.km1mer
                    dstlab = dst.km1mer
                    dotFh.write("  %s -> %s [label=\"\"] ;\n" % (srclab, dstlab))
        dotFh.write("}\n")
