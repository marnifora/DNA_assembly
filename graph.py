from statistics import mean
import copy


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

    class Contig:
        """ Contig object is the set of nodes which create the contig
            sequence. Attribute weight is the average of all weights
            of nodes in the path. """

        def __init__(self, path):
            self.path = path
            self.weight = 0
            self.seq = ''

        def __str__(self):
            return self.seq

        def set_params(self, k):
            w = self.path[0].weights
            self.seq = self.path[0].km1mer
            for n in self.path[1:]:
                self.seq += n.km1mer[k - 2:]
                w += n.weights
            if not w:
                w = [0]
            self.weight = mean(w)

    def __init__(self, strIter, k, wrong_kmers, thresh, name, output=None, verbose=False):
        """ Build de Bruijn multigraph given string iterator and k-mer
            length k """

        self.name = name
        self.verbose = verbose
        if output is not None:
            self.output = output
        else:
            self.output = '%s_contigs.fasta' % name
        self.thresh = max(thresh, 1)
        self.k = k
        self.nodes = {}  # maps k-1-mers to Node objects
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
        # cutting edges with weight below thresh
        cutedges = self.cut_graph()
        # removing not connected nodes
        notconnected = 0
        for n in list(self.nodes.values()).copy():
            if len(n.children) == 0 and len(n.parents) == 0:
                del (self.nodes[n.km1mer])
                notconnected += 1
        # establishing head-nodes
        self.head = [n for n in self.nodes.values() if len(n.parents) == 0 and len(n.children) > 0]
        # merging linear nodes
        self.done = []
        for h in self.head:
            self.merge(h, True)

        if self.verbose:
            print('Threshold = %d' % thresh)
            print('Number of wrong kmers = %d' % len(wrong_kmers))
            print('Number of removed rare kmers = %d' % removed)
            print('Number of cut edges = %d' % cutedges)
            print('Number of removed, not-connected nodes = %d' % notconnected)
            print('Number of nodes = %d' % len(self.nodes))
            print('Real heads = %d' % len(self.head))
            for n in self.nodes.values():
                print('%s\t%d\t%d\t%s' % (n.km1mer, len(n.parents), len(n.children), n in self.done))
            print('Number of nodes after merging: %d' % len(self.nodes))

        # establishing heads after merging
        self.head = [n for n in self.nodes.values() if len(n.parents) == 0 and len(n.children) > 0]
        # looking for not-connected nodes
        self.solo = [n for n in self.nodes.values() if len(n.parents) == 0 and len(n.children) == 0]

        # looking for contigs based on not-connected, merged nodes
        nodes = copy.deepcopy(self.nodes)
        self.contigs = []
        for n in self.solo:
            if len(n.km1mer) >= 300:
                contig = self.Contig([n])
                contig.set_params(self.k)
                self.contigs.append(contig)
                del (nodes[n.km1mer])

        if self.verbose:
            print('Real heads after merging = %d' % len(self.head))
            print('Solo nodes after merging = %d' % len(self.solo))
            if self.contigs:
                print('Number of contigs based on solo nodes = %d' % len(self.contigs))
            else:
                print('No contigs based on solo nodes found!')

        # looking for contigs based on more than one node
        self.cut_contigs(nodes)

        n = set()
        for c in self.contigs:
            for cc in c.path:
                n.add(cc)
        self.used = len(n)/len(self.nodes)

        if self.contigs and self.verbose:
            print('Number of contigs = %d' % len(self.contigs))
            print('Mean number of nodes in one contig = %.2f' % mean([len(l.path) for l in self.contigs]))

    def cut_graph(self):
        """ Delete edges which weight is lower than threshold. """

        cut = 0
        for n in self.nodes.values():
            chs = list(n.children.keys())
            for ch in chs:
                if n.children[ch] <= self.thresh:
                    cut += 1
                    del (n.children[ch])
                    del (ch.parents[n])
        return cut

    def contigs_to_file(self):
        """ Write found contigs into fasta file. """

        o = open(self.output, 'w')
        for i, contig in enumerate(self.contigs):
            o.write('>contig%d\n%s\n' % (i, contig))
        o.close()

    def remove_rare_kmers(self, wrong_kmers):
        """ Remove kmers from input list wrong_kmers,
            but only if given kmer is not head (no input
            edges) or tail (no output edges). """

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

    def merge(self, node, direction):
        """ Merge nodes which are connected only with each other.
            If direction is True then function goes down the graph
            (from parents to children), conversely if direction
            is False. """

        if node not in self.done:
            self.done.append(node)
            if direction:
                downstream = node.children
                upstream = node.parents
            else:
                downstream = node.parents
                upstream = node.children
            if len(upstream) > 1:
                for n in upstream.keys():
                    self.merge(n, not direction)
                for n in downstream.keys():
                    self.merge(n, direction)
                return 0
            elif len(downstream) == 1:
                nnode = next(iter(downstream.keys()))
                if direction:
                    if len(nnode.parents) != 1:
                        for n in [i for i in nnode.parents if i != node]:
                            self.merge(n, not direction)
                        self.merge(nnode, direction)
                        return 0
                    else:
                        del (self.nodes[nnode.km1mer])
                        nnode.km1mer = node.km1mer + nnode.km1mer[self.k - 2:]
                        nnode.parents = node.parents
                        for p in node.parents:
                            p.children[nnode] = p.children[node]
                            del (p.children[node])
                else:
                    if len(nnode.children) != 1:
                        for n in [i for i in nnode.children if i != node]:
                            self.merge(n, not direction)
                        self.merge(nnode, direction)
                        return 0
                    else:
                        del (self.nodes[nnode.km1mer])
                        nnode.km1mer = nnode.km1mer[:-self.k + 2] + node.km1mer
                        nnode.children = node.children
                        for p in node.children:
                            p.parents[nnode] = p.parents[node]
                            del (p.parents[node])

                nnode.weights = node.weights + [downstream[nnode]]
                if nnode.km1mer in self.nodes:
                    unode = self.nodes[nnode.km1mer]
                    unode.weights.append(nnode.weights)
                    for el in nnode.children:
                        unode.children[el] = unode.children.setdefault(el, 0) + nnode.children[el]
                    for el in nnode.parents:
                        unode.parents[el] = unode.parents.setdefault(el, 0) + nnode.parents[el]
                    self.merge(unode, direction)
                else:
                    del (self.nodes[node.km1mer])
                    self.nodes[nnode.km1mer] = nnode
                    self.merge(nnode, direction)

            else:
                for ch in downstream.keys():
                    self.merge(ch, direction)
                return 0
        else:
            return 0

    def all_contigs(self, node, contig):
        """ Return all possible contigs started from the given node. """

        if node not in self.done:
            self.done.append(node)
            contig.append(self.nodes[node.km1mer])
            if len(node.children) == 0:
                self.clist.append(contig)
                return 0
            elif len(node.children) == 1:
                self.all_contigs(next(iter(node.children.keys())), contig)
            else:
                for ch in node.children.keys():
                    self.all_contigs(ch, contig)
        elif contig:
            self.clist.append(contig)
            return 0

    def find_contigs(self, heads):
        """ Return the best contig for every given head-node. """

        contigs = {}
        for h in heads:
            self.done = []
            self.clist = []
            self.all_contigs(h, [])
            best = None
            max_weight = 0
            for c in self.clist:
                cc = self.Contig(c)
                cc.set_params(self.k)
                if len(cc.seq) >= 300 and cc.weight > max_weight:
                    max_weight = cc.weight
                    best = cc
            if best is not None:
                contigs[h] = best
        return contigs

    def cut_contigs(self, nodes):
        """ Searching for the best not-overlapping contigs. """

        altcontigs = False
        contigs = False
        cc = True
        i = 1
        while cc != contigs:
            cc = contigs
            heads = [n for n in nodes.values() if len(n.parents) == 0]
            if len(heads) == 0 or altcontigs:
                if nodes:
                    maks = max([len(n.children)-len(n.parents) for n in nodes.values()])
                    heads = [n for n in nodes.values() if len(n.children)-len(n.parents) == maks]
                else:
                    break
            contigs = self.find_contigs(heads)
            if not contigs:
                altcontigs = True
                i += 1
                continue
            maks = max([len(el.seq)*el.weight for el in contigs.values()])
            k = [kk for kk in contigs.keys() if len(contigs[kk].seq)*contigs[kk].weight == maks][0]
            self.contigs.append(contigs[k])
            for n in contigs[k].path:
                nn = nodes[n.km1mer]
                for ch in nn.children:
                    del (nodes[ch.km1mer].parents[nn])
                for pa in nn.parents:
                    del (nodes[pa.km1mer].children[nn])
                del (nodes[nn.km1mer])
            i += 1

    def to_dot(self):
        """ Write dot representation to given filehandle.  If 'weights'
            is true, label edges corresponding to distinct k-1-mers
            with weights, instead of writing a separate edge for each
            copy of a k-1-mer. """

        dotFh = open('%s_graph.dot' % self.name, 'w')
        dotFh.write("digraph \"Graph\" {\n")
        dotFh.write("  bgcolor=\"transparent\";\n")
        for node in self.nodes.values():
            lab = node.km1mer
            dotFh.write("  %s [label=\"%s\"] ;\n" % (lab, lab))
        for src in self.nodes.values():
            srclab = src.km1mer
            for dst, v in src.children.items():
                dstlab = dst.km1mer
                dotFh.write("  %s -> %s [label=\"%d\"] ;\n" % (srclab, dstlab, v))
        dotFh.write("}\n")
        dotFh.close()

    def to_csv(self):

        file = open('%s_table.csv' % self.name, 'w')
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
