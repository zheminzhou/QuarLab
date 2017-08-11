import sys, json, re
from multiprocessing import Pool
import dendropy as dp, numpy as np

def get_splits(tre) :
    solo, duo, quartet = {}, {}, {}
    taxa = {}
    for node in tre.postorder_node_iter() :
        taxa[node] = list(set(node.annotations.get_value('taxa', []) + ([node.taxon.label] if node.taxon is not None else [])))
        for i1, t1 in enumerate(taxa[node]) :
            solo[t1] = 0.0
            for t2 in taxa[node][(i1+1):] :
                duo[(t1, t2)] = duo[(t2, t1)] = [0.0, 0.0, 0.0]
        nodes = [node] + node.child_nodes()
        for i1, n1 in enumerate(nodes) :
            for n2 in nodes[(i1+1):] :
                for t1 in taxa[n1] :
                    for t2 in taxa[n2] :
                        if t1 < t2 :
                            duo[(t1, t2)] = [solo[t1] + solo[t2], solo[t1], solo[t2]]
                        else :
                            duo[(t2, t1)] = [solo[t1] + solo[t2], solo[t2], solo[t1]]
        taxa[node] += [t for n in node.child_nodes() for t in taxa[n]]
        if node.edge_length is not None :
            for t in taxa[node] :
                solo[t] += node.edge_length
    n_str, solo = len(solo), sorted(solo)

    for i1, a in enumerate(solo) :
        for i2 in xrange(i1+1, n_str) :
            b = solo[i2]
            for i3 in xrange(i2+1, n_str) :
                c = solo[i3]
                for i4 in xrange(i3+1, n_str) :
                    d = solo[i4]

                    dists = [[duo[(a,b)][0] + duo[(c,d)][0], 1, duo[(a,b)][1:] + duo[(c,d)][1:]], 
                             [duo[(a,c)][0] + duo[(b,d)][0], 2, duo[(a,c)][1:] + duo[(b,d)][1:]], 
                             [duo[(a,d)][0] + duo[(b,c)][0], 3, duo[(a,d)][1:] + duo[(b,c)][1:]]]
                    dists.sort(key=lambda x:x[0])

                    if dists[0][0] +2e-8 >= dists[1][0] :
                        dists[0] = [duo[(a,b)][0] + duo[(c,d)][0], 0, duo[(a,b)][1:] + duo[(c,d)][1:] + [0.0]]
                    else :
                        dists[0][2].append((dists[1][0]-dists[0][0])/2.0)
                    quartet[(a,b,c,d)] = dists[0]
    return quartet

if __name__ == '__main__' :
    trees = []
    for fname in sys.argv[1:] :
        try:
            trees.extend(dp.TreeList.get_from_path(fname, 'nexus'))
        except :
            trees.extend(dp.TreeList.get_from_path(fname, 'newick'))
    for id, t in enumerate(trees) :
        if not t.label :
            t.label = str(id)
    pool = Pool(5)
    quartets = pool.map(get_splits, trees)
    
    print '#Tree\t#Split\t#Internal\t#Tip_1\t#Tip_2\t#Tip_3\t#Tip_4'
    for tre, quartet in zip(trees, quartets) :
        for tips, topology in quartet.iteritems() :
            if topology[1] == 1 :
                print '{0}\t{1},{2}:{3},{4}'.format(tre.label, *tips) + '\t{4}\t{0}\t{1}\t{2}\t{3}'.format(*topology[2])
            elif topology[1] == 2 :
                print '{0}\t{1},{3}:{2},{4}'.format(tre.label, *tips) + '\t{4}\t{0}\t{1}\t{2}\t{3}'.format(*topology[2])
            elif topology[1] == 3 :
                print '{0}\t{1},{4}:{2},{3}'.format(tre.label, *tips) + '\t{4}\t{0}\t{1}\t{2}\t{3}'.format(*topology[2])
            elif topology[1] == 0 :
                print '{0}\t{1},{2},{3},{4}'.format(tre.label, *tips) + '\t{4}\t{0}\t{1}\t{2}\t{3}'.format(*topology[2])
