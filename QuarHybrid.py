import sys, os, re

if __name__ == '__main__' :
    q1 = {}
    
    for fname in sys.argv[1:] :
        with open(fname) as fin :
            for line in fin :
                if line.startswith('#') : continue
                part = line.strip().split('\t')
                if part[1] in q1 :
                    q1[part[1]].append(fname+':'+part[0])
                else :
                    q1[part[1]] = [fname+':'+part[0]]
    
    for q, trees in q1.iteritems() :
        print '{0}\t{1}'.format(q, ','.join(trees))