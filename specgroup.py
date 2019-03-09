##### CONFIGURATION ######

'''
REQUIRES: https://pypi.org/project/python-Levenshtein/ AND matprotlib
INPUT: input fasta filepath
OUTPUT: output fasta filepath
PROGRESS: show when steps are completed
GRAPH: show graph of cluster size frequencies
MAXSIZE: maximum allowed cluster size 
DISTTYPE: distance function to use: either 1 or 2
    1: edit distance
    2: edit distance divided by maximum length of the two strings
MAGIC: magic parameter to tune the distance function
    if None, it will automatically select parameter for the chosen dist function
'''

INPUT = 'human_no_red.fasta'
OUTPUT = 'human_grouped.fasta'
MAXSIZE = 20
DISTTYPE = 1
PROGRESS = True
GRAPH = False
MAGIC = 1.2


##########################


import Levenshtein as edit
try: 
    import matplotlib.pyplot as plt
    can_graph = True
except ImportError:
    can_graph = False
import random
import time

#process params to determine cutoff & distance vars
if MAGIC == None: 
    MAGIC = {1: 1.2, 2: 0.86}[DISTTYPE]
if DISTTYPE == 1: 
    CUTOFF = lambda x: max([2, len(x.pep)/MAGIC]) #Peptide -> number
    #the original MAGIC param here was 10, but the clusters are best when it's much lower
    DIST = lambda a, b: edit.distance(a, b) #string, string -> number
elif DISTTYPE == 2: 
    CUTOFF = lambda x: MAGIC
    DIST = lambda a, b: edit.distance(a, b)/float(max([len(a), len(b)]))
GRAPH = GRAPH and can_graph

class ClusterCollection: 
    '''A collection of Cluster classes, with writing and summarizing methods'''
    def __init__(self, clusters): 
        self.clusters = clusters
    def summarize(self): 
        #makes a graph of cluster sizes vs frequences
        L = [0 for _ in range(1, self.clusters[0].maxsize+2)]
        for i in self.clusters: 
            L[len(i.peptides)] += 1 
        plt.plot(range(1, self.clusters[0].maxsize+2), L, 'g', linewidth = 3)
        plt.title("Cluster size distribution")
        plt.xlabel("Cluster size")
        plt.ylabel("Frequency")
        plt.show()
    def write(self, filepath): 
        #writes centers, contents of clusters to file
        f = open(filepath, 'w')
        for i in range(len(self.clusters)): 
            f.write(">Group " + str(i+1) + " center\n")
            f.write(str(self.clusters[i].center) + '\n')
            for j in self.clusters[i].peptides: 
                f.write(">Group " + str(i+1) + "\n")
                f.write(str(j) + '\n')
        f.close()
    def merge(self): 
        pass
    
class Cluster: 
    '''A collection of similar Peptide classes, supporting addition and printing'''
    def __init__(self, center, peptides = []): 
        self.peptides = peptides
        self.center = center
        self.cutoff = CUTOFF(self.center)
        self.maxsize = MAXSIZE
    def add(self, pep): 
        if DIST(self.center.pep, pep.pep) > self.cutoff or len(self.peptides) == self.maxsize: 
            return False
        else: 
            self.peptides.append(pep)
            return True
    def __str__(self): 
        return '\n'.join(map(str, [self.center] + self.peptides))

class Peptide: 
    '''Pairs a peptide with a description, in case we need that at some point '''
    def __init__(self, peptide, description = None): 
        self.pep = peptide
        self.desc = description
    def __str__(self): 
        return self.pep
        
class PeptideCollection: 
    '''A collection of peptides, supporting clustering'''
    def __init__(self, peptides): 
        self.peptides = peptides
    def cluster(self): 
        #returns a ClusterCollection of the detected clusters
        clusters = []
        current = None
        i = 0
        makenew = True
        while i < len(self.peptides):
            if makenew: 
                del current
                current = Cluster(self.peptides[i])
                current.peptides = [] #fixes pass-by-reference errors 
                makenew = False
            elif not(current.add(self.peptides[i])): 
                clusters.append(current)
                makenew = True
                i -= 1
            i += 1
        clusters.append(current)
        return ClusterCollection(clusters)
            
def readFastaFile(filepath): 
    #returns PeptideCollection from filepath
    f = open(filepath).readlines() + ['>']
    peptides = []
    desc = ''
    pep = ''
    maxlen = 0
    for i in f: 
        if i[0] == '>': 
            if pep != '': 
                peptides.append(Peptide(pep, desc))
            desc = i.strip('\r\n')
            if len(pep) > maxlen: 
                maxlen = len(pep)
            pep = ''
        else: 
            pep += i.strip('\r\n')
    peptides = sorted(peptides, key = lambda x: (len(x.pep), x.pep))
    _peptides = [peptides[0]]
    for i in range(1, len(peptides)): 
        if peptides[i].pep != peptides[i-1].pep: 
            _peptides.append(peptides[i])
    #for i in peptides: 
    #    i.pep += (maxlen - len(i.pep)) * '*'
    return PeptideCollection(_peptides)

def cluster(infile, outfile, progress = True, summarize = True): 
    if progress: 
        print ('Reading input from ' + infile + '...')
        ti = time.time()
    P = readFastaFile(infile)
    if progress: 
        print ('File read (took ' + str(round(time.time()-ti, 4)) + ' s). Grouping...')
        ti = time.time()
    Q = P.cluster()
    if progress: 
        print ('Grouped (took ' + str(round(time.time()-ti, 4)) + ' s). Writing output to ' + outfile + '...')
        ti = time.time()
    Q.write(outfile)
    if progress: 
        print ('Writing complete! (took ' + str(round(time.time()-ti, 4)) + ' s)')
    if summarize: 
        Q.summarize()

cluster(INPUT, OUTPUT, PROGRESS, GRAPH)
