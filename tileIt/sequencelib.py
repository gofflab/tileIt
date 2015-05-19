#/usr/bin/env python
import string,operator,random,math
from itertools import product

######
#Parsers
######
def FastaIterator(handle):
    """
    Generator function to iterate over fasta records in <handle>:
    Use in a loop to apply to each Seq record contained in a .fasta file
    Input: record handle as obtained by handle = open(<file>,'r')
    Returns an iterator across Sequences in file
    """
    #Skip any header text
    while True:
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line [0] == ">":
            break
    
    while True:
        if line[0] <> ">":
            raise ValueError("Records in Fasta files should start with a '>' character")
        name = line[1:].rstrip()
        lines = []
        line = handle.readline()
        while True:
            if not line : break
            if line[0] == ">" : break
            lines.append(line.rstrip().replace(" ",""))
            line = handle.readline()
        #Return record then continue
        newSeq = {'name':name,'sequence':"".join(lines)}
        yield newSeq
        
        if not line : return #StopIteration
    assert False, "Should not reach this line"

def dbSNPFastaIterator(handle):
    """
    Generator function to iterate over fasta records FROM dbSNP in <handle>:
    Use in a loop to apply to each Seq record contained in a .fasta file
    Input: record handle as obtained by handle = open(<file>,'r')
    Returns an iterator across Sequences in file
    """
    #Skip any header text
    while True:
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line [0] == ">":
            break
    
    while True:
        if line[0] <> ">":
            raise ValueError("Records in Fasta files should start with a '>' character")
        header = line[1:].rstrip()
        fields = header.split("|")
        newSeq = {'name' : fields[2].split(" ")[0]}
        for field in fields[3:]:
            k,v = field.split("=")
            v = v.strip('"')
            newSeq[k] = v
        lines = []
        line = handle.readline()
        while True:
            if not line : break
            if line[0] == ">" : break
            lines.append(line.rstrip().replace(" ",""))
            line = handle.readline()
        #Return record then continue
        newSeq['sequence'] = "".join(lines)
        try:
            newSeq['pos'] = int(newSeq['pos'])
        except KeyError:
            pass

        yield newSeq
        
        if not line : return #StopIteration
    assert False, "Should not reach this line"
    

###
#Generic Sequence tools
###

def complement(s):
    comp = {'A':'T',
            'T':'A',
            'U':'A',
            'G':'C',
            'C':'G',
            'Y':'R',
            'R':'Y',
            'S':'S',
            'W':'W',
            'K':'M',
            'M':'K',
            'B':'V',
            'D':'H',
            'H':'D',
            'V':'B',
            'N':'N',
            }
    complseq = [comp[base] for base in s]
    return complseq

def reverse_complement(s):
    seq = list(s)
    seq.reverse()
    return ''.join(complement(seq))

def rcomp(s):
    """Does same thing as reverse_complement only cooler"""
    return s.translate(string.maketrans("ATUGCYRSWKMBDHVN","TAACGRYSWMKVHDBN"))[::-1]

def getTm(seq):
    Tm = 79.8 + 18.5*math.log10(0.05) + (58.4 * getGC(seq)) + (11.8 * getGC(seq)**2) - (820/len(seq))
    return Tm

def getGC(seq):
    return (seq.count('C')+seq.count('G')+seq.count('c')+seq.count('g'))/float(len(seq))

def gc_content(seq):
    gc = mcount(seq, 'GCgc')
    at = mcount(seq, 'ATUatu')
    return 100*gc/float((gc+at))

def mcount(s, chars):
    # sums the counts of appearances of each char in chars
    count = 0
    for char in chars:
        count = count+string.count(s,char)
    return count

def prob_seq(seq, pGC=.5):
    # given a GC content, what is the probability  
    # of getting the particular sequence
        
    assert(0<=pGC<=1)
    # the probability of obtaining sequence seq
    # given a background gc probability of .5
    ps = []
    for char in seq:
        if char in 'CG': ps.append(pGC/2)
        elif char in 'AT': ps.append((1-pGC)/2)
        else: raise "Unexpected char: ",char
    return reduce(operator.mul, ps, 1)

def transcribe(seq):
    RNA = seq.replace('T', 'U')  
    return RNA

def GenRandomSeq(length, type='DNA'):
    if type == 'DNA':
        chars = ['A','T','G','C']
    if type == 'RNA':
        chars = ['A','U','G','C']
    return ''.join([random.choice(chars) for i in range(length)])

def seed():
    random.seed()
    
def draw(distribution):
    sum=0
    r = random.random()
    for i in range(0,len(distribution)):
        sum += distribution[i]
        if r< sum:
            return i

def makeDistFromFreqs(freqs):
    res = []
    chars = ['A','T','C','G']
    cum = 0
    res.append(cum)
    for i in chars:
        cum += freqs[i]
        res.append(cum)
    return res

def genRandomFromDist(length,freqs):
    """Generates a random sequence of length 'length' drawing from a distribution of
    base frequencies in a dictionary"""
    myDist = makeDistFromFreqs(freqs)
    chars = ['A','T','C','G']
    return ''.join([chars[draw(myDist)] for i in range(length)])

#############
# IUPAC tools
##############
# iupacdict = {'A':'A',
#     'C':'C',
#     'G':'G',
#     'T':'T',
#     'M':'[AC]',
#     'R':'[AG]',
#     'W':'[AT]',
#     'S':'[CG]',
#     'Y':'[CT]',
#     'K':'[GT]',
#     'V':'[ACG]',
#     'H':'[ACT]',
#     'D':'[AGT]',
#     'B':'[CGT]',
#     'X':'[ACGT]',
#     'N':'[ACGT]'}

iupacdict = {'A':'A',
    'C':'C',
    'G':'G',
    'T':'T',
    'a':'A',
    'c':'C',
    'g':'G',
    't':'T',
    'M':'AC',
    'R':'AG',
    'W':'AT',
    'S':'CG',
    'Y':'CT',
    'K':'GT',
    'V':'ACG',
    'H':'ACT',
    'D':'AGT',
    'B':'CGT',
    'X':'ACGT',
    'N':'ACGT',
    'n':'ACGT'}

def disambiguateIUPAC(seq,d=iupacdict):
    return list(map("".join, product(*map(d.get, seq))))

###########
#Motif Tools
###########
def allindices(string, sub, listindex=[], offset=0):
    i = string.find(sub, offset)
    while i >= 0:
        listindex.append(i)
        i = string.find(sub, i + 1)
    return listindex

def find_all(seq, sub):
    #print "Looking for %s in %s"%(sub,seq)
    found = []
    next = string.find(seq,sub)
    while next != -1:
        found.append(next)
        next = string.find(seq,sub,next+1)
    return found

def kmer_dictionary_counts(seq,k,dic={}):
    """Returns a dictionary of k,v = kmer:'count of kmer in seq'"""
    for i in range(0, len(seq)-k):
        subseq = seq[i:][:k]
        #if not dic.has_key(subseq): dic[subseq] = 1
        #else: dic[subseq] = dic[subseq] + 1
        #OR
        dic[subseq] = 1 + dic.get(subseq,0)
    return dic

def kmer_dictionary(seq,k,dic={},offset=0):
    """Returns dictionary of k,v = kmer:'list of kmer start positions in seq' """
    for i in range(0,len(seq)-k):    
        subseq = seq[i:][:k]
        dic.setdefault(subseq,[]).append(i+1)
    return dic

def get_seeds(iter,seeds={}):
    counter = 0
    for i in iter:
        counter+=1
        if counter%10000==0:
            print "%d" % counter
        i.CSToDNA()
        seed = i.sequence[1:8]
        seeds[seed] = 1 + seeds.get(seed,0)
    return seeds