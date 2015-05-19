#!/usr/bin/env python
'''
Created on Jun 25, 2009

@author: lgoff
'''
# import genomelib
import copy
import numpy as np
import algorithms
import os,sys,random,string,commands
import urllib2
from sequencelib import dbSNPFastaIterator


#Common
RNAFOLD = 'RNAfold -noPS'

#This is very human-specific at this point

class Interval:
    """Basic interval class, try to use ChipInterval or SeqInterval if possible...
        At this point, the Interval class is rather human specific so avoid calls to self.fetchSequence() or self.getChrNum(), etc...
    """
    def __init__(self, chr, start, end, strand="*", score=0.0, readcount = -1,name="",sequence = "",data={},genome="hg18"):
        
        #Check if creating new instance from old instance as 1st arg
        if isinstance(chr,Interval):
            interval = chr
            
            #Copy information from other instance
            self.chr = interval.chr
            self.start = interval.start
            self.end = interval.end
            self.strand = interval.strand
            self.score = interval.score
            self.readcount = interval.readcount
            self.sequence = interval.sequence
            self.data = copy.copy(interval.data)
            self.genome = interval.genome
            self.TSS = interval.TSS
   
        else:
            #default settings for new init
            self.chr=chr
            self.start = int(start)
            self.end = int(end)
            self.strand = strand
            if self.strand == "+":
                self.TSS = self.start
            elif self.strand == "-":
                self.TSS = self.end
            self.score = float(score) #can be used as a proxy for count, but I would prefer the readcount attribute
            self.readcount = int(readcount)
            if name == "":
                self.name = "%s:%d-%d:%s" % (self.chr,self.start,self.end,self.strand)
            else:
                self.name = name
            self.sequence = sequence
            self.structure = ''
            self.data = copy.copy(data)
            self.genome = genome
            self.startIndex = -1
            self.endIndex = -1
                   
    def getTSS(self):
        if self.strand == "+":
            self.TSS = self.start
        elif self.strand == "-":
            self.TSS = self.end
        return self.TSS
    
    def addChild(self, child):
        """Adds child node to self.children"""
        #assert child not in self.children
        #if child not in self.children:
        child.parents.append(self)
        self.children.append(child)
    
    def removeChild(self, child):
        """Removes child node from self.children (not sure how or if this works. Don't trust it yet)"""
        child.parents.remove(self)
        self.children.remove(child)
    
    def childScores(self):
        """Returns list of scores for each interval in self.children"""
        return [x.score for x in self.children]
    
    def makeValMap(self,value = 'readcount'):
        """Check these two to see which one is right..."""
        self.valMap = np.zeros(len(self))
        self.valMap = self.valMap-1
        myTmp = []
        for x in range(0,len(self)):
            myTmp.append([])
        for i in self.children:
            for j in range(i.start,i.end+1):
                myTmp[j-self.start].append(i.__dict__[value])
        for nt in range(0,len(myTmp)):
            if len(myTmp[nt])>0:
                self.valMap[nt]=sum(myTmp[nt])/len(myTmp[nt])
     
    def __iter__(self):
        return iter(self.sequence)
    
    def __getitem__(self,key):
        return self.sequence[key]
    
    def __repr__(self):
        if self.name == "":
            return "%s:%d-%d:%s" % (self.chr,self.start,self.end,self.strand)
        else:
            return self.name

    def __neg__(self):
        strandLookup = {"+":"-","-":"+"}
        newStrand = strandLookup[self.strand]
        return Interval(self.chr,self.start,self.end,newStrand,self.score,self.readcount)
    
    def __len__(self):
        return self.end-self.start+1
    
    def __str__(self):
        if self.sequence != "":
            return self.sequence
        else:
            return self.name
    
    def __cmp__(self,b):
        if self.equals(b):return 0
        chrTest = cmp(self.getChrNum(),b.getChrNum())
        if chrTest==0:
            mid1 = (self.start+self.end)/2 
            mid2 = (b.start+b.end)/2
            return cmp(mid1,mid2)
        else:
            return chrTest 
    
    def windows(self,windowSize):
        """Generator that yields windows across the interval self.start -- self.end"""
        for i in range(0,len(self)-windowSize):
            yield (i,i+windowSize)
    
    def toBed(self,value = 'score'):
        """Change value to readcount to return number of reads within interval"""
        return "%s\t%d\t%d\t%s\t%.2f\t%s" %(self.chr,self.start,self.end,self.name,self.__dict__[value],self.strand)
    
    def toUCSC(self):
        return "%s:%d-%d" % (self.chr,self.start,self.end)
    
    def toStringNumIGV(self):
        return "%s\t%d" % (self.chr.replace("chr",""),self.start)
    
    def toFasta(self):
        return ">%s\n%s" % (self.name,self.sequence)
    
    def getString(self):
        return "%s:%d-%d:%s" % (self.chr,self.start,self.end,self.strand)
    
    def getScore(self):
        return self.score
    
    def getStrand(self):
        return self.strand
    
    def mature(self,start,end):
        """Can be used to treat self as a microRNA Precursor.  By using matureStart and matureEnd you can define miRNA boundaries."""
        self.matureStart = start
        self.matureEnd = end
    
#    def overlaps_old(self,b):
#        """Return true if b overlaps self"""
#        if b.chr != self.chr :return False
#        if (self.end>b.end and self.end<b.end) or (b.end>=self.start and b.end<=self.end) or (b.start==self.start) or (b.end == self.end):
#            return True
#        else:
#            return False
    
    def overlaps(self,b):
        """Return true if b overlaps self"""
        if b.chr != self.chr :return False
        if (self.start <= b.start and b.start <=self.end) or (self.start >= b.start and self.start <= b.end):
            return True
        else:
            return False
    
    def distance(self,b,enforceStrand=False):
        """Returns absolute distance between self and another interval start positions.  
        Returns -1 if they are on different chromosome. If enforceStrand=True, then this function requires that both intervals
        be on the same strand. If they aren't, -1 is returned.
        """
        if b.chr != self.chr:
            return -1
        if enforceStrand:
            if self.strand!=b.strand:
                return -1
        else:
            return abs(self.start-b.start)
    
    def distanceBetweenTSS(self,b):
        """
        Returns the distance between two interval TSSs.
        """ 
        if self.chr != b.chr:
            return False
        if self.strand == "+":
            return b.TSS-self.TSS
        elif self.strand == "-":
            return self.TSS-b.TSS
        else:
            return False
    
    def findDist(self,b):
        """
        """
        if self.strand == "+" and b.strand == "+":
            return b.start-self.TSS
        elif self.strand == "+" and b.strand == "-":
            return b.end-self.TSS
        elif self.strand == "-" and b.strand =="+":
            return self.TSS-b.start
        elif self.strand == "-" and b.strand == "-":
            return self.TSS-b.end
        
    def isFullyContained(self,b):
        """Returns True if b is fully contained within self"""
        if b.chr != self.chr: return False
        if(b.start>=self.start and b.end<=self.end):return True
        else:
            return False
    
    def equals(self,b):
        """Returns True if b has the same start and end as self"""
        if (self.chr != b.chr): return False
        if (self.start==b.start and self.end == b.end):return True
        else:
            return False
    
    def getChrNum(self):
        """Assumes human (hg18) but fetches a chromosome 'number' to be used for sorting"""
        chrLookup = {"X":23,"x":23,"Y":24,"y":24}
        if self.chr.startswith("chr"):
            num = self.chr[3:]
            if num in ("X","x","Y","y"):
                num = chrLookup[num]
            return int(num)
        else: return self.chr
    
    def fetchSequence(self):
        if self.genome != "":
            genome = genomelib.pygrConnect(self.genome)
            seq = genome[self.chr][self.start-1:self.end]
            if self.strand =="-":
                seq = -seq
            self.sequence = str(seq)
        else:
            self.sequence = ''
        return self.sequence
    
    def fetchSequence2(self,contig = None):
        """Trying to be faster than fetchSequence by providing the pygr chromosome as an argument ('contig').  This should help avoid having to make multiple calls and speed
        up the sequence retrieval.
        """
        if contig == None:
            connection = genomelib.pygrConnect(self.genome)
            seq = connection[self.chr][self.start-1:self.end]
        else:
            seq = contig[self.start-1:self.end]
        if self.strand == "-":
            seq = -seq
        self.sequence = str(seq)
        return self.sequence

    def transcribe(self):
        """Makes sequence into RNA"""
        self.sequence = self.sequence.replace("T","U")
        return
    
    def getGC(self):
        """Returns GC fraction of self.sequence"""
        numGC = self.sequence.upper().count("G") + self.sequence.upper().count("C")
        self.gc = float(numGC)/len(self.sequence)
        return self.gc 
    
    def getPromoter(self,promUp=2000,promDown=0):
        if self.strand == "+":
            align = Interval(self.chr,self.start-promUp,self.start+promDown,self.strand,score=self.score,name=self.name+"_promoter")
        elif self.strand == "-":
            align = Interval(self.chr,self.end-promDown,self.end+promUp,self.strand,score=self.score,name = self.name+"_promoter")
        return align  
    
    def fold(self):
        command = "echo '%s' | %s" % (self.sequence,RNAFOLD)
        output = commands.getoutput(command)
        if len(output.split())>2:
            self.structure,self.mfe = output.split()[1:]
            self.mfe = float(self.mfe.strip("(").rstrip(")"))
        else:
            self.structure,self.mfe = ("nan","nan")
        return

    def getStructureFasta(self):
        return ">%s\n%s" % (self.name,self.structure)

    def isPlus(self):
        if self.strand=="+":
            return True
        else:
            return False
        
    def isMinus(self):
        if self.strand=="-":
            return True
        else:
            return False
    
    def nmer_dictionary(self,n,dic={}):
        """
        Returns nmer_dictionary from self.sequence
        """
        if self.sequence == "":
            self.fetchSequence()
        self.sequence = self.sequence.upper()
        for i in range(0,len(self.sequence)-n):
            subseq = self.sequence[i:][:n]
            dic[subseq]=1+dic.get(subseq,0)
            del subseq
        return dic

    def intersects(self,b,start='start',end='end',offset=0):
        if self.chr == b.chr and self.strand==b.strand:
            return not(self.start>b.end+offset or b.start>self.end+offset)
        else:
            return False
    
    def grow5_prime(self,length):
        if self.strand == "+":
            self.start = self.start-length
        elif self.strand == "-":
            self.end = self.end+length
    
    def grow3_prime(self,length):
        if self.strand == "+":
            self.end = self.end+length
        elif self.strand == "-":
            self.start = self.start-length

class SNP(Interval):
    def __init__(self,chr,start,end,strand="*",score=0.0,name="",refAllele="",nonrefAllele="",genome="hg19"):
        Interval.__init__(self,chr,start,end,strand,score=score,name=name,sequence=sequence,genome=genome)
        self.refAllele = refAllele
        self.nonrefAllele = nonrefAllele

    def fetchSequence(self,halfWidth=0):
        #chrom = 'chr1'
        flankStart = self.start-halfWidth
        flandEnd = self.end+halfWidth
        urlbase = "http://genome.ucsc.edu/cgi-bin/das/%s/dna?segment=%s:%d,%d"
        request = urllib2.Request(urlbase % (self.chr,self.chr,flankStart,flankEnd))
        u = urllib2.urlopen(request)
        tree = ElementTree.parse(u)
        rootElem = tree.getroot()
        sequenceElem = rootElem.find("SEQUENCE")
        DNAElem = sequenceElem.find("DNA")
        return DNAElem.text.replace("\n","")

class dbSNP():
    def __init__(self,name="",snpPos=-1,score=0.0):
        self.name = name
        self.score = score
        self.snpPos = snpPos
        self.fetchSequence()

    def fetchSequence(self):
        baseurl = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.cgi?db=snp&id=%s&rettype=fasta&retmode=text"
        request = urllib2.Request(baseurl % self.name)
        response = urllib2.urlopen(request)
        #parse fasta record
        fastaIter = dbSNPFastaIterator(response)
        record = fastaIter.next()
        #Update chr, start, end, sequence snpPos, snpClass
        self.sequence = record['sequence']
        self.alleles = record['alleles']
        self.taxid = record['taxid']
        self.snpPos = record['pos']-1
        self.GMAF = record['GMAF']
        self.varClass = record['class']

    def __len__(self):
        return len(self.sequence)


class SplicedInterval(Interval):
    """Extends Interval and Adds/overwrites methods to incorporate spliced elements"""
    def __init__(self, chr, start, end, strand="*",exonLengths=[],exonOffsets=[],score=0.0, readcount = -1,name="",sequence = "",data={},genome="hg18"):
        Interval.__init__(self,chr,start,end,strand,score=score, readcount = readcount,name=name,sequence = sequence,data=data,genome=genome)
        self.exonLengths = [int(x) for x in exonLengths.rstrip(",").split(",")]
        self.exonOffsets = [int(x) for x in exonOffsets.rstrip(",").split(",")]
        self.exonStarts = [self.start+self.exonOffsets[i] for i in xrange(0,len(self.exonOffsets))]
        self.exonEnds = [self.start+self.exonOffsets[i]+self.exonLengths[i] for i in xrange(0,len(self.exonStarts))]
        self.numExons = len(self.exonStarts)
    
    def __len__(self):
        return self.CDSlen()
    
    def intervalLen(self):
        """Length of genomic footprint for self (ie. end-start+1)"""
        return self.end-self.start+1
    
    def CDSlen(self):
        """Returns length of the exons"""
        return sum(self.exonLengths)
    
    def getExons(self):
        """Returns list of intervals corresponding to exonic sequences for self"""
        rtrn = []
        for i in range(0,len(self.exonStarts)):
            rtrn.append(Interval(self.chr,self.exonStarts[i],self.exonEnds[i],self.strand,name = self.name+"_exon_"+str(i+1)))
        return rtrn
    
    def getIntrons(self):
        """Returns list of intervals corresponding to intronic sequences for self"""
        rtrn = []
        for i in range(0,len(self.exonStarts)-1):
            rtrn.append(Interval(self.chr,self.exonEnds[i]+1,self.exonStarts[i+1]-1))
        return rtrn
    
    def fetchSplicedSequence(self):
        """Self explanatory"""
        connection = genomelib.pygrConnect(self.genome)
        components = []
        for i in xrange(0,len(self.exonStarts)):
            components.append(connection[self.chr][self.exonStarts[i]:self.exonStarts[i]+self.exonLengths[i]])
        if self.strand =="-":
            components = [-x for x in components]
            components = components[::-1]
        self.splicedSequence = "".join([str(x) for x in components])
        self.sequence = self.splicedSequence
    
    def toFasta(self):
        """Return fasta format"""
        return ">%s\n%s" % (self.name,self.splicedSequence)
    
    def toBed(self,value = 'score',rgb='0,0,0'):
        """Change value to readcount to return number of reads within interval"""
        return "%s\t%d\t%d\t%s\t%.2f\t%s\t%d\t%d\t%s\t%d\t%s\t%s" %(self.chr,self.start,self.end,self.name,self.__dict__[value],self.strand,self.start,self.end,rgb,len(self.exonStarts),",".join([str(x) for x in self.exonLengths]),",".join([str(x) for x in self.exonOffsets]))
    
    def makePNG(self,outDir=os.getcwd(),tmpFname='temp.R'):
        """
        Draws transcript structure of the interval to the file 'self.name'.png
        """
        rscript = """
name<-'%s'
contig<-'%s'
start<-%d
end<-%d
strand<-'%s'
exonLengths<-c(%s)
exonOffsets<-c(%s)
myLen<-end-start+1

png(filename=paste('%s/',name,'.png',sep=''),width=900,height=300)
plot.new()
plot.window(xlim=c(start,end),ylim=c(0,3))
axis(1)
title(xlab=contig)
title(main=name)
lines(seq(start,end+1),rep(1,myLen+1),col='blue',lwd=2,lend='butt')

segments(start+exonOffsets,rep(1,length(exonOffsets)),start+exonOffsets+exonLengths,rep(1,length(exonOffsets)),col='blue',lwd=20,lend='butt')
if (strand=='+'){
    arrows(start,1.5,(start+(myLen*0.05)),1.5,length=0.125,lwd=1.5,angle=30,col='black')
} else if (strand=='-') {
    arrows(end,0.5,(end-(myLen*0.05)),0.5,length=0.125,lwd=1.5,angle=30,col='black')
}


dev.off()""" % (self.name,self.chr,self.start,self.end,self.strand,",".join([str(x) for x in self.exonLengths]),",".join([str(x) for x in self.exonOffsets]),outDir)
        tmpHandle = open(tmpFname,'w')
        print >>tmpHandle, rscript
        tmpHandle.close()
        commands.getoutput('R CMD BATCH --vanilla %s' % tmpFname)
        os.remove(tmpFname)
        return
                    
                
########
#Generic interval operations
########

def findIntervalPos(intervals,pos):
    """Find the first interval that starts after 'pos' in a sorted list of 'Intervals'"""
    low,top = algorithms.binsearch(intervals,pos-1,lambda a,b: cmp(a.start,b))
    return (low,top)

def findInterval(intervals,interval):
    """Find an interval in a sorted list of 'intervals'"""
    low,ind = algorithms.binsearch(intervals,interval.start-1,lambda a,b: cmp(a.start,b))
    return (low,ind)
    
def iterChrom(intervals,start,end,index = None):
    """An iterator that walks down a sorted list of intervals"""
    
    nintervals = len(intervals)
    #find index
    if index == None:
        #find starting index by binary search
        index = findIntervalPos(intervals,start)
        if index == None:
            return
    
    #walk down chromosome
    while index < nintervals and intervals[index].start < end:
        yield intervals[index]
        index+=1

def overlaps(interval,intervals):
    """Find the intervals in list 'intervals' that overlap interval"""
    return [x for x in intervals if interval.overlaps(x)]

def intervalLookup(intervals,key = "ID"):
    """
    Returns a dict lookup of regions based on a key (default = "ID")
    """
    lookup = {}
    
    for interval in intervals:
        ikey = None
        
        if key in interval.data:
            ikey = interval.data[key]
        else:
            ikey = key(interval)
        
        if ikey is not None:
            assert ikey not in lookup, Exception("duplicate key '%s'" % ikey)
            lookup[ikey] = interval
    
    return lookup

def joinIntervalsSum(myIntervals,start='start',end='end',score='readcount',sampleName=".",offset=0):
    """This will return a list of non-overlapping intervals and sum their scores (score)"""
    
    if not myIntervals: return myIntervals
    non_overlapping = []
    sep = {'+':[],'-':[]}
    
    print "Splitting intervals by strand"
    for i in myIntervals:
        sep[i.strand].append(i)
    
    print "Joining intervals..."
    for strand in sep.keys():
        print strand
        intervals = sep[strand]
        intervals.sort()
        
        
        current = copy.copy(intervals[0])
        for x in intervals[1:]:
            next = copy.copy(x)
            if current.intersects(next, start=start, end=end,offset=offset):
                current.end = max(current.end,next.end)
                current.__dict__[score] = current.__dict__[score]+next.__dict__[score]
            else:
                current.name = sampleName
                non_overlapping.append(current)
                current = copy.copy(next)
        current.name=sampleName
        non_overlapping.append(current)
    print "Sorting intervals"
    non_overlapping.sort()
    print "Done"
    return non_overlapping

def intervals2wig(iter,sampleName="",outDir=os.getcwd(),scratchDir=os.getcwd()):
    """
    Slow method to take iterator over intervals and convert to a .wig track.
    It's messy, but it should work...
    """
    seqs = {}
    count = 0
    print "Preparing Dictionary of alignments\nEach '.' is 10000 alignments"
    for interval in iter:
        count = count+1
        if count % 10000 == 0:
            sys.stdout.write(".")
        if count % 100000 == 0:
            print "\n%d" % (count)
        if not seqs.has_key(interval.chr):
            seqs[interval.chr]={'+':scratchDir+"/"+GenRandom(),'-':scratchDir+"/"+GenRandom()}
        FILE = open(seqs[interval.chr][interval.strand],'a')
        for i in range(interval.start,len(interval)+1):
            print >>FILE, "%d\t%d" % (i,interval.readcount)
    print "Done preparing dictionary, Begin sort and write"
    chrKeys = seqs.keys()
    chrKeys.sort()
    for chr in chrKeys:
        print "Printing " + chr
        strands = seqs[chr].keys()
        for strand in strands:
            INPUT = open(seqs[chr][strand],'r')
            filename = outDir + "/%s_%s_%s.wig" % (sampleName,chr,strand)
            OUTPUT = open(filename,'w')
            OUTPUT.write("track type=wiggle_0 name='%s_%s_%s' description='Wiggle Track for read alignment of %s sample to %s'\n" % (sampleName,chr,strand,sampleName,chr))
            print strand
            positions = {}
            while True:
                line = INPUT.readline()
                if not line: break
                pos,obs = line.split("\t")
                pos,obs = int(pos),int(obs)
                try: positions[pos]=positions[pos]+obs
                except KeyError: positions[pos]=obs
            posKeys = positions.keys()
            posKeys.sort()
            for pos in posKeys:
                wigLine = "%s\t%d\t%d\t%d" % (chr,int(pos),int(pos)+1,positions[pos])
                print >>OUTPUT, wigLine
            os.remove(seqs[chr][strand])
    return

##################
#Interval utilities and parsers
##################
bed_fields = ['chr','start','end','label','score','strand']

def parseBed(fname):
    """
    Generator that returns an iterator over spliced or unspliced BED entries.
    Iterates as Interval or SplicedInterval objects.
    """
    
    handle=open(fname,'r')
    for line in handle:
        if line.startswith("#"):
            continue
        if line.startswith("track") or line.startswith("browser"):
            continue
        vals=line.rstrip().split("\t")
        chr = vals[0]
        start = int(vals[1])
        end = int(vals[2])
        if len(vals)>=3:
            strand = vals[5]
            score = float(vals[4])
            name = vals[3]
        res = Interval(chr,start,end)
        if len(vals)>3:
            res.strand = strand
            res.score = score
            res.name = name
            res = Interval(chr,start,end,strand=strand,score=score,name=name)
        if len(vals)>6:
            res = SplicedInterval(res.chr,res.start,res.end,res.strand,score=res.score,name=res.name,exonLengths=vals[10],exonOffsets=vals[11])
        #res=dict(zip(bed_fields,vals))
        #res['start'],res['end'],res['score'] = int(res['start']),int(res['end']),int(res['score'])
        yield res

def parseBedSNP(fname):
    """
    Generator that returns an iterator over spliced or unspliced BED entries.
    Iterates as Interval or SplicedInterval objects.
    """
    
    handle=open(fname,'r')
    for line in handle:
        if line.startswith("#"):
            continue
        if line.startswith("track") or line.startswith("browser"):
            continue
        vals=line.rstrip().split("\t")
        chr = vals[0]
        start = int(vals[1])
        end = int(vals[2])
        if len(vals)>=3:
            strand = vals[5]
            score = float(vals[4])
            name = vals[3]
        res = Interval(chr,start,end)
        if len(vals)>3:
            res.strand = strand
            res.score = score
            res.name = name
            res = Interval(chr,start,end,strand=strand,score=score,name=name)
        if len(vals)>6:
            res = SNP(res.chr,res.start,res.end,res.strand,score=res.score,name=res.name,exonLengths=vals[10],exonOffsets=vals[11])
        #res=dict(zip(bed_fields,vals))
        #res['start'],res['end'],res['score'] = int(res['start']),int(res['end']),int(res['score'])
        yield res

def preprocessBed(fname):
    """
    Returns a list of Intervals from a bed file.
    """
    res = {}
    iter = parseBed(fname)
    for i in iter:
        res.setdefault(i.chr,[])
        res[i.chr].append(i)
    for k in res.keys():
        res[k].sort()
    return res

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
        if line[0] <>">":
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


def makeTSSMap(TSSBedfile,compareBedFile,flankSize=1000):
    """
    Makes a 'sense' and 'antisense' map of compareBedFile start positions within the flanking interval defined by flankSize around the start positions of TSSBedFile.
    Only increments when there is a start, does not add expression value (score).
    """
    compareDict = preprocessBed(compareBedFile)
    sys.stderr.write("Processing file: %s\n" ) % (compareBedFile)
    sense = np.zeros(2*flankSize+1)
    antisense = np.zeros(2*flankSize+1)
    
    iter = parseBed(TSSBedfile)
    sys.stderr.write("Iterating over TSSs from %s\n") % TSSBedfile
    count = 0
    for i in iter:
        if count % 100 == 0: sys.stderr.write("%d\n") % count
        count +=1
        for j in compareDict[i.chr]:
            myDist = i.distanceBetweenTSS(j)
            if abs(myDist)<=flankSize:
                if i.strand == j.strand:
                    sense[myDist+flankSize]+=1
                elif i.strand != j.strand:
                    antisense[myDist+flankSize]+=1
    return sense,antisense
        
def fetchRefSeqDict(RefSeqBed="/fg/compbio-t/lgoff/magda/references/human/transcriptome/hg18/hg18_RefSeq.bed"):
    """
    Returns a dictionary of RefSeq intervals using default hg18 RefSeq file...
    ***Depricated*** Use dbConn.fetchRefSeq instead to grab refSeq genes from UCSC
    """
    res = {}
    iter = parseBed(RefSeqBed)
    for i in iter:
        res[i.name] = i
    return res

def fetchRefSeqByChrom(RefSeqBed="/fg/compbio-t/lgoff/magda/references/human/transcriptome/hg18/hg18_RefSeq.bed"):
    """
    Returns a dictionary of RefSeq Intervals in sub-dictionaries by chromosome and strand
    ***Depricated*** Use dbConn.fetchRefSeq instead to grab refSeq genes from UCSC
    """
    res = {}
    iter = parseBed(RefSeqBed)
    for i in iter:
        res.setdefault(i.chr,{})
        res[i.chr].setdefault(i.strand,[])
        res[i.chr][i.strand].append(i)
    return res

def makeTSSBed(fname,outFname):
    iter = parseBed(fname)
    outHandle = open(outFname,'w')
    for i in iter:
        myInterval = copy.copy(i)
        if myInterval.strand == "+":
            myInterval.end = myInterval.start
        elif myInterval.strand == "-":
            myInterval.start = myInterval.end
        print >>outHandle, myInterval.toBed()        

def parseGalaxyCons(fname):
    """Parses bed-like output of conservation fetch from Galaxy webserver"""
    handle=open(fname,'r')
    for line in handle:
        if line.startswith("#"):
            continue
        if line.startswith("track") or line.startswith("browser"):
            continue
        vals=line.rstrip().split("\t")
        chr = vals[0]
        start = int(vals[1])
        end = int(vals[2])
        strand = vals[5]
        #Field[6] contains the average phastCons score
        score = float(vals[6])
        name = vals[3]
        res = Interval(chr,start,end,strand=strand,score=score,name=name)
        #res=dict(zip(bed_fields,vals))
        #res['start'],res['end'],res['score'] = int(res['start']),int(res['end']),int(res['score'])
        yield res

def findNearest(myInterval,IntervalList):
    """It would be nice to write some sort of binary search for Intervals"""
    
    myDist = 9999999999999999999
    res = 0
    for i in IntervalList:
        distance = myInterval.distance(i)
        if distance > 0 and distance < myDist:
            myDist = distance
            res = i
    return res 

def GenRandom(length = 10, chars=string.letters+string.digits):
    """
    Generates random string (by default, length=10)
    """
    return ''.join([random.choice(chars) for i in range(length)])