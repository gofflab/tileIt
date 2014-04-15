#!/usr/bin/env python
import sequencelib,getopt,sys,re
from itertools import tee,izip

# Base class for a guide RNA
# Takes a 23mer 20 Guide + 3 PAM and processes it accordingly
class TileError(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class TileRNA:
	def __init__(self,sequence,seqName,startPos,prefix='',suffix=''):
		self.sequence = sequence
		self.start = startPos
		self.seqName = seqName
		self.name = "%s:%d" % (self.seqName,self.start)
		self.prefix = prefix
		self.suffix = suffix
		#Validate guide
		self.getGC()

	def validate(self):
		pass

	def __repr__(self):
		return "%s:%s  %s  %s" % (self.name,self.prefix,self.sequence.upper(),self.suffix)

	def toFasta(self):
		return ">%s\n%s" % (self.name,self.sequence)
	
	def toBed(self):
		pass

	def getGC(self):
		self.GC = sequencelib.gc_content(self.sequence)

	def __hash__(self):
		return hash(self.sequence.upper())

	def __eq__(self,other):
		if self.sequence.upper() == other.sequence.upper():
			return True
		else:
			return False

	def __len__(self):
		return len(self.prefix) + len(self.sequence) + len(self.suffix)



def usage():
	sys.stderr.write("""
TileIt: simply python utility for oligo sequence generation from fasta files
-----------------------------------------------------------------------------
Usage:
	tileIt.py [options] <filename.fasta>

General Options:
	-h/--help               Print this helpful help page     [ default: NA ]
	-o/--output             Output filename                  [ default: stdout ]
	-v/--verbose            Doesn't do much yet..            [ default: False ]

Oligo Options:
	-p/--prefix             String added to each designed oligo 5' end (length is subtracted from max-tile-length)   [ default: T7 Universal ]
	-s/--suffix             String added to each designed oligo 3' end (length is subtracted from max-tile-length)   [ default: M13 reverse (-24) ]
	-l/--max-tile-length    Integer. Maximum length of designed oligos                                               [ default: 150 ]
	-w/--tile-step          Integer.  Window step size for designed oligos                                           [ default: 1 ]

""")

#######################
# Helper functions
#######################
def onlyNucleic(seq,set=['a','c','g','t','u','A','C','G','T','U','n','N']):
	for c in seq:
		if c not in set: 
			return False
	return True

# def findUnique(tiles):
#	return list(set(tiles))

def findUnique(tiles):
	seen = set()
	res = []
	for tile in tiles:
		if tile not in seen:
			seen.add(tile)
			res.append(tile)
	return res

#######################
# Scan input sequence #
#######################
def scanSequence(sequence,seqName,tileStep=1,prefix='',suffix='',maxTileSize=150):
	tileSize = maxTileSize-len(prefix)-len(suffix)
	tiles = []

	#Pre-compute number of chunks to emit
	numOfChunks = ((len(sequence)-tileSize)/tileStep) + 1
    
    #Tile across sequence
	for i in xrange(0,numOfChunks*tileStep,tileStep):
		tiles.append(TileRNA(sequence=sequence[i:i+tileSize],seqName=seqName,startPos=i+1,prefix=prefix,suffix=suffix))
	return tiles

def main():
	#######################
	# Variables
	#######################
	universalPrimers = {
		'M13 forward sequencing primer (-20)': 'GTAAAACGACGGCCAGT',
		'M13 forward sequencing primer (-40)': 'GTTTTCCCAGTCACGAC',
		'M13 forward sequencing primer (-47)': 'CGCCAGGGTTTTCCCAGTCACGAC',
		'M13 reverse sequencing primer (-24)': 'AACAGCTATGACCATG',
		'M13 reverse sequencing primer (-48)': 'AGCGGATAACAATTTCACACAGGA',
		'T3 promoter': 'ATTAACCCTCACTAAAGGGA',
		'T7 universal primer': 'TAATACGACTCACTATAGGG',
		'T7 Terminator': 'GCTAGTTATTGCTCAGCGG',
		'SP6 promoter': 'CATACGATTTAGGTGACACTATAG',
		'SP6 universal primer': 'ATTTAGGTGACACTATAG',
		'VF2': 'tgccacctgacgtctaagaa',
		'VR': 'attaccgcctttgagtgagc',
	}

	maxTileSize = 150
	tileStep = 1
	prefix = universalPrimers['T7 universal primer']
	suffix = universalPrimers['M13 reverse sequencing primer (-24)']

	try:
		opts,args = getopt.getopt(sys.argv[1:],"ho:vp:s:l:w:",["help","output=","verbose","prefix=","suffix=","max-tile-length=","tile-step="])
	except getopt.GetoptError as err:
		print(err)
		usage()
		sys.exit(2)
	output = None
	verbose = False
	for o,a in opts:
		if o == "-v":
			verbose = True
		elif o in ("-h","--help"):
			usage()
			sys.exit()
		elif o in ("-o","--output"):
			output = a
		elif o in ("-p","--prefix"):
			if onlyNucleic(a):
				prefix = a
			else:
				raise(TileError,"Prefix is not a valid nucleic acid sequence.")
		elif o in ("-s","--suffix"):
			if onlyNucleic(a):
				suffix = a
			else:
				raise(TileError,"Suffix is not a valid nucleic acid sequence.")
		elif o in ("-l","--max-tile-length"):
			maxTileSize = int(a)
		elif o in ("-w","--tile-step"):
			tileStep = int(a)
		else:
			assert False, "Unhandled option"
	# Grab fname as remainder argument
	try:
		fname = str(args[0])
		handle = open(fname,'r')
	except:
		usage()
		sys.exit(2)
		#print fname
	#Find window size
	tileSize = maxTileSize-len(prefix)-len(suffix)
	
	#Main workflow
	fastaIter = sequencelib.FastaIterator(handle)
	tiles = []
	for mySeq in fastaIter:
		tmpTiles = scanSequence(mySeq['sequence'],mySeq['name'],tileStep=tileStep,maxTileSize=maxTileSize,prefix=prefix,suffix=suffix)
		tiles += tmpTiles

	# Remove duplicate tile sequences
	tiles = findUnique(tiles)

	for tile in tiles:
		print "%s\t%0.2f\t%d" % (tile,tile.GC,len(tile))

	print len(tiles)

def test():
	fname = "test/Firre.fasta"
	handle = open(fname,'r')
	fastaIter = sequencelib.FastaIterator(handle)

	mySeq = fastaIter.next()
	tiles = scanSequence(mySeq['sequence'],mySeq['name'])

	print len(tiles)
	# Remove duplicate tile sequences
	tiles = findUnique(tiles)

	# for tile in tiles:
	# 	print "%s\t%0.2f\t%d" % (tile,tile.GC,len(tile))

	print len(tiles)

if __name__ == "__main__":
	main()