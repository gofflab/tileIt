#!/usr/bin/env python

###########################################
# Flavor of tileIt.py designed to generate
# MPRA oligos of a given size, centered on
# a particular SNP, + or - the particular
# risk allele
#
# usage: python SNP_tiles.py [options] snp_info.txt
#
#
###########################################

import urllib2
from xml.etree import ElementTree
import sequencelib,getopt,sys,re,intervallib
from Bio import Restriction
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from itertools import tee,izip
import copy
from utils import pp
import math

# Base class for a guide RNA
# Takes a 23mer 20 Guide + 3 PAM and processes it accordingly
class TileError(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class Tile:
	def __init__(self,sequence,seqName,startPos,prefix='',suffix='',tag=''):
		self.sequence = sequence
		self.startPos = int(startPos)
		self.start = int(startPos)
		self.seqName = seqName
		self.name = "%s:%d-%d" % (self.seqName,self.start,self.start+len(self.sequence))
		self.prefix = prefix
		self.suffix = suffix
		self.tag = tag
		#Validate guide

	def validate(self):
		self.getGC()

	def compiledPrefix(self):
		""" Check prefix for '@' indicating position to add tag'"""
		tagPos = self.prefix.find('@')
		if tagPos == -1:
			return self.prefix
		else:
			return self.prefix[:tagPos]+self.tag+self.prefix[tagPos:]

	def compiledSuffix(self):
		""" Check suffix for '@' indicating position to add tag'"""
		tagPos = self.suffix.find('@')
		if tagPos == -1:
			return self.suffix
		else:
			return self.suffix[:tagPos]+self.tag+self.suffix[tagPos+1:]

	def __repr__(self):
		return "%s:%s" % (self.name,self.oligoSequence())

	def __str__(self):
		#return "%s\t%0.2f\t%d" % (self.__repr__(),self.GC,len(self))
		return "%s" % (self.__repr__())

	def toFasta(self):
		return ">%s\n%s" % (self.name,self.sequence)

	def toBed(self):
		pass

	def GC(self):
		return sequencelib.gc_content(self.oligoSequence())

	def oligoSequence(self):
		return self.compiledPrefix()+self.sequence+self.compiledSuffix()

	def __hash__(self):
		return hash(self.oligoSequence())

	def __eq__(self,other):
		#if self.sequence.upper() == other.sequence.upper():
		if self.oligoSequence() == other.oligoSequence():
			return True
		else:
			return False

	def __len__(self):
		return len(self.oligoSequence())

	def __cmp__(self,other):
		return cmp((self.seqName, self.startPos, self.name),(other.seqName, other.startPos, other.name))

	def toFasta(self):
		return ">%s\n%s" % (self.name,self.oligoSequence())

	def tileFasta(self):
		"""Only write tile sequence to fasta"""
		return ">%s\n%s" % (self.name,self.sequence)

class SNPTile(Tile):
	def __init__(self,sequence,seqName,snpPos,alleles,snpClass,startPos=1,prefix='',suffix='',tag='',GMAF=''):
		#super(SNPTile, self).__init__(*args, **kwargs)
		Tile.__init__(self,sequence=sequence,seqName=seqName,startPos=int(startPos),prefix=prefix,suffix=suffix,tag=tag)
		self.snpPos=int(snpPos)
		self.GMAF=GMAF
		#self.snpIndexPos=len(self.sequence)/2
		self.snpClass=snpClass
		self.alleles=alleles
		self.numVars=len(self.alleles.split("/"))

	def oligoSequence(self,ambiguous=True):
		pass

	def splitTiles(self):
		"""Returns a list of Tiles with unambiguous bases."""
		pass

	def __hash__(self):
		return hash(self.sequence)

	def __cmp__(self,other):
		return cmp(hash(self.sequence),hash(other.sequence))


def usage():
	sys.stderr.write("""
TileIt: Python utility for oligo sequence generation from SNP 'rs' identifiers in input file <filename.txt>
-----------------------------------------------------------------------------
Usage:
	SNP_tiles.py [options] <filename.txt>

General Options:
	-h/--help               Print this helpful help page     [ default: NA ]
	-o/--output             Output filename                  [ default: stdout ]
	-v/--verbose            Doesn't do much yet..            [ default: False ]

Oligo Options:
	-p/--prefix             String added to each designed oligo 5' end (length is subtracted from max-tile-length)   [ default: T7 Universal ]
	-s/--suffix             String added to each designed oligo 3' end (length is subtracted from max-tile-length)   [ default: M13 reverse (-24) ]
	-l/--max-tile-length    Integer. Maximum length of designed oligos                                              [ default: 150 ]
	-w/--tile-step          Integer. Window step size for designed oligos                                           [ default: 1 ]
	-t/--tag				Logical argument whether or not to add tag to prefix or suffix position marked with '@'.	[ default: F ]
	-e/--tag-length			Integer. Length of barcode 'tag' to add at position indicated by '@'.					[ default: 10 ]
	-n/--num-tags-per-tile	Integer. Number of unique tags to make per tile											[ default: 5 ]
	-r/--restriction-sites 	String. Comma-separated list of restriction enzymes whose cutting sites are to be avoided.
""")

#######################
# Helper functions
#######################
def onlyNucleic(seq,set=['a','c','g','t','u','A','C','G','T','U','n','N','@']):
	for c in seq:
		if c not in set:
			return False
	return True

def findUnique(tiles):
	return list(set(tiles))

# def findUnique(tiles):
# 	seen = set()
# 	res = []
# 	for tile in tiles:
# 		if tile not in seen:
# 			seen.add(tile)
# 			res.append(tile)
# 	return res

def estimateAffixLength(sequence,tagLength):
	tagHits = sequencelib.mcount(sequence, '@')
	if tagHits == 0:
		return len(sequence)
	elif tagHits > 1:
		raise TileError("""You can only have one instance of 'tag' per tile""")
	elif tagHits == 1:
		return len(sequence) + tagLength - 1 #-1 is required upone removal of '@' tag


def buildTags(numTags,tagLength,sites=None):
	tmpTags = set()
	while len(tmpTags)<numTags:
		tmpTag = sequencelib.GenRandomSeq(tagLength,type="DNA")
		if sites != None:
			if hasRestrictionSites(tmpTag,sites):
				continue
		tmpTags.add(tmpTag)
	return list(tmpTags)

def hasRestrictionSites(sequence,sites):
	#Parse sites
	sites = sites.split(",")
	rb = Restriction.RestrictionBatch(sites)

	#Get Bio.Seq object
	amb = IUPACAmbiguousDNA()
	tmpSeq = Seq(sequence,amb)

	#Search for sites
	res = rb.search(tmpSeq)

	#Sum hits
	totalSites = 0
	for v in res.values():
		totalSites += len(v)

	if totalSites > 0:
		return True
	else:
		return False

def warnRestrictionSites(sequence,name,sites):
	sites = sites.split(",")
	rb = Restriction.RestrictionBatch(sites)

	#Get Bio.Seq object
	amb = IUPACAmbiguousDNA()
	tmpSeq = Seq(sequence,amb)

	#Search for sites
	res = rb.search(tmpSeq)

	#Sum hits
	totalSites = 0
	for v in res.values():
		totalSites += len(v)

	if totalSites > 0:
		print >>sys.stderr, "Warning: The following positions in '%s' will be masked from tiles due to incompatible restictions sites:" % (name)
		pp(res)
	else:
		pass


#######################
# Scan input sequence #
#######################
#TODO: Modify this so that it only gets the window of appropriate size.  We will add prefix and suffix afterwards.

def scanSequence(sequence,seqName,tileStep=1,tileSize=150):
	tiles = []

	#Pre-compute number of chunks to emit
	numOfChunks = ((len(sequence)-tileSize)/tileStep) + 1

    #Tile across sequence
	for i in xrange(0,numOfChunks*tileStep,tileStep):
		tiles.append(Tile(sequence=sequence[i:i+tileSize],seqName=seqName,startPos=i+1))
	return tiles

###################
# Make tiles from SNP flanks
###################

def makeTileFromSnp(snp,halfWidth=50):
	flankStart = int(snp.snpPos-halfWidth)
	flankEnd = int(snp.snpPos+halfWidth)
	#print >>sys.stderr, "%s\t%s" % (flankStart,flankEnd)
	subSequence = snp.sequence[flankStart:flankEnd]
	tile = SNPTile(sequence=subSequence,seqName=snp.name,alleles=snp.alleles,snpClass=snp.varClass,startPos=1,snpPos=1+int(halfWidth))
	return tile

###################
# Reporting
###################

def outputTable(tiles,outHandle=sys.stdout):
	outputKeys=["name","seqName","startPos","oligoSequence","prefix","suffix","tag","GC"]
	print >>outHandle, "\t".join(outputKeys)

	for tile in tiles:
		vals = []
		for k in outputKeys:
			v = getattr(tile,k)
			if callable(v):
				v = v()
			vals.append(str(v))
		print >>outHandle, "\t".join(vals)

# def main():
# 	#Argument handling
# 	try:
# 		opts,args = getopt.getopt(sys.argv[1:],"ho:vp:s:l:w:te:n:r:",["help","output=","verbose","prefix=","suffix=","max-tile-length=","tile-step=","tag","tag-length=","num-tags-per-tile=","restriction-sites="])
# 	except getopt.GetoptError as err:
# 		print(err)
# 		usage()
# 		sys.exit(2)
# 	output = None
# 	verbose = False
# 	for o,a in opts:
# 		if o == "-v":
# 			verbose = True
# 		elif o in ("-h","--help"):
# 			usage()
# 			sys.exit()
# 		elif o in ("-o","--output"):
# 			output = a
# 		elif o in ("-p","--prefix"):
# 			if onlyNucleic(a):
# 				prefix = a
# 			else:
# 				raise(TileError,"Prefix is not a valid nucleic acid sequence.")
# 		elif o in ("-s","--suffix"):
# 			if onlyNucleic(a):
# 				suffix = a
# 			else:
# 				raise(TileError,"Suffix is not a valid nucleic acid sequence.")
# 		elif o in ("-l","--max-tile-length"):
# 			maxTileSize = int(a)
# 		elif o in ("-w","--tile-step"):
# 			tileStep = int(a)
# 		elif o in ("-t","--tag"):
# 			addTag = True
# 		elif o in ("-e","--tag-length"):
# 			tagLength = int(a)
# 		elif o in ("-n","--num-tags-per-tile"):
# 			numTagsPerTile = int(a)
# 		elif o in ("-r","--restriction-sites"):
# 			sites = a.strip()
# 			if sites == "":
# 				sites = None
# 		else:
# 			assert False, "Unhandled option"

# 	# Grab fname as remainder argument
# 	try:
# 		fname = str(args[0])
# 		handle = open(fname,'r')
# 	except:
# 		usage()
# 		sys.exit(2)
# 		#print fname

# 	#Estimate prefix and suffix length
# 	prefixLength = estimateAffixLength(prefix,tagLength)
# 	suffixLength = estimateAffixLength(suffix,tagLength)

# 	#Find window size
# 	tileSize = maxTileSize-prefixLength-suffixLength
# 	halfWidth = floor(tileSize/2)

# 	#Read in input file

# 	#Fetch all tile sequences
# 	fastaIter = sequencelib.FastaIterator(handle)
# 	tiles = []
# 	for mySeq in fastaIter:
# 		#Warn about masked regions
# 		if sites != None:
# 			warnRestrictionSites(mySeq['sequence'],mySeq['name'],sites)

# 		#Get tiles from sequence
# 		tmpTiles = scanSequence(mySeq['sequence'],mySeq['name'],tileStep=tileStep,tileSize=tileSize)
# 		tiles += tmpTiles

# 	# Remove duplicate tile sequences
# 	tiles = findUnique(tiles)

# 	# Check tiles for restriction sites
# 	if sites != None:
# 		cleanTiles = set()
# 		for tile in tiles:
# 			if hasRestrictionSites(tile.sequence, sites):
# 				continue
# 			else:
# 				cleanTiles.add(tile)
# 		tiles = list(cleanTiles)

# 	#determine number of tags needed
# 	numTagsReq = len(tiles) * numTagsPerTile

# 	tags = buildTags(numTagsReq,tagLength,sites=sites)

# 	assert len(tags) == len(tiles) * numTagsPerTile

# 	#Create numTagsPerTile tiles for each sequence
# 	tmpTiles = set()
# 	#Add prefix, suffix, and tag
# 	for i in xrange(len(tiles)):
# 		for j in xrange(numTagsPerTile):
# 			tmpTile = copy.copy(tiles[i])
# 			tmpTile.name = "%s:%d" % (tmpTile.name,j)
# 			tmpTile.prefix = prefix
# 			tmpTile.suffix = suffix
# 			tmpTile.tag = tags.pop()
# 			tmpTiles.add(tmpTile)

# 	tiles = list(tmpTiles)
# 	tiles.sort()

# 	#Just for QC
# 	#for i in xrange(10):
# 	outputTable(tiles)

# 	print >>sys.stderr, "There are a total of %d unique tiles" % len(tiles)

def test():
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
		'VR': 'attaccgcctttgagtgagc'
	}

	MPRATags = {
		'prefix': 'ACTGGCCGCTTCACTG',
		'suffix': 'GGTACCTCTAGA@AGATCGGAAGAGCGTCG'
	}

	maxTileSize = 150
	tileStep = 1
	tagLength = 10
	numTagsPerTile = 5
	prefix = MPRATags['prefix']
	suffix = MPRATags['suffix']
	sites = "KpnI,XbaI,SfiI"


	#Argument handling
	try:
		opts,args = getopt.getopt(sys.argv[1:],"ho:vp:s:l:w:te:n:r:",["help","output=","verbose","prefix=","suffix=","max-tile-length=","tile-step=","tag","tag-length=","num-tags-per-tile=","restriction-sites="])
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
		elif o in ("-t","--tag"):
			addTag = True
		elif o in ("-e","--tag-length"):
			tagLength = int(a)
		elif o in ("-n","--num-tags-per-tile"):
			numTagsPerTile = int(a)
		elif o in ("-r","--restriction-sites"):
			sites = a.strip()
			if sites == "":
				sites = None
		else:
			assert False, "Unhandled option"
	# Grab fname as remainder argument
	try:
		fname = str(args[0])
		#handle = open(fname,'r')
	except:
		usage()
		sys.exit(2)
		#print fname

	#Estimate prefix and suffix length
	prefixLength = estimateAffixLength(prefix,tagLength)
	suffixLength = estimateAffixLength(suffix,tagLength)

	#Find window size
	tileSize = maxTileSize-prefixLength-suffixLength
	halfWidth = math.floor(tileSize/2)

	# Read in input file (list of 'rs' ids)
	rsIds = set()
	inputHandle = open(fname,'r')
	for line in inputHandle:
		rsIds.add(line.rstrip())
	rsIds = list(rsIds)

	# Get dbSNP objects from rsIds
	snps = []
	snpCount = 0
	for rsid in rsIds[:50]:
		snpCount += 1
		if snpCount%10 == 0:
			print >>sys.stderr, "."
		snps.append(intervallib.dbSNP(name=rsid))

	#Test that allele position is correct
	for snp in snps:
		print >>sys.stderr, "%s\t%s" % (snp.sequence[snp.snpPos],snp.alleles)


	SNPtiles = []
	for snp in snps:
		snpTile = makeTileFromSnp(snp,halfWidth=halfWidth)
		SNPtiles.append(snpTile)

	#Just some debug reporting. can be removed.
	for t in SNPtiles:
		print >>sys.stderr, "%s\n\t%s\n\t%s:%d:%s:%d" %(t,t.sequence,t.alleles,t.snpPos,t.sequence[t.snpPos-1],t.numVars)

	#Warn if any restriction sites might complicate matters
	for snpTile in SNPtiles:
		if sites != None:
			warnRestrictionSites(snpTile.sequence,snpTile.name,sites)

	# Check tiles for restriction sites and remove those that have matches
	if sites != None:
		cleanTiles = set()
		for tile in SNPtiles:
			if hasRestrictionSites(tile.sequence, sites):
				continue
			else:
				cleanTiles.add(tile)
		SNPtiles = list(cleanTiles)

	#determine number of tags needed
	tileCount = 0
	for tile in SNPtiles:
		tileCount += tile.numVars
	numTagsReq = tileCount * numTagsPerTile

	print >>sys.stderr, "%d total tags requested" % numTagsReq

	# tags = buildTags(numTagsReq,tagLength,sites=sites)

	# assert len(tags) == len(tiles) * numTagsPerTile

	# #Create numTagsPerTile tiles for each sequence
	# tmpTiles = set()

	# #Add prefix, suffix, and tag
	# for i in xrange(len(tiles)):
	# 	for j in xrange(numTagsPerTile):
	# 		tmpTile = copy.copy(tiles[i])
	# 		tmpTile.name = "%s:%d" % (tmpTile.name,j)
	# 		tmpTile.prefix = prefix
	# 		tmpTile.suffix = suffix
	# 		tmpTile.tag = tags.pop()
	# 		tmpTiles.add(tmpTile)

	# tiles = list(tmpTiles)
	# tiles.sort()

	# #Just for QC
	# #for i in xrange(10):
	# outputTable(tiles)


if __name__ == "__main__":
	test()
