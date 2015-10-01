############
# Barcode Utils
##############

from tileIt import Tile


# Utility functions
def buildTags(numTags,tagLength,sites=None):
	tmpTags = set()
	while len(tmpTags)<numTags:
		tmpTag = sequencelib.GenRandomSeq(tagLength,type="DNA")
		if tmpTag[:2] == "TC": 	# This is specific to XbaI sites...the first two bases after XbaI cannot be 'TC' because the 'GAtc' motif is a target for dam methylation
				continue 		# Going to keep this in here because it shouldn't mess up things to badly anyways.
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
		#pp(res)
	else:
		pass

###################
# Reporting
###################

def outputTable(tiles,outHandle=sys.stdout):
	outputKeys=["name","seqName","startPos","oligoSequenceWithBarcode"]
	print >>outHandle, "\t".join(outputKeys)

	for tile in tiles:
		vals = []
		for k in outputKeys:
			v = getattr(tile,k)
			if callable(v):
				v = v()
			vals.append(str(v))
		print >>outHandle, "\t".join(vals)

################
# Main
################
def main():
	# Variables
	numTagsPerTile = 5
	tagLength = 7
	badSites = "KpnI,XbaI,SfiI"


	# Import merged tile list (from SNP_tile output and tileIt output)

	# Create Tiles for each line
	rawTiles = []
	for line in handle:
		rawTiles.append(Tile(sequence=,name=,startPos=1,prefix=,suffix=,))

	# Determine total number of barcodes needed
	tagsNeeded = len(rawTiles)*numTagsPerTile


	# Generate unique barcodes
	buildTags(tagsNeeded,tagLength,sites=)

	# Add tags to Tiles
	# you will need to create new Tiles (numTagsPerTile for each) and append an integer 0-numTagsPerTile to the Tile name

	# Output Tiles (merge prefix, sequence, and suffix)


if __name__ == "__main__":
	main()
