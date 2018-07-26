import sys


if len(sys.argv) < 2:
	print "Usage: python BedToAVInput.py <in.erds> <out.avinput>"
	sys.exit()

infile = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')


# The general gyst of this is as follows:
# We have:
# chrom	start	end	size	type	score	precise_boundary	ref_cn	infer_cn
# 1	13320001	13373600	53600	DEL	1469.68	0	2	1
# 1	17676185	17678109	1925	DEL	NA	1	2	1
# Where: precise_boundary is a boolean 
# reference_cn is the number of copies in the ref genome
# inferred_cn is copies in sample



# We want to match example AV input, which for simplicity I'm just going to mirror a deletion:
# 13	20797176	21105944	0	-	comments: a 342kb deletion encompassing GJB6, associated with hearing loss

# So I will produce:
# 1	13320001	13373600	0	-	comments: type=DEL|qual=1469.68|imprecise|refCopy=2|altCopy=1
# 1	13320001	13373600	53600	DEL	1469.68	0	2	1

# 1	17676185	17678109	0	-	comments: type=DEL|qual=NA|precise|refCopy=2|altCopy=1
# 1	17676185	17678109	1925	DEL	NA	1	2	1

# For Dups:
# 1	724001	727000	3000	DUP	31.88	0	2	12
# Same thing


for line in infile:
	line = line.strip('\n')
	cols = line.split('\t')
	chrom = cols[0]
	start = cols[1]
	end = cols[2]
	vartype = cols[4]
	score = cols[5]
	if cols[6]=='0':
		confidence='imprecise'
	elif cols[6]=='1':
		confidence='precise'
	refcopy = cols[7]
	altcopy = cols[8]
	
	OtherInfo = 'type=%s|qual=%s|%s|refCopy=%s|altCopy=%s'%(vartype,score,confidence,refcopy,altcopy)
	outfile.write('%s\t%s\t%s\t0\t-\tcomments: %s\n'%(chrom,start,end,OtherInfo))





