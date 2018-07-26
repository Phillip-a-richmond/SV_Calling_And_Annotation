import sys
import re

#Fixed the reciprocal overlap function January 04, 2018
#This files takes the output from annovar SV annotation and DGV + 1000G overlap (e.g. <sample>.metaSV.annovar.sorted.hg19_multianno.denovo.withsummaryHPO.DGV.1000Goverlap.bed, and calculates the overlap fraction between the sample SV, and any overlapping DECIPHER SV
#This is necessary because while annovar annovar can specify the fraction of overlap of the sample SV (in our case, 70%), it does not consider the fraction of overlap with respect to the database SV.

#usage: python annovar_reciprocaloverlap_DECIPHER.py <sample>.metaSV.annovar.sorted.hg19_multianno.<inheritance>.withsummaryHPO.DGV.1000Goverlap.bed <sample>.metaSV.annovar.sorted.hg19_multianno.<inheritance>.withsummaryHPO.DGV.1000G.DEC.overlap.bed

infile = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')


for line in infile:
	cols=line.strip('\n').split('\t')
	chr=cols[0]
	start=int(cols[1])
	end=int(cols[2])
	len=end-start
	DECIPHERfreq=cols[23]
	DECIPHERcoord=cols[24].split(",")
	recoverlaplist=[]
	for record in DECIPHERcoord:
		if record != '.':
			DECIPHERrecord=record.strip('"').split(":")
			DECIPHERstart=int(DECIPHERrecord[1])
			DECIPHERend=int(DECIPHERrecord[2])
			DECIPHERlen=float(DECIPHERend-DECIPHERstart)
			startMax=max(start,DECIPHERstart)
			endMin=min(end,DECIPHERend)
			recoverlap=str((endMin-startMax)/DECIPHERlen)
			recoverlaplist.append(recoverlap)
		else:	
			recoverlap="."
	recoverlaplist=",".join(recoverlaplist)
	colsnew=cols[0:24]
	colsnew.append(recoverlaplist)
	colsnew= colsnew + cols[24:]
	newline="\t".join(colsnew)
	outfile.write("%s\n"%newline)
