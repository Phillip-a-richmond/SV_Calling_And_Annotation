import sys
import re

#Fixed the reciprocal overlap function January 04, 2018
#This files takes the output from annovar SV annotation and DGV + 1000G + DECIPHER overlap (e.g. <sample>.metaSV.annovar.sorted.hg19_multianno.denovo.withsummaryHPO.DGV.1000G.DEC.overlap.bed, and calculates the overlap fraction between the sample SV, and any overlapping ISCA SV
#This is necessary because while annovar annovar can specify the fraction of overlap of the sample SV (in our case, 70%), it does not consider the fraction of overlap with respect to the database SV.

#usage: python annovar_reciprocaloverlap_ISCA.py <sample>.metaSV.annovar.sorted.hg19_multianno.<inheritance>.withsummaryHPO.DGV.1000G.DEC.overlap.bed <sample>.metaSV.annovar.sorted.hg19_multianno.<inheritance>.withsummaryHPO.DGV.1000G.DEC.ISCA.overlap.bed

infile = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')


for line in infile:
	cols=line.strip('\n').split('\t')
	chr=cols[0]
	start=int(cols[1])
	end=int(cols[2])
	len=end-start
	ISCAPathcoord=cols[35].split(",")
	recoverlaplist=[]
	for record in ISCAPathcoord:
		if record != '.':
			ISCArecord=record.strip('"').split(":")
			ISCAstart=int(ISCArecord[1])
			ISCAend=int(ISCArecord[2])
			ISCAlen=float(ISCAend-ISCAstart)
			startMax=max(start,ISCAstart)
			endMin=min(end,ISCAend)
			recoverlap=str((endMin-startMax)/ISCAlen)
			recoverlaplist.append(recoverlap)
		else:	
			recoverlap="."
	recoverlaplist=",".join(recoverlaplist)
	colsnew=cols[0:35]
	colsnew.append(recoverlaplist)
	colsnew= colsnew + cols[35:]
	newline="\t".join(colsnew)
	outfile.write("%s\n"%newline)
