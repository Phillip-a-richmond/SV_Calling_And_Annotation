import sys
import re

#Fixed reciprocal overlap function January 04, 2018
#This files takes the output from annovar SV annotation and DGV overlap (e.g. <sample>.metaSV.annovar.sorted.hg19_multianno.denovo.withsummaryHPO.DGVoverlap.bed, and calculates the overlap fraction between the sample SV, and any overlapping 1000G SV
#This is necessary because while annovar annovar can specify the fraction of overlap of the sample SV (in our case, 70%), it does not consider the fraction of overlap with respect to the database SV.

#usage: python annovar_reciprocaloverlap_1000G.py <sample>.metaSV.annovar.sorted.hg19_multianno.<inheritance>.withsummaryHPO.DGV.overlap.bed <sample>.metaSV.annovar.sorted.hg19_multianno.<inheritance>.withsummaryHPO.DGV.1000Goverlap.bed



infile = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')


for line in infile:
	cols=line.strip('\n').split('\t')
	chr=cols[0]
	start=int(cols[1])
	end=int(cols[2])
	len=end-start
	KGfreq=cols[20]
	KGcoord=cols[21].split(",")
	recoverlaplist=[]
	for record in KGcoord:
		if record != '.':
			KGrecord=record.strip('"').split(":")
			KGstart=int(KGrecord[1])
			KGend=int(KGrecord[2])
			KGlen=float(KGend-KGstart)
                        startMax=max(start,KGstart)
                        endMin=min(end,KGend)
                        recoverlap=str((endMin-startMax)/KGlen)
                        recoverlaplist.append(recoverlap)
		else:	
			recoverlap="."
	recoverlaplist=",".join(recoverlaplist)
	colsnew=cols[0:21]
	colsnew.append(recoverlaplist)
	colsnew= colsnew + cols[21:]
	newline="\t".join(colsnew)
	outfile.write("%s\n"%newline)
