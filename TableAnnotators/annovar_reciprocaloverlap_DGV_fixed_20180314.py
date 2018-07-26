import sys
import re

#This files takes the output from annovar SV annotation (e.g. <sample>.metaSV.annovar.sorted.hg19_multianno.denovo.withsummaryHPO.bed, and calculates the overlap fraction between the sample SV, and any overlapping DGV SV
#This is necessary because while annovar annovar can specify the fraction of overlap of the sample SV (in our case, 70%), it does not consider the fraction of overlap with respect to the database SV. 

#usage: python annovar_reciprocaloverlap_DGV.py <sample>.metaSV.annovar.sorted.hg19_multianno.<inheritance>.withsummaryHPO.bed <sample>.metaSV.annovar.sorted.hg19_multianno.<inheritance>.withsummaryHPO.DGVoverlap.bed

infile = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')


for line in infile:
	cols=line.strip('\n').split('\t')
	chr=cols[0]
	start=int(cols[1])
	end=int(cols[2])
	len=end-start
	DGVfreq=cols[17]
	DGVcoord=cols[18].split(",")
	recoverlaplist=[]
	for record in DGVcoord:
		if record != '.':
			DGVrecord=record.strip('"').split(":")
			DGVstart=int(DGVrecord[1])
			DGVend=int(DGVrecord[2])
			DGVlen=float(DGVend-DGVstart)
                        startMax=max(start,DGVstart)
                        endMin=min(end,DGVend)
                        recoverlap=str((endMin-startMax)/DGVlen)
                        recoverlaplist.append(recoverlap)
		else:	
			recoverlap="."
	recoverlaplist=",".join(recoverlaplist)
	colsnew=cols[0:18]
	colsnew.append(recoverlaplist)
	colsnew= colsnew + cols[18:]
	newline="\t".join(colsnew)
	outfile.write("%s\n"%newline)
