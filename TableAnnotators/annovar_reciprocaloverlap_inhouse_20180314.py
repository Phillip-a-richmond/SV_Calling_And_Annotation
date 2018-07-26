import sys
import re


#This files takes the output from annovar SV annotation and intersection with the inhouse SV database
#and adds the fraction of the inhouse SV that is overlapped by the sample SV


#usage: python annovar_reciprocaloverlap_1000G.py metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.Inhouse.bed metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.Inhouse.overlap.bed



infile = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')


for line in infile:
        cols=line.strip('\n').split('\t')
        chr=cols[0]
        start=int(cols[1])
        end=int(cols[2])
        len=end-start
        InHousefreq=cols[11].split(",")
        InHousecoord=cols[12].split(",")
        recoverlaplist=[]
        for record in InHousecoord:
                if record != '.':
                        InHouserecord=record.strip('"').split(":")
                        InHousestart=int(InHouserecord[1])
                        InHouseend=int(InHouserecord[2])
                        InHouselen=float(InHouseend-InHousestart)
                        startMax=max(start,InHousestart)
                        endMin=min(end,InHouseend)
                        recoverlap=str((endMin-startMax)/InHouselen)
                        recoverlaplist.append(recoverlap)
                else:
                        recoverlap="."
        recoverlaplist=",".join(recoverlaplist)
        colsnew=cols[0:12]
        colsnew.append(recoverlaplist)
        colsnew= colsnew + cols[12:]
        newline="\t".join(colsnew)
        outfile.write("%s\n"%newline)





	









	
