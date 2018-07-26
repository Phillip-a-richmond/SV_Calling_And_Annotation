import sys
import re


#This files takes the output from annovar SV annotation and intersection with the inhouse SV database
#and adds the fraction of the inhouse SV that is overlapped by the sample SV.
#It also adds a column containing the allele frequency of the inhouse SV with the greatest fraction of overlap with the sample SV to facilitate filtering. 


#usage: python annovar_reciprocaloverlap_1000G.py metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.Inhouse.bed metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.Inhouse.overlap.bed



infile = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')


for line in infile:
        cols=line.strip('\n').split('\t')
        chr=cols[0]
        start=int(cols[1])
        end=int(cols[2])
        length=end-start
        type=cols[44]
        print type
        InHousefreq=cols[11].split(",")
        InHousecoord=cols[12].split(",")
        recoverlaplist=[]
        AClist=[]
        for record, AC in zip(InHousecoord, InHousefreq):
                if record != '.':
                        InHouserecord=record.strip('"').split(":")
                        InHousestart=int(InHouserecord[1])
                        InHouseend=int(InHouserecord[2])
                        if type == '<INS>':
                                InHouselen = 1
                        else:
                                InHouselen=float(InHouseend-InHousestart)
                        startMax=max(start,InHousestart)
                        endMin=min(end,InHouseend)
                        recoverlap=str((endMin-startMax)/InHouselen)
                        recoverlaplist.append(recoverlap)
                        if recoverlap >= 0.9:
                                AC=AC.split(":")
                                AC=AC[1]
                                AClist.append(AC)
                else:
                        recoverlap="."
        recoverlaplist=",".join(recoverlaplist)
        if len(AClist) == 0:
                ACmax = "."
        else:
                ACmax=max(AClist)
        colsnew=cols[0:3]
        length=str(length)
        colsnew.append(length)
        colsnew=colsnew + cols[3:12]
        colsnew.append(recoverlaplist)
        colsnew.append(ACmax)
        colsnew= colsnew + cols[12:]
        newline="\t".join(colsnew)
        outfile.write("%s\n"%newline)





	









	
