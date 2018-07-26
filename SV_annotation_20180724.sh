#!/bin/bash
#PBS -N sample_annovar_OMIMall
#PBS -V
#PBS -o /workdir/sample_annovar.o
#PBS -e /workdir/sample_annovar.e
#PBS -m bea
#PBS -M mcouse@bcchr.ca
#PBS -l walltime=200:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=3gb


SAMPLE='sample'

WORKDIR='workdir'
SAMPLESV=$WORKDIR/${SAMPLE}Out/variants.vcf.gz
SAMPLEERDS=$WORKDIR/
#Get variants in proband
zgrep DEL $SAMPLESV > $WORKDIR/${SAMPLE}.DEL.preformat
zgrep DUP $SAMPLESV > $WORKDIR/${SAMPLE}.DUP.preformat
zgrep INV $SAMPLESV > $WORKDIR/${SAMPLE}.INV.preformat
zgrep INS $SAMPLESV > $WORKDIR/${SAMPLE}.INS.preformat

#Convert to bed using annovar
perl /mnt/causes-data01/data/Databases/annovar/convert2annovar.pl $WORKDIR/${SAMPLE}.DEL.preformat \
	-format vcf4 \
	-outfile  $WORKDIR/${SAMPLE}.DEL \
	-includeinfo 

perl /mnt/causes-data01/data/Databases/annovar/convert2annovar.pl $WORKDIR/${SAMPLE}.DUP.preformat \
        -format vcf4 \
        -outfile  $WORKDIR/${SAMPLE}.DUP \
        -includeinfo

perl /mnt/causes-data01/data/Databases/annovar/convert2annovar.pl $WORKDIR/${SAMPLE}.INV.preformat \
        -format vcf4 \
        -outfile  $WORKDIR/${SAMPLE}.INV \
        -includeinfo

perl /mnt/causes-data01/data/Databases/annovar/convert2annovar.pl $WORKDIR/${SAMPLE}.INS.preformat \
        -format vcf4 \
        -outfile  $WORKDIR/${SAMPLE}.INS \
        -includeinfo

# rm  $WORKDIR/${SAMPLE}*preformat

#Add 1 bp to INS to faciliate comparison to inhouse INS

python /mnt/causes-data01/new/COUSE/SVscripts/add_1bp_to_INS.py $WORKDIR/${SAMPLE}.INS $WORKDIR/${SAMPLE}.INS.plus1

mv $WORKDIR/${SAMPLE}.INS.plus1 $WORKDIR/${SAMPLE}.INS


#Annotate DEL with RefSeq via Annovar
perl /mnt/causes-data01/data/Databases/annovar/table_annovar.pl $WORKDIR/${SAMPLE}.DEL /mnt/causes-data01/data/Databases/annovar/humandb \
-buildver hg19 -out $WORKDIR/${SAMPLE}.metaSV.DEL.annovar -nastring . -remove -otherinfo -protocol refGene,genomicSuperDups,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed \
-bedfile InHouseDB_SV_50_20180507_DEL_name.bed,\
InHouseDB_SV_50_20180507_DEL_name.bed,\
OMIM_phenotypes_with_coordinate_20180312.txt,\
ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype_coordinates.txt,\
RefSeqGene_Summaries_20171012_coordinates.txt,\
GeneImprintMetaImprintRefSeqCombinedNoIndexHGNCCoordinates.txt,\
DGV.GS.March2016.50percent.LossSep.Final.hg19.freq.appended.bed,\
DGV.GS.March2016.50percent.LossSep.Final.hg19.freq.appended.bed,\
1000Gphase3DELfreqadded.bed,\
1000Gphase3DELfreqadded.bed,\
DECIPHER_population_DEL_freqadded.bed,\
DECIPHER_population_DEL_freqadded.bed,\
hg19_HI_Predictions_Version3.bed,\
exac-final-cnv.gene.scores071316.bed,\
exac-final-cnv.gene.scores071316.bed,\
exac-final-cnv.gene.scores071316.bed,\
exac-final-cnv.gene.scores071316.bed,\
iscaCuratedBenignLoss.bed,\
iscaCuratedBenignLoss.bed,\
iscaCuratedPathogenicLoss.bed,\
iscaCuratedPathogenicLoss.bed,\
RLCRs_no_Repeat_Masker.txt,\
GeneHancer_hg19.sorted.bed,\
total.combined.domain.named \
-operation g,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r \
-arg ',-colsWanted 5,-colsWanted 9 -minqueryfrac 0.9,-colsWanted 8 -minqueryfrac 0.9,-colsWanted 5,-colsWanted 5,-colsWanted 5,-colsWanted 5,\
-colsWanted 5 -minqueryfrac 0.5,\
-colsWanted 6 -minqueryfrac 0.5,\
-colsWanted 10 -minqueryfrac 0.5,\
-colsWanted 11 -minqueryfrac 0.5,\
-colsWanted 14 -minqueryfrac 0.5,\
-colsWanted 15 -minqueryfrac 0.5,\
-colsWanted 4,\
-colsWanted 19,\
-colsWanted 20,\
-colsWanted 21,\
-colsWanted 22,\
-colsWanted 11 -minqueryfrac 0.5,\
-colsWanted 12 -minqueryfrac 0.5,\
-colsWanted 11 -minqueryfrac 0.5,\
-colsWanted 12 -minqueryfrac 0.5,\
-colsWanted 4 -minqueryfrac 0.9,\
-colsWanted 5,\
-colsWanted 4'


#Annotate DUP with RefSeq via Annovar
perl /mnt/causes-data01/data/Databases/annovar/table_annovar.pl $WORKDIR/${SAMPLE}.DUP /mnt/causes-data01/data/Databases/annovar/humandb \
-buildver hg19 -out $WORKDIR/${SAMPLE}.metaSV.DUP.annovar -nastring . -remove -otherinfo -protocol refGene,genomicSuperDups,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed \
-bedfile InHouseDB_SV_50_20180507_DUP_name.bed,\
InHouseDB_SV_50_20180507_DUP_name.bed,\
OMIM_phenotypes_with_coordinate_20180312.txt,\
ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype_coordinates.txt,\
RefSeqGene_Summaries_20171012_coordinates.txt,\
GeneImprintMetaImprintRefSeqCombinedNoIndexHGNCCoordinates.txt,\
DGV.GS.March2016.50percent.GainSep.Final.hg19.freq.appended.bed,\
DGV.GS.March2016.50percent.GainSep.Final.hg19.freq.appended.bed,\
1000Gphase3DUPfreqadded.bed,\
1000Gphase3DUPfreqadded.bed,\
DECIPHER_population_DUP_freqadded.bed,\
DECIPHER_population_DUP_freqadded.bed,\
hg19_HI_Predictions_Version3.bed,\
exac-final-cnv.gene.scores071316.bed,\
exac-final-cnv.gene.scores071316.bed,\
exac-final-cnv.gene.scores071316.bed,\
exac-final-cnv.gene.scores071316.bed,\
iscaCuratedBenignGain.bed,\
iscaCuratedBenignGain.bed,\
iscaCuratedPathogenicGain.bed,\
iscaCuratedPathogenicGain.bed,\
RLCRs_no_Repeat_Masker.txt,\
GeneHancer_hg19.sorted.bed,\
total.combined.domain.named \
-operation g,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r \
-arg ',-colsWanted 5,-colsWanted 9 -minqueryfrac 0.9,-colsWanted 8 -minqueryfrac 0.9,-colsWanted 5,-colsWanted 5,-colsWanted 5,-colsWanted 5,\
-colsWanted 5 -minqueryfrac 0.5,\
-colsWanted 6 -minqueryfrac 0.5,\
-colsWanted 10 -minqueryfrac 0.5,\
-colsWanted 11 -minqueryfrac 0.5,\
-colsWanted 14 -minqueryfrac 0.5,\
-colsWanted 15 -minqueryfrac 0.5,\
-colsWanted 4,\
-colsWanted 19,\
-colsWanted 20,\
-colsWanted 21,\
-colsWanted 22,\
-colsWanted 11 -minqueryfrac 0.5,\
-colsWanted 12 -minqueryfrac 0.5,\
-colsWanted 11 -minqueryfrac 0.5,\
-colsWanted 12 -minqueryfrac 0.5,\
-colsWanted 4 -minqueryfrac 0.9,\
-colsWanted 5,\
-colsWanted 4'



#Annotate INV with RefSeq via Annovar
perl /mnt/causes-data01/data/Databases/annovar/table_annovar.pl $WORKDIR/${SAMPLE}.INV /mnt/causes-data01/data/Databases/annovar/humandb \
-buildver hg19 -out $WORKDIR/${SAMPLE}.metaSV.INV.annovar -nastring . -remove -otherinfo -protocol refGene,genomicSuperDups,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed \
-bedfile InHouseDB_SV_50_20180507_INV_name.bed,\
InHouseDB_SV_50_20180507_INV_name.bed,\
OMIM_phenotypes_with_coordinate_20180312.txt,\
ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype_coordinates.txt,\
RefSeqGene_Summaries_20171012_coordinates.txt,\
GeneImprintMetaImprintRefSeqCombinedNoIndexHGNCCoordinates.txt,dummy.bed,\
dummy.bed,\
1000Gphase3INVfreqadded.bed,\
1000Gphase3INVfreqadded.bed,\
dummy.bed,\
dummy.bed,\
hg19_HI_Predictions_Version3.bed,\
exac-final-cnv.gene.scores071316.bed,\
exac-final-cnv.gene.scores071316.bed,\
exac-final-cnv.gene.scores071316.bed,\
exac-final-cnv.gene.scores071316.bed,\
dummy.bed,\
dummy.bed,\
dummy.bed,\
dummy.bed,\
RLCRs_no_Repeat_Masker.txt,\
GeneHancer_hg19.sorted.bed,\
total.combined.domain.named \
-operation g,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r \
-arg ',-colsWanted 5,-colsWanted 9 -minqueryfrac 0.9,-colsWanted 8 -minqueryfrac 0.9,-colsWanted 5,-colsWanted 5,-colsWanted 5,-colsWanted 5,,,\
-colsWanted 10 -minqueryfrac 0.5,\
-colsWanted 11 -minqueryfrac 0.5,\
-colsWanted 4 -minqueryfrac 0.5,\
-colsWanted 4 -minqueryfrac 0.5,\
-colsWanted 4,\
-colsWanted 19,\
-colsWanted 20,\
-colsWanted 21,\
-colsWanted 22,\
,,,,\
-colsWanted 4 -minqueryfrac 0.9,\
-colsWanted 5,\
-colsWanted 4'




#annotate INS with RefSeq via Annovar
perl /mnt/causes-data01/data/Databases/annovar/table_annovar.pl $WORKDIR/${SAMPLE}.INS /mnt/causes-data01/data/Databases/annovar/humandb \
-buildver hg19 -out $WORKDIR/${SAMPLE}.metaSV.INS.annovar -nastring . -remove -otherinfo -protocol refGene,genomicSuperDups,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed \
-bedfile InHouseDB_SV_50_20180507_INS_name.bed,\
InHouseDB_SV_50_20180507_INS_name.bed,\
OMIM_phenotypes_with_coordinate_20180312.txt,\
ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype_coordinates.txt,\
RefSeqGene_Summaries_20171012_coordinates.txt,\
GeneImprintMetaImprintRefSeqCombinedNoIndexHGNCCoordinates.txt,\
dummy.bed,\
dummy.bed,\
1000Gphase3INSfreqadded.bed,\
1000Gphase3INSfreqadded.bed,\
dummy.bed,\
dummy.bed,\
hg19_HI_Predictions_Version3.bed,\
exac-final-cnv.gene.scores071316.bed,\
exac-final-cnv.gene.scores071316.bed,\
exac-final-cnv.gene.scores071316.bed,\
exac-final-cnv.gene.scores071316.bed,\
dummy.bed,\
dummy.bed,\
dummy.bed,\
dummy.bed,\
RLCRs_no_Repeat_Masker.txt,\
GeneHancer_hg19.sorted.bed,\
total.combined.domain.named \
-operation g,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r,r \
-arg ',-colsWanted 5,-colsWanted 9,-colsWanted 8,-colsWanted 5,-colsWanted 5,-colsWanted 5,-colsWanted 5,,,\
-colsWanted 10 -minqueryfrac 0.5,\
-colsWanted 11 -minqueryfrac 0.5,\
-colsWanted 4 -minqueryfrac 0.5,\
-colsWanted 4 -minqueryfrac 0.5,\
-colsWanted 4,\
-colsWanted 19,\
-colsWanted 20,\
-colsWanted 21,\
-colsWanted 22,\
,,,,\
-colsWanted 4 -minqueryfrac 0.9,\
-colsWanted 5,\
-colsWanted 4'


#remove header
tail -n +2 $WORKDIR/${SAMPLE}.metaSV.DEL.annovar.hg19_multianno.txt > $WORKDIR/${SAMPLE}.metaSV.DEL.annovar.hg19_multianno.bed
tail -n +2 $WORKDIR/${SAMPLE}.metaSV.DUP.annovar.hg19_multianno.txt > $WORKDIR/${SAMPLE}.metaSV.DUP.annovar.hg19_multianno.bed
tail -n +2 $WORKDIR/${SAMPLE}.metaSV.INV.annovar.hg19_multianno.txt > $WORKDIR/${SAMPLE}.metaSV.INV.annovar.hg19_multianno.bed
tail -n +2 $WORKDIR/${SAMPLE}.metaSV.INS.annovar.hg19_multianno.txt > $WORKDIR/${SAMPLE}.metaSV.INS.annovar.hg19_multianno.bed


cat $WORKDIR/${SAMPLE}.metaSV.DEL.annovar.hg19_multianno.bed \
$WORKDIR/${SAMPLE}.metaSV.DUP.annovar.hg19_multianno.bed \
$WORKDIR/${SAMPLE}.metaSV.INV.annovar.hg19_multianno.bed \
$WORKDIR/${SAMPLE}.metaSV.INS.annovar.hg19_multianno.bed \
> $WORKDIR/${SAMPLE}.metaSV.annovar.hg19_multianno.bed


/opt/tools/bedtools/bin/sortBed -i $WORKDIR/${SAMPLE}.metaSV.annovar.hg19_multianno.bed > $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.bed


#Add overlap calculations for SVs
python /mnt/causes-data01/new/COUSE/SVscripts/annovar_reciprocaloverlap_DGV_fixed_20180314.py $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.bed \
	$WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.overlap.bed
python /mnt/causes-data01/new/COUSE/SVscripts/annovar_reciprocaloverlap_1000G_fixed_20180314.py $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.overlap.bed \
	$WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.overlap.bed
python /mnt/causes-data01/new/COUSE/SVscripts/annovar_reciprocaloverlap_DECIPHER_fixed_20180314.py $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.overlap.bed \
        $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.overlap.bed
python /mnt/causes-data01/new/COUSE/SVscripts/annovar_reciprocaloverlap_ISCA_fixed_20180314.py $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.overlap.bed \
        $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.ben.overlap.bed
python /mnt/causes-data01/new/COUSE/SVscripts/annovar_reciprocaloverlap_ISCA_pathogenic_fixed_20180314.py $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.ben.overlap.bed \
        $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.bed

rm $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.overlap.bed $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.overlap.bed \
$WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.overlap.bed $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.ben.overlap.bed


#Remove 'NAME=' that annovar adds to annotation
awk '{ gsub(/Name=/,""); print }' $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.bed \
        > $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.bed.NoName

mv $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.bed.NoName \
        $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.bed

#Calculate fraction overlap with inhouse SV

python /mnt/causes-data01/new/COUSE/SVscripts/annovar_reciprocaloverlap_inhouse_20180507.py $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.bed \
	$WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.InHouse.overlap.bed

#Add header
HEADER="Chr\tStart\tEnd\tLength\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tSegDup\tInHouseAC\tInHouse.Overlap\tInHouse.AC.LargestOverlap\tInHouse.Coordinates\tOMIMphenotype\tHPO\tEntrezSummary\tImprint\tDGV.GS.Freq\tDGV.GS.Overlap\tDGV.GS.Coordinates\t1000g.AF\t1000g.Overlap\t1000g.Coordinates\tDECIPHER.popCNV.freq\tDECIPHER.popCNV.Overlap\tDECIPHER.popCNV.Coordinates\tHIscore\tExACdel.score\tExACdup.score\tExACcnv.intolerance.score\tExACflag\tISCABenign\tISCABenign.Overlap\tISCABenign.Coordinates\tISCAPathogenic\tISCAPathogenic.Overlap\tISCAPathogenic.Coordinates\tLCR\tGeneHancer\tTAD\tChr\tStart\tID\tRef\tAlt\tFilter\tQuality\tOtherInfo\tGT\tGenotype"
NEWHEADER=$(echo -e $HEADER)

sed  "1s/.*/$NEWHEADER/"  $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.InHouse.overlap.bed > \
$WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.Inhouse.header.bed

#Remove 'Score=' that annovar adds to annotations

awk '{ gsub(/Score=/,""); print }' $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.Inhouse.header.bed \
        > $WORKDIR/${SAMPLE}.metaSV.annotated.SV.50.70.parentgt.txt

rm $WORKDIR/${SAMPLE}.metaSV.annotated.SV.NoName.txt
rm $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.Inhouse.header.bed
rm $WORKDIR/${SAMPLE}.metaSV.*.annovar.hg19_multianno.bed 
rm $WORKDIR/${SAMPLE}.metaSV.*.annovar.hg19_multianno.txt
rm $WORKDIR/${SAMPLE}*overlap
rm $WORKDIR/${SAMPLE}*bed.*
rm $WORKDIR/${SAMPLE}*avinput
rm $WORKDIR/${SAMPLE}*.overlap
rm $WORKDIR/${SAMPLE}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.Inhouse.test.bed#
rm $WORKDIR/${SAMPLE}.metaSV.annovar.hg19_multianno.bed $WORKDIR/${SAMPLE}.metaSV.annovar.sorted*























