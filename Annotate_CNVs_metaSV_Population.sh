#!/bin/bash
#PBS -N Population_STR_ExpansionHunter
#PBS -V
#PBS -o /mnt/causes-vnx1/CAUSES/CNVS/Annotation.o
#PBS -e /mnt/causes-vnx1/CAUSES/CNVS/Annotation.e
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
## Set the total memory for the job
#PBS -l mem=24gb
## Set the max walltime for the job
#PBS -l walltime=200:00:00
## Set the total number of processors for the job
#PBS -l nodes=1:ppn=6
NSLOTS=$PBS_NUM_PPN
umask 0002
source /opt/tools/hpcenv.sh





WORKING_DIR='/mnt/causes-vnx1/CAUSES/CNVS/ANNOTATED/'
GENOME_FASTA='/mnt/causes-data01/data/GENOMES/hg19/FASTA/hg19.fa'
# G249-1.metaSV 



for WORKDIR in /mnt/causes-vnx1/CAUSES/CNVS/*Out
do
	echo $WORKDIR
        IFS='/' read -a array <<< "${WORKDIR}"
        metaSVDIR=${array[-1]}
        IFS='O' read -a array2 <<< "${metaSVDIR}"
        SAMPLE_ID=${array2[0]}
        echo $SAMPLE_ID
	VCF=${WORKDIR}/variants.vcf.gz
	ls $VCF

	# Get variants into DEL and DUP sub files
	zgrep DEL $VCF > ${WORKDIR}/${SAMPLE_ID}.DEL.preformat
	zgrep DUP $VCF > ${WORKDIR}/${SAMPLE_ID}.DUP.preformat

	#Convert to bed using annovar
	perl /mnt/causes-data01/data/Databases/annovar/convert2annovar.pl $WORKDIR/${SAMPLE_ID}.DEL.preformat \
        	-format vcf4 \
        	-outfile  $WORKDIR/${SAMPLE_ID}.DEL \
        	-includeinfo

	perl /mnt/causes-data01/data/Databases/annovar/convert2annovar.pl $WORKDIR/${SAMPLE_ID}.DUP.preformat \
        	-format vcf4 \
        	-outfile  $WORKDIR/${SAMPLE_ID}.DUP \
        	-includeinfo

	
	#Annotate DEL with RefSeq via Annovar
	
	perl /mnt/causes-data01/data/Databases/annovar/table_annovar.pl $WORKDIR/${SAMPLE_ID}.DEL /mnt/causes-data01/data/Databases/annovar/humandb \
	--buildver hg19 --verbose --out $WORKDIR/${SAMPLE_ID}.metaSV.DEL.annovar --nastring . --remove --otherinfo --protocol refGene,genomicSuperDups,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed \
	--bedfile InHouseDB_SV_100_20180724_DEL_name.bed,\
InHouseDB_SV_100_20180724_DEL_name.bed,\
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


	# Annotate DUPs	

	perl /mnt/causes-data01/data/Databases/annovar/table_annovar.pl $WORKDIR/${SAMPLE_ID}.DUP /mnt/causes-data01/data/Databases/annovar/humandb \
	-buildver hg19 -out $WORKDIR/${SAMPLE_ID}.metaSV.DUP.annovar -nastring . -remove -otherinfo -protocol refGene,genomicSuperDups,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed \
	-bedfile InHouseDB_SV_100_20180724_DUP_name.bed,\
InHouseDB_SV_100_20180724_DUP_name.bed,\
OMIM_phenotypes_with_coordinate_20180312.txt,\
ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype_coordinates.txt,\
RefSeqGene_Summaries_20171012_coordinates.txt,\
GeneImprintMetaImprintRefSeqCombinedNoIndexHGNCCoordinates.txt,\
DGV.GS.March2016.50percent.LossSep.Final.hg19.freq.appended.bed,\
DGV.GS.March2016.50percent.LossSep.Final.hg19.freq.appended.bed,\
1000Gphase3DUPfreqadded.bed,\
1000Gphase3DUPfreqadded.bed,\
DECIPHER_population_DUP_freqadded.bed,\
DECIPHER_population_DUP_freqadded.bed,\
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

	echo "finished annotation"


	#remove header
	tail -n +2 $WORKDIR/${SAMPLE_ID}.metaSV.DEL.annovar.hg19_multianno.txt > $WORKDIR/${SAMPLE_ID}.metaSV.DEL.annovar.hg19_multianno.bed
	tail -n +2 $WORKDIR/${SAMPLE_ID}.metaSV.DUP.annovar.hg19_multianno.txt > $WORKDIR/${SAMPLE_ID}.metaSV.DUP.annovar.hg19_multianno.bed

	# concatenate dups and dels together
	cat $WORKDIR/${SAMPLE_ID}.metaSV.DEL.annovar.hg19_multianno.bed \
	$WORKDIR/${SAMPLE_ID}.metaSV.DUP.annovar.hg19_multianno.bed \
	> $WORKDIR/${SAMPLE_ID}.metaSV.annovar.hg19_multianno.bed

	/opt/tools/bedtools/bin/sortBed -i $WORKDIR/${SAMPLE_ID}.metaSV.annovar.hg19_multianno.bed > $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.bed

	#Add overlap calculations for SVs
	python /mnt/causes-data01/new/COUSE/SVscripts/annovar_reciprocaloverlap_DGV_fixed_20180314.py $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.bed \
	        $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.overlap.bed
	python /mnt/causes-data01/new/COUSE/SVscripts/annovar_reciprocaloverlap_1000G_fixed_20180314.py $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.overlap.bed \
	        $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.overlap.bed
	python /mnt/causes-data01/new/COUSE/SVscripts/annovar_reciprocaloverlap_DECIPHER_fixed_20180314.py $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.overlap.bed \
	        $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.overlap.bed
	python /mnt/causes-data01/new/COUSE/SVscripts/annovar_reciprocaloverlap_ISCA_fixed_20180314.py $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.overlap.bed \
	        $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.ben.overlap.bed
	python /mnt/causes-data01/new/COUSE/SVscripts/annovar_reciprocaloverlap_ISCA_pathogenic_fixed_20180314.py $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.ben.overlap.bed \
	        $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.bed
	
	rm $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.overlap.bed $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.overlap.bed \
	$WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.overlap.bed $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.ben.overlap.bed


	#Remove 'NAME=' that annovar adds to annotation
	awk '{ gsub(/Name=/,""); print }' $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.bed \
        > $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.bed.NoName

	mv $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.bed.NoName \
        $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.bed

	#Calculate fraction overlap with inhouse SV
	python /mnt/causes-data01/new/COUSE/SVscripts/annovar_reciprocaloverlap_inhouse_20180507.py $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.bed \
        $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.InHouse.overlap.bed


	#Add header
	HEADER="Chr\tStart\tEnd\tLength\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tSegDup\tInHouseAC\tInHouse.Overlap\tInHouse.AC.LargestOverlap\tInHouse.Coordinates\tOMIMphenotype\tHPO\tEntrezSummary\tImprint\tDGV.GS.Freq\tDGV.GS.Overlap\tDGV.GS.Coordinates\t1000g.AF\t1000g.Overlap\t1000g.Coordinates\tDECIPHER.popCNV.freq\tDECIPHER.popCNV.Overlap\tDECIPHER.popCNV.Coordinates\tHIscore\tExACdel.score\tExACdup.score\tExACcnv.intolerance.score\tExACflag\tISCABenign\tISCABenign.Overlap\tISCABenign.Coordinates\tISCAPathogenic\tISCAPathogenic.Overlap\tISCAPathogenic.Coordinates\tLCR\tGeneHancer\tTAD\tChr\tStart\tID\tRef\tAlt\tFilter\tQuality\tOtherInfo\tGT\tGenotype"
	NEWHEADER=$(echo -e $HEADER)

	sed  "1s/.*/$NEWHEADER/"  $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.InHouse.overlap.bed > \
	$WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.Inhouse.header.bed

	#Remove 'Score=' that annovar adds to annotations
	awk '{ gsub(/Score=/,""); print }' $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.Inhouse.header.bed \
        > $WORKDIR/${SAMPLE_ID}.metaSV.annotated.SV.50.70.parentgt.txt

	# Clean up
	rm $WORKDIR/${SAMPLE_ID}.metaSV.annotated.SV.NoName.txt
	rm $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.Inhouse.header.bed
	rm $WORKDIR/${SAMPLE_ID}.metaSV.*.annovar.hg19_multianno.bed
	rm $WORKDIR/${SAMPLE_ID}.metaSV.*.annovar.hg19_multianno.txt
	rm $WORKDIR/${SAMPLE_ID}*overlap
	rm $WORKDIR/${SAMPLE_ID}*bed.*
	rm $WORKDIR/${SAMPLE_ID}*avinput
	rm $WORKDIR/${SAMPLE_ID}*.overlap
	rm $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted.hg19_multianno.OMIMall.withsummaryHPO.OMIM.DGV.1000G.DEC.ISCA.overlap.Inhouse.test.bed#
	rm $WORKDIR/${SAMPLE_ID}.metaSV.annovar.hg19_multianno.bed $WORKDIR/${SAMPLE_ID}.metaSV.annovar.sorted*

	cp $WORKDIR/${SAMPLE_ID}.metaSV.annotated.SV.50.70.parentgt.txt /mnt/causes-vnx1/CAUSES/CNVS/ANNOTATED/


done


