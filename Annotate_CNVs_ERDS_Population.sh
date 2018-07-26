#!/bin/bash
#PBS -N Population_CNV_Annotation
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
# G249-1.ERDS 



for WORKDIR in /mnt/causes-vnx1/CAUSES/CNVS/*ERDS
do
	echo $WORKDIR
        IFS='/' read -a array <<< "${WORKDIR}"
        ERDSDIR=${array[-1]}
        IFS='.' read -a array2 <<< "${ERDSDIR}"
        SAMPLE_ID=${array2[0]}
        echo $SAMPLE_ID
	DUPEVENTS=${WORKDIR}/${SAMPLE_ID}.dup.events
	DELEVENTS=${WORKDIR}/${SAMPLE_ID}.del.events

	#Convert to avinput 
	python /mnt/causes-vnx1/PIPELINES/SV_Calling_And_Annotation/ERDS2AVInput.py $DUPEVENTS $WORKDIR/${SAMPLE_ID}.DUP
	python /mnt/causes-vnx1/PIPELINES/SV_Calling_And_Annotation/ERDS2AVInput.py $DELEVENTS $WORKDIR/${SAMPLE_ID}.DEL

	#Annotate DEL with RefSeq via Annovar
	perl /mnt/causes-data01/data/Databases/annovar/table_annovar.pl --verbose $WORKDIR/${SAMPLE_ID}.DEL /mnt/causes-data01/data/Databases/annovar/humandb \
	--buildver hg19 --out $WORKDIR/${SAMPLE_ID}.ERDS.DEL.annovar --nastring . --remove --otherinfo --protocol refGene,genomicSuperDups,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed \
	--bedfile InHouseDB_SV_100_20180724_DEL_name.bed,\
InHouseDB_SV_500_20180724_DEL_name.bed,\
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
	-arg ',-colsWanted 4,-colsWanted 9 -minqueryfrac 0.9,-colsWanted 9 -minqueryfrac 0.9,-colsWanted 5,-colsWanted 5,-colsWanted 5,-colsWanted 5,\
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
	-buildver hg19 -out $WORKDIR/${SAMPLE_ID}.ERDS.DUP.annovar -nastring . -remove --otherinfo -protocol refGene,genomicSuperDups,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed,bed \
	-bedfile InHouseDB_SV_100_20180724_DUP_name.bed,\
InHouseDB_SV_500_20180724_DUP_name.bed,\
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
	-arg ',-colsWanted 4,-colsWanted 9 -minqueryfrac 0.9,-colsWanted 9 -minqueryfrac 0.9,-colsWanted 5,-colsWanted 5,-colsWanted 5,-colsWanted 5,\
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


	#remove header
	tail -n +2 $WORKDIR/${SAMPLE_ID}.ERDS.DEL.annovar.hg19_multianno.txt > $WORKDIR/${SAMPLE_ID}.ERDS.DEL.annovar.hg19_multianno.bed
	tail -n +2 $WORKDIR/${SAMPLE_ID}.ERDS.DUP.annovar.hg19_multianno.txt > $WORKDIR/${SAMPLE_ID}.ERDS.DUP.annovar.hg19_multianno.bed

	# concatenate dups and dels together
	cat $WORKDIR/${SAMPLE_ID}.ERDS.DEL.annovar.hg19_multianno.bed \
	$WORKDIR/${SAMPLE_ID}.ERDS.DUP.annovar.hg19_multianno.bed \
	> $WORKDIR/${SAMPLE_ID}.ERDS.annovar.hg19_multianno.bed

	/opt/tools/bedtools/bin/sortBed -i $WORKDIR/${SAMPLE_ID}.ERDS.annovar.hg19_multianno.bed > $WORKDIR/${SAMPLE_ID}.ERDS.annovar.sorted.hg19_multianno.sorted.bed

	# Add header, and get rid of stupid Name= from Annovar, and Score=
	cat /mnt/causes-vnx1/PIPELINES/SV_Calling_And_Annotation/TableAnnotators/Annotation_CNV_ERDS_Header.txt \
	$WORKDIR/${SAMPLE_ID}.ERDS.annovar.sorted.hg19_multianno.sorted.bed | \
	sed -e 's/Score=//g' | \
	sed -e 's/Name=//g' > $WORKDIR/${SAMPLE_ID}.ERDS.annovar.sorted.hg19_multianno.sorted.header.bed

	python /mnt/causes-vnx1/PIPELINES/SV_Calling_And_Annotation/TableAnnotators/CNVTable2Final.py \
	-i $WORKDIR/${SAMPLE_ID}.ERDS.annovar.sorted.hg19_multianno.sorted.header.bed \
	-o $WORKDIR/${SAMPLE_ID}.ERDS.annovar.sorted.hg19_multianno.sorted.header.final.tsv 
	
	cp $WORKDIR/${SAMPLE_ID}.ERDS.annovar.sorted.hg19_multianno.sorted.header.final.tsv /mnt/causes-vnx1/CAUSES/CNVS/ANNOTATED/


done


