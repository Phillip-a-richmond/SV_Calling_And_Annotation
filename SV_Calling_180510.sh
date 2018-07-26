#!/bin/bash
#PBS -N G439-1_SVpipeline
#PBS -V
#PBS -o /mnt/causes-data02/new/Process/G439/SV/G439-1_SVpipeline.o
#PBS -e /mnt/causes-data02/new/Process/G439/SV/G439-1_SVpipeline.e
#PBS -m bea
#PBS -M jmwenifumbo@bcchr.ca
#PBS -l walltime=250:00:00
#PBS -l nodes=1:ppn=12
#PBS -l mem=45gb

source /opt/tools/hpcenv.sh

REF='/mnt/causes-data01/data/GENOMES/GSC/GRCh37-lite.fa'
CHROM='/mnt/causes-data01/data/GENOMES/GSC/hg19a_per_chr_fastas/'
BAM='A80151_36_lanes_dupsFlagged.bam'
BAMDIR='/mnt/causes-data02/new/Process/G439/'
SAMPLE='G439-1'
WORKDIR='/mnt/causes-data02/new/Process/G439/SV/'
GENOME_FASTA='/mnt/causes-data01/data/GENOMES/GSC/GRCh37-lite.fa'

#Index bam file
#/opt/tools/samtools-1.2/samtools index $BAM

#Extract discordant reads
/opt/tools/samtools-1.2/samtools view \
	-b \
	-F 1294 \
$BAMDIR/$BAM > $WORKDIR/${SAMPLE}.discordants.unsorted.bam

#Extract split reads
/opt/tools/samtools-1.2/samtools view -h $BAMDIR/$BAM \
	| /opt/tools/lumpy-0.2.11/scripts/extractSplitReads_BwaMem -i stdin \
	| /opt/tools/samtools-1.2/samtools view -Sb - > $WORKDIR/${SAMPLE}.splitters.unsorted.bam

##Sort alignments
/opt/tools/samtools-1.2/samtools sort $WORKDIR/${SAMPLE}.discordants.unsorted.bam $WORKDIR/${SAMPLE}.discordants
/opt/tools/samtools-1.2/samtools sort $WORKDIR/${SAMPLE}.splitters.unsorted.bam $WORKDIR/${SAMPLE}.splitters 

#Obtain insert size metrics
/opt/tools/samtools-1.2/samtools view $BAMDIR/$BAM \
	| tail -n+100000 \
	| /opt/tools/lumpy-0.2.11/scripts/pairend_distro.py \
	-r 150 \
	-X 4 \
	-N 10000 \
	-o $WORKDIR/${SAMPLE}.histo > $WORKDIR/${SAMPLE}.insertsize

INSMEAN=$(awk '{print $1}' $WORKDIR/${SAMPLE}.insertsize | cut -d ':' -f2)
INSSTD=$(awk '{print $2}' $WORKDIR/${SAMPLE}.insertsize | cut -d ':' -f2)

#Lumpy
/opt/tools/lumpy-0.2.11/lumpy \
-mw 4 \
-tt 0 \
-pe id:$WORKDIR/$SAMPLE,bam_file:$WORKDIR/${SAMPLE}.discordants.bam,histo_file:$WORKDIR/${SAMPLE}.histo,mean:$INSMEAN,stdev:$INSSTD,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,min_mapping_threshold:20 \
-sr id:$WORKDIR/$SAMPLE,bam_file:$WORKDIR/${SAMPLE}.splitters.bam,back_distance:20,weight:1,min_mapping_threshold:20 > $WORKDIR/${SAMPLE}.LUMPY.vcf

##CNVnator
/opt/tools/CNVnator/src/cnvnator -root $WORKDIR/${SAMPLE}.root -genome $REF -tree $BAMDIR/$BAM
/opt/tools/CNVnator/src/cnvnator -root $WORKDIR/${SAMPLE}.root -genome $REF -his 100 -d $CHROM
/opt/tools/CNVnator/src/cnvnator -root $WORKDIR/${SAMPLE}.root -genome $REF -stat 100
/opt/tools/CNVnator/src/cnvnator -root $WORKDIR/${SAMPLE}.root -genome $REF -partition 100
/opt/tools/CNVnator/src/cnvnator -root $WORKDIR/${SAMPLE}.root -genome $REF -call 100 > $WORKDIR/${SAMPLE}_CNVcall_100

##Pindel
python /opt/tools/pindel-0.2.5b6/pindel_config.py $BAMDIR/$BAM $INSMEAN $WORKDIR/$SAMPLE
/opt/tools/pindel-0.2.5b6/pindel \
        -f $REF \
        -i $BAMDIR/${BAM}_config.txt \
        -c ALL \
        -o $WORKDIR/$SAMPLE \
        --number_of_threads 12 \
        -M 4 \
        -N \
        -r false \
        -t false \
        -I false \
        -x 2 \
        -J /mnt/causes-data01/data/Databases/hg19_centromeres_telomeres.bed 

# CNV Calling with ERDS
perl /opt/tools/erds1.1/erds_pipeline.pl \
        -b $BAMDIR/$BAM \
        -v $BAMDIR/${SAMPLE}_haplotypecaller_split_norm_filter.vcf.gz \
        -o $WORKDIR/${SAMPLE}.ERDS/ \
        -r $GENOME_FASTA

##metaSV
run_metasv.py --bam $BAMDIR/$BAM --reference $REF \
        --lumpy_vcf $WORKDIR/${SAMPLE}.LUMPY.vcf \
        --cnvnator_native $WORKDIR/${SAMPLE}_CNVcall_100 \
        --pindel_native $WORKDIR/${SAMPLE}_D $WORKDIR/${SAMPLE}_SI \
        --filter_gaps \
        --sample $WORKDIR/$SAMPLE --num_threads 12 \
        --outdir $WORKDIR/${SAMPLE}Out \
        --workdir $WORKDIR/${SAMPLE}Work \
        --spades /opt/tools/SPAdes-3.10.1/bin/spades.py \
        --age /opt/tools/AGE/age_align \
        --svs_to_assemble INV \
        --minsvlen 15

