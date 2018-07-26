#########################################
# Author: Phillip Richmond              #
# Contact: phillip.a.richmond@gmail.com #
# License: open source GNU              #
#########################################

# The purpose of this python script is to take in a gemini table file, and add some annotations. 
# The annotations are all gene-based at this point, and are simple to add from text files. The heavy lifting of finding the gene has already been done
# This script also assumes only 1 gene being in the annotation column. For cases of overlapping genes, or intergenic variants, this script will not handle those annotations well (as of April 23rd 2018)

###########
# Imports #
###########

import sys, argparse
import re


########################
# Function Definitions #
########################

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--InFile",help="This is a file that has come as an output from the GEMINI query/built-in functions (e.g. de novo)",type=str,required=True)
    parser.add_argument("-o","--OutFile",help="Output Gemini Table file with Additional Annotations",type=str,required=True)
    parser.add_argument("-R","--Remove",help='Remove variants that fall below reciprocal overlap calculations',action='store_true',default=False)
    args = parser.parse_args()
    infilename=args.InFile
    outfilename=args.OutFile
    remove=args.Remove
    return infilename,outfilename,remove


# This script relies on a few different files in order to add some manual annotations. Those are defined in the main function at the bottom.


## Get a dictionary for the RefSeq Gene Summaries, gene is key, value is summary
def GetSummaryDict(SummaryFileName):
    SummaryFile = open(SummaryFileName,'r')
    SummaryForGenes = {}
    for line in SummaryFile:
        Number,Symbol,Summary = line.strip('\n').split("\t")
        SummaryForGenes[Symbol]=Summary
    return SummaryForGenes

## Get a dictionary for the OMIM gene map, gene is key, value is Mim#
def GetOMIM_Gene2MimDict(OMIM_Gene2MimFileName):
    OMIMfile = open(OMIM_Gene2MimFileName,'r')
    GenesToMim = {}
    for line in OMIMfile:
        if line[0]=='#':
            continue
        cols=line.strip('\n').split('\t')
        mimNumber= cols[0]
        if cols[1]!='gene':
            continue
        GeneSymbol=cols[3]
        GenesToMim[GeneSymbol]=mimNumber
    return GenesToMim

## Get a dictionary for the OMIM gene to phenotypes, gene is key, phenotypes are the values, with multiple phenos split by '|'
def GetOMIM_Gene2PhenoDict(OMIM_Gene2PhenoFileName):
	OMIMfile = open(OMIM_Gene2PhenoFileName,'r')
	GeneToPheno = {}
	for line in OMIMfile:
		if line[0]=='#':
			continue
		cols = line.strip('\n').split('\t')
		gene = cols[0]
		phenos = "|".join(cols[1:])
		GeneToPheno[gene]=phenos	
	return GeneToPheno

#### 
# RVIS file looks like this:
#GENE	ALL_0.01%	%ALL_0.01%
#A1BG	-0.35	29.49

def GetRVISDict(RVISFILE):
	infile = open(RVISFILE,'r')
	Gene2RVIS = {}
	headerline = infile.readline()
	for line in infile:
		gene,rvis_score,rvis_pct = line.strip('\n').split('\t')
		Gene2RVIS[gene] = [rvis_score,rvis_pct]
	#print Gene2RVIS
	return Gene2RVIS


# PLI file looks like this:
# gene	pLI	pRec
#AGRN	0.17335234048116	0.826647657926747

def GetPLIDict(PLIFILE):
	infile = open(PLIFILE,'r')
	Gene2PLI = {}
	headerline = infile.readline()
	for line in infile:
		gene,pli_score,pli_pct = line.strip('\n').split('\t')
		Gene2PLI[gene] = [pli_score,pli_pct]
	#print Gene2PLI
	return Gene2PLI

def GetHPODict(HPOFILE):
	infile = open(HPOFILE,'r')
	Gene2HPO={}
	headerline = infile.readline()
	for line in infile:
		entrez,gene,pheno,hpo = line.strip('\n').split('\t')
		if Gene2HPO.has_key(gene):
			Gene2HPO[gene].append("%s|%s"%(pheno,hpo))
		else:
			Gene2HPO[gene]=["%s|%s"%(pheno,hpo)]
	#print Gene2HPO
	return Gene2HPO

def GetMESHOPDict(MESHOPFILE):
        infile = open(MESHOPFILE,'r')
        Gene2MESHOP={}
        for line in infile:
                gene,meshop = line.strip('\n').split('\t')
                if Gene2MESHOP.has_key(gene):
                        Gene2MESHOP[gene].append(meshop)
                else:
                        Gene2MESHOP[gene]=[meshop]
        #print Gene2MESHOP
        return Gene2MESHOP	




# This function takes in the dictionaries desribed above, reads in the gemini infile, and outputs the gemini outfile
def AddColumnsToTable(GeminiInFileName,GeminiOutFileName,Gene2Pheno,Gene2Mim,GeneSummary,Gene2PLI,Gene2RVIS,Gene2HPO,Gene2MESHOP):
	infile = open(GeminiInFileName,'r')
	outfile = open(GeminiOutFileName,'w')
	header = infile.readline()
	header = header.strip('\n')
	outfile.write("%s\tOMIM_Entry\tOMIM_Phenotypes\tRVIS_Score\tRVIS_Pct\tpLI_Score\tpLI_Pct\tHPO\tMeSHOP\tGeneSummary\n"%header)
	for line in infile:
		line = line.strip('\n')
		cols = line.split("\t")
		# initialize a list of genes
		geneList = cols[6].split(',')
		# I will store all of the annotations for each of these genes within arrays. They will by default be listed in the same manner as the genes are, and they will be ';' separated, just like the genes
		pLIScores = []
		pLIPcts = []
		omimPhenos = []
		omimKeys = []
		hpoTerms = []
		RVISScores = []
		RVISPcts = []
		MeSHOPs = []
		geneSummaries = []
		for gene in geneList:
			# Check for Omim phenotype, add if there, if not make it '.'
			if Gene2Pheno.has_key(gene):
				omim_pheno=Gene2Pheno[gene]
			else:
				omim_pheno='.'
			# Check for Omim gene, add hyperlink for excel if there, if not make it '.'
			if Gene2Mim.has_key(gene):
				omim_key=Gene2Mim[gene]
				# Hyperlink won't work here like it does for GEMINI
				omim_hyperlink='=HYPERLINK(\"http://www.omim.org/entry/%s\")'%omim_key
			else:
				omim_hyperlink='.'
				omim_key='.'

			# HPO		
			if Gene2HPO.has_key(gene):
				hpo="|".join(Gene2HPO[gene])
			else:
				hpo = '.'	
			
			# PLI
			if Gene2PLI.has_key(gene):
				pLI_score = Gene2PLI[gene][0]
				pLI_pct = Gene2PLI[gene][1]
			else:
				pLI_score = '.'
				pLI_pct = '.'		
			
			# RVIS	
			if Gene2RVIS.has_key(gene):
				rvis_score = Gene2RVIS[gene][0]
				rvis_pct = Gene2RVIS[gene][1]
			else:
				rvis_score = '.'
				rvis_pct = '.'			

			# MESHOP
			if Gene2MESHOP.has_key(gene):
                	        meshop="|".join(Gene2MESHOP[gene])
                	else:
                	        meshop = '.'	

 			# Check for gene summary, if there, add it, if not, make it '.' 
			if GeneSummary.has_key(gene):
				gene_summary=GeneSummary[gene]
			else:
				gene_summary='.'
				
			pLIScores.append(pLI_score)
			pLIPcts.append(pLI_pct)
	                omimPhenos.append(omim_pheno)
	                omimKeys.append(omim_key)
	                hpoTerms.append(hpo)
	                RVISScores.append(rvis_score)
	                RVISPcts.append(rvis_pct)
	                MeSHOPs.append(meshop)
	                geneSummaries.append(gene_summary)
			

		# Now I'll collapse them all
		omim_pheno=";".join(omimPhenos)
		omim_key=';'.join(omimKeys)
		rvis_score=";".join(RVISScores)
		rvis_pct=";".join(RVISPcts)
		pLI_score=";".join(pLIScores)
		pLI_pct=";".join(pLIPcts)
		hpo=";".join(hpoTerms)
		meshop=";".join(MeSHOPs)
		gene_summary=";".join(geneSummaries)
	

		# join the cols
		newline = "\t".join(cols)

		outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(newline,omim_pheno,omim_key,rvis_score,rvis_pct,pLI_score,pLI_pct,hpo,meshop,gene_summary))

		#reset variables
		gene_summary = '.'
		omim_pheno = '.'
		omim_key = '.'
		rvis_score = '.'
		rvis_pct = '.'
		pLI_score = '.'
		pLI_pct = '.'			
		hpo = '.'
		meshop = '.'

# This function returns the number of nucleotide overlap between A and B, divided by the length of B. Since we already know that nucleotides in A are at least 50% covered by B, we want to now know how many nucleotides in B are covered by A
# This is important because several database variants may completely overlap a small single observed variant, but they are inherently different 
def Reciprocal(Astart,Aend,Bstart,Bend):
	Alen=int(Aend)-int(Astart)
	Blen=int(Bend)-int(Bstart)
	array = [Astart,Aend,Bstart,Bend]
	sortedarray=sorted(array)
	Inner=int(sortedarray[1])-int(sortedarray[2])
	overlap=(float(Inner)/Blen)
	return overlap


# This function is designed to take in the expanded table intermediate file (after additions have been made that draw upon the gene column from AddColumnsToTable)
# It will read in the full table, and run reciprocal overlap on each column that is in the column list argument
# Remove is a boolean set by an argument given to the program. if true, then it will remove variants that don't pass reciprocal overlap
def RunReciprocalOverlap(infilename,remove,outfilename):
	infile = open(infilename,'r')
	outfile = open(outfilename,'w')
	header=infile.readline()
	headerCols = header.strip('\n').split('\t')
	outfile.write(header)
	for line in infile:
		cols=line.strip('\n').split('\t')
		VarStart=int(cols[1])
		VarEnd=int(cols[2])
		for i in range(len(cols)):
			# Initialize a new column value
			newCol = cols[i]
			# this is standardized in the header
			if "_Coordinate" in headerCols[i]:
				# check to see if there is anything there, cuz there may not be
				if len(cols[i]) < 4:
					continue
				# Define a newRegions array
				newRegions = []
				# split the regions that are appearing in that column by commas
				regions=cols[i].split(',')
				# Go through the regions one by one
				for region in regions:
					# Split based on the colon
					regionCols=region.split(':')
					# if there is an AC value in the region Cols, then it goes AC:#:chrom:start:end
					if regionCols[0]=='AC':
						RegionStart=regionCols[3]
						RegionEnd=regionCols[4]
						# Generate the overlap of B with A as a proportion of B
						# Here, overlap a given region with the called variant, and report as a fraction of that region size
						overlap=Reciprocal(VarStart,VarEnd,RegionStart,RegionEnd)
						# Then add that right into the region cols at the end, so now looks like: AC:#:chrom:start:end:overlap
						regionCols.append('%f'%overlap)
						# Join it together
						newRegionCols=":".join(regionCols)
						newRegions.append(newRegionCols)
					# Otherwise, it goes: chrom:start:end, rest of the logic is same as above
					else:
						print regions
						print region
						print headerCols[i]
						RegionStart=regionCols[1]
                                                RegionEnd=regionCols[2]
                                                overlap=Reciprocal(VarStart,VarEnd,RegionStart,RegionEnd)
                                                regionCols.append('(%f)'%overlap)
                                                newRegionCols=":".join(regionCols)
                                                newRegions.append(newRegionCols)
				newColumn = ','.join(newRegions)
				cols[i] = newColumn
			else:
				continue

		newline = '\t'.join(cols)
		outfile.write("%s\n"%newline)
	return


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

# These are hard coded locations for database files. These files are small, so they are all just text files.
# They contain gene-based information
    SummaryFileName = '/mnt/causes-data01/data/Databases/RefSeqGene_Summaries_270316.txt'
    Gene2MimFileName = '/mnt/causes-data01/data/Databases/OMIM_mim2gene'
    Gene2Disease = '/mnt/causes-data01/data/Databases/OMIM_phenotype_genelist'
    PLI = '/mnt/causes-data01/data/Databases/TOLERANCE/PLI_March2016.txt'
    RVIS = '/mnt/causes-data01/data/Databases/TOLERANCE/RVIS_March2016.txt'
    MESHOP = '/mnt/causes-data01/data/Databases/gene2pubmedBG-hum-gene2pubmed-gene-mesh-p_ONLYMESHDISEAS_P-valuecorrected_withGeneSymbols.txt'
    HPO = '/mnt/causes-data01/data/Databases/ALL_SOURCES_FREQUENT_FEATURES_genes_to_phenotype.txt'

    # Read in the annotations
    GeneSummaries = GetSummaryDict(SummaryFileName)
    Gene2Mim = GetOMIM_Gene2MimDict(Gene2MimFileName)
    Gene2Pheno = GetOMIM_Gene2PhenoDict(Gene2Disease)
    Gene2PLI = GetPLIDict(PLI)
    Gene2RVIS = GetRVISDict(RVIS)
    Gene2HPO = GetHPODict(HPO)
    Gene2MESHOP = GetMESHOPDict(MESHOP)
    
    # Read in the options from argparse	
    CNVInfile,CNVOutfile,Remove=GetOptions()
    # Add columns based on gene identifier column
    AddColumnsToTable(CNVInfile,'%s_allCols'%CNVInfile,Gene2Pheno,Gene2Mim,GeneSummaries,Gene2PLI,Gene2RVIS,Gene2HPO,Gene2MESHOP)

    RunReciprocalOverlap('%s_allCols'%CNVInfile,Remove,CNVOutfile)
    

    sys.exit()



