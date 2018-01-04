# Before running the script load the modules present in MODULES file

# This is a script for combined MRNA and MiRNA analysis
# The MRNA and MiRNA analyses are combined for easy use of intermediate files between both analyses

# Needs gene counts matrix files for both MRNA and MiRNA RNA-seq samples  
# and GTF annotation for reference and RNA assembly data, and samples information (Samples.txt) file,
# in the current directory or specified by their complete path

# Start fresh
rm(list=ls())

# Working directory
# Set the working directory, which has all function files
WD <- "/home/smalkaram/Projects/alwayIntegratedMiRNAmRNADeAnalysis"
setwd(WD)

# Start Log
sink()
sink(file=paste(format(Sys.Date(), format="%m-%d-%Y"), format(Sys.time(), format="%H-%M-%S"), "log", sep="."), split=TRUE)

# Load Libraries & Functions
source("PACKAGES")
source("FUNCTIONS.R")

# Set the necessary Inputs
# MRNA_COUNTS_FILE = Counts file for all samples of mRNA generated using HISAT2 and STRINGTIE pipeline (See PRE_MRNA_ANALYSIS subdirectory)
# MiRNA_COUNTS_FILE = Counts file for all samples of miRNA generated using miRDeep2 pipeline (See PRE_MIRNA_ANALYSIS subdirectory)
# MiRNA_MAP_FILE = A mapping file that maps the sample names and the identifiers output from miRDeep2
# PVAL = Adjusted pvalue cutoff
# TYPE = can be either gene or transcript depending on your counts files
# GTF = file generated from the HISAT2 and STRINGTIE pipeline
# refGTF = reference GTF file for the genome assembly used for alignments
# MAX_TARGETS = maximum number of miRNA targets to output in DEresult tables
# CUTOFF_READ = cutoff # for minimum reads for any gene in any sample
# TYPE = "gene" # necessary for compatibility with pipeline scripts

TYPE <- "gene"
PVAL <- 0.05
refGTF <- "gencode.vM14.primary_assembly.annotation.gtf"
GTF <- "Merged.gtf"
MRNA_COUNTS_FILE <- "stringtie.gene.count.matrix"
MiRNA_COUNTS_FILE <- "miRNAs_expressed_all_samples_20170823142906.csv"
MiRNA_MAP_FILE <- "mapping_file2.txt" 
Samples_FILE<- "Samples.txt"
CUTOFF_READ <- 10
MAX_TARGETS<-10

PC(paste("Working directory", WD, sep=":"))
#--------------------------------------------------------------------------------------------------------------
# DESEQ analysis for mRNA 
#--------------------------------------------------------------------------------------------------------------

PC("Analyzing mRNA")

# Reads MRNA Counts file
PC("Reading Counts file", 1)
CountsOrig <- read.delim(MRNA_COUNTS_FILE, sep = "\t", header = TRUE, 
           row.names = 1)

# Cleanup gene names
PC("Cleanup gene names", 1)
Counts <- MergeGeneNamesInExpMat(GTF, refGTF, CountsOrig)
colnames(Counts) <- gsub("[.]", "_", colnames(Counts))
Counts <- Counts[rowSums(Counts) > CUTOFF_READ, ]

# Read sample information file
PC("Reading sample information file", 1)
ColData<-read.table(Samples_FILE, header=TRUE, sep="\t", row.names=1)
genotype<-as.factor(ColData$genotype)
time<-as.factor(ColData$time)
diet<-as.factor(ColData$diet)
newvar<-as.factor(ColData$newvar)

# Compute DE for all samples (for comparisons between 80W and 20W)
PC("Computing DE between 20W and 80W", 1)
Counts <- Counts[, rownames(ColData)]
PC("Samples", 2)
print(colnames(Counts))
dds <- DESeqDataSetFromMatrix(countData = Counts, colData = ColData, 
		design = formula(~ newvar))
dds$genotype<-relevel(dds$genotype, ref="WT")
dds$diet<-relevel(dds$diet, ref="CR")
dds <- DESeq(dds)


# Compute DE for 20W samples
PC("Computing DE within 20W", 1)
Counts20W<-Counts[,grep("20W", colnames(Counts))]
PC("Samples", 2)
print(colnames(Counts20W))
ColData20W<-ColData[colnames(Counts20W),]
Counts20W <- Counts20W[, rownames(ColData20W)]
dds20W <- DESeqDataSetFromMatrix(countData = Counts20W, colData = ColData20W, 
		design = formula(~ genotype + diet + genotype:diet))
dds20W$genotype<-relevel(dds20W$genotype, ref="WT")
dds20W$diet<-relevel(dds20W$diet, ref="CR")
dds20W <- DESeq(dds20W)

# Compute DE for 80W samples
PC("Computing DE within 80W", 1)
Counts80W<-Counts[,grep("80W", colnames(Counts))]
PC("Samples", 2)
print(colnames(Counts80W))
ColData80W<-ColData[colnames(Counts80W),]
Counts80W <- Counts80W[, rownames(ColData80W)]
dds80W <- DESeqDataSetFromMatrix(countData = Counts80W, colData = ColData80W, 
		design = formula(~ genotype + diet + genotype:diet))
dds80W$genotype<-relevel(dds80W$genotype, ref="WT")
dds80W$diet<-relevel(dds80W$diet, ref="CR")
dds80W <- DESeq(dds80W)


PC("Extracting results", 1) 
PC("20W", 2)
PC("genotype", 3)
KO_WT_CR_20W<-results(dds20W, contrast=c("genotype", "KO", "WT"))
OE_WT_CR_20W<-results(dds20W, contrast=c("genotype", "OE", "WT"))
PC("diet", 3)
AL_CR_WT_20W<-results(dds20W, contrast=c("diet", "AL", "CR"))
PC("interactions", 3)
KO_WT_int_AL_CR_20W<-results(dds20W, name="genotypeKO.dietAL")
OE_WT_int_AL_CR_20W<-results(dds20W, name="genotypeOE.dietAL")
KO_WT_AL_20W<-results(dds20W, contrast=list(c("diet_AL_vs_CR", "genotypeKO.dietAL")))
OE_WT_AL_20W<-results(dds20W, contrast=list(c("diet_AL_vs_CR", "genotypeOE.dietAL")))

PC("80W", 2)
PC("genotype", 3)
KO_WT_CR_80W<-results(dds80W, contrast=c("genotype", "KO", "WT"))
OE_WT_CR_80W<-results(dds80W, contrast=c("genotype", "OE", "WT"))
PC("diet", 3)
AL_CR_WT_80W<-results(dds80W, contrast=c("diet", "AL", "CR"))
PC("interactions", 3)
KO_WT_int_AL_CR_80W<-results(dds80W, name="genotypeKO.dietAL")
OE_WT_int_AL_CR_80W<-results(dds80W, name="genotypeOE.dietAL")
KO_WT_AL_80W<-results(dds80W, contrast=list(c("diet_AL_vs_CR", "genotypeKO.dietAL")))
OE_WT_AL_80W<-results(dds80W, contrast=list(c("diet_AL_vs_CR", "genotypeOE.dietAL")))

# Comparison between 20W and 80W 
PC("time", 2)
KO80WCR_KO20WCR<-results(dds, contrast=c("newvar", "KO80WCR", "KO20WCR"))
KO80WAL_KO20WAL<-results(dds, contrast=c("newvar", "KO80WAL", "KO20WAL"))
OE80WCR_OE20WCR<-results(dds, contrast=c("newvar", "OE80WCR", "OE20WCR"))
OE80WAL_OE20WAL<-results(dds, contrast=c("newvar", "OE80WAL", "OE20WAL"))
WT80WCR_WT20WCR<-results(dds, contrast=c("newvar", "WT80WCR", "WT20WCR"))
WT80WAL_WT20WAL<-results(dds, contrast=c("newvar", "WT80WAL", "WT20WAL"))

#Collect all above results into this list
PC("MRNA result variables", 1)
RESULTS <-c('KO_WT_CR_20W', 'OE_WT_CR_20W', 'AL_CR_WT_20W', 'KO_WT_int_AL_CR_20W', 'OE_WT_int_AL_CR_20W',
 'KO_WT_AL_20W', 'OE_WT_AL_20W', 'KO_WT_CR_80W', 'OE_WT_CR_80W', 'AL_CR_WT_80W', 'KO_WT_int_AL_CR_80W',
 'OE_WT_int_AL_CR_80W', 'KO_WT_AL_80W', 'OE_WT_AL_80W', "KO80WCR_KO20WCR", "KO80WAL_KO20WAL", 
 "OE80WCR_OE20WCR", "OE80WAL_OE20WAL", "WT80WCR_WT20WCR", "WT80WAL_WT20WAL")
PC(paste(RESULTS, collapse= " "), 2)

# Get Ensembl-MGI map
PC("Creating Ensembl-MGI map", 1)
library(biomaRt)
MART <- ""
while (class(MART) != "Mart") {
	try(MART <- useDataset(as.character("mmusculus_gene_ensembl"), mart = useMart("ENSEMBL_MART_ENSEMBL")))
}
GENE_MGI <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "mgi_symbol"), 
		mart = MART)
GENE_MGI <- ReplaceEmptyMappings(GENE_MGI, "hgnc_symbol", "mgi_symbol")
GENE_MGI <- ReplaceEmptyMappings(GENE_MGI, "ensembl_gene_id", "mgi_symbol")
GENE_MGI <- unlist(lapply(Map2List(GENE_MGI, "ensembl_gene_id", "mgi_symbol"), 
				function(x) paste(x, collapse = ",")))


# Clean up results tables
# Replace NAs in pvalue and padj column with 1. Add foldchange and absolute foldchange column
# Add MGI code column and prefix the results with "MRNA"
PC("Cleaning up results", 1)
for (VAR in RESULTS) {
	PC(VAR, 2)
	assign(VAR, replaceDFColNAs(get(VAR), "pvalue", 1))
	assign(VAR, replaceDFColNAs(get(VAR), "padj", 1))
	assign(VAR, addMGI4GenCode(get(VAR), GENE_MGI))
	assign(VAR, addFC(get(VAR), logFCCOL = "log2FoldChange"))
	assign(VAR, addAbsFC(get(VAR), FCCOL = "foldChange"))
	assign(paste("MRNA", VAR, sep="_"), addUpDown(get(VAR), log2FCCOL = "log2FoldChange"))
}


#--------------------------------------------------------------------------------------------------------------
# DESEQ analysis for miRNA 
#--------------------------------------------------------------------------------------------------------------
PC("Analyzing miRNA")


# Reads MRNA Counts file
PC("Reading Counts file", 1)
CountsOrig <- read.delim(MiRNA_COUNTS_FILE, sep = "\t", header = TRUE)
CountsOrig <- CountsOrig[CountsOrig$total != 0, ]

# Clean up
PC("Cleanup gene names", 1)
CountsOrig$refname <- gsub("chr.*-", "", CountsOrig[, 1])
CountsOrig$reftotal <- paste(CountsOrig$refname, CountsOrig$total, 
		sep = "_")
CountsOrig <- takeFirstCommonRow(CountsOrig, "reftotal")
CountsOrig <- CountsOrig[, -grep("reftotal", colnames(CountsOrig))]
CountsOrig <- CountsOrig[, -grep("read_count|total|norm", colnames(CountsOrig), 
				perl = TRUE, ignore.case = TRUE)]
CountsOrig <- addByColFactor(CountsOrig, "refname")

# Read MiRNA map file
PC("Read miRDeep2 id to samples map", 1)
MAP <- read.table(MiRNA_MAP_FILE, header = FALSE)
NAMES <- MAP[, 2]
MAP <- gsub("_trimmed.fq.gz", "", MAP[, 1])
names(MAP) <- NAMES

# Change column names to that in map
PC("Change miRDeep2 id to sample names", 1)
colnames(CountsOrig) <- as.vector(MAP[colnames(CountsOrig)])
SAMPLES <- unique(gsub("_S[0-9]+_L[0-9]_[0-9]$", "", colnames(CountsOrig)))
Counts <- sapply(SAMPLES, function(x) GETSUM(x, CountsOrig))

Counts <- Counts[rowSums(Counts) > CUTOFF_READ, ]

# Read sample information file (No Need to run again as it is done for MRNA)
#PC("Reading sample information file", 1)
#ColData<-read.table(Samples_FILE, header=TRUE, sep="\t", row.names=1)
#genotype<-as.factor(ColData$genotype)
#time<-as.factor(ColData$time)
#diet<-as.factor(ColData$diet)
#newvar<-as.factor(ColData$newvar)

PC("Compute DE", 1)
# Compute DE for all samples (for comparisons between 80W and 20W)
PC("Computing DE between 20W and 80W", 2)
Counts <- Counts[, rownames(ColData)]
PC("Samples", 2)
print(colnames(Counts))
dds <- DESeqDataSetFromMatrix(countData = Counts, colData = ColData, 
		design = formula(~ newvar))
dds$genotype<-relevel(dds$genotype, ref="WT")
dds$diet<-relevel(dds$diet, ref="CR")
dds <- DESeq(dds)


# Compute DE for 20W samples
PC("Computing DE within 20W", 2)
Counts20W<-Counts[,grep("20W", colnames(Counts))]
PC("Samples", 2)
print(colnames(Counts20W))
ColData20W<-ColData[colnames(Counts20W),]
Counts20W <- Counts20W[, rownames(ColData20W)]
dds20W <- DESeqDataSetFromMatrix(countData = Counts20W, colData = ColData20W, 
		design = formula(~ genotype + diet + genotype:diet))
dds20W$genotype<-relevel(dds20W$genotype, ref="WT")
dds20W$diet<-relevel(dds20W$diet, ref="CR")
dds20W <- DESeq(dds20W)

# Compute DE for 80W samples
PC("Computing DE within 80W", 2)
Counts80W<-Counts[,grep("80W", colnames(Counts))]
PC("Samples", 2)
print(colnames(Counts80W))
ColData80W<-ColData[colnames(Counts80W),]
Counts80W <- Counts80W[, rownames(ColData80W)]
dds80W <- DESeqDataSetFromMatrix(countData = Counts80W, colData = ColData80W, 
		design = formula(~ genotype + diet + genotype:diet))
dds80W$genotype<-relevel(dds80W$genotype, ref="WT")
dds80W$diet<-relevel(dds80W$diet, ref="CR")
dds80W <- DESeq(dds80W)

PC("Extracting results", 1) 
PC("20W", 2)
PC("genotype", 3)
KO_WT_CR_20W<-results(dds20W, contrast=c("genotype", "KO", "WT"))
OE_WT_CR_20W<-results(dds20W, contrast=c("genotype", "OE", "WT"))
PC("diet", 3)
AL_CR_WT_20W<-results(dds20W, contrast=c("diet", "AL", "CR"))
PC("interactions", 3)
KO_WT_int_AL_CR_20W<-results(dds20W, name="genotypeKO.dietAL")
OE_WT_int_AL_CR_20W<-results(dds20W, name="genotypeOE.dietAL")
KO_WT_AL_20W<-results(dds20W, contrast=list(c("diet_AL_vs_CR", "genotypeKO.dietAL")))
OE_WT_AL_20W<-results(dds20W, contrast=list(c("diet_AL_vs_CR", "genotypeOE.dietAL")))

PC("80W", 2)
PC("genotype", 3)
KO_WT_CR_80W<-results(dds80W, contrast=c("genotype", "KO", "WT"))
OE_WT_CR_80W<-results(dds80W, contrast=c("genotype", "OE", "WT"))
PC("diet", 3)
AL_CR_WT_80W<-results(dds80W, contrast=c("diet", "AL", "CR"))
PC("interactions", 3)
KO_WT_int_AL_CR_80W<-results(dds80W, name="genotypeKO.dietAL")
OE_WT_int_AL_CR_80W<-results(dds80W, name="genotypeOE.dietAL")
KO_WT_AL_80W<-results(dds80W, contrast=list(c("diet_AL_vs_CR", "genotypeKO.dietAL")))
OE_WT_AL_80W<-results(dds80W, contrast=list(c("diet_AL_vs_CR", "genotypeOE.dietAL")))

# Comparison between 20W and 80W 
PC("time", 2)
KO80WCR_KO20WCR<-results(dds, contrast=c("newvar", "KO80WCR", "KO20WCR"))
KO80WAL_KO20WAL<-results(dds, contrast=c("newvar", "KO80WAL", "KO20WAL"))
OE80WCR_OE20WCR<-results(dds, contrast=c("newvar", "OE80WCR", "OE20WCR"))
OE80WAL_OE20WAL<-results(dds, contrast=c("newvar", "OE80WAL", "OE20WAL"))
WT80WCR_WT20WCR<-results(dds, contrast=c("newvar", "WT80WCR", "WT20WCR"))
WT80WAL_WT20WAL<-results(dds, contrast=c("newvar", "WT80WAL", "WT20WAL"))

#Collect all above results into this list
RESULTS <-c('KO_WT_CR_20W', 'OE_WT_CR_20W', 'AL_CR_WT_20W', 'KO_WT_int_AL_CR_20W', 'OE_WT_int_AL_CR_20W',
		'KO_WT_AL_20W', 'OE_WT_AL_20W', 'KO_WT_CR_80W', 'OE_WT_CR_80W', 'AL_CR_WT_80W', 'KO_WT_int_AL_CR_80W',
		'OE_WT_int_AL_CR_80W', 'KO_WT_AL_80W', 'OE_WT_AL_80W', "KO80WCR_KO20WCR", "KO80WAL_KO20WAL", 
		"OE80WCR_OE20WCR", "OE80WAL_OE20WAL", "WT80WCR_WT20WCR", "WT80WAL_WT20WAL")
PC("MiRNA result variables", 1)
PC(paste(RESULTS, collapse= " "), 2)

# Clean up results tables
# Replace NAs in pvalue and padj column with 1. Add foldchange and absolute foldchange column
# prefix the results with "MiRNA"
PC("Cleaning up results", 1)
for (VAR in RESULTS) {
	PC(VAR, 2)
	assign(VAR, replaceDFColNAs(get(VAR), "pvalue", 1))
	assign(VAR, replaceDFColNAs(get(VAR), "padj", 1))
	assign(VAR, addFC(get(VAR), logFCCOL = "log2FoldChange"))
	assign(VAR, addAbsFC(get(VAR), FCCOL = "foldChange"))
	assign(paste("MiRNA", VAR, sep="_"), addUpDown(get(VAR), log2FCCOL = "log2FoldChange"))
}

#--------------------------------------------------------------------------------------------------------------
# Combine the results and find miRNA targets
PC("Combine DE results")

PC("Get a common list of MRNA present in all results", 1)
# Get the common list of all MRNA names from all DE results. Remove novel genes identified by stringtie
MRNA_QUERY <- unique(unlist(lapply(RESULTS, function(X) { rownames(get(paste("MRNA", X, sep = "_"))) })))
MRNA_QUERY <- gsub("[.][0-9]*$", "", MRNA_QUERY)
MRNA_QUERY <- MRNA_QUERY[-grep("^MSTRG", MRNA_QUERY)]
#Get targets for this list
PC("Get targets", 2)
MRNA_TARGET_LIST <- get.multimir(org = "mmu", target = MRNA_QUERY, summary = FALSE, table = "all", add.link = FALSE, 
    use.tibble = FALSE, predicted.site = "conserved", predicted.cutoff = 10000, predicted.cutoff.type = "n")

PC("Get a common list of MiRNA present in all results", 1)
# Get the common list of all MiRNA names from all DE results. 
MiRNA_QUERY <- unique(unlist(lapply(RESULTS, function(X) { rownames(get(paste("MiRNA", X, sep = "_"))) })))
#Get targets for this list
PC("Get targets", 2)
MiRNA_TARGET_LIST <- get.multimir(org = "mmu", mirna = MiRNA_QUERY, summary = FALSE, table = "all", add.link = FALSE, 
    use.tibble = FALSE, predicted.site = "conserved", predicted.cutoff = 10000, predicted.cutoff.type = "n")

#--------------------------------------------------------------------------------------------------------------
PC("Writing final result tables")
# Analyze the results of DESeq2 results for both MRNA and MiRNA
# And the targets 
# Provide a variable 'VAR', which has DESeq2 results
# Writes DE tables with extra columns as for both MRNA and MiRNA 
# at 0.05, 0.1 and 1(all) adjusted pvalue cutoffs sorting on absolute FoldChange 

# If this is on compute cluster with torque resource manager, 
# enable PARALLEL=TRUE only if its requirements, as below are met, otherwise use PARALLEL=FALSE

#PARALLEL=TRUE
PARALLEL=FALSE

if (PARALLEL == TRUE)
{
	PC("Writing results in parallel mode", 1)
	# PARALLEL PROCESSING USING CLUSTER-------------------
	# Results are stored in respective directories with same name as the result variable name
	# This needs the R script file Alway_mRNA_microRNA_analysis.Sub1.R and a bash script file run_Alway_mRNA_microRNA_analysis.Sub1.sh that actually submits the R script to the torque queue
	# Need two files: 'MODULES' : having the required modules and 'PACKAGES': having the libraries to be loaded
	# Need two files: 'MODULES' : having the required modules and 'PACKAGES': having the libraries to be loaded
	PC("Saving current environment", 1)
	save(list=ls(all.names=TRUE, .GlobalEnv), ascii=FALSE, file=paste("Alway_mRNA_microRNA_analysis.RData", FILE, sep="/"))

	for (VAR in RESULTS) {
		PC(VAR, 2)
		system(paste("run_Alway_mRNA_microRNA_analysis.Sub1.sh  Alway_mRNA_microRNA_analysis.Sub1.R", VAR, paste=" "))
	}
}
# ------------------------------------


if (PARALLEL == FALSE)
{
# SERIAL PROCESSING EACH RESULT-------------------
# All results are stored in the current directory
PC("Writing results in serial mode", 1)

for (VAR in RESULTS) {
	PC(VAR, 1)
	MiRNA_VAR <- paste("MiRNA", VAR, sep = "_")
	MiRNA_DE <- get(MiRNA_VAR)
	MRNA_VAR <- paste("MRNA", VAR, sep = "_")
	MRNA_DE <- get(MRNA_VAR)
	MRNA_DE <- MRNA_DE[-grep("^MSTRG", rownames(MRNA_DE)), ]
	rownames(MRNA_DE) <- gsub("[.][0-9]*$", "", rownames(MRNA_DE))
	MRNA_QUERY <- rownames(MRNA_DE)
	MiRNA_QUERY <- rownames(MiRNA_DE)

	#MiRNA
	PC('MiRNA', 2)
	PREV_COLUMNS <- colnames(MiRNA_DE)
	ADD_COLUMNS <- c("val.sig.up", "val.sig.down", "val.nsig.up", "val.nsig.down",
			"pred.sig.up", "pred.sig.down", "pred.nsig.up", "pred.nsig.down",
			"val.sig.up.freq", "val.sig.down.freq", "val.nsig.up.freq", "val.nsig.down.freq",
			"pred.sig.up.freq", "pred.sig.down.freq", "pred.nsig.up.freq",
			"pred.nsig.down.freq")
	TEMP <- cbind(as.data.frame(MiRNA_DE), t(sapply(MiRNA_QUERY, function(QUERY) {
								EXTRACT_TARGETS(QUERY, MiRNA_TARGET_LIST, "mature_mirna_id", "target_ensembl",
										MAX_TARGETS, MRNA_DE) })))
	colnames(TEMP) <- c(PREV_COLUMNS, ADD_COLUMNS)
	WRITE_SIG_TABLES(DF = TEMP, NAME = MiRNA_VAR, PVAL = "padj", ALPHA = 0.05,
			FC = "absFoldChange")
	WRITE_SIG_TABLES(DF = TEMP, NAME = MiRNA_VAR, PVAL = "padj", ALPHA = 0.1,
			FC = "absFoldChange")
	WRITE_SIG_TABLES(DF = TEMP, NAME = MiRNA_VAR, PVAL = "padj", ALPHA = 1,
			FC = "absFoldChange")

	# MRNA
	PC('MRNA', 2)
	PREV_COLUMNS <- colnames(MRNA_DE)
	TEMP <- cbind(as.data.frame(MRNA_DE), t(sapply(MRNA_QUERY, function(QUERY) {
								EXTRACT_TARGETS(QUERY, MRNA_TARGET_LIST, "target_ensembl", "mature_mirna_id",
										MAX_TARGETS, MiRNA_DE)})))
	colnames(TEMP) <- c(PREV_COLUMNS, ADD_COLUMNS)
	WRITE_SIG_TABLES(DF = TEMP, NAME = MRNA_VAR, PVAL = "padj", ALPHA = 0.05,
			FC = "absFoldChange")
	WRITE_SIG_TABLES(DF = TEMP, NAME = MRNA_VAR, PVAL = "padj", ALPHA = 0.1,
			FC = "absFoldChange")
	WRITE_SIG_TABLES(DF = TEMP, NAME = MRNA_VAR, PVAL = "padj", ALPHA = 1,
			FC = "absFoldChange")
	}
}
# ------------------------------------

sink()
