#! /usr/bin/env R
# This is a sub script and needs an R data file with saved objects from the main analysis file

# Load the arguments
args<-commandArgs(TRUE)
THISVAR <- as.character(args[1])

# Load saved objects
load("Alway_mRNA_microRNA_analysis.RData", .GlobalEnv)

# Load the functions/packages
source("PACKAGES")
source("FUNCTIONS.R")

# Get the variable to run for this instance
VAR<-THISVAR

	PC(VAR, 1)

	# This is the MiRNA variable name and the DE result object
	MiRNA_VAR <- paste("MiRNA", VAR, sep = "_")
	MiRNA_DE <- get(MiRNA_VAR)

	# This is the MRNA variable name and the DE result object
	MRNA_VAR <- paste("MRNA", VAR, sep = "_")
	MRNA_DE <- get(MRNA_VAR)

	# Clean up MRNA DE to not have novel genes
	MRNA_DE <- MRNA_DE[-grep("^MSTRG", rownames(MRNA_DE)), ]
	rownames(MRNA_DE) <- gsub("[.][0-9]*$", "", rownames(MRNA_DE))

	# Get the names for MRNA and MiRNA genes
	MRNA_QUERY <- rownames(MRNA_DE)
	MiRNA_QUERY <- rownames(MiRNA_DE)
	
	#MiRNA
	PREV_COLUMNS <- colnames(MiRNA_DE)
	# These are additional columns to be added to DE tables
	ADD_COLUMNS <- c("val.sig.up", "val.sig.down", "val.nsig.up", "val.nsig.down",
			"pred.sig.up", "pred.sig.down", "pred.nsig.up", "pred.nsig.down",
			"val.sig.up.freq", "val.sig.down.freq", "val.nsig.up.freq", "val.nsig.down.freq",
			"pred.sig.up.freq", "pred.sig.down.freq", "pred.nsig.up.freq",
			"pred.nsig.down.freq")

	# This code extracts MRNA targets for MiRNA genes
	TEMP <- cbind(as.data.frame(MiRNA_DE), t(sapply(MiRNA_QUERY, function(QUERY) {
								EXTRACT_TARGETS(QUERY, MiRNA_TARGET_LIST, "mature_mirna_id", "target_ensembl",
										MAX_TARGETS, MRNA_DE) })))
	colnames(TEMP) <- c(PREV_COLUMNS, ADD_COLUMNS)

	# Write the DE tables
	WRITE_SIG_TABLES(DF = TEMP, NAME = MiRNA_VAR, PVAL = "padj", ALPHA = 0.05,
			FC = "absFoldChange")
	WRITE_SIG_TABLES(DF = TEMP, NAME = MiRNA_VAR, PVAL = "padj", ALPHA = 0.1,
			FC = "absFoldChange")
	WRITE_SIG_TABLES(DF = TEMP, NAME = MiRNA_VAR, PVAL = "padj", ALPHA = 1,
			FC = "absFoldChange")
	
	# MRNA
	PC('MRNA', 2)
	PREV_COLUMNS <- colnames(MRNA_DE)
	# This code extracts MiRNA targets for MRNA genes
	TEMP <- cbind(as.data.frame(MRNA_DE), t(sapply(MRNA_QUERY, function(QUERY) {
								EXTRACT_TARGETS(QUERY, MRNA_TARGET_LIST, "target_ensembl", "mature_mirna_id",
										MAX_TARGETS, MiRNA_DE)})))
	colnames(TEMP) <- c(PREV_COLUMNS, ADD_COLUMNS)

	# Write the DE tables
	WRITE_SIG_TABLES(DF = TEMP, NAME = MRNA_VAR, PVAL = "padj", ALPHA = 0.05,
			FC = "absFoldChange")
	WRITE_SIG_TABLES(DF = TEMP, NAME = MRNA_VAR, PVAL = "padj", ALPHA = 0.1,
			FC = "absFoldChange")
	WRITE_SIG_TABLES(DF = TEMP, NAME = MRNA_VAR, PVAL = "padj", ALPHA = 1,
			FC = "absFoldChange")
