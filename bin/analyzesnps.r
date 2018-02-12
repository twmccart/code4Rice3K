#!/usr/bin/env Rscript

## R script to calculate synonymous and nonsynonymous SNPs
## in the genome of Rice

## Chromosomes and and their coordinates are as follows:
  ## chr01      1 - 43270923
  ## chr02      1 - 35937250
  ## chr03      1 - 36413819
  ## chr04      1 - 35502694
  ## chr05      1 - 29958434
  ## chr06      1 - 31248787
  ## chr07      1 - 29697621
  ## chr08      1 - 28443022
  ## chr09      1 - 23012720
  ## chr10      1 - 23207287
  ## chr11      1 - 29021106
  ## chr12      1 - 27531856

#---------------------------------------------------------------
## Parsing command line arguments:

args <- commandArgs(TRUE)
vcf.file <- as.character(args[1])
ref.file <- as.character(args[2])
gff.file <- as.character(args[3])
out.file <- as.character(args[4])
num.cols <- as.numeric(args[5])
t.id <- as.character(args[6])
from.pos <- as.numeric(args[7])
to.pos <- as.numeric(args[8])

#---------------------------------------------------------------
# Code for popGenome package

library(PopGenome)

GENOME.class <- readVCF(vcf.file, numcols=num.cols, tid=t.id, frompos=from.pos, topos=to.pos,
		       	include.unknown=TRUE, approx=FALSE, gffpath=gff.file)

GENOME.class <- set.synnonsyn(GENOME.class, ref.chr=ref.file, save.codons=TRUE)

genes <- splitting.data(GENOME.class, subsites = "gene")
genesSize <- length(genes@region.names)
genesLOC <- genes@region.names

syn <- rep(NA, genesSize)
nonsyn <- rep(NA, genesSize)
ka_ks <- rep(NA, genesSize)
for (i in 1:genesSize) {
  syn[i] <- sum(genes@region.data@synonymous[[i]]==1,na.rm=TRUE)
  nonsyn[i] <- sum(genes@region.data@synonymous[[i]]==0,na.rm=TRUE)
  ka_ks[i] <- nonsyn[i]/syn[i]
}
ka.ks <- round(ka_ks, digits=2)

# Obtaining gene coordinates:
start <- sapply(genesLOC, function(x)
{ return(as.numeric(strsplit(x, " ") [[1]][1]))})

end <- sapply(genesLOC, function(x)
{ return(as.numeric(strsplit(x, " ") [[1]][3]))})

# Gene lengths:
length <- end - start + 1
posData <- cbind(start, end)

# gff <- read.delim("gff/IRGSP-1.0.gff3", header=F, comment.char="#",
#                   row.names=NULL)
# gff.info <- subset(gff2, V1 == "chr01" & V3 == "gene")

# Extracting info from the gff3 file:
gff.attr <- get.feature.names(genes, gff.file=gff.file, chr=t.id)
gff.attr.split <- strsplit(gff.attr, "[=;]+")
feature <- sapply(gff.attr.split, "[[", 12)
gene.ID <- sapply(gff.attr.split, "[[", 6)
#

#---------------------------------------------------------------
## Compiling results:
row.num <- rep(1:genesSize)
res <- data.frame(gene.ID, posData, length, syn, nonsyn, ka.ks, feature, row.names=row.num)
write.table(res, file=out.file, sep="\t", row.names=FALSE)

#-----------------------------------------------------------Finn

