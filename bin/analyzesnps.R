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
library(GenomicRanges)
library(seqinr)

GENOME.class <- readVCF(vcf.file, numcols=num.cols, tid=t.id, frompos=from.pos, topos=to.pos,
		       	include.unknown=TRUE, approx=FALSE, gffpath=gff.file)

GENOME.class <- set.synnonsyn(GENOME.class, ref.chr=ref.file, save.codons=TRUE)

# Splitting data into "genes" regions:
genes <- splitting.data(GENOME.class, subsites = "gene")
genes.size <- length(genes@region.names)
genes.loc <- genes@region.names

# Obtaining gene coordinates:
start <- sapply(genes.loc, function(x)
{ return(as.numeric(strsplit(x, " ") [[1]][1]))})

end <- sapply(genes.loc, function(x)
{ return(as.numeric(strsplit(x, " ") [[1]][3]))})

# Gene lengths:
length <- end - start + 1
pos.data <- cbind(start, end)

# Calculate syn, non-syn, and Mu for each gene: 
syn <- rep(NA, genes.size)
nonsyn <- rep(NA, genes.size)
ka_ks <- rep(NA, genes.size)
m_u <- rep(NA, genes.size)
for (i in 1:genes.size) {
  syn[i] <- sum(genes@region.data@synonymous[[i]]==1,na.rm=TRUE)
  nonsyn[i] <- sum(genes@region.data@synonymous[[i]]==0,na.rm=TRUE)
  ka_ks[i] <- nonsyn[i]/syn[i]
  m_u[i] <- ((syn[i]+nonsyn[i]) * 1000) / length[i]
}
ka.ks <- round(ka_ks, digits=2)
mu <- round(m_u, digits=2)

# Extracting info from the gff3 file:
gff.attr <- get.feature.names(genes, gff.file=gff.file, chr=t.id)
gff.attr.split <- strsplit(gff.attr, "[=;]+")
feature <- sapply(gff.attr.split, "[[", 12)
gene.ID <- sapply(gff.attr.split, "[[", 6)

# Extracting GC contents for the targeted genes:
gr <- GRanges(
	      seqnames = rep(t.id, genes.size),
	      ranges = IRanges(start=start, end=end)
	      )

chrom.list <- c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06",
		"chr07", "chr08", "chr09", "chr10", "chr11", "chr12")
idx <- match(t.id, chrom.list)
fas <- read.fasta(file=ref.file)[[idx]]

extract.seqs <- lapply(gr, function(x) {
    if (as.character(strand(x)) == '-') {
        comp(fas[end(x):start(x)])
      } else {
	 fas[start(x):end(x)]}})

gc <- lapply(extract.seqs, function(x)
{ round(GC(x), digits=2) })
GC <- unlist(gc) 

#---------------------------------------------------------------
## Compiling results:
row.num <- rep(1:genes.size)
res <- data.frame(gene.ID, pos.data, length, syn, nonsyn, ka.ks, GC, mu, feature, row.names=row.num)
write.table(res, file=out.file, sep="\t", row.names=FALSE)

#---------------------------------------------------------------
#
