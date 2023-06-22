#!/usr/bin/env Rscript
refGene <- read.table("gencode.vM25.annotation.bed")
#colnames(refGene) <- c("TXCHROM", "TXSTART", "TXEND", "GENEID", "TXNAME", "TXSTRAND")
#refGene <- refGene[grep("random|Un|hap", refGene$V1, invert=TRUE),]
refGene.excluded <- refGene[(refGene$V3-refGene$V2)<=1000,]
refGene <- refGene[(refGene$V3-refGene$V2)>500,]

#prepare regions of gene body
 
bodies <- refGene
bodies$V2[bodies$V6 == "+"] <- bodies$V2[bodies$V6 == "+"]+250
bodies$V3[bodies$V6 == "-"] <- bodies$V3[bodies$V6 == "-"]-250
bodies <-unique(bodies)
write.table(bodies, "Ref_bodies.bed", row.names=FALSE, quote=FALSE)
