#Plotting Manhattan plot, Dot Plot, and Violin Plot for FPKM values for four stages, namely mouseESCs, Lepto/Zygotene, Pachytene, and Diplotene
#Lastly plotting Gene density plot for all four stages.

setwd("/workdir/jl3285/GRO-seq_mESCs/ManhattanPlot_all")
#install.packages("readxl")
library(dplyr)
library(readxl)
library(ggplot2) 
#Sys.setenv(https_proxy="http://cbsulsrv13.biohpc.cornell.edu:3128")
library(ggpubr)
library(dplyr)
library(tidyverse)

Ref_bodies <- read.table("Ref_bodies.bed")
Ref_bodies <- Ref_bodies[,c(1,2,3)] 

###################### starting off with mouse Embryonic stem cell data from GRO-seq (adelman et al. 2015) ##################

upload_excel <- function(file_path) {
  
  data <- read.table(file_path, header = FALSE, col.names = c("chrom","start", "end",	"accession",	"mRNA_size",	"gene_strand",	"Frag_count",	"FPM",	"FPKM"))
  print(data[2,2])
  data <- data[,c(1,2,3,4,9)]
  return(data)
}

Merge_Replicates <- function(rep1,rep2) {
  #duplicated <- D_R1[duplicated(D_R1[,c("chrom", "start", "end", "accession")]),] #zero
  #df <-anti_join( D_R1, D_R2, by = c("chrom", "start", "end")) #zero
  df <- merge(rep1, rep2, by = c("chrom", "start", "end", "accession"))
  df$FPKM <- rowMeans(df[,5:6]) #calculating the average
  df <- df[,c(1,2,3,4,7)]
  return(df)
}

Plot_PrepareData <- function(data) {
  
  data$BP <- (data$end + data$start)/2 #midpoint of gene 
  data <- data[,c(4,1,6,5)]
  colnames(data) <- c("SNP", "CHR", "BP", "P") #P=FPKM, BP=midpoint position, CHR=chromosome, SNP=geneName
  data <- data %>%filter(CHR != "chrM")
  data <- data[-c(49663, 49664),] #this row has very high FPKM so I'm excluding this outlier
  
  #Arrange the dataset by chromosome number
  data$Number <- as.numeric(gsub("[^0-9]+", "", data$CHR)) 
  data$Number[data$CHR=="chrX"] <- 20
  data$Number[data$CHR=="chrY"] <- 21
  #data %>%mutate(NUMBERS = ifelse(CHR == "chrX", 20, NUMBERS))
  data <- data[order(data$Number, data$BP), ] #order via chromosome Number and the position
  #data$P <- log2(data$P)
  return(data)
}


myManhattan <- function(df, graph.title = "", 
                        col = c("lightblue", "navy"), even.facet = TRUE, 
                        font.size = 12, axis.size = 0.5, significance = NULL, report = FALSE,
                        inf.corr = 0.95, y.step = 10, point.size = 1){
  #myMin <- min(df$P[df$P != 0]) * inf.corr
  #df$P[df$P == 0] <- myMin
  require(ggplot2)
  require(stats)
  y.title <- expression((FPKM))
  
  if (length(col) > length(unique(df$CHR))){
    chrom.col <- col[1:length(unique(df$CHR))]
  } else if (!(length(col) > length(unique(df$CHR)))){
    chrom.col <- rep(col, length(unique(df$CHR))/length(col))
    if (length(chrom.col) < length(unique(df$CHR))){
      dif <- length(unique(df$CHR)) - length(chrom.col)
      #chrom.col <- c(chrom.col, col[1:dif])
      chrom.col <- c(chrom.col, "lightblue")
    }
  }
  y.max <- floor(max((df$P))) + 1
  
  if (y.max %% y.step != 0){ #checking the remainder
    y.diff <- floor(y.max/10)
    y.max <- y.diff*10 + 10 #since the y.step is 2
  }
  
  chromosome_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                        "chrX", "chrY")
  
  df$CHR <- factor(df$CHR, levels = chromosome_order)
  df <- df[order(df$CHR), ]
  #All_Density <- All_Density[order(All_Density$M_Frequency, decreasing = TRUE), ]
  #chrom <- sub("chr", "", All_Density$CHR)
  df$CHR <- factor(df$CHR,  labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19','X','Y'))
 
   g <- ggplot(df) +
    geom_point(aes(BP, (P), colour = as.factor(CHR)), size = point.size)
  
  g <- g + facet_grid(.~CHR, scale = "free_x", space = "free_x", switch = "x")
  
  
  g <- g + scale_colour_manual(values = chrom.col) +
    scale_y_continuous(expand = c(0, 0), limit = c(0, 40),
                       breaks = seq(from = 0, to = 40, by = 10)) +
    scale_x_continuous() +
    theme(strip.background = element_blank(), legend.position = "none",
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.spacing.x=unit(0.1, "lines"),
          axis.line.y = element_line(size = axis.size, color = "black"),
          axis.ticks.y = element_line(size = axis.size, color = "black"),
          axis.ticks.length = unit(axis.size * 10, "points"),
          plot.title = element_text(hjust = (0.5), size = font.size + 8),
          axis.title.y = element_text(size = font.size + 5),
          axis.title.x = element_text(size = font.size + 5),
          axis.text = element_text(size = font.size),
          strip.text.x = element_text(size = font.size))+
    labs(title = graph.title, x = "", y = y.title) 
    #geom_hline(yintercept = 1.242006278, linetype = "dashed", color = "black") #AUTOSOMAL MEAN +
   #geom_hline(yintercept = 0.581725095, linetype = "dashed", color = "red") #X mean
  
  return(g)
}

stage_list <- c("adelman_Rep1", "adelman_Rep2")

for(stage in stage_list) {
  print(stage)
  #stage="adelman_Rep1"
  pattern2=".FPKM.xls"
  object_name <- paste0(stage) #adelman_Rep1
  values <- paste0(stage, pattern2) #ChRO_D_R1.FPKM.FPKM.xls
  print(object_name)
  print(values)
  assign(object_name, upload_excel(values))  #calls the function upload_excel() to upload the files 
}

merge_data <- Merge_Replicates(adelman_Rep1, adelman_Rep1)

M <- Plot_PrepareData(merge_data)

save(M, file="mESCs_FPKM_ReplicateMerged.Rdata")

#Manhattan Plot:
mESC <- myManhattan(M, graph.title="Mouse Embryonic Stem cell")

########################## Spermatogenic stem cell data (Leptotene/Zygotene, Pachytene, Diplotene) from PRO-seq data (Alexander et al. 2023) ####################################################

Merge_Replicates <- function(rep1,rep2,rep3,rep4) {
  #duplicated <- D_R1[duplicated(D_R1[,c("chrom", "start", "end", "accession")]),] #zero
  #df <-anti_join( D_R1, D_R2, by = c("chrom", "start", "end")) #zero
  
  df <- merge(merge(merge(rep1, rep2, by = c("chrom", "start", "end", "accession"), all = TRUE), rep3, by = c("chrom", "start", "end", "accession"), all = TRUE), rep4, by = c("chrom", "start", "end", "accession"), all = TRUE)
  df$FPKM <- rowMeans(df[, 5:8]) #calculating the average
  df <- df[,c(1,2,3,4,9)]
  
  return(df)
}

Plot_PrepareData <- function(data, filename) {
  
  data$BP <- (data$end + data$start)/2 #midpoint of gene 
  
  if (filename == "P") {
    data$FPKM <- data$FPKM * 2.7153 # Normalization: multiplying by radioactive run-on scale factor
  } else if (filename == "D") {
    data$FPKM <- data$FPKM * 1.2474
  } else if (filename == "LZ") {
    data$FPKM <- data$FPKM * 1
  }
  
  data <- data[,c(4,1,6,5)]
  colnames(data) <- c("SNP", "CHR", "BP", "P") #P=FPKM, BP=midpoint position, CHR=chromosome, SNP=geneName
  data <- data %>%filter(CHR != "chrM")
  #data <- data[-113543,] #this row has very high FPKM so I'm excluding this outlier
  
  #Arrange the dataset by chromosome number
  data$Number <- as.numeric(gsub("[^0-9]+", "", data$CHR)) 
  data$Number[data$CHR=="chrX"] <- 20
  data$Number[data$CHR=="chrY"] <- 21
  #data %>%mutate(NUMBERS = ifelse(CHR == "chrX", 20, NUMBERS))
  data <- data[order(data$Number, data$BP), ] #order via chromosome Number and the position
  #data$P <- log2(data$P)
  return(data)
}

myManhattan <- function(df, graph.title = "", 
                        col = c("lightblue", "navy"), even.facet = TRUE, 
                        font.size = 12, axis.size = 0.5, significance = NULL, report = FALSE,
                        inf.corr = 0.95, y.step = 10, point.size = 1){
  #myMin <- min(df$P[df$P != 0]) * inf.corr
  #df$P[df$P == 0] <- myMin
  require(ggplot2)
  require(stats)
  y.title <- expression((FPKM))
  
  if (length(col) > length(unique(df$CHR))){
    chrom.col <- col[1:length(unique(df$CHR))]
  } else if (!(length(col) > length(unique(df$CHR)))){
    chrom.col <- rep(col, length(unique(df$CHR))/length(col))
    if (length(chrom.col) < length(unique(df$CHR))){
      dif <- length(unique(df$CHR)) - length(chrom.col)
      #chrom.col <- c(chrom.col, col[1:dif])
      chrom.col <- c(chrom.col, "lightblue")
    }
  }
  y.max <- floor(max((df$P))) + 1
  
  if (y.max %% y.step != 0){ #checking the remainder
    y.diff <- floor(y.max/10)
    y.max <- y.diff*10 + 10 #since the y.step is 10
  }
  
  chromosome_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                       "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                       "chrX", "chrY")
  
  #D$chrom <- sub("chr", "", D$CHR)
  #D$chrom <- factor(D$chrom, levels = unique(D$chrom),labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19','X','Y'))
  
  df$CHR <- factor(df$CHR, levels = chromosome_order)
  df <- df[order(df$CHR), ]
  #All_Density <- All_Density[order(All_Density$M_Frequency, decreasing = TRUE), ]
  #chrom <- sub("chr", "", All_Density$CHR)
  df$CHR <- factor(df$CHR,  labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19','X','Y'))
  
  g <- ggplot(df) +
    geom_point(aes(BP, (P), colour = as.factor(CHR)), size = point.size)
  
  g <- g + facet_grid(.~CHR, scale = "free_x", space = "free_x", switch = "x")
  
  g <- g + scale_colour_manual(values = chrom.col) +
    scale_y_continuous(expand = c(0, 0), limit = c(0, 30),
                       breaks = seq(from = 0, to = 30, by = y.step)) +
    scale_x_continuous() +
    theme(strip.background = element_blank(), legend.position = "none",
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.spacing.x=unit(0.1, "lines"),
          axis.line.y = element_line(size = axis.size, color = "black"),
          axis.ticks.y = element_line(size = axis.size, color = "black"),
          axis.ticks.length = unit(axis.size * 10, "points"),
          plot.title = element_text(hjust = (0.5), size = font.size + 8),
          axis.title.y = element_text(size = font.size + 5),
          axis.title.x = element_text(size = font.size + 5),
          axis.text = element_text(size = font.size),
          strip.text.x = element_text(size = font.size))+
    labs(title = graph.title, x = "Chromosome", y = y.title) 
    #geom_hline(yintercept = autosome_mean, linetype = "dashed", color = "black") #AUTOSOMAL MEAN +
    #geom_hline(yintercept = X_mean, linetype = "dashed", color = "red") #X mean
  
  return(g)
}

stage_list <- c("ChRO_D", "ChRO_P", "ChRO_LZ")

for(stage in stage_list) {
  print(stage)
  #stage="ChRO_D"
  pattern2=".FPKM.FPKM.xls"
  
  for(i in 1:4) { 
    object_name <- paste0(paste(stage, "R", sep = "_"), i)
    values <- paste0(paste(stage, "R", sep = "_"), i, pattern2) #ChRO_D_R1.FPKM.FPKM.xls
    print(object_name)
    assign(object_name, upload_excel(values))  #calls the function upload_excel() to upload the files 
  }
}

D <- Merge_Replicates(ChRO_D_R1,ChRO_D_R2,ChRO_D_R3,ChRO_D_R4)
P <- Merge_Replicates(ChRO_P_R1,ChRO_P_R2,ChRO_P_R3,ChRO_P_R4)
LZ <- Merge_Replicates(ChRO_LZ_R1,ChRO_LZ_R2,ChRO_LZ_R3,ChRO_LZ_R4)

D <- Plot_PrepareData(D, "D")
P <- Plot_PrepareData(P, "P")
LZ <- Plot_PrepareData(LZ, "LZ")

#this row has very high FPKM so I'm excluding this gene for convineience 
D <- D[-113543,] #gene is ERF1, ENSMUST00000205406.1
P <- P[-113543,]
NUMBER <- which(P$SNP == "ENSMUST00000198477.2")
P <- P[-NUMBER,]
LZ <- LZ[-113543,]

save(list = c("LZ", "P", "D", "Ref_bodies"), file = "SCs_FPKM_ReplicateMerged.Rdata")

#Manhattan Plot
Diplo <- myManhattan(D, graph.title="Diplotene")
Pachy <- myManhattan(P, graph.title="Pachytene")
LeptoZygo <- myManhattan(LZ,  graph.title="Lepto/Zygotene")

#Arrange all the graph in a single plot.
Plot_all <- ggarrange(mESC, LeptoZygo, Pachy, Diplo, ncol = 2, nrow = 2)

####################################################################################
#Removing zeros and NAs from the data.
load(file ="mESCs_FPKM_ReplicateMerged.Rdata")
load(file ="SCs_FPKM_ReplicateMerged.Rdata")
D$P[D$P==0] <- NA 
D <- D[!is.na(D$P),] #106420
D <- D[D$P>=1,]  #52043
Ddot <- D

P$P[P$P==0] <- NA 
P <- P[!is.na(P$P),] #106420
P <- P[P$P>=1,]  #52043
Pdot <- P

LZ$P[LZ$P==0] <- NA 
LZ <- LZ[!is.na(LZ$P),] #106420
LZ <- LZ[LZ$P>=1,]  #52043
LZdot <- LZ

M$P[M$P==0] <- NA 
M <- M[!is.na(M$P),] #106420
M <- M[M$P>=1,]  #52043
Mdot <- M

### Dot plot and bar plot of median and quartiles################################################
#D_summary <- Ddot %>%  group_by(CHR) %>% summarize(med = median(P), Q1 = quantile(P, 0.25), Q3 = quantile(P, 0.75))
#P_summary <- Pdot %>%  group_by(CHR) %>% summarize(med = median(P), Q1 = quantile(P, 0.25), Q3 = quantile(P, 0.75))
#LZ_summary <-LZdot %>%  group_by(CHR) %>% summarize(med = median(P), Q1 = quantile(P, 0.25), Q3 = quantile(P, 0.75))
#M_summary <-Mdot %>%  group_by(CHR) %>% summarize(med = median(P), Q1 = quantile(P, 0.25), Q3 = quantile(P, 0.75))

D_summary <- Ddot %>%  group_by(CHR) %>% summarize(med = median(P), Q1 = quantile(P, 0.25), Q3 = quantile(P, 0.90))
P_summary <- Pdot %>%  group_by(CHR) %>% summarize(med = median(P), Q1 = quantile(P, 0.25), Q3 = quantile(P, 0.90))
LZ_summary <-LZdot %>%  group_by(CHR) %>% summarize(med = median(P), Q1 = quantile(P, 0.25), Q3 = quantile(P, 0.90))
M_summary <-Mdot %>%  group_by(CHR) %>% summarize(med = median(P), Q1 = quantile(P, 0.25), Q3 = quantile(P, 0.90))

M_Y <- c("chrY", NA, NA, NA)
M_summary <- rbind(M_summary, M_Y)

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", size = (25), hjust = 0.5), 
                      legend.title = element_blank(), 
                      legend.text = element_blank(), 
                      axis.title = element_text(family = "Helvetica", size = (20), colour = "black"),
                      axis.text = element_text(family = "Helvetica", colour = "black", size = (16)))

# Define a custom color palette with alternating colors
col <- c("dodgerblue4", "deepskyblue2")
unique_chrs <- unique(D_summary$CHR)
chrom.col <- rep(col, length.out = length(unique_chrs))

#Diplotene data
D_summary$chrom <- sub("chr", "", D_summary$CHR)
chromosome_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                      "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                      "chrX", "chrY")

D_summary$CHR <- factor(D_summary$CHR, levels = chromosome_order)
D_summary <- D_summary[order(D_summary$CHR), ]
D_summary$CHR <- factor(D_summary$CHR,  labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19','X','Y'))

#Generating a dot plot from median, Quartile 25, and Quartile 75
D_dot <- ggplot(D_summary, aes(x = CHR, colour = as.factor(CHR))) +
  geom_pointrange(aes(y = med, ymin = Q1, ymax = Q3), show.legend = FALSE) +
  geom_point(aes(y = med), shape = 16, size = 2) +
  labs(y = 'FPKM', title = 'Diplotene', x = '') +
  theme_bw() + mynamestheme +
  scale_colour_manual(values = chrom.col) + theme(legend.position = "none") 
  #theme(panel.grid.major = element_blank())

#####################
#Pachytene data
P_summary$chrom <- sub("chr", "", P_summary$CHR)
P_summary$CHR <- factor(P_summary$CHR, levels = chromosome_order)
P_summary <- P_summary[order(P_summary$CHR), ]
P_summary$CHR <- factor(P_summary$CHR,  labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19','X','Y'))

#Generating a dot plot from median, Quartile 25, and Quartile 75
P_dot <- ggplot(P_summary, aes(x = CHR, colour = as.factor(CHR))) +
  geom_pointrange(aes(y = med, ymin = Q1, ymax = Q3), show.legend = FALSE) +
  geom_point(aes(y = med), shape = 16, size = 2) +
  labs(y = 'FPKM', title = 'Pachytene', x = '') +
  theme_bw() + mynamestheme +
  scale_colour_manual(values = chrom.col) + theme(legend.position = "none")

#####################
#Leptotene/Zygotene data
LZ_summary$chrom <- sub("chr", "", LZ_summary$CHR)
LZ_summary$CHR <- factor(LZ_summary$CHR, levels = chromosome_order)
LZ_summary <- LZ_summary[order(LZ_summary$CHR), ]
LZ_summary$CHR <- factor(LZ_summary$CHR,  labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19','X','Y'))

#Generating a dot plot from median, Quartile 25, and Quartile 75
LZ_dot <- ggplot(LZ_summary, aes(x = CHR, colour = as.factor(CHR))) +
  geom_pointrange(aes(y = med, ymin = Q1, ymax = Q3), show.legend = FALSE) +
  geom_point(aes(y = med), shape = 16, size = 2) +
  labs(y = 'FPKM', title = 'Leptotene/Zygotene', x = '') +
  theme_bw() + mynamestheme +
  scale_colour_manual(values = chrom.col) + theme(legend.position = "none")

#####################
#Mouse Embryonic stem cell data

M_summary$chrom <- sub("chr", "", M_summary$CHR)
M_summary$CHR <- factor(M_summary$CHR, levels = chromosome_order)
M_summary <- M_summary[order(M_summary$CHR), ]
M_summary$CHR <- factor(M_summary$CHR,  labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19','X','Y'))
range(M_summary$Q1)
M_summary$med <- as.numeric(M_summary$med)
M_summary$Q1 <- as.numeric(M_summary$Q1)
M_summary$Q3 <- as.numeric(M_summary$Q3)

#Generating a dot plot from median, Quartile 25, and Quartile 75
M_dot <- ggplot(M_summary, aes(x = CHR, colour = as.factor(CHR))) +
  geom_pointrange(aes(y = med, ymin = Q1, ymax = Q3), show.legend = FALSE) +
  geom_point(aes(y = med), shape = 16, size = 2) +
  labs(y = 'FPKM', title = 'mouseESCs', x = '') + scale_y_continuous(breaks=seq(1,11,by=2)) +
  theme_bw() + mynamestheme +
  scale_colour_manual(values = chrom.col) + theme(legend.position = "none")

#Arranging all graph in a single plot.
Plot_dot <- ggarrange(M_dot, LZ_dot, P_dot, D_dot, ncol = 2, nrow = 2)
Plot_dot

######################################################################################
#For gene density of transcribed genes (FPKM >= 1)

D_Density <- table(D$CHR) %>% data.frame()
LZ_Density <- table(LZ$CHR) %>% data.frame()
P_Density <- table(P$CHR) %>% data.frame()
M_Density <- table(M$CHR) %>% data.frame()
# Add a new factor level "chrY" to the data frame
M_Density$Var1 <- factor(M_Density$Var1, levels = unique(c(levels(M_Density$Var1), "chrY")))
# Create a new row for "chrY" and add it to the data frame
gene <- data.frame(Var1 = "chrY", Freq = 0)
M_Density <- rbind(M_Density, gene)

colnames(D_Density) <- c("CHR", "D_Frequency")
colnames(LZ_Density)  <- c("CHR", "LZ_Frequency")
colnames(P_Density)  <- c("CHR", "P_Frequency")
colnames(M_Density)  <- c("CHR", "M_Frequency")

All_Density <- merge(M_Density, LZ_Density, by = "CHR")
All_Density <- merge(All_Density, P_Density, by = "CHR")
All_Density <- merge(All_Density, D_Density, by = "CHR")

col <- c("lightblue", "navy")

unique_chrs <- unique(D_summary$CHR)
chrom.col <- rep(col, length.out = length(unique_chrs))

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", size = (25), hjust = 0.5), 
                      legend.title = element_blank(), 
                      legend.text = element_blank(), 
                      axis.title = element_text(family = "Helvetica", size = (20), colour = "black"),
                      axis.text = element_text(family = "Helvetica", colour = "black", size = (16)))

chromosome_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                      "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                      "chrX", "chrY")
All_Density$chrom <- sub("chr", "", All_Density$CHR)
All_Density$CHR <- factor(All_Density$CHR, levels = chromosome_order)
All_Density <- All_Density[order(All_Density$CHR), ]
All_Density$CHR <- factor(All_Density$CHR,  labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19','X','Y'))

#Density plot
MDensityPlot <- ggplot(All_Density, aes(x = CHR, y = M_Frequency, fill = as.factor(CHR))) +
  geom_bar(stat = "identity") +
  labs(y = "Gene Density", x = "", title = "mESCs") +
  theme_bw() + mynamestheme + scale_y_continuous(breaks=seq(0,4000,by=1000)) +
  scale_fill_manual(values = chrom.col) + theme(legend.position = "none")

LZDensityPlot <- ggplot(All_Density, aes(x = CHR, y = LZ_Frequency, fill = as.factor(CHR))) +
  geom_bar(stat = "identity") +
  labs(y = "Gene Density", x = "", title = "Lepto/Zygotene") +
  theme_bw() + mynamestheme + scale_y_continuous(breaks=seq(0,4000,by=1000)) +
  scale_fill_manual(values = chrom.col) + theme(legend.position = "none")

PDensityPlot <- ggplot(All_Density, aes(x = CHR, y = P_Frequency, fill = as.factor(CHR))) +
  geom_bar(stat = "identity") +
  labs(y = "Gene Density", x = "", title = "Pachytene") +
  theme_bw() + mynamestheme + scale_y_continuous(breaks=seq(0,6200,by=1000)) +
  scale_fill_manual(values = chrom.col) + theme(legend.position = "none")

DDensityPlot <- ggplot(All_Density, aes(x = CHR, y = D_Frequency, fill = as.factor(CHR))) +
  geom_bar(stat = "identity") +
  labs(y = "Gene Density", x = "", title = "Diplotene") +
  theme_bw() + mynamestheme + scale_y_continuous(breaks=seq(0,4000,by=1000)) +
  scale_fill_manual(values = chrom.col) + theme(legend.position = "none")

Plot_density <- ggarrange(MDensityPlot, LZDensityPlot, PDensityPlot, DDensityPlot, ncol = 2, nrow = 2)
Plot_density

#For gene density of all genes 
colnames(Ref_bodies)<- c("chrom","start", "end")
Ref_bodies_Density <- table(Ref_bodies$chrom) %>% data.frame()
colnames(Ref_bodies_Density) <- c("CHR", "Freq")
chromosome_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                      "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                      "chrX", "chrY", "chrM")

Ref_bodies_Density$CHR <- factor(Ref_bodies_Density$CHR, levels = chromosome_order)
Ref_bodies_Density <- Ref_bodies_Density[order(Ref_bodies_Density$CHR), ]

Ref_bodies_Density$CHR <- factor(Ref_bodies_Density$CHR, labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19','X','Y','M'))
unique_chrs <- unique(Ref_bodies_Density$CHR)
chrom.col <- rep(col, length.out = length(unique_chrs))

Ref_DensityPlot <- ggplot(Ref_bodies_Density, aes(x = CHR, y = Freq, fill = as.factor(CHR))) +
  geom_bar(stat = "identity") +
  labs(y = "Gene Density", x = "", title = "All Gene Density") +
  theme_bw() + mynamestheme  +
  scale_fill_manual(values = chrom.col) + theme(legend.position = "none")

#####################################################################

# Subsetting ChrX and Chr10 (having similar gene density) and comparing FPKM:
col <- c("dodgerblue4", "deepskyblue2")
M_subset <- subset(M_summary, CHR == "X"|  CHR == 18)

M_subset_plot <- ggplot(M_subset, aes(x = CHR, colour = as.factor(CHR))) +
  geom_pointrange(aes(y = med, ymin = Q1, ymax = Q3), show.legend = FALSE) +
  geom_point(aes(y = med), shape = 16, size = 2) +  scale_y_continuous(breaks=seq(0,5,by=1))+
  labs(y = 'FPKM', title = 'mouseESCs', x = '')  +
  theme_bw() + mynamestheme +
  scale_colour_manual(values = col) + theme(legend.position = "none")

# Subsetting ChrX and Chr10 (having similar gene density):
LZ_subset <- subset(LZ_summary, CHR == "X"|  CHR == 18)

LZ_subset_plot <- ggplot(LZ_subset, aes(x = CHR, colour = as.factor(CHR))) +
  geom_pointrange(aes(y = med, ymin = Q1, ymax = Q3), show.legend = FALSE) +
  geom_point(aes(y = med), shape = 16, size = 2) +
  labs(y = 'FPKM', title = 'Lepto/Zygotene', x = '')  +
  theme_bw() + mynamestheme +
  scale_colour_manual(values = col) + theme(legend.position = "none")

# Subsetting ChrX and Chr10 (having similar gene density):
P_subset <- subset(P_summary, CHR == "X"|  CHR == 18)

P_subset_plot <- ggplot(P_subset, aes(x = CHR, colour = as.factor(CHR))) +
  geom_pointrange(aes(y = med, ymin = Q1, ymax = Q3), show.legend = FALSE) +
  geom_point(aes(y = med), shape = 16, size = 2) +
  labs(y = 'FPKM', title = 'Pachytene', x = '')  +
  theme_bw() + mynamestheme +
  scale_colour_manual(values = col) + theme(legend.position = "none")

# Subsetting ChrX and Chr10 (having similar gene density):
D_subset <- subset(D_summary, CHR == "X"|  CHR == 18)

D_subset_plot <- ggplot(D_subset, aes(x = CHR, colour = as.factor(CHR))) +
  geom_pointrange(aes(y = med, ymin = Q1, ymax = Q3), show.legend = FALSE) +
  geom_point(aes(y = med), shape = 16, size = 2) +  scale_y_continuous(breaks=seq(0,5,by=1))+
  labs(y = 'FPKM', title = 'Diplotene', x = '')  +
  theme_bw() + mynamestheme +
  scale_colour_manual(values = col) + theme(legend.position = "none")

Plot_subset <- ggarrange(M_subset_plot, LZ_subset_plot, P_subset_plot, D_subset_plot, ncol = 2, nrow = 2)
Plot_subset

####################################################################################
#Empirical cumulative density function (ECDF)
mynamestheme <- theme(plot.title = element_text(family = "Helvetica", size = (25), hjust = 0.5), 
                      legend.title = element_blank(), 
                      legend.text = element_text(family = "Helvetica", size = (15), colour = "black"), 
                      axis.title = element_text(family = "Helvetica", size = (20), colour = "black"),
                      axis.text = element_text(family = "Helvetica", colour = "black", size = (16)))
my_colors = c("green", "darkgreen", "forestgreen", "seagreen", "mediumseagreen", "limegreen", "springgreen", "yellowgreen", "chartreuse", "darkorchid2",
              "blue", "darkblue", "mediumblue", "navyblue", "royalblue", "dodgerblue", "steelblue", "palegreen", "lightblue", "red", "gold1")

ECDF_Diplo <- ggplot(D, aes(x=P, col=CHR)) + labs(y = '\n', title = 'Diplotene', x = '\nFPKM')  +
  # stat_ecdf() function is used to plot ECDF plot
  stat_ecdf() + theme_bw() + mynamestheme +scale_color_manual(values = my_colors)+ xlim(-1,25)

ECDF_Pachy <-ggplot(P, aes(x=P, col=CHR)) + labs(y = 'emperical\n cumulative density function\n (ecdf)\n', title = 'Pachytene', x = '\nFPKM')  +
  # stat_ecdf() function is used to plot ECDF plot
  stat_ecdf() + theme_bw() + mynamestheme +scale_color_manual(values = my_colors)+ xlim(-1,50) 

ECDF_LZ <-ggplot(LZ, aes(x=P, col=CHR)) + labs(y = '\n', title = 'Lepto/Zygotene', x = '')  +
  # stat_ecdf() function is used to plot ECDF plot
  stat_ecdf() + theme_bw() + mynamestheme +scale_color_manual(values = my_colors)+ xlim(-1,25) 

ECDF_M <- ggplot(M, aes(x=P, col=CHR)) + labs(y = 'emperical\n cumulative density function\n (ecdf)\n', title = 'mESCs', x = '')  +
  # stat_ecdf() function is used to plot ECDF plot
  stat_ecdf() + theme_bw() + mynamestheme +scale_color_manual(values = my_colors)+ xlim(-1,25)

ggarrange(ECDF_M, ECDF_LZ, ECDF_Pachy, ECDF_Diplo, nrow=2, ncol=2, common.legend = TRUE, legend = "right") 
#facet_grid(~chr,scales="free")

############################Violin Plot ########################################
# Create a violin plot using ggplot2
mynamestheme <- theme(plot.title = element_text(family = "Helvetica", size = (25), hjust = 0.5), 
                      legend.title = element_blank(), 
                      legend.text = element_blank(), 
                      axis.title = element_text(family = "Helvetica", size = (20), colour = "black"),
                      axis.text = element_text(family = "Helvetica", colour = "black", size = (16)))

# Define a custom color palette with alternating colors
col <- c("dodgerblue4", "deepskyblue2")
unique_chrs <- unique(D$CHR)
chrom.col <- rep(col, length.out = length(unique_chrs))

D$chrom <- sub("chr", "", D$CHR)
D$chrom <- factor(D$chrom, levels = unique(D$chrom),labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19','X','Y'))

Violin_D<- ggplot(D, aes(x = chrom, y = P, fill = as.factor(chrom))) +
  geom_violin() +
  labs(title = "Diplotene", x = "", y = "FPKM")  +
  theme_bw() + mynamestheme + ylim(0,70) +
  scale_colour_manual(values = chrom.col) + theme(legend.position = "none") 

LZ$chrom <- sub("chr", "", LZ$CHR)
LZ$chrom <- factor(LZ$chrom, levels = unique(LZ$chrom),labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19','X','Y'))

Violin_LZ<- ggplot(LZ, aes(x = chrom, y = P, fill = as.factor(chrom))) +
  geom_violin() +
  labs(title = "Lepto/Zygotene", x = "", y = "FPKM")  + ylim(0,70) +
  theme_bw() + mynamestheme +
  scale_colour_manual(values = chrom.col) + theme(legend.position = "none") 

P$chrom <- sub("chr", "", P$CHR)
P$chrom <- factor(P$chrom, levels = unique(P$chrom),labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19','X','Y'))

Violin_P<- ggplot(P, aes(x = chrom, y = P, fill = as.factor(chrom))) +
  geom_violin() +
  labs(title = "Pachytene", x = "", y = "FPKM") + ylim(0,70) +
  theme_bw() + mynamestheme +
  scale_colour_manual(values = chrom.col) + theme(legend.position = "none") 

#M <- rbind(M, c(NA,"chrY",NA,NA,NA,"Y"))
M$chrom <- sub("chr", "", M$CHR)
M$chrom <- factor(M$chrom, levels = unique(M$chrom),labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19','X'))

Violin_M <- ggplot(M, aes(x = chrom, y = P, fill = as.factor(chrom))) +
  geom_violin() +
  labs(title = "mESCs", x = "", y = "FPKM") +
  theme_bw() + mynamestheme +  scale_y_continuous(breaks=seq(0,253,by=50)) + ylim(0,70) +
  scale_colour_manual(values = chrom.col) + theme(legend.position = "none") 

ViolinAll <- ggarrange(Violin_M, Violin_LZ, Violin_P, Violin_D, nrow = 2, ncol = 2)

#################################################################################################
#
