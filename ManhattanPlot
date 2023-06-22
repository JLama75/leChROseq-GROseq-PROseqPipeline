#Generating Manhattan plot from ChRO-seq data using FPKM to visualize the gene expression in each chromosome.
#Plot to compare genome-wide gene expression between three different sub-stages of Prophase I during spermatogenesis in mice (Leptotene/Zygotene, Pachytene, and Diplotene).
#Reference genome used- mm10

setwd("/local/workdir/jl3285/ChRO_seq/mm10geneAnnotation/FPKM")
#install.packages("readxl")
library(dplyr)
library(readxl)

Ref_bodies <- read.table("Ref_bodies.bed") #Bed file containing the read counts overlapping gene bodies 
Ref_bodies <- Ref_bodies[,c(1,2,3)] 
colnames(Ref_bodies)<- c("chrom","start", "end")

upload_excel <- function(file_path) {
  
  data <- read.table(file_path, header = FALSE, col.names = c("chrom","start", "end",	"accession",	"mRNA_size",	"gene_strand",	"Frag_count",	"FPM",	"FPKM"))
  print(data[2,2])
  data <- data[,c(1,2,3,4,9)]
  return(data)
}


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
  
  df$CHR <- factor(df$CHR, levels = unique(df$CHR), labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10','11', '12', '13', '14', '15', '16', '17', '18', '19', 'X', 'Y'))
  
  g <- ggplot(df) +
    geom_point(aes(BP, (P), colour = as.factor(CHR)), size = point.size)
  
  g <- g + facet_grid(.~CHR, scale = "free_x", space = "free_x", switch = "x")

  
  g <- g + scale_colour_manual(values = chrom.col) +
    scale_y_continuous(expand = c(0, 0), limit = c(0, y.max),
                       breaks = seq(from = 0, to = y.max, by = y.step)) +
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
LZ <- LZ[-113543,]

dev.off()
myManhattan(D, graph.title="Diplotene")
dev.off()

# Specify the file path and name for the PDF output
output_file <- "Manhattan_Diplotene_FPKM.pdf"
# Save the plot to PDF using ggsave
ggsave(output_file, width = 9, height = 4, dpi = 300)

dev.off()
myManhattan(P, graph.title="Pachytene", y.step = 50)
# Specify the file path and name for the PDF output
output_file <- "Manhattan_Pachytenee_FPKM.pdf"
# Save the plot to PDF using ggsave
ggsave(output_file, width = 9, height = 4, dpi = 300)

dev.off()
myManhattan(LZ,  graph.title="Lepto/Zygotene")
# Specify the file path and name for the PDF output
output_file <- "Manhattan_LeptoZygotene_FPKM.pdf"

# Save the plot to PDF using ggsave
ggsave(output_file, width = 9, height = 4, dpi = 300)

#final plot dimension: 
#plot dimention: 12.74X4.08 for R1
