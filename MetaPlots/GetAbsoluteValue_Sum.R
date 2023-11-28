#!/usr/bin/env Rscript
cat("Running R to generate absolute value!\n")

df<- read.table("tmp.merge.bedGraph")
df$V4 <- abs(df$V4)
df$V6 <- df$V4 + df$V5
df <- df[,c(1,2,3,6)]

rownames(df) <- NULL
colnames(df) <- NULL
write.table(df, "Output.bedGraph")
