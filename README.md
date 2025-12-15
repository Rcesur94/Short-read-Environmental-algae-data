#Project overview 







# Environmental setup 



# R script 
#load required libraries
library(ggplot2)

#Read the table into R 
df <- read.table("marker_hits.tsv",
                 header = TRUE,
                 sep = "\t",
stringsAsFactors = FALSE)

# check table to see this is what you think it should look like 
head(df)

#extract contig length 
df$Contig_length <- as.numeric(
    sub(".*length_([0-9]+)_.*", "\\1", df$Contig)
 )

#Convert E-Value to base log 10
df$log10_Evalue <- -log10(as.numeric(df$E_value))

#Plot graph
ggplot(df, aes(x = Contig_length, y = log10_Evalue, color = Hit)) +
   geom_point(size = 2, alpha = 0.8) +
   theme_minimal() +
   labs(
         x = "Contig length (bp)",
        y = expression(-log[10](E-value)),
title = "HMM Marker vs Viral Contig Length"
    )

<img width="1103" height="373" alt="image" src="https://github.com/user-attachments/assets/a6d83475-16f1-4a06-a748-16ca97ffc9c6" />

