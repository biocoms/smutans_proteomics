# This script details how the reduced protein set was obtained from the raw Proteome Discoverer output.

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(proteinDiscover)
library(DBI)
library(vsn)

##### Violin plot function #####
violin <- function(dat, title) {
  violin_dat <- pivot_longer(dat, 
                             cols = starts_with("Abundance"), 
                             names_to = "Sample") %>%
    mutate(Sample = gsub("Abundance", "", Sample)) 
  violin_dat$Sample <- factor(violin_dat$Sample, 
                              levels = unique(violin_dat$Sample))
  
  ggplot(violin_dat, aes(x = Sample, y = value, fill = Sample)) +
    geom_violin() +
    ggtitle(title) + 
    labs(x = "Sample", y = "Abundance") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = "none")
}

###############################

ProtDiscProcess <- function(filename, sample) {
  temppd <- dbOpen(filename)
  finalpd <- dbGetProteinTable(temppd, columns = c("Accession",
                                                   "Abundances")) %>%
    dfTransformRaws()
  colnames(finalpd) <- c("Accession", paste0("Abundance", sample))
  finalpd
}

# Connect to raw Proteome Discoverer output files
super1pd <- ProtDiscProcess("P1387-01_Supernatant_1.pdResult", "Supernatant1")
super2pd <- ProtDiscProcess("P1387-02_Supernatant_2.pdResult", "Supernatant2")
super3pd <- ProtDiscProcess("P1387-03_Supernatant_3.pdResult", "Supernatant3")
intra1pd <- ProtDiscProcess("P1387-04_Intra_1.pdResult", "Intrasample1")
intra2pd <- ProtDiscProcess("P1387-05_Intra_2.pdResult", "Intrasample2")
intra3pd <- ProtDiscProcess("P1387-06_Intra_3.pdResult", "Intrasample3")

# Get full list of unique accessions from all samples
supernatantpd <- merge(super1pd, super2pd, by = "Accession", all = TRUE) %>%
  merge(super3pd, by = "Accession", all = TRUE)

intrapd <- merge(intra1pd, intra2pd, by ="Accession", all = TRUE) %>%
  merge(intra3pd, by = "Accession", all = TRUE)

write.table(supernatantpd, "Supernatant.txt", quote = FALSE, row.names = FALSE, sep = " \t")
write.table(intrapd, "Intrasample.txt", quote = FALSE, row.names = FALSE, sep = " \t")

##### Intrasample Results #####

intra_raw <- as.matrix(intrapd[,2:4])
dimnames(intra_raw) <- list(intrapd$Accession, colnames(intrapd)[2:4])

# Perform vsn normalization
meansdbef <- meanSdPlot(intra_raw, plot = FALSE)

intra_norm <- justvsn(intra_raw)
dimnames(intra_norm) <- list(intrapd$Accession, paste0(colnames(intra_raw), "Normalized"))

meansdaft <- meanSdPlot(intra_norm, plot = FALSE)

meansdbef$gg + ggtitle("Mean-SD Plot Before VSN Normalization")
meansdaft$gg + ggtitle("Mean-SD Plot After VSN Normalization")

# Create violin plots
violin(as.data.frame(intra_raw), "Violin Plot of Raw Data")

violin(as.data.frame(intra_norm), "Violin Plot of Normalized Data")

# Create reduced dataset - each accession can only have one NA value
intra_red <- intra_norm[rowSums(is.na(intra_norm)) < 2,] %>%
  as.data.frame()

# Remove PRTC
intra_red <- intra_red[!grepl("PRTC", rownames(intra_red)),]

# Take row means and order by magnitude
intra_red$MeanNormalizedAbundance <- rowMeans(intra_red, na.rm = TRUE)
intra_red <- intra_red[order(intra_red$MeanNormalizedAbundance, decreasing = TRUE),]

# merge with original abundance
intra_red <- intra_red %>%
  as.data.frame() %>%
  mutate(Accession = rownames(intra_red))
intra_red <- merge(intra_red, intrapd, by = "Accession", 
                   all.y = FALSE, sort = FALSE)

write.csv(intra_red, "Intrasample2+_MutualProteins.csv", row.names = FALSE)

##### Supernatant Results #####

super_raw <- as.matrix(supernatantpd[,2:4])
dimnames(super_raw) <- list(supernatantpd$Accession, colnames(supernatantpd)[2:4])

# Perform vsn normalization
meansdbef <- meanSdPlot(super_raw, plot = FALSE)

super_norm <- justvsn(super_raw, minDataPointsPerStratum = 20)
dimnames(super_norm) <- list(supernatantpd$Accession, paste0(colnames(super_raw), "Normalized"))

meansdaft <- meanSdPlot(super_norm, plot = FALSE)

meansdbef$gg + ggtitle("Mean-SD Plot Before VSN Normalization")
meansdaft$gg + ggtitle("Mean-SD Plot After VSN Normalization")

# Create violin plots
violin(as.data.frame(super_raw), "Violin Plot of Raw Data")

violin(as.data.frame(super_norm), "Violin Plot of Normalized Data")

# Create reduced dataset - each accession can only have one NA value
super_red <- super_norm[rowSums(is.na(super_norm)) < 2,] %>% 
  as.data.frame()

# Remove PRTC
super_red <- super_red[!grepl("PRTC", rownames(super_red)),]

# Take row means and order by magnitude
super_red$MeanNormalizedAbundance <- rowMeans(super_red, na.rm = TRUE)
super_red <- super_red[order(super_red$MeanNormalizedAbundance, decreasing = TRUE),]

# Merge with original abundance
super_red <- super_red %>%
  mutate(Accession = rownames(super_red))
super_red <- merge(super_red, supernatantpd, by = "Accession",
                               all.y = FALSE, sort = FALSE)

write.csv(super_red, "Supernatant2+_MutualProteins.csv", row.names = FALSE)
