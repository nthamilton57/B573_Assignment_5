# Install Packages --------------------------------------------------------
#run each of these lines individually. If there is an option to update, click yes.
install.packages("BiocManager")
install.packages("forcats")
install.packages("stringr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readr")
install.packages("tidyr")
install.packages("survminer")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("org.Hs.eg.db")


# a. Download ‘GSE33126’ from Gene Expression Omnibus ------------------------------------------------------
library(GEOquery)     #loading up the GEOquery library
my_id <- "GSE33126"     #this id is from the GEO website that corresponds to a microarray experiment. this is the link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33126. each sample listed on the website is a person. the RAW data is the fasta files. teh non-normalized data file is convereted to featurecounts (slide 7 in the ppt)
gse <- getGEO(my_id)    #this brings the entire study into our local environment as a large list of lists
gse <- gse[[1]]     #this simplifies the data
subj = pData(gse)   #this gets the phenotypic data for each sample (subjects) as a table (and shows a bunch of information about each sample).
expr = exprs(gse)   #this is the expression data for each sample (18 samples) for each gene (the 48803 number). the row names are the illumina gene names
annot = fData(gse)  #gives information about the genes

# b. Print the expression levels present in the data
print(expr)


# c. Log-normalize the data -----------------------------------------------------------
#summary(exprs(gse))   #looks mostly standardized as the max/mins/quartiles are mostly the same. Since the values are so big, we can log transform it to make it smaller/more manageable
exprs(gse) <- log2(exprs(gse))    #transform the data to smaller values/more manageable

# d. Create a box-plot illustrating the log-normalized data
boxplot(exprs(gse),outline=FALSE) #visualize that the data is normal
summary(exprs(gse))

# e. Load the corresponding phenotype data and print it -------------------------------------------------------------
library(dplyr)
sampleInfo <- pData(gse)  #get phenotypic data again in new variable

# Use select and rename from dplyr to transform the data
sampleInfo <- sampleInfo %>%
  dplyr::select(source_name_ch1, characteristics_ch1.1) %>%
  dplyr::rename(group = source_name_ch1, patient = characteristics_ch1.1)

print(sampleInfo)

# f. Perform Clustering Between your samples and display it as a heat map ----------------------------------------------------------
library(pheatmap)

## argument use="c" stops an error if there are any missing data points
corMatrix <- cor(exprs(gse),use="c")    #get correlation matrix of the expression profiles
pheatmap(corMatrix) #visualizes it with a heat map. Make sure to check the range. this correlation range is so high because the samples are so similar (all humans, all have some cancer, etc)


# g. Recreate the heat map from part e, including annotations for the patient id as well as patient group
rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo)   #this shows that the big separation in correlation is due to the group which makes sense cause cancer vs non cancer


# h. Create a scatter plot of your data utilizing Principle Component Analysis (PCA) vectors--------------------------------------------
#PCA can reduce the number of dimensions (48000 genes to a lot less)
library(ggrepel)
library(ggplot2)

pca <- prcomp(t(exprs(gse)))
## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>%
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste("Patient", patient))) + geom_point() + geom_text_repel()

# i. Differential Expression -------------------------------------------------
library(limma)
design <- model.matrix(~0+sampleInfo$group)   #the ~0+sampleInfo$group is like the equivalent of the linear regression formula where 0 is the intercept for each subject (baseline for each subject)

## the column names are a bit ugly, so we will rename
colnames(design) <- c("Normal","Tumour")

## calculate median expression level
cutoff <- median(exprs(gse))

## TRUE or FALSE for whether each gene is "expressed" in each sample. this will get rid of genes that are expressed lower than the median
is_expressed <- exprs(gse) > cutoff

## Identify genes expressed in more than 2 samples
keep <- rowSums(is_expressed) > 2

## subset to just those expressed genes
gse <- gse[keep,]

#here we are fitting the model to see if the group affects the gene expression
fit <- lmFit(exprs(gse), design)

#n order to perform the differential analysis, we have to define the contrast that we are interested in. In our case we only have two groups and one contrast of interest. Multiple contrasts can be defined in the makeContrasts function.
#Here we are looking at differences between Tumour and Normal
contrasts <- makeContrasts(Tumour - Normal, levels=design)    #shows that -1 is Normal and 1 is the Tumour cells

fit2 <- contrasts.fit(fit, contrasts) #runs the differential expression with the contrasts now
fit2 <- eBayes(fit2)  #applies empirical bayes formula. gives a more accurate p-value based on the actual data from this dataset
    #this shows the genes that are the most significant in differentiating between normal tissue and tissue
    #positive logFC means that it is more present in tumourous genes where negative logFC means that it is present in less tumourous genes

#for each of these genes we found, we would want to find out what they actually are and give more information
anno <- fData(gse)
anno <- dplyr::select(anno, Symbol, Entrez_Gene_ID, Chromosome, Cytoband)  #selects the information we want to show for each gene
fit2$genes <- anno
top_10 <- topTable(fit2)  #shows the topTable with the additional information
print(top_10)

# j. Volcano Plot ------------------------------------------------------------
full_results <- topTable(fit2, number=Inf)    #include all results
full_results <- tibble::rownames_to_column(full_results,"ID")

p_cutoff <- 0.05
fc_cutoff <- 1
topN <- 10

full_results %>%
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>%
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, Symbol,"")) %>%
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")


# k. Create an MA plot -------------------------------------------------------

# Extract the log-fold changes (logFC) and average expression (AveExpr) values
logFC <- fit2$coefficients[, "Tumour - Normal"]
aveExpr <- rowMeans(exprs(gse))

# Create a variable indicating DEGs (e.g., based on adjusted p-value cutoff)
cutoff <- 0.05
DEGs <- full_results$adj.P.Val < cutoff

# Create a color vector based on DEG status
colors <- ifelse(DEGs, "red", "blue")

ma_plot <- ggplot(data = NULL, aes(x = aveExpr, y = logFC, col = colors)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "MA Plot",
       x = "Average Expression (AveExpr)",
       y = "Log-Fold Change (logFC)",
       color = "DEG Status") +
  scale_color_manual(values = c("red", "blue"), labels = c("DEG", "Non-DEG")) +
  theme_minimal()

# Display the MA plot
print(ma_plot)

# l. Use the identified DEGs to determine over-represented Reactome pathways --------------------------------------------
#https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html#pathway-visualization

BiocManager::install('ReactomePA')
library(ReactomePA)

de <- (full_results$Entrez_Gene_ID)[abs(full_results$logFC) > 1.5]

x <- enrichPathway(gene=de, pvalueCutoff = 0.05, readable=TRUE)
x <- as.data.frame(x)

print(x)

# m. Perform Gene Set Enrichment analysis on your identified DEGs ---------

genelist <- as.vector(full_results$logFC)
genelist <- sort(genelist, decreasing = TRUE)
genelist <- setNames(genelist, full_results$Entrez_Gene_ID)

y <- gsePathway(genelist,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE)

print(y)


# n. Visualize your pathway as a network ----------------------------------
viewPathway("Nitric oxide stimulates guanylate cyclase",
            readable = TRUE)


# o. Disease Enrichment Analviysis of Top 20 Genes -----------------------------
#https://yulab-smu.top/biomedical-knowledge-mining-book/dose-enrichment.html

library(DOSE)
y <- gseDO(genelist,
           minGSSize     = 120,
           pvalueCutoff  = 0.2,
           pAdjustMethod = "BH",
           verbose       = FALSE)
print(y)
