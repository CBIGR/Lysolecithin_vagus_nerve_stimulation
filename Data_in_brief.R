#!/usr/bin/Rscript

library(data.table)
library(limma)
library(stringr)
library(ggplot2)
library(ggfortify)
library(ggpubr)
#library(ggrepel)
library(EnhancedVolcano)
library(msigdbr)
#library(enrichplot)
#library(ReactomePA)
library(org.Rn.eg.db)
#library(clusterProfiler)
library(KEGGREST)

setwd('/home/boris/Documents/PhD/gut_brain/VNS_Helen/DiB_rebuttal')

##########################################################################################################################
#### Read and explore data
##########################################################################################################################

dat_DE <- fread('../RE-ANALYSIS/MS_results_PRC-6146_DE_imp_quant_protein.csv', header =TRUE, data.table=FALSE, fill=TRUE, sep ='\t')
dat_DE <- dat_DE[, c(3,9:18)]
colnames(dat_DE) <- c('Protein_names', 'DE_cVNS_1', 'DE_cVNS_2', 'DE_cVNS_3', 'DE_cVNS_4', 'DE_cVNS_5', 'DE_sham_1', 'DE_sham_2', 'DE_sham_3', 'DE_sham_4', 'DE_sham_5')

dat_RE <- fread('../RE-ANALYSIS/MS_results_PRC-6146_RE_imp_quant_protein.csv', header =TRUE, data.table=FALSE, fill=TRUE, sep ='\t')
dat_RE <- dat_RE[, c(3,9:18)]
colnames(dat_RE) <- c('Protein_names', 'RE_cVNS_1', 'RE_cVNS_2', 'RE_cVNS_3', 'RE_cVNS_4', 'RE_cVNS_5', 'RE_sham_1', 'RE_sham_2', 'RE_sham_3', 'RE_sham_4', 'RE_sham_5')

# Merge data on protein names
dat <- merge(dat_DE, dat_RE, by = 'Protein_names')
write.csv(dat, file ='./All_samples_imp_quant_protein.csv', row.names = FALSE)

#### Exploratory PCA with imputed data

# Set row names
rownames(dat) <- dat$Protein_names; dat$Protein_names <- NULL

dat_t <- data.frame(t(dat))
dat_t$condition <- substr(rownames(dat_t), 1,2)
dat_t$treatment <- substr(rownames(dat_t), 4,7)
# Rearrange data
dat_t <- dat_t[, c(ncol(dat_t)-1,ncol(dat_t), 1:(ncol(dat_t)-2))]

dat_t[is.na(dat_t)] = 0

pca <- prcomp(dat_t[,-c(1,2)], center = TRUE, scale. = TRUE)
pca_plot <- autoplot(pca, data = dat_t, colour = 'condition', shape ='treatment')
pca_plot <- pca_plot + theme_light() + ggtitle(' cVNS and SHAM') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_label_repel(aes(label = rownames(dat_t)), size = 2, box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', segment.size = 0.25, nudge_x = 0.15, nudge_y = 0.15)


pca_plot <- autoplot(pca, data = dat_t, colour = 'condition', shape = 'treatment', size = 3)
pca_plot + theme_light() + 
  ggtitle('PCA projection of cVNS and sham samples during demyelination and remyelination') + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, margin = margin(b = 20)),  # Adjust b value for bottom margin
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)  # Adjust these values for overall plot margins
  )

ggsave('./PCA_all_samples.pdf')





##########################################################################################################################
#### Try differential analysis using limma
##########################################################################################################################

conditions <- c('DE', 'RE')
treatments <- c('cVNS', 'sham')

# Create metadata
metadata <- data.frame(row.names = rownames(dat_t))
metadata$sample <- rownames(dat_t)
metadata$condition <- factor(substr(rownames(dat_t), 1, 2), levels = conditions)
metadata$treatment <- factor(substr(rownames(dat_t), 4, 7), levels = treatments)
metadata$treatment <- relevel(metadata$treatment, ref = 'sham')

TS <- paste(metadata$condition, metadata$treatment, sep = '.')
TS

TS <- factor(TS, levels = c('DE.cVNS', 'DE.sham', 'RE.cVNS', 'RE.sham'))
design2 <- model.matrix(~0+TS)
colnames(design2) <- levels(TS)
fit <- lmFit(dat, design2)

cont.matrix <- makeContrasts(
  DE_cVNS_vs_DE_sham = DE.cVNS - DE.sham, # The effect of cVNS in DE
  RE_cVNS_vs_RE_sham = RE.cVNS - RE.sham, # The effect of cVNS in RE
  DE_cVNS_vs_RE_cVNS = DE.cVNS - RE.cVNS, # The difference between DE and RE in cVNS
  DE_sham_vs_RE_sham = DE.sham - RE.sham, # The difference between DE and RE in sham (control)
  levels=design2
)

# Make a directory for each of the contrasts
dir.create('./limma/DE_cVNS_vs_DE_sham', recursive = T)
dir.create('./limma/RE_cVNS_vs_RE_sham', recursive = T)
dir.create('./limma/DE_cVNS_vs_RE_cVNS', recursive = T)
dir.create('./limma/DE_sham_vs_RE_sham', recursive = T)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

DE_cVNS_vs_DE_sham <- topTable(fit2, coef="DE_cVNS_vs_DE_sham", number = Inf)
RE_cVNS_vs_RE_sham <- topTable(fit2, coef="RE_cVNS_vs_RE_sham", number = Inf)
DE_cVNS_vs_RE_cVNS <- topTable(fit2, coef="DE_cVNS_vs_RE_cVNS", number = Inf)
DE_sham_vs_RE_sham <- topTable(fit2, coef="DE_sham_vs_RE_sham", number = Inf)

# Write tables
write.table(DE_cVNS_vs_DE_sham, file = './limma/DE_cVNS_vs_DE_sham/DE_cVNS_vs_DE_sham.txt', sep = '\t', quote = FALSE)
write.table(RE_cVNS_vs_RE_sham, file = './limma/RE_cVNS_vs_RE_sham/RE_cVNS_vs_RE_sham.txt', sep = '\t', quote = FALSE)
write.table(DE_cVNS_vs_RE_cVNS, file = './limma/DE_cVNS_vs_RE_cVNS/DE_cVNS_vs_RE_cVNS.txt', sep = '\t', quote = FALSE)
write.table(DE_sham_vs_RE_sham, file = './limma/DE_sham_vs_RE_sham/DE_sham_vs_RE_sham.txt', sep = '\t', quote = FALSE)

# Write lists of up- and down-regulated genes
#write.table(DE_cVNS_vs_DE_sham[DE_cVNS_vs_DE_sham$logFC >= 1 & DE_cVNS_vs_DE_sham$adj.P.Val <= 0.05,], file = './limma/DE_cVNS_vs_DE_sham/DE_cVNS_vs_DE_sham_up_padj_0.05_LFC_1.txt', sep = '\t', quote = FALSE)
#write.table(DE_cVNS_vs_DE_sham[DE_cVNS_vs_DE_sham$logFC <= -1 & DE_cVNS_vs_DE_sham$adj.P.Val <= 0.05,], file = './limma/DE_cVNS_vs_DE_sham/DE_cVNS_vs_DE_sham_down_padj_0.05_LFC_-1.txt', sep = '\t', quote = FALSE)
#write.table(RE_cVNS_vs_RE_sham[RE_cVNS_vs_RE_sham$logFC >= 1 & RE_cVNS_vs_RE_sham$adj.P.Val <= 0.05,], file = './limma/RE_cVNS_vs_RE_sham/RE_cVNS_vs_RE_sham_up_padj_0.05_LFC_1.txt.txt', sep = '\t', quote = FALSE)
#write.table(RE_cVNS_vs_RE_sham[RE_cVNS_vs_RE_sham$logFC <= -1 & RE_cVNS_vs_RE_sham$adj.P.Val <= 0.05,], file = './limma/RE_cVNS_vs_RE_sham/RE_cVNS_vs_RE_sham_down_padj_0.05_LFC_-1.txt.txt', sep = '\t', quote = FALSE)
write.table(DE_cVNS_vs_RE_cVNS[DE_cVNS_vs_RE_cVNS$logFC >= 1 & DE_cVNS_vs_RE_cVNS$adj.P.Val <= 0.05,], file = './limma/DE_cVNS_vs_RE_cVNS/DE_cVNS_vs_RE_cVNS_up_padj_0.05_LFC_1.txt.txt', sep = '\t', quote = FALSE)
write.table(DE_cVNS_vs_RE_cVNS[DE_cVNS_vs_RE_cVNS$logFC <= -1 & DE_cVNS_vs_RE_cVNS$adj.P.Val <= 0.05,], file = './limma/DE_cVNS_vs_RE_cVNS/DE_cVNS_vs_RE_cVNS_down_padj_0.05_LFC_-1.txt.txt', sep = '\t', quote = FALSE)
write.table(DE_sham_vs_RE_sham[DE_sham_vs_RE_sham$logFC >= 1 & DE_sham_vs_RE_sham$adj.P.Val <= 0.05,], file = './limma/DE_sham_vs_RE_sham/DE_sham_vs_RE_sham_up_padj_0.05_LFC_1.txt.txt', sep = '\t', quote = FALSE)
write.table(DE_sham_vs_RE_sham[DE_sham_vs_RE_sham$logFC <= -1 & DE_sham_vs_RE_sham$adj.P.Val <= 0.05,], file = './limma/DE_sham_vs_RE_sham/DE_sham_vs_RE_sham_down_padj_0.05_LFC_-1.txt.txt', sep = '\t', quote = FALSE)
 

length(rownames(DE_sham_vs_RE_sham[DE_sham_vs_RE_sham$logFC >= 1 & DE_sham_vs_RE_sham$adj.P.Val <= 0.05,]))
length(rownames(DE_sham_vs_RE_sham[DE_sham_vs_RE_sham$logFC <= -1 & DE_sham_vs_RE_sham$adj.P.Val <= 0.05,]))

length(rownames(DE_cVNS_vs_RE_cVNS[DE_cVNS_vs_RE_cVNS$logFC >= 1 & DE_cVNS_vs_RE_cVNS$adj.P.Val <= 0.05,]))
length(rownames(DE_cVNS_vs_RE_cVNS[DE_cVNS_vs_RE_cVNS$logFC <= -1 & DE_cVNS_vs_RE_cVNS$adj.P.Val <= 0.05,]))

# WIth LFC cutoff >= 2

write.table(DE_cVNS_vs_RE_cVNS[DE_cVNS_vs_RE_cVNS$logFC >= 2 & DE_cVNS_vs_RE_cVNS$adj.P.Val <= 0.05,], file = './limma/DE_cVNS_vs_RE_cVNS/DE_cVNS_vs_RE_cVNS_up_padj_0.05_LFC_2.txt.txt', sep = '\t', quote = FALSE)
write.table(DE_cVNS_vs_RE_cVNS[DE_cVNS_vs_RE_cVNS$logFC <= -2 & DE_cVNS_vs_RE_cVNS$adj.P.Val <= 0.05,], file = './limma/DE_cVNS_vs_RE_cVNS/DE_cVNS_vs_RE_cVNS_down_padj_0.05_LFC_-2.txt.txt', sep = '\t', quote = FALSE)
write.table(DE_sham_vs_RE_sham[DE_sham_vs_RE_sham$logFC >= 2 & DE_sham_vs_RE_sham$adj.P.Val <= 0.05,], file = './limma/DE_sham_vs_RE_sham/DE_sham_vs_RE_sham_up_padj_0.05_LFC_2.txt.txt', sep = '\t', quote = FALSE)
write.table(DE_sham_vs_RE_sham[DE_sham_vs_RE_sham$logFC <= -2 & DE_sham_vs_RE_sham$adj.P.Val <= 0.05,], file = './limma/DE_sham_vs_RE_sham/DE_sham_vs_RE_sham_down_padj_0.05_LFC_-2.txt.txt', sep = '\t', quote = FALSE)


length(rownames(DE_sham_vs_RE_sham[DE_sham_vs_RE_sham$logFC >= 2 & DE_sham_vs_RE_sham$adj.P.Val <= 0.05,]))
length(rownames(DE_sham_vs_RE_sham[DE_sham_vs_RE_sham$logFC <= -2 & DE_sham_vs_RE_sham$adj.P.Val <= 0.05,]))

length(rownames(DE_cVNS_vs_RE_cVNS[DE_cVNS_vs_RE_cVNS$logFC >= 2 & DE_cVNS_vs_RE_cVNS$adj.P.Val <= 0.05,]))
length(rownames(DE_cVNS_vs_RE_cVNS[DE_cVNS_vs_RE_cVNS$logFC <= -2 & DE_cVNS_vs_RE_cVNS$adj.P.Val <= 0.05,]))


#### Read genes that are interesting in the context of neuroinflammation
neuroinflammation_genes <- fread('./Proteins_neuroinflammation.csv', sep = '\t')
# If `Gene names` is an empty string, put in the `Protein IDs`
neuroinflammation_genes$`Gene names`[neuroinflammation_genes$`Gene names` == ''] <- neuroinflammation_genes$`Protein IDs`[neuroinflammation_genes$`Gene names` == '']
inflammation_genes <- neuroinflammation_genes[neuroinflammation_genes$Inflammation == 'x',]$`Gene names`
neuro_genes <- neuroinflammation_genes[neuroinflammation_genes$`Neurological revelant` == 'x',]$`Gene names`
all_genes <- unique(c(inflammation_genes, neuro_genes))
# Remove empty strings from gene sets
all_genes <- all_genes[all_genes != '']
inflammation_genes <- inflammation_genes[inflammation_genes != '']
neuro_genes <- neuro_genes[neuro_genes != '']


protein_gene <- fread('../RE-ANALYSIS/MS_results_PRC-6146_DE_imp_quant_protein.csv', header =TRUE, data.table=FALSE, fill=TRUE, sep ='\t')[,c(3,4,5)]
colnames(protein_gene) <- c('Protein_names', 'Gene_names', 'Description')
rownames(protein_gene) <- protein_gene$Protein_names

# Add column Gene_name to toptables
DE_cVNS_vs_RE_cVNS$Gene_names <- protein_gene[rownames(DE_cVNS_vs_RE_cVNS),]$Gene_names
DE_sham_vs_RE_sham$Gene_names <- protein_gene[rownames(DE_sham_vs_RE_sham),]$Gene_names
DE_cVNS_vs_DE_sham$Gene_names <- protein_gene[rownames(DE_cVNS_vs_DE_sham),]$Gene_names
RE_cVNS_vs_RE_sham$Gene_names <- protein_gene[rownames(RE_cVNS_vs_RE_sham),]$Gene_names



# Define differential expression criteria
pval_threshold <- 0.05
log2fc_threshold <- 1

# Define function to create volcano plot
create_volcano_plot <- function(df, title, genes_to_highlight) {
  #Split title on space and get first element
  folder <- unlist(strsplit(title, ' '))[1]
  
  #df <- DE_cVNS_vs_RE_cVNS
  EnhancedVolcano(
    df,
    lab = df$Gene_names,
    selectLab = genes_to_highlight,
    x = 'logFC',
    y = 'adj.P.Val',
    pCutoff = pval_threshold,
    FCcutoff = log2fc_threshold,
    #xlim = c(-5, 5),
    title = title,
    subtitle = '',
    pointSize = 3,
    labSize = 5.5) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0('./limma/', folder,'/', title,'_volcano.png'), width = 10, height = 10, dpi = 300)
}


# Create volcano plots
#pdf('./results/limma/DE_cVNS_vs_DE_sham/DE_cVNS_vs_DE_sham_volcano.pdf')
create_volcano_plot(DE_cVNS_vs_DE_sham, 'DE_cVNS_vs_DE_sham - Neuroinflammation', inflammation_genes)
create_volcano_plot(DE_cVNS_vs_DE_sham, 'DE_cVNS_vs_DE_sham - Neurologically relevant', neuro_genes)
#dev.off()

#pdf('./results/limma/RE_cVNS_vs_RE_sham/RE_cVNS_vs_RE_sham_volcano.pdf')
create_volcano_plot(RE_cVNS_vs_RE_sham, 'RE_cVNS_vs_RE_sham - Neuroinflammation', inflammation_genes)
create_volcano_plot(RE_cVNS_vs_RE_sham, 'RE_cVNS_vs_RE_sham - Neurologically relevant', neuro_genes)
#dev.off()

#pdf('./results/limma/DE_cVNS_vs_RE_cVNS/DE_cVNS_vs_RE_cVNS_volcano.pdf')
create_volcano_plot(DE_cVNS_vs_RE_cVNS, 'DE_cVNS_vs_RE_cVNS - Neuroinflammation', inflammation_genes)
create_volcano_plot(DE_cVNS_vs_RE_cVNS, 'DE_cVNS_vs_RE_cVNS - Neurologically relevant', neuro_genes)
#dev.off()

#pdf('./results/limma/DE_sham_vs_RE_sham/DE_sham_vs_RE_sham_volcano.pdf')
create_volcano_plot(DE_sham_vs_RE_sham, 'DE_sham_vs_RE_sham - Neuroinflammation', inflammation_genes)
create_volcano_plot(DE_sham_vs_RE_sham, 'DE_sham_vs_RE_sham - Neurologically relevant', neuro_genes)
#dev.off()


# Define differential expression criteria
pval_threshold <- 0.05
log2fc_threshold <- 2


create_volcano_plot(DE_cVNS_vs_DE_sham, 'DE_cVNS_vs_DE_sham - Neuroinflammation 2', inflammation_genes)
create_volcano_plot(DE_cVNS_vs_DE_sham, 'DE_cVNS_vs_DE_sham - Neurologically relevant 2', neuro_genes)
#dev.off()

#pdf('./results/limma/RE_cVNS_vs_RE_sham/RE_cVNS_vs_RE_sham_volcano.pdf')
create_volcano_plot(RE_cVNS_vs_RE_sham, 'RE_cVNS_vs_RE_sham - Neuroinflammation 2', inflammation_genes)
create_volcano_plot(RE_cVNS_vs_RE_sham, 'RE_cVNS_vs_RE_sham - Neurologically relevant 2', neuro_genes)
#dev.off()

#pdf('./results/limma/DE_cVNS_vs_RE_cVNS/DE_cVNS_vs_RE_cVNS_volcano.pdf')
create_volcano_plot(DE_cVNS_vs_RE_cVNS, 'DE_cVNS_vs_RE_cVNS - Neuroinflammation 2', inflammation_genes)
create_volcano_plot(DE_cVNS_vs_RE_cVNS, 'DE_cVNS_vs_RE_cVNS - Neurologically relevant 2', neuro_genes)
#dev.off()

#pdf('./results/limma/DE_sham_vs_RE_sham/DE_sham_vs_RE_sham_volcano.pdf')
create_volcano_plot(DE_sham_vs_RE_sham, 'DE_sham_vs_RE_sham - Neuroinflammation 2', inflammation_genes)
create_volcano_plot(DE_sham_vs_RE_sham, 'DE_sham_vs_RE_sham - Neurologically relevant 2', neuro_genes)
#dev.off()

pval_threshold <- 0.05
log2fc_threshold <- 1


create_volcano_plot(DE_cVNS_vs_DE_sham, 'DE_cVNS_vs_DE_sham', all_genes)
create_volcano_plot(RE_cVNS_vs_RE_sham, 'RE_cVNS_vs_RE_sham', all_genes)
create_volcano_plot(DE_cVNS_vs_RE_cVNS, 'DE_cVNS_vs_RE_cVNS', all_genes)
create_volcano_plot(DE_sham_vs_RE_sham, 'DE_sham_vs_RE_sham', all_genes)

library(pheatmap)
# Create heatmap that shows Z-scored protein expression for proteins related to neuroinflammation or neurologically relevant proteins

# inflammation_proteins <- neuroinflammation_genes[neuroinflammation_genes$Inflammation == 'x',]$`Protein IDs`
# neuro_proteins <- neuroinflammation_genes[neuroinflammation_genes$`Neurological revelant` == 'x',]$`Protein IDs`


create_heatmap <- function(df, title_full, protein_set, samples, folder, plot_title) {
  
  # df <- dat
  # title_full <- 'DE_cVNS_vs_RE_cVNS_neurologically_relevant'
  # protein_set <- all_genes
  # samples <- c('DE_cVNS', 'RE_cVNS', 'gene')
  # # 
  annotation_df <- data.frame(
    row.names = colnames(df),
    group = factor(substr(colnames(df), 1, 7))
  )
  
  df$gene <- protein_gene[rownames(df),]$Gene_names
  df <- df[df$gene %in% protein_set, (substr(colnames(df),1,7) %in% samples)]
  
 
  #remove rows with duplicate gene
  df <- df[!duplicated(df$gene),]
  
  rownames(df) <- df$gene; df$gene <- NULL

  # Z-score the data
  significant_genes_zscore <- t(scale(t(df)))
  # Define color palette
  my_palette <- colorRampPalette(c("blue", "white", "red"))(50)
  # Create heatmap
  pheatmap(
    significant_genes_zscore, 
    annotation_col = annotation_df, 
    show_rownames = TRUE, 
    show_colnames = FALSE, 
    fontsize = 8, 
    fontsize_row = 9, 
    fontsize_col = 8, 
    color = my_palette, 
    cluster_rows = TRUE, 
    cluster_cols = TRUE,
    main = paste0(plot_title, '\n'), 
    filename = paste0('./limma/', folder, '/', title_full, '_heatmap.png')
  ) 
  
  #ggsave(paste0('./limma/', folder,'/', title,'_heatmap.png'))
}

samples <- c('DE_cVNS', 'RE_cVNS', 'gene')
create_heatmap(dat, 'DE_cVNS_vs_RE_cVNS_neurologically_relevant', neuro_genes, samples, 'DE_cVNS_vs_RE_cVNS', plot_title = 'DE_cVNS_vs_RE_cVNS - Neurologically relevant')
create_heatmap(dat, 'DE_cVNS_vs_RE_cVNS_neuroinflammation', inflammation_genes, samples, 'DE_cVNS_vs_RE_cVNS', plot_title = 'DE_cVNS_vs_RE_cVNS - Neuroinflammation')
create_heatmap(dat, 'DE_cVNS_vs_RE_cVNS_all', all_genes, samples, 'DE_cVNS_vs_RE_cVNS', plot_title = 'DE_cVNS_vs_RE_cVNS')

samples <- c('DE_sham', 'RE_sham', 'gene')
create_heatmap(dat, 'DE_sham_vs_RE_sham_neurologically_relevant', neuro_genes, samples, 'DE_sham_vs_RE_sham', plot_title = 'DE_sham_vs_RE_sham - Neurologically relevant')
create_heatmap(dat, 'DE_sham_vs_RE_sham_neuroinflammation', inflammation_genes, samples, 'DE_sham_vs_RE_sham' , plot_title = 'DE_sham_vs_RE_sham - Neuroinflammation')
create_heatmap(dat, 'DE_sham_vs_RE_sham_all', all_genes, samples, 'DE_sham_vs_RE_sham', plot_title = 'DE_sham_vs_RE_sham')
