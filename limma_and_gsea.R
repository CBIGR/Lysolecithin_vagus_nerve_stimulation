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
library(enrichplot)
library(ReactomePA)
library(org.Rn.eg.db)
library(clusterProfiler)
library(KEGGREST)

setwd('/home/boris/Documents/PhD/gut_brain/VNS_Helen/RE-ANALYSIS/')

##########################################################################################################################
#### Read and explore data
##########################################################################################################################

dat_DE <- fread('MS_results_PRC-6146_DE_imp_quant_protein.csv', header =TRUE, data.table=FALSE, fill=TRUE, sep ='\t')
dat_DE <- dat_DE[, c(3,9:18)]
colnames(dat_DE) <- c('Protein_names', 'DE_cVNS_1', 'DE_cVNS_2', 'DE_cVNS_3', 'DE_cVNS_4', 'DE_cVNS_5', 'DE_sham_1', 'DE_sham_2', 'DE_sham_3', 'DE_sham_4', 'DE_sham_5')

dat_RE <- fread('MS_results_PRC-6146_RE_imp_quant_protein.csv', header =TRUE, data.table=FALSE, fill=TRUE, sep ='\t')
dat_RE <- dat_RE[, c(3,9:18)]
colnames(dat_RE) <- c('Protein_names', 'RE_cVNS_1', 'RE_cVNS_2', 'RE_cVNS_3', 'RE_cVNS_4', 'RE_cVNS_5', 'RE_sham_1', 'RE_sham_2', 'RE_sham_3', 'RE_sham_4', 'RE_sham_5')

# Merge data on protein names
dat <- merge(dat_DE, dat_RE, by = 'Protein_names')
write.csv(dat, file ='./results/All_samples_imp_quant_protein.csv', row.names = FALSE)

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
pca_plot <- pca_plot + theme_light() + ggtitle('Remyelination cVNS and SHAM') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_label_repel(aes(label = rownames(dat_t)), size = 2, box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', segment.size = 0.25, nudge_x = 0.15, nudge_y = 0.15)

pdf('./results/PCA_all_samples.pdf')
pca_plot
dev.off()



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
dir.create('./results/limma/DE_cVNS_vs_DE_sham')
dir.create('./results/limma/RE_cVNS_vs_RE_sham')
dir.create('./results/limma/DE_cVNS_vs_RE_cVNS')
dir.create('./results/limma/DE_sham_vs_RE_sham')

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

DE_cVNS_vs_DE_sham <- topTable(fit2, coef="DE_cVNS_vs_DE_sham", number = Inf)
RE_cVNS_vs_RE_sham <- topTable(fit2, coef="RE_cVNS_vs_RE_sham", number = Inf)
DE_cVNS_vs_RE_cVNS <- topTable(fit2, coef="DE_cVNS_vs_RE_cVNS", number = Inf)
DE_sham_vs_RE_sham <- topTable(fit2, coef="DE_sham_vs_RE_sham", number = Inf)

# Write tables
write.table(DE_cVNS_vs_DE_sham, file = './results/limma/DE_cVNS_vs_DE_sham/DE_cVNS_vs_DE_sham.txt', sep = '\t', quote = FALSE)
write.table(RE_cVNS_vs_RE_sham, file = './results/limma/RE_cVNS_vs_RE_sham/RE_cVNS_vs_RE_sham.txt', sep = '\t', quote = FALSE)
write.table(DE_cVNS_vs_RE_cVNS, file = './results/limma/DE_cVNS_vs_RE_cVNS/DE_cVNS_vs_RE_cVNS.txt', sep = '\t', quote = FALSE)
write.table(DE_sham_vs_RE_sham, file = './results/limma/DE_sham_vs_RE_sham/DE_sham_vs_RE_sham.txt', sep = '\t', quote = FALSE)

# Write lists of up- and down-regulated genes
write.table(DE_cVNS_vs_DE_sham[DE_cVNS_vs_DE_sham$logFC > 1 & DE_cVNS_vs_DE_sham$adj.P.Val < 0.05,], file = './results/limma/DE_cVNS_vs_DE_sham_up.txt', sep = '\t', quote = FALSE)
write.table(DE_cVNS_vs_DE_sham[DE_cVNS_vs_DE_sham$logFC < -1 & DE_cVNS_vs_DE_sham$adj.P.Val < 0.05,], file = './results/limma/DE_cVNS_vs_DE_sham_down.txt', sep = '\t', quote = FALSE)
write.table(RE_cVNS_vs_RE_sham[RE_cVNS_vs_RE_sham$logFC > 1 & RE_cVNS_vs_RE_sham$adj.P.Val < 0.05,], file = './results/limma/RE_cVNS_vs_RE_sham_up.txt', sep = '\t', quote = FALSE)
write.table(RE_cVNS_vs_RE_sham[RE_cVNS_vs_RE_sham$logFC < -1 & RE_cVNS_vs_RE_sham$adj.P.Val < 0.05,], file = './results/limma/RE_cVNS_vs_RE_sham_down.txt', sep = '\t', quote = FALSE)
write.table(DE_cVNS_vs_RE_cVNS[DE_cVNS_vs_RE_cVNS$logFC > 1 & DE_cVNS_vs_RE_cVNS$adj.P.Val < 0.05,], file = './results/limma/DE_cVNS_vs_RE_cVNS_up.txt', sep = '\t', quote = FALSE)
write.table(DE_cVNS_vs_RE_cVNS[DE_cVNS_vs_RE_cVNS$logFC < -1 & DE_cVNS_vs_RE_cVNS$adj.P.Val < 0.05,], file = './results/limma/DE_cVNS_vs_RE_cVNS_down.txt', sep = '\t', quote = FALSE)
write.table(DE_sham_vs_RE_sham[DE_sham_vs_RE_sham$logFC > 1 & DE_sham_vs_RE_sham$adj.P.Val < 0.05,], file = './results/limma/DE_sham_vs_RE_sham_up.txt', sep = '\t', quote = FALSE)
write.table(DE_sham_vs_RE_sham[DE_sham_vs_RE_sham$logFC < -1 & DE_sham_vs_RE_sham$adj.P.Val < 0.05,], file = './results/limma/DE_sham_vs_RE_sham_down.txt', sep = '\t', quote = FALSE)
 

length(rownames(DE_sham_vs_RE_sham[DE_sham_vs_RE_sham$logFC > 1 & DE_sham_vs_RE_sham$adj.P.Val < 0.05,]))
length(rownames(DE_sham_vs_RE_sham[DE_sham_vs_RE_sham$logFC < -1 & DE_sham_vs_RE_sham$adj.P.Val < 0.05,]))

length(rownames(DE_cVNS_vs_RE_cVNS[DE_cVNS_vs_RE_cVNS$logFC > 1 & DE_cVNS_vs_RE_cVNS$adj.P.Val < 0.05,]))
length(rownames(DE_cVNS_vs_RE_cVNS[DE_cVNS_vs_RE_cVNS$logFC < -1 & DE_cVNS_vs_RE_cVNS$adj.P.Val < 0.05,]))


pdf('./results/limma/DE_cVNS_vs_DE_sham/DE_cVNS_vs_DE_sham_volcano.pdf')
EnhancedVolcano(DE_cVNS_vs_DE_sham, lab = rownames(DE_cVNS_vs_DE_sham), x = 'logFC', y = 'adj.P.Val', pCutoff = 0.05, FCcutoff = 1, xlim = c(-5,5), title = 'DE_cVNS_vs_DE_sham')
dev.off()

pdf('./results/limma/RE_cVNS_vs_RE_sham/RE_cVNS_vs_RE_sham_volcano.pdf')
EnhancedVolcano(RE_cVNS_vs_RE_sham, lab = rownames(RE_cVNS_vs_RE_sham), x = 'logFC', y = 'adj.P.Val', pCutoff = 0.05, FCcutoff = 1, xlim = c(-5,5), title = 'RE_cVNS_vs_RE_sham')
dev.off()

pdf('./results/limma/DE_cVNS_vs_RE_cVNS/DE_cVNS_vs_RE_cVNS_volcano.pdf')
EnhancedVolcano(DE_cVNS_vs_RE_cVNS, lab = rownames(DE_cVNS_vs_RE_cVNS), x = 'logFC', y = 'adj.P.Val', pCutoff = 0.05, FCcutoff = 1, xlim = c(-5,5), title = 'DE_cVNS_vs_RE_cVNS')
dev.off()

pdf('./results/limma/DE_sham_vs_RE_sham/DE_sham_vs_RE_sham_volcano.pdf')
EnhancedVolcano(DE_sham_vs_RE_sham, lab = rownames(DE_sham_vs_RE_sham), x = 'logFC', y = 'adj.P.Val', pCutoff = 0.05, FCcutoff = 1, xlim = c(-5,5), title = 'DE_sham_vs_RE_sham')
dev.off()


protein_gene <- fread('MS_results_PRC-6146_DE_imp_quant_protein.csv', header =TRUE, data.table=FALSE, fill=TRUE, sep ='\t')[,c(3,4,5)]
colnames(protein_gene) <- c('Protein_names', 'Gene_names', 'Description')
rownames(protein_gene) <- protein_gene$Protein_names; protein_gene$Protein_names <- NULL

DE_cVNS_vs_DE_sham <- merge(DE_cVNS_vs_DE_sham, protein_gene, by=0)
RE_cVNS_vs_RE_sham <- merge(RE_cVNS_vs_RE_sham, protein_gene, by=0)
DE_sham_vs_RE_sham <- merge(DE_sham_vs_RE_sham, protein_gene, by=0)
DE_cVNS_vs_RE_cVNS <- merge(DE_cVNS_vs_RE_cVNS, protein_gene, by=0)

##########################################################################################################################
#### Make gene sets
##########################################################################################################################
 
## Load some gene sets from MSigDB as well
gene_sets_df <- msigdbr(species = "Rattus norvegicus")
C2_gene_set <- msigdbr(species = "Rattus norvegicus", category = "C2") %>% dplyr::select(gs_name, gene_symbol)
C3_gene_set <- msigdbr(species = "Rattus norvegicus", category = "C3") %>% dplyr::select(gs_name, gene_symbol)
C7_gene_set <- msigdbr(species = "Rattus norvegicus", category = "C7") %>% dplyr::select(gs_name, gene_symbol)
C8_gene_set <- msigdbr(species = "Rattus norvegicus", category = "C8") %>% dplyr::select(gs_name, gene_symbol)
# 
# # # Select gene sets from C2 with KEGG in the name
C2_KEGG <- C2_gene_set %>% filter(grepl("KEGG", gs_name))

### KEGG gene sets - download broke down last time I tried
#!/usr/bin/Rscript


pathways.list.rno <- keggList("pathway", organism = "rno")
pathway.codes.rno <- sub("path:", "", names(pathways.list.rno))


# genes.by.pathway <- sapply(pathway.codes.rno,
#                            function(pwid){
#                              pw <- keggGet(pwid)
#                              if (is.null(pw[[1]]$GENE)) return(NA)
#                              pw2 <- pw[[1]]$GENE[c(FALSE,TRUE)] # may need to modify this to c(FALSE, TRUE) for other organisms
#                              pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
# 
#                              #names(pw2) <- pw[[1]]$NAME
#                              return(pw2)
#                            }
# )
# 
# pathway.to.name <- sapply(pathway.codes.rno,
#                           function(pwid){
#                             pw <- keggGet(pwid)
#                             if (is.null(pw[[1]]$NAME)) return(NA)
#                             name <- pw[[1]]$NAME
#                             # Remove ' - Rattus norvegicus (rat)' from the end of the name
#                             name <- gsub(" - Rattus norvegicus \\(rat\\)", "", name)
#                             return(name)
#                           }
# )

# Save genes.by.pathway and pathway.to.name objects
# save(genes.by.pathway, file = "pathway.to.gene.rda")
# save(pathway.to.name, file = "pathway.to.name.rda")

# Load genes.by.pathway and pathway.to.name objects
load("pathway.to.gene.rda")
load("pathway.to.name.rda")

head(genes.by.pathway)

# Convert the list to a dataframe
pathways_df <- data.frame(
  pathway = rep(names(genes.by.pathway), sapply(genes.by.pathway, length)),
  gene = unlist(genes.by.pathway)
)

pathway_names_df <- data.frame(
  pathway_id = names(pathway.to.name),
  pathway_name = unlist(pathway.to.name)
)


pathways_df_KEGG <- merge(pathways_df, pathway_names_df, by.x = 'pathway', by.y = 'pathway_id')

pathways_df_KEGG$pathway <- NULL
colnames(pathways_df_KEGG) <- c('gene_symbol', 'gs_name'); pathways_df_KEGG <- pathways_df_KEGG[,c(2,1)]

##### Gene sets from literature

setwd('/home/boris/Documents/PhD/gut_brain/VNS_Helen/RE-ANALYSIS/')

# Astrocyte secreted proteins from Dowel 2010 Journal of Proteome research
Dowell_astrocyte_secreted <- fread('./results/Literature_gene_sets/Dowell2010_IDs_mapped.tsv')
intersect(Dowell_astrocyte_secreted$To, protein_gene$`Gene_names`) # 0 -> accession IDs to gene names using UniProt
Dowell_astrocyte_secreted$gs_name <- 'Dowell2010_Astrocyte_secreted'; colnames(Dowell_astrocyte_secreted) <- c('acc', 'gene_symbol', 'gs_name')
Dowell_astrocyte_secreted <- Dowell_astrocyte_secreted[,c(2,3)]

# Oligodendrocyte proteome from Dumont 2007
Dumont_oligodendrocytes <- fread('./results/Literature_gene_sets/Dumont2007_oligodendrocytes.csv')
# Remove rows without Gene_symbol, they have an empty string
Dumont_oligodendrocytes <- Dumont_oligodendrocytes[Dumont_oligodendrocytes$Gene_symbol != '',]
#Dumont_oligodendrocytes <- Dumont_oligodendrocytes[!is.na(Dumont_oligodendrocytes$Gene_symbol),]
Dumont_oligodendrocytes$gs_name <- 'Dumont2007_oligodendrocytes'; 
Dumont_oligodendrocytes <- Dumont_oligodendrocytes[,c('Gene_symbol','gs_name')]; 
colnames(Dumont_oligodendrocytes) <- c('gene_symbol', 'gs_name')

# We need to map human genes to mouse orthologs, I will use the table provided by the mouse metabolic atlas for this https://github.com/SysBioChalmers/Mouse-GEM
human2mouse <- fread('/home/boris/Documents/PhD/resources/models/mouse-GEM/human2MouseOrthologs.tsv')[, c(2,4)]

# Santiago 2023 microglia proteomics (multiple cell states)
Santiago <- fread('./results/Literature_gene_sets/Santiago2023_cell_proteomes.csv')

Santiago_microglia_protective_up <- Santiago[(Santiago$`IL10-Control` <= 0.05 & Santiago$`diff IL10-Control` > 0), 'V1']
Santiago_microglia_protective_down <- Santiago[(Santiago$`IL10-Control` <= 0.05 & Santiago$`diff IL10-Control` < 0), 'V1']
Santiago_microglia_inflammatory_up <- Santiago[(Santiago$`LPS-Control` <= 0.05 & Santiago$`diff LPS-Control` > 0), 'V1']
Santiago_microglia_inflammatory_down <- Santiago[(Santiago$`LPS-Control` <= 0.05 & Santiago$`diff LPS-Control` < 0), 'V1']
Santiago_microglia_homeostatic_up <- Santiago[(Santiago$`TGFB-Control` <= 0.05 & Santiago$`diff TGFB-Control` > 0), 'V1']
Santiago_microglia_homeostatic_down <- Santiago[(Santiago$`TGFB-Control` <= 0.05 & Santiago$`diff TGFB-Control` < 0), 'V1']

Santiago_microglia_protective_up <- data.table(gene_symbol = Santiago_microglia_protective_up, gs_name = 'Santiago_microglia_protective_up')
Santiago_microglia_protective_down <- data.table(gene_symbol = Santiago_microglia_protective_down, gs_name = 'Santiago_microglia_protective_down')
Santiago_microglia_inflammatory_up <- data.table(gene_symbol = Santiago_microglia_inflammatory_up, gs_name = 'Santiago_microglia_inflammatory_up')
Santiago_microglia_inflammatory_down <- data.table(gene_symbol = Santiago_microglia_inflammatory_down, gs_name = 'Santiago_microglia_inflammatory_down')
Santiago_microglia_homeostatic_up <- data.table(gene_symbol = Santiago_microglia_homeostatic_up, gs_name = 'Santiago_microglia_homeostatic_up')
Santiago_microglia_homeostatic_down <- data.table(gene_symbol = Santiago_microglia_homeostatic_down, gs_name = 'Santiago_microglia_homeostatic_down')

colnames(Santiago_microglia_protective_up) <- c('gene_symbol', 'gs_name')
colnames(Santiago_microglia_protective_down) <- c('gene_symbol', 'gs_name')
colnames(Santiago_microglia_inflammatory_up) <- c('gene_symbol', 'gs_name')
colnames(Santiago_microglia_inflammatory_down) <- c('gene_symbol', 'gs_name')
colnames(Santiago_microglia_homeostatic_up) <- c('gene_symbol', 'gs_name')
colnames(Santiago_microglia_homeostatic_down) <- c('gene_symbol', 'gs_name')


Microglia_Zhong <- fread('./results/Literature_gene_sets/ZHONG_PFC_C3_MICROGLIA.v2023.2.Hs.tsv')
genes_Zhong <- as.character(Microglia_Zhong[17,2])
genes_Zhong <- strsplit(genes_Zhong, split = ',')
# Create df with columns gene symbol and gs_name
Zhong_microglia <- data.table(gene_symbol = unlist(genes_Zhong), gs_name = 'Zhong_microglia')
# Map to Human proteins
Zhong_microglia <- merge(Zhong_microglia, human2mouse, by.x = 'gene_symbol', by.y = 'fromSymbol')
Zhong_microglia <- Zhong_microglia[,c(3,2)]; colnames(Zhong_microglia) <- c('gene_symbol', 'gs_name')


gene_sets_list <- list(Zhong_microglia,
                        Santiago_microglia_protective_up, Santiago_microglia_protective_down, 
                        Santiago_microglia_inflammatory_up, Santiago_microglia_inflammatory_down, 
                        Santiago_microglia_homeostatic_up, Santiago_microglia_homeostatic_down,
                        Dowell_astrocyte_secreted,
                        Dumont_oligodendrocytes)

gene_sets_literature <- rbindlist(gene_sets_list)
gene_sets_literature <- gene_sets_literature[,c(2,1)] # Rearrange columns

##########################################################################################################################
#### GSEA
##########################################################################################################################

DE_gene_list <- DE_cVNS_vs_DE_sham$t
names(DE_gene_list) <- DE_cVNS_vs_DE_sham$Gene_names
RE_gene_list <- RE_cVNS_vs_RE_sham$t
names(RE_gene_list) <- RE_cVNS_vs_RE_sham$Gene_names
sham_gene_list <- DE_sham_vs_RE_sham$t
names(sham_gene_list) <- DE_sham_vs_RE_sham$Gene_names
cVNS_gene_list <- DE_cVNS_vs_RE_cVNS$t
names(cVNS_gene_list) <- DE_cVNS_vs_RE_cVNS$Gene_names

# Sort in decreasing order
DE_gene_list = sort(DE_gene_list, decreasing = TRUE)
RE_gene_list = sort(RE_gene_list, decreasing = TRUE)
sham_gene_list = sort(sham_gene_list, decreasing = TRUE)
cVNS_gene_list = sort(cVNS_gene_list, decreasing = TRUE)

# Convert gene names to Entrez IDs
DE_gene_list_mapping <- bitr(names(DE_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)
RE_gene_list_mapping <- bitr(names(RE_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)
sham_gene_list_mapping <- bitr(names(sham_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)
cVNS_mapping_mapping <- bitr(names(cVNS_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)

DE_cVNS_vs_DE_sham_extended <- merge(DE_cVNS_vs_DE_sham, DE_gene_list_mapping, by.x = 'Gene_names', by.y = 'SYMBOL')
RE_cVNS_vs_RE_sham_extended <- merge(RE_cVNS_vs_RE_sham, RE_gene_list_mapping, by.x = 'Gene_names', by.y = 'SYMBOL')
DE_sham_vs_RE_sham_extended <- merge(DE_sham_vs_RE_sham, sham_gene_list_mapping, by.x = 'Gene_names', by.y = 'SYMBOL')
DE_cVNS_vs_RE_cVNS_extended <- merge(DE_cVNS_vs_RE_cVNS, cVNS_mapping_mapping, by.x = 'Gene_names', by.y = 'SYMBOL')

DE_gene_list_extended <- DE_cVNS_vs_DE_sham_extended$t
names(DE_gene_list_extended) <- DE_cVNS_vs_DE_sham_extended$ENTREZID
RE_gene_list_extended <- RE_cVNS_vs_RE_sham_extended$t
names(RE_gene_list_extended) <- RE_cVNS_vs_RE_sham_extended$ENTREZID
sham_gene_list_extended <- DE_sham_vs_RE_sham_extended$t
names(sham_gene_list_extended) <- DE_sham_vs_RE_sham_extended$ENTREZID
cVNS_gene_list_extended <- DE_cVNS_vs_RE_cVNS_extended$t
names(cVNS_gene_list_extended) <- DE_cVNS_vs_RE_cVNS_extended$ENTREZID

# Sort in decreasing order
DE_gene_list_extended = sort(DE_gene_list_extended, decreasing = TRUE)
RE_gene_list_extended = sort(RE_gene_list_extended, decreasing = TRUE)
sham_gene_list_extended = sort(sham_gene_list_extended, decreasing = TRUE)
cVNS_gene_list_extended = sort(cVNS_gene_list_extended, decreasing = TRUE)


# Get the position of the specified entry
# position_DE <- which(names(DE_gene_list) == target_entry)
# position_RE <- which(names(RE_gene_list) == target_entry)
# position_sham <- which(names(sham_gene_list) == target_entry)
# position_cVNS <- which(names(cVNS_gene_list) == target_entry)

# keyType This is the source of the annotation (gene ids). The options vary for each annotation. 

organism<-org.Rn.eg.db
#keytypes(organism)
#head(keys(organism, keytype="SYMBOL"))

setwd('./results/limma/')

##############################################################
#### Demyelination

DE_gse_all <- gseGO(geneList=DE_gene_list, 
                    ont ="ALL", # Which ontology: cellular component, molecular function, biological process
                    keyType = "SYMBOL", #genesymbol
                    #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                    nPermSimple = 1000,
                    minGSSize = 3, # Minimum number of genes for enrichment
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = organism, 
                    pAdjustMethod = "BH")

pdf('./DE_cVNS_vs_DE_sham/DE_cVNS_vs_DE_sham_GO_all_full_ranked_geneset.pdf')
dotplot(DE_gse_all, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=9)) + theme(axis.text.x = element_text(size=9))
dev.off()


DE_gse_all_BP <- gseGO(geneList=DE_gene_list, 
                       ont ="BP", # Which ontology: cellular component, molecular function, biological process
                       keyType = "SYMBOL", #genesymbol
                       #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                       nPermSimple = 1000,
                       minGSSize = 3, # Minimum number of genes for enrichment
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

pdf('./DE_cVNS_vs_DE_sham/DE_GO_BP_full_ranked_geneset.pdf')
dotplot(DE_gse_all_BP, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

DE_gse_all_CC <- gseGO(geneList=DE_gene_list,
                       ont ="CC", # Which ontology: cellular component, molecular function, biological process
                       keyType = "SYMBOL", #genesymbol
                       #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                       nPermSimple = 1000,
                       minGSSize = 3, # Minimum number of genes for enrichment
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

pdf('./DE_cVNS_vs_DE_sham/DE_GO_CC_full_ranked_geneset.pdf')
dotplot(DE_gse_all_CC, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

DE_gse_all_MF <- gseGO(geneList=DE_gene_list,
                       ont ="MF", # Which ontology: cellular component, molecular function, biological process
                       keyType = "SYMBOL", #genesymbol
                       #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                       nPermSimple = 1000,
                       minGSSize = 3, # Minimum number of genes for enrichment
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

pdf('./DE_cVNS_vs_DE_sham/DE_GO_MF_full_ranked_geneset.pdf')
dotplot(DE_gse_all_MF, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()


DE_C2 <- GSEA(geneList = DE_gene_list,
              TERM2GENE = C2_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
              )

DE_C3 <- GSEA(geneList = DE_gene_list,
              TERM2GENE = C3_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
              )
pdf('./DE_cVNS_vs_DE_sham/DE_C3_full_ranked_geneset.pdf')
dotplot(DE_C3, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

DE_C7 <- GSEA(geneList = DE_gene_list,
              TERM2GENE = C7_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "none",
              )

# pdf('./DE_C7_full_ranked_geneset.pdf')
# dotplot(DE_C7, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
# dev.off()

DE_C8 <- GSEA(geneList = DE_gene_list,
              TERM2GENE = C8_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
)

# DE_C2_KEGG <- GSEA(geneList = DE_gene_list,
#               TERM2GENE = C2_KEGG,
#               nPerm = 1000,
#               minGSSize = 3,
#               maxGSSize = 800,
#               pvalueCutoff = 0.05,
#               verbose = TRUE,
#               seed = 1234,
#               pAdjustMethod = "none",
# )
# 
# pdf('./DE_cVNS_vs_DE_sham/DE_C2_KEGG_full_ranked_geneset_Pnotadjusted.pdf')
# dotplot(DE_C2_KEGG, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
# dev.off()


# DE_KEGG <- GSEA(geneList = DE_gene_list,
#                 TERM2GENE = pathways_df_KEGG,
#                 nPerm = 1000,
#                 minGSSize = 3,
#                 maxGSSize = 800,
#                 pvalueCutoff = 0.05,
#                 verbose = TRUE,
#                 seed = 1234,
#                 pAdjustMethod = "none",
# )
# 
# library(fgsea)
# 
# pdf('./DE_cVNS_vs_DE_sham/DE_KEGG_full_ranked_geneset_Pnotadjusted.pdf')
# dotplot(DE_KEGG, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
# dev.off()

# pdf('./DE_C8_full_ranked_geneset.pdf')
# dotplot(DE_C8, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
# dev.off()

DE_literature_gene_sets <- GSEA(geneList = DE_gene_list,
                                TERM2GENE = gene_sets_literature,
                                nPerm = 1000,
                                minGSSize = 3,
                                maxGSSize = 800,
                                pvalueCutoff = 0.05,
                                verbose = TRUE,
                                seed = 1234,
                                pAdjustMethod = "BH")

DE_literature_gene_sets_table <- DE_literature_gene_sets@result
# write table
write.table(DE_literature_gene_sets_table, file = "./DE_cVNS_vs_DE_sham/DE_literature_gene_sets_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

pdf('./DE_cVNS_vs_DE_sham/DE_literature_full_ranked_geneset.pdf')
dotplot(DE_literature_gene_sets, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()


pdf('./DE_cVNS_vs_DE_sham/gseaplot_literature_gene_sets.pdf')
gseaplot2(DE_literature_gene_sets, geneSetID = 1:nrow(DE_literature_gene_sets), pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot", base_size = 8, rel_heights = c(2,0.5, 1)) 
dev.off()

pdf('./DE_cVNS_vs_DE_sham/gseaplot_literature_gene_sets_noTable.pdf')
gseaplot2(DE_literature_gene_sets, geneSetID = 1:nrow(DE_literature_gene_sets), pvalue_table = FALSE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot", base_size = 8, rel_heights = c(2,0.5, 1)) 
dev.off()

# Reactome
y <- gsePathway(DE_gene_list_extended, 
                organism = "rat",
                
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)

pdf('./DE_cVNS_vs_DE_sham/DE_cVNS_vs_DE_sham_Reactome_enrichment.pdf')
dotplot(y, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

##############################################################
#### Remyelination

RE_gse_all <- gseGO(geneList=RE_gene_list,
                    ont ="ALL", # Which ontology: cellular component, molecular function, biological process
                    keyType = "SYMBOL", #genesymbol
                    #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                    nPermSimple = 1000,
                    minGSSize = 3, # Minimum number of genes for enrichment
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = organism, 
                    pAdjustMethod = "BH")

pdf('./RE_cVNS_vs_RE_sham/RRE_cVNS_vs_RE_sham_GO_all_full_ranked_geneset.pdf')
dotplot(RE_gse_all, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=9))
dev.off()

RE_gse_all_CC <- gseGO(geneList=RE_gene_list,
                       ont ="CC", # Which ontology: cellular component, molecular function, biological process
                       keyType = "SYMBOL", #genesymbol
                       #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                       nPermSimple = 1000,
                       minGSSize = 3, # Minimum number of genes for enrichment
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

pdf('./RE_cVNS_vs_RE_sham/RE_GO_CC_full_ranked_geneset.pdf')
dotplot(RE_gse_all_CC, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

RE_gse_all_BP <- gseGO(geneList=RE_gene_list,
                       ont ="BP", # Which ontology: cellular component, molecular function, biological process
                       keyType = "SYMBOL", #genesymbol
                       #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                       nPermSimple = 1000,
                       minGSSize = 3, # Minimum number of genes for enrichment
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

pdf('./RE_cVNS_vs_RE_sham/RE_GO_BP_full_ranked_geneset.pdf')
dotplot(RE_gse_all_BP, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()


RE_gse_all_MF <- gseGO(geneList=RE_gene_list,
                       ont ="MF", # Which ontology: cellular component, molecular function, biological process
                       keyType = "SYMBOL", #genesymbol
                       #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                       nPermSimple = 1000,
                       minGSSize = 3, # Minimum number of genes for enrichment
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

pdf('./RE_cVNS_vs_RE_sham/RE_GO_MF_full_ranked_geneset.pdf')
dotplot(RE_gse_all_MF, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()


RE_C2 <- GSEA(geneList = RE_gene_list,
              TERM2GENE = C2_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
)

RE_C3 <- GSEA(geneList = RE_gene_list,
              TERM2GENE = C3_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
)
# pdf('./RE_C3_full_ranked_geneset.pdf')
# dotplot(RE_C3, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
# dev.off()

RE_C7 <- GSEA(geneList = RE_gene_list,
              TERM2GENE = C7_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
)

# pdf('./RE_C7_full_ranked_geneset.pdf')
# dotplot(RE_C7, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
# dev.off()

RE_C8 <- GSEA(geneList = RE_gene_list,
              TERM2GENE = C8_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
)

# pdf('./RE_C8_full_ranked_geneset.pdf')
# dotplot(RE_C8, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
# dev.off()

RE_C2_KEGG <- GSEA(geneList = RE_gene_list,
                   TERM2GENE = C2_KEGG,
                   nPerm = 1000,
                   minGSSize = 3,
                   maxGSSize = 800,
                   pvalueCutoff = 0.05,
                   verbose = TRUE,
                   seed = 1234,
                   pAdjustMethod = "none",
)

pdf('./RE_cVNS_vs_RE_sham/RE_C2_KEGG_full_ranked_geneset_Pnotadjusted.pdf')
dotplot(RE_C2_KEGG, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

RE_KEGG <- GSEA(geneList = RE_gene_list,
                TERM2GENE = pathways_df_KEGG,
                nPerm = 1000,
                minGSSize = 3,
                maxGSSize = 800,
                pvalueCutoff = 0.05,
                verbose = TRUE,
                seed = 1234,
                pAdjustMethod = "none",
)

pdf('./RE_cVNS_vs_RE_sham/RE_KEGG_full_ranked_geneset_Pnotadjusted.pdf')
dotplot(RE_KEGG, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()



RE_literature_gene_sets <- GSEA(geneList = RE_gene_list,
                                TERM2GENE = gene_sets_literature,
                                nPerm = 1000,
                                minGSSize = 3,
                                maxGSSize = 800,
                                pvalueCutoff = 0.05,
                                verbose = TRUE,
                                seed = 1234,
                                pAdjustMethod = "BH")

RE_literature_gene_sets_table <- RE_literature_gene_sets@result
# werite table
write.table(RE_literature_gene_sets_table, file = "./RE_cVNS_vs_RE_sham/RE_literature_gene_sets_table.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

pdf('./RE_cVNS_vs_RE_sham/RE_literature_full_ranked_geneset.pdf')
dotplot(RE_literature_gene_sets, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

#tmp <- data.frame(RE_literature_gene_sets@result)

pdf('./RE_cVNS_vs_RE_sham/gseaplot_literature_gene_sets.pdf')
gseaplot2(RE_literature_gene_sets, geneSetID = 1:nrow(RE_literature_gene_sets_table), pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")
dev.off()


y <- gsePathway(RE_gene_list_extended, 
                organism = "rat",
                
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)

pdf('./RE_cVNS_vs_RE_sham/RE_cVNS_vs_RE_sham_Reactome_enrichment.pdf')
dotplot(y, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()


##############################################################
#### sham


DE_gse_all <- gseGO(geneList=sham_gene_list, 
                    ont ="ALL", # Which ontology: cellular component, molecular function, biological process
                    keyType = "SYMBOL", #genesymbol
                    #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                    nPermSimple = 1000,
                    minGSSize = 3, # Minimum number of genes for enrichment
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = organism, 
                    pAdjustMethod = "BH")

pdf('./DE_sham_vs_RE_sham/sham_GO_all_full_ranked_geneset.pdf')
dotplot(DE_gse_all, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

DE_gse_all_BP <- gseGO(geneList=sham_gene_list, 
                       ont ="BP", # Which ontology: cellular component, molecular function, biological process
                       keyType = "SYMBOL", #genesymbol
                       #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                       nPermSimple = 1000,
                       minGSSize = 3, # Minimum number of genes for enrichment
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

pdf('./DE_sham_vs_RE_sham/sham_GO_BP_full_ranked_geneset.pdf')
dotplot(DE_gse_all_BP, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

DE_gse_all_CC <- gseGO(geneList=sham_gene_list,
                       ont ="CC", # Which ontology: cellular component, molecular function, biological process
                       keyType = "SYMBOL", #genesymbol
                       #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                       nPermSimple = 1000,
                       minGSSize = 3, # Minimum number of genes for enrichment
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

pdf('./DE_sham_vs_RE_sham/sham_GO_CC_full_ranked_geneset.pdf')
dotplot(DE_gse_all_CC, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

DE_gse_all_MF <- gseGO(geneList=sham_gene_list,
                       ont ="MF", # Which ontology: cellular component, molecular function, biological process
                       keyType = "SYMBOL", #genesymbol
                       #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                       nPermSimple = 1000,
                       minGSSize = 3, # Minimum number of genes for enrichment
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

pdf('./DE_sham_vs_RE_sham/sham_GO_MF_full_ranked_geneset.pdf')
dotplot(DE_gse_all_MF, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()


sham_C2 <- GSEA(geneList = sham_gene_list,
              TERM2GENE = C2_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
)

sham_C3 <- GSEA(geneList = sham_gene_list,
              TERM2GENE = C3_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
)
# pdf('./RE_C3_full_ranked_geneset.pdf')
# dotplot(RE_C3, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
# dev.off()

sham_C7 <- GSEA(geneList = sham_gene_list,
              TERM2GENE = C7_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
)

# pdf('./RE_C7_full_ranked_geneset.pdf')
# dotplot(sham_C7, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
# dev.off()

sham_C8 <- GSEA(geneList = sham_gene_list,
              TERM2GENE = C8_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
)


sham_C2_KEGG <- GSEA(geneList = sham_gene_list,
                   TERM2GENE = C2_KEGG,
                   nPerm = 1000,
                   minGSSize = 3,
                   maxGSSize = 800,
                   pvalueCutoff = 0.05,
                   verbose = TRUE,
                   seed = 1234,
                   pAdjustMethod = "none",
)

pdf('./DE_sham_vs_RE_sham/sham_C2_KEGG_full_ranked_geneset_Pnotadjusted.pdf')
dotplot(sham_C2_KEGG, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()


sham_KEGG <- GSEA(geneList = sham_gene_list,
                TERM2GENE = pathways_df_KEGG,
                nPerm = 1000,
                minGSSize = 3,
                maxGSSize = 800,
                pvalueCutoff = 0.05,
                verbose = TRUE,
                seed = 1234,
                pAdjustMethod = "BH",
)

# pdf('./DE_sham_vs_RE_sham/sham_KEGG_full_ranked_geneset_Pnotadjusted.pdf')
# dotplot(sham_KEGG, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
# dev.off()



sham_literature_gene_sets <- GSEA(geneList = sham_gene_list,
                                TERM2GENE = gene_sets_literature,
                                nPerm = 1000,
                                minGSSize = 3,
                                maxGSSize = 800,
                                pvalueCutoff = 0.05,
                                verbose = TRUE,
                                seed = 1234,
                                pAdjustMethod = "BH")


#sham_literature_gene_sets_table <- data.frame(sham_literature_gene_sets@result)

# pdf('./DE_sham_vs_RE_sham/RE_literature_full_ranked_geneset.pdf')
# dotplot(sham_literature_gene_sets, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
# dev.off()

# tmp <- data.frame(sham_literature_gene_sets@result)
# 
# pdf('./DE_sham_vs_RE_sham/gseaplot_literature_gene_sets.pdf')
# gseaplot2(sham_literature_gene_sets, geneSetID = 1:nrow(sham_literature_gene_sets_table), pvalue_table = TRUE,
#           color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")
# dev.off()

y <- gsePathway(sham_gene_list_extended, 
                organism = "rat",
                
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)

pdf('./DE_sham_vs_RE_sham/Reactome_enrichment.pdf')
dotplot(y, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()



##################################################################################################
### cVNS

DE_gse_all <- gseGO(geneList=cVNS_gene_list, 
                    ont ="ALL", # Which ontology: cellular component, molecular function, biological process
                    keyType = "SYMBOL", #genesymbol
                    #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                    nPermSimple = 1000,
                    minGSSize = 3, # Minimum number of genes for enrichment
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = organism, 
                    pAdjustMethod = "BH")

pdf('./DE_cVNS_vs_RE_cVNS/cVNS_GO_all_full_ranked_geneset.pdf')
dotplot(DE_gse_all, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

DE_gse_all_BP <- gseGO(geneList=cVNS_gene_list, 
                       ont ="BP", # Which ontology: cellular component, molecular function, biological process
                       keyType = "SYMBOL", #genesymbol
                       #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                       nPermSimple = 1000,
                       minGSSize = 3, # Minimum number of genes for enrichment
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

pdf('./DE_cVNS_vs_RE_cVNS/cVNS_GO_BP_full_ranked_geneset.pdf')
dotplot(DE_gse_all_BP, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

DE_gse_all_CC <- gseGO(geneList=cVNS_gene_list,
                       ont ="CC", # Which ontology: cellular component, molecular function, biological process
                       keyType = "SYMBOL", #genesymbol
                       #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                       nPermSimple = 1000,
                       minGSSize = 3, # Minimum number of genes for enrichment
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

pdf('./DE_cVNS_vs_RE_cVNS/cVNS_GO_CC_full_ranked_geneset.pdf')
dotplot(DE_gse_all_CC, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

DE_gse_all_MF <- gseGO(geneList=cVNS_gene_list,
                       ont ="MF", # Which ontology: cellular component, molecular function, biological process
                       keyType = "SYMBOL", #genesymbol
                       #nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                       nPermSimple = 1000,
                       minGSSize = 3, # Minimum number of genes for enrichment
                       maxGSSize = 800, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

pdf('./DE_cVNS_vs_RE_cVNS/cVNS_GO_MF_full_ranked_geneset.pdf')
dotplot(DE_gse_all_MF, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()


DE_C2 <- GSEA(geneList = cVNS_gene_list,
              TERM2GENE = C2_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
)

DE_C3 <- GSEA(geneList = cVNS_gene_list,
              TERM2GENE = C3_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
)
# pdf('./cVNS_C3_full_ranked_geneset.pdf')
# dotplot(DE_C3, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
# dev.off()

DE_C7 <- GSEA(geneList = cVNS_gene_list,
              TERM2GENE = C7_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
)

# pdf('./DE_C7_full_ranked_geneset.pdf')
# dotplot(DE_C7, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
# dev.off()

DE_C8 <- GSEA(geneList = cVNS_gene_list,
              TERM2GENE = C8_gene_set,
              nPerm = 1000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              seed = 1234,
              pAdjustMethod = "BH",
)

# pdf('./DE_C8_full_ranked_geneset.pdf')
# dotplot(DE_C8, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
# dev.off()

cVNS_C2_KEGG <- GSEA(geneList = cVNS_gene_list,
                   TERM2GENE = C2_KEGG,
                   nPerm = 1000,
                   minGSSize = 3,
                   maxGSSize = 800,
                   pvalueCutoff = 0.05,
                   verbose = TRUE,
                   seed = 1234,
                   pAdjustMethod = "none",
)

pdf('./DE_cVNS_vs_RE_cVNS/cVNS_C2_KEGG_full_ranked_geneset_Pnotadjusted.pdf')
dotplot(cVNS_C2_KEGG, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()


cVNS_KEGG <- GSEA(geneList = cVNS_gene_list,
                  TERM2GENE = pathways_df_KEGG,
                  nPerm = 1000,
                  minGSSize = 3,
                  maxGSSize = 800,
                  pvalueCutoff = 0.05,
                  verbose = TRUE,
                  seed = 1234,
                  pAdjustMethod = "none",
)

pdf('./DE_cVNS_vs_RE_cVNS/cVNS_KEGG_full_ranked_geneset_Pnotadjusted.pdf')
dotplot(cVNS_KEGG, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()


# pdf('./DE_cVNS_vs_RE_cVNS/Astrocytes_cVNS_enrichment_plot.pdf')
# gseaplot(cVNS_KEGG, geneSetID = 25, by = "runningScore", title = DE_KEGG$Description[25])
# dev.off()


cVNS_literature_gene_sets <- GSEA(geneList = cVNS_gene_list,
                                  TERM2GENE = gene_sets_literature,
                                  nPerm = 1000,
                                  minGSSize = 3,
                                  maxGSSize = 800,
                                  pvalueCutoff = 0.05,
                                  verbose = TRUE,
                                  seed = 1234,
                                  pAdjustMethod = "BH")

cVNS_literature_gene_sets_table <- data.frame(cVNS_literature_gene_sets@result)

pdf('./DE_cVNS_vs_RE_cVNS/literature_full_ranked_geneset.pdf')
dotplot(cVNS_literature_gene_sets, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()


pdf('./DE_cVNS_vs_RE_cVNS/gseaplot_literature_gene_sets.pdf')
gseaplot2(cVNS_literature_gene_sets, geneSetID = 1:nrow(cVNS_literature_gene_sets_table), pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")
dev.off()


cVNS_mapping <- bitr(names(cVNS_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)
# Replace gene symbols with Entrez IDs as names for cVNS_gene_list


y <- gsePathway(cVNS_gene_list_extended, 
                organism = "rat",
                
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)

pdf('./DE_cVNS_vs_RE_cVNS/Reactome_enrichment.pdf')
dotplot(y, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()



