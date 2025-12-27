library(GEOquery)
library(limma)
library(edgeR)
library(Biobase)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

########################################################################

gse <- getGEO("GSE50223", GSEMatrix = TRUE, getGPL = TRUE)
gse_data <- gse[[1]]
pheno <- pData(gse_data)
diagnosisv <- factor(pheno$`diagnosis:ch1`, levels = c("H", "P"))
subjectv   <- factor(pheno$`subject:ch1`)
designv <- model.matrix(~ diagnosisv)
corfit <- duplicateCorrelation(exprs(gse_data), designv, block = subjectv)
fit <- lmFit(exprs(gse_data), designv, block = subjectv, correlation = corfit$consensus)
fitfinal <- eBayes(fit)
results <- topTable(fitfinal, coef = "diagnosisvP", number = Inf)
deg_up <- results[results$P.Value < 0.05 & results$logFC > 0, ]
deg_down <- results[results$P.Value < 0.05 & results$logFC < 0, ]
annot <- fData(gse_data)[, c("ID", "Symbol")]
deg_up$ID <- rownames(deg_up)
deg_up_annot <- merge(deg_up, annot, by = "ID")
deg_down$ID <- rownames(deg_down)
deg_down_annot <- merge(deg_down, annot, by = "ID")
hadb <- read.csv("D:/HADB.csv", stringsAsFactors = FALSE)
overlap_up <- deg_up_annot$Symbol[deg_up_annot$Symbol %in% hadb$GENE]
overlap_down <- deg_down_annot$Symbol[deg_down_annot$Symbol %in% hadb$GENE]
deg_all <- rbind(deg_up_annot, deg_down_annot)
article_genes <- c("PIK3C3","BECN1","EDEM1","SIRT2","UVRAG",
                   "ULK2","EIF2AK3","ZFYVE1","CASP4","PELP1",
                   "PEX14","CALCOCO2","CDKN1B","TM9SF1","NPC1",
                   "BNIP3","HGS","GABARAPL2")
overlap_all <- data.frame(Symbol = c(overlap_up, overlap_down))
overlap_all <- merge(overlap_all,
                     deg_all[, c("ID", "Symbol", "logFC", "P.Value", "adj.P.Val")],
                     by = "Symbol",
                     all.x = TRUE)
overlap_match <- overlap_all[overlap_all$Symbol %in% article_genes, ]

##########################################################################

ggplot(overlap_match, aes(x = logFC, y = -log10(P.Value),
                    color = ifelse(logFC > 0, "Up", "Down"))) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text(aes(label = Symbol), vjust = -0.5, size = 3, check_overlap = TRUE) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue")) +
  labs(title = "Volcano Plot of DEGs",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  theme_minimal()

##########################################################################

top_up   <- overlap_match[order(overlap_match$P.Value), ][1:50, ]
top_down <- overlap_match[order(overlap_match$P.Value), ][1:50, ]
top_ids  <- c(top_up$ID, top_down$ID)
valid_ids <- intersect(top_ids, rownames(exprs(gse_data)))
exprs_top <- exprs(gse_data)[valid_ids, ]
exprs_top <- exprs_top[!duplicated(rownames(exprs_top)), ]
id_to_symbol <- overlap_match[, c("ID", "Symbol")]
rownames(exprs_top) <- id_to_symbol$Symbol[match(rownames(exprs_top), rownames(exprs_top))]
exprs_top <- exprs_top[complete.cases(rownames(exprs_top)), ]
annotation_col <- data.frame(
  Type = factor(pheno$`diagnosis:ch1`, levels = c("H","P"))
)
rownames(annotation_col) <- rownames(pheno)
sample_order <- rownames(annotation_col)[order(annotation_col$Type)]
pheatmap(exprs_top[, sample_order],
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         main = "Up- and Down-regulated DEGs")

###########################################################################

genesv <- overlap_match$Symbol
entrez_ids <- bitr(genesv, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
entrez_ids <- unique(entrez_ids$ENTREZID)
ego <- enrichGO(gene          = entrez_ids,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05)
barplot(ego, showCategory = 20, title = "GO Enrichment", x = "Count")
dotplot(ego, showCategory = 20, title = "GO Enrichment", x = "GeneRatio")

##########################################################################

ekegg <- enrichKEGG(gene         = entrez_ids,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)
barplot(ekegg, showCategory = 20, title = "KEGG Pathway Enrichment", x = "Count")
dotplot(ekegg, showCategory = 20, title = "KEGG Pathway Enrichment", x = "GeneRatio")

#########################################################################
rid2symbol <- annot
rownames(rid2symbol) <- rid2symbol$ID
r2 <- results
r2$Symbol <- rid2symbol$Symbol[match(rownames(r2), rid2symbol$ID)]
r2 <- r2[r2$Symbol %in% hadb$GENE, ]
genes_df <- bitr(r2$Symbol,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db)
genes_df <- genes_df[!is.na(genes_df$ENTREZID), ]
genes_df <- genes_df[!duplicated(genes_df$ENTREZID), ]
r2 <- merge(r2, genes_df, by.x = "Symbol", by.y = "SYMBOL", all.x = TRUE)
r2$ENTREZID <- genes_df$ENTREZID[match(r2$Symbol, genes_df$SYMBOL)]
r2_unique <- r2[!is.na(r2$ENTREZID), ]
r2_unique <- r2_unique[order(abs(r2_unique$logFC), decreasing = TRUE), ]
r2_unique <- r2_unique[!duplicated(r2_unique$ENTREZID), ]
gene_list <- r2_unique$logFC
names(gene_list) <- r2_unique$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
gsea_kegg <- gseKEGG(geneList     = gene_list,
                     organism     = "hsa",
                     minGSSize    = 1,
                     pvalueCutoff = 100,
                     scoreType = "std",
                     verbose      = TRUE)
print(gsea_kegg@result$ID)
gseaplot2(gsea_kegg,
          geneSetID = c("hsa04217",  # Autophagy – animal
                        "hsa04142",  # Mitophagy – animal
                        "hsa04659",  # Shigellosis
                        "hsa05022",  # Pathways of neurodegeneration – multiple diseases
                        "hsa04136",  # Autophagy – other
                        "hsa05014",  # Amyotrophic lateral sclerosis (ALS)
                        "hsa05010",  # Alzheimer disease
                        "hsa04068",  # FoxO signaling pathway
                        "hsa05017",  # Spinocerebellar ataxia
                        "hsa05167",  # Kaposi sarcoma-associated herpesvirus infection
                        "hsa04140"), # Autophagy – animal (duplicate entry, same as hsa04217)
          title = "GSEA Enrichment Plot")

###########################################################################

