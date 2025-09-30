# ===============================
# 1. Load libraries
# ===============================
library(clusterProfiler)
library(org.Hs.eg.db)    # ???? ?????????????? ?????? human
library(enrichplot)
library(DOSE)

# ===============================
# 2. Define your gene list
# ===============================
genes <- c("PTGES", "PLA2G4A", "ALOX5", "PTGS2", "PTGIS", "ALOX15", "PTGES2")

# ===============================
# 3. Convert gene symbols ??? Entrez IDs
# ===============================
entrez <- bitr(genes, fromType = "SYMBOL", 
               toType = "ENTREZID", 
               OrgDb = org.Hs.eg.db)

entrez_ids <- entrez$ENTREZID

# ===============================
# 4. Run GO enrichment
# ===============================
ego <- enrichGO(gene          = entrez_ids,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",   # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2)

# ===============================
# 5. Run KEGG enrichment
# ===============================
ekegg <- enrichKEGG(gene         = entrez_ids,
                    organism     = "hsa",   # hsa = human
                    pvalueCutoff = 0.05)

# ===============================
# 6. Export results to CSV
# ===============================
write.csv(as.data.frame(ego),   "GO_Enrichment.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg), "KEGG_Enrichment.csv", row.names = FALSE)

# ===============================
# 7. Add similarity for network plots
# ===============================
ego_sim   <- pairwise_termsim(ego)
ekegg_sim <- pairwise_termsim(ekegg)

# ===============================
# 8. Visualization
# ===============================

# Dotplots
png("GO_dotplot.png", width=900, height=600)
dotplot(ego, showCategory=15, title="GO BP enrichment (PTGS2 related genes)")
dev.off()

png("KEGG_dotplot.png", width=900, height=600)
dotplot(ekegg, showCategory=15, title="KEGG pathway enrichment")
dev.off()

# Enrichment map plots
png("GO_emapplot.png", width=900, height=600)
emapplot(ego_sim, showCategory=15, layout="kk")
dev.off()

png("KEGG_emapplot.png", width=900, height=600)
emapplot(ekegg_sim, showCategory=15, layout="kk")
dev.off()

# Cnet plots (gene–pathway networks)
png("GO_cnetplot.png", width=900, height=600)
cnetplot(ego_sim, showCategory=10)
dev.off()

png("KEGG_cnetplot.png", width=900, height=600)
cnetplot(ekegg_sim, showCategory=10)
dev.off()

# ===============================
# 9. Preview results in console
# ===============================
head(as.data.frame(ego))
head(as.data.frame(ekegg))
