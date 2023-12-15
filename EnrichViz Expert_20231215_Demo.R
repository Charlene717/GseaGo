# Load required libraries
if(!require("DOSE")) install.packages("DOSE"); library(DOSE)
if(!require("clusterProfiler")) install.packages("clusterProfiler"); library(clusterProfiler)
if(!require("org.Hs.eg.db")) install.packages("org.Hs.eg.db"); library(org.Hs.eg.db)
if(!require("enrichplot")) install.packages("enrichplot"); library(enrichplot)

# Provided gene list
geneList <- c("APH1A", "ARRB1", "CCND1", "CUL1", "DLL1", "DTX1", "DTX2", "DTX4",
              "FBXW11", "FZD1", "FZD5", "FZD7", "HES1", "HEYL", "JAG1", "KAT2A",
              "LFNG", "MAML2", "NOTCH1", "NOTCH2", "NOTCH3", "PPARD", "PRKCA",
              "PSEN2", "PSENEN", "RBX1", "SAP30", "SKP1", "ST3GAL6", "TCF7L2",
              "WNT2", "WNT5A", "TOP2A", "TP53", "BRCA")


# Convert gene symbols to Entrez IDs
entrez_ids <- unlist(mget(geneList, org.Hs.egSYMBOL2EG, ifnotfound=NA))
entrez_ids <- na.omit(entrez_ids) # Remove NAs

# Check if there are any mapped genes
if (length(entrez_ids) == 0) {
  stop("No genes were mapped. Please check your gene symbols.")
}

# Perform enrichment analysis
edo <- enrichDGN(entrez_ids)

## Bar Plot
library(enrichplot)
barplot(edo, showCategory=20)
mutate(edo, qscore = -log(p.adjust, base=10)) %>%
  barplot(x="qscore")

###############################################################################
## Demo: https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#bar-plot

library(DOSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de)

##
# 将Entrez ID转换为符号
gene_symbols <- mapIds(org.Hs.eg.db, keys = de, column = "SYMBOL", keytype = "ENTREZID")

# 输出符号
gene_symbols <- as.character(gene_symbols)


###############################################################################

## Bar Plot
library(enrichplot)
barplot(edo, showCategory=20)
mutate(edo, qscore = -log(p.adjust, base=10)) %>%
  barplot(x="qscore")

## Dot plot
edo2 <- gseDO(entrez_ids)
# edo2 <- gseDO(geneList)
dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")


## Gene-Concept Network
# convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

# 1
p1 <- cnetplot(edox, foldChange=geneList)
# categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

# 2
p1 <- cnetplot(edox, node_label="category",
               cex_label_category = 1.2)
p2 <- cnetplot(edox, node_label="gene",
               cex_label_gene = 0.8)
p3 <- cnetplot(edox, node_label="all")
p4 <- cnetplot(edox, node_label="none",
               color_category='firebrick',
               color_gene='steelblue')
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

# 3
set.seed(123)
x <- list(A = letters[1:10], B=letters[5:12], C=letters[sample(1:26, 15)])
p1 <- cnetplot(x)

set.seed(123)
d <- setNames(rnorm(26), letters)
# p2 <- cnetplot(x, foldChange=d) +
#   scale_color_gradient2(name='associated data', low='darkgreen', high='firebrick')
p2 <- cnetplot(x, foldChange=d)

cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:2])


## Heatmap-like functional classification
p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

## Tree plot
edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

## Enrichment Map
edo <- pairwise_termsim(edo)
p1 <- emapplot(edo)
p2 <- emapplot(edo, cex_category=1.5)
p3 <- emapplot(edo, layout="kk")
p4 <- emapplot(edo, cex_category=1.5,layout="kk")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])


## Biological theme comparison
library(clusterProfiler)
data(gcSample)
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
xx <- pairwise_termsim(xx)
p1 <- emapplot(xx)
p2 <- emapplot(xx, legend_n=2)
p3 <- emapplot(xx, pie="count")
p4 <- emapplot(xx, pie="count", cex_category=1.5, layout="kk")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

## UpSet Plot
upsetplot(edo)
upsetplot(kk2)

## ridgeline plot for expression distribution of GSEA result
ridgeplot(edo2)

## running score and preranked list of GSEA result
p1 <- gseaplot(edo2, geneSetID = 1, by = "runningScore", title = edo2$Description[1])
p2 <- gseaplot(edo2, geneSetID = 1, by = "preranked", title = edo2$Description[1])
p3 <- gseaplot(edo2, geneSetID = 1, title = edo2$Description[1])
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
gseaplot2(edo2, geneSetID = 1, title = edo2$Description[1])
