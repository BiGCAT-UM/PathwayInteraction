load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
  print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  status
}
lung.expr <- read.csv(system.file("extdata","data-lung-cancer.csv", package="rWikiPathways"),stringsAsFactors = FALSE)
nrow(lung.expr)
head(lung.expr)
#
up.genes <- lung.expr[lung.expr$log2FC > 1 & lung.expr$adj.P.Value < 0.05, 1] 
dn.genes <- lung.expr[lung.expr$log2FC < -1 & lung.expr$adj.P.Value < 0.05, 1]
bkgd.genes <- lung.expr[,1]
# find out up.genes
up.genes.entrez <- clusterProfiler::bitr(up.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
cat("\n\nWhich column contains my new Entrez IDs?\n")
head(up.genes.entrez)
dn.genes.entrez <- bitr(dn.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez <- bitr(bkgd.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


#
wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt")
wp2gene <- readPathwayGMT(wp.hs.gmt)
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
wpid2gene
wpid2name


ewp.up <- clusterProfiler::enricher(
  gene = up.genes.entrez[[2]],
  universe = bkgd.genes.entrez[[2]],
  pAdjustMethod = "fdr",
  pvalueCutoff = 1, #p.adjust cutoff; relaxed for demo purposes
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)
ewp.up <- DOSE::setReadable(ewp.up, org.Hs.eg.db, keyType = "ENTREZID")
head(ewp.up)

barplot(ewp.up, showCategory = 20)
dotplot(ewp.up, showCategory = 20)
emapplot(ewp.up, showCategory = 20)


############## GSEA ##################
lung.expr$fcsign <- sign(lung.expr$log2FC)
lung.expr$logfdr <- -log10(lung.expr$P.Value)
# Below is used.. signed logfdr
lung.expr$sig <- lung.expr$logfdr/lung.expr$fcsign
sig.lung.expr.entrez<-merge(lung.expr, bkgd.genes.entrez, by.x = "GeneID", by.y = "ENSEMBL")
gsea.sig.lung.expr <- sig.lung.expr.entrez[,8]
names(gsea.sig.lung.expr) <- as.character(sig.lung.expr.entrez[,9])
gsea.sig.lung.expr <- sort(gsea.sig.lung.expr,decreasing = TRUE)

gwp.sig.lung.expr <- clusterProfiler::GSEA(
  gsea.sig.lung.expr,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

gwp.sig.lung.expr.df = data.frame(ID=gwp.sig.lung.expr$ID,
                                  Description=gwp.sig.lung.expr$Description,
                                  enrichmentScore=gwp.sig.lung.expr$enrichmentScore,
                                  NES=gwp.sig.lung.expr$NES,
                                  pvalue=gwp.sig.lung.expr$pvalue,
                                  p.adjust=gwp.sig.lung.expr$p.adjust,
                                  rank=gwp.sig.lung.expr$rank,
                                  leading_edge=gwp.sig.lung.expr$leading_edge
)
nrow(gwp.sig.lung.expr.df[which(gwp.sig.lung.expr.df$NES > 1),]) #pathways enriched for upregulated lung cancer genes
nrow(gwp.sig.lung.expr.df[which(gwp.sig.lung.expr.df$NES < -1),]) #pathways enriched for downregulated lung cancer genes






