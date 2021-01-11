# getReactome.R
# download reactome pathway information and return 

library(dplyr)
library(SPARQL)
ENDPOINT <- "http://sparql.wikipathways.org/sparql"
OPTIONS <- NULL
SPARQL_PREFIX <- paste("PREFIX dc: <http://purl.org/dc/elements/1.1/> 
                  PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
                  PREFIX cur: <http://vocabularies.wikipathways.org/wp#Curation:>
                  PREFIX rdfs:    <http://www.w3.org/2000/01/rdf-schema#>
                  PREFIX dcterms: <http://purl.org/dc/terms/> 
                  ")

reactome_query <- paste(SPARQL_PREFIX,
                        "SELECT DISTINCT ?pathwayRef (str(?titleLit) as ?title) (str(?wpidIdentifier) as ?wpid)
  WHERE {
  ?pathwayRef wp:ontologyTag cur:Reactome_Approved ;
           dc:title ?titleLit ;
           dcterms:identifier ?wpidIdentifier .
  }")

reactome_pathways <- SPARQL(ENDPOINT,
                            reactome_query,
                            ns=SPARQL_PREFIX,
                            extra=OPTIONS
)$results

reactome_pathway_list <- list()
no_node_pathway <- NULL
# load Protein/geneProduct/RNA for each pathway
# if no data retrieved: add to a list
for (pathway_index in 1:length(reactome_pathways$wpid)){
  if (pathway_index %% 100 == 0){print(pathway_index)}
  pathway_wpid <- reactome_pathways$wpid[pathway_index]
  one_pathway_query <- paste(SPARQL_PREFIX,
                             "SELECT DISTINCT 
    (str(?wpidIdentifier) as ?wpid) 
    (fn:substring(?ncbiGeneId,33) as ?entrez)
    (fn:substring(?ensId,32) as ?ensembl)
    WHERE {
      VALUES ?type {wp:Protein wp:GeneProduct wp:Rna}
      ?geneProduct a ?type . 
      ?geneProduct dcterms:isPartOf ?pathwayRef .
      ?geneProduct wp:bdbEnsembl ?ensId .
      ?geneProduct wp:bdbEntrezGene ?ncbiGeneId .
      ?pathwayRef a wp:Pathway .
      ?pathwayRef dcterms:identifier",
                             paste("\"", pathway_wpid, "\"", sep=""),
                             "^^xsd:string . 
      ?pathwayRef dcterms:identifier ?wpidIdentifier .}")
  
  one_pathway <- SPARQL(ENDPOINT,
                        one_pathway_query,
                        ns=SPARQL_PREFIX,
                        extra=OPTIONS
  )$results
  reactome_pathway_list[[pathway_index]] <- one_pathway
  if (nrow(one_pathway)==0){
    no_node_pathway <- c(no_node_pathway, pathway_wpid)
  }
}
reactome_pathway_df <- do.call(rbind, reactome_pathway_list)
reactome_pathway_extended <- merge(x=reactome_pathway_df,
                                   y=reactome_pathways,
                                   by="wpid",
                                   all.x=TRUE)
reactome2gene <- unique(reactome_pathway_extended %>% dplyr::select(wpid, entrez))
colnames(reactome2gene) <- c("wpid", "gene")
reactome2name <- unique(reactome_pathway_extended %>% dplyr::select(wpid, title))
colnames(reactome2name) <- c("wpid", "name")
print(head(reactome2gene))
print(head(reactome2name))


