#Install rdf package; 
more info here: http://www.r-bloggers.com/sparql-with-r-in-less-than-5-minutes/
  install.packages('SPARQL')
install.packages('ggplot2')

#Get library to perform SPARQL queries(library(SPARQL) 
# SPARQL querying package
library(ggplot2)
##Perform query and save results
##Endpoint Wikidata:
endpointwd <- "https://query.wikidata.org/sparql"
##Test query:
querywd <-"SELECT ?item ?itemLabel WHERE { 
?item wdt:P31 wd:Q146. SERVICE wikibase:label 
{ bd:serviceParam wikibase:language 
'[AUTO_LANGUAGE],en'. }}"

qd <- SPARQL(endpointwd,querywd,curl_args=list(useragent=R.version.string))
df <- qd$results
## Cytogenic location!

