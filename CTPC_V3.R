# CTPC V3
library(tidyr)
setwd("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3")
#now include more samples
original=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/CTPC_origianl_meta.csv")
pcta=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/pan_cell_line/PCTA_annotiation.csv")
pcta=subset(pcta,(pcta$DepmapModelType=="PRAD")|(pcta$DepmapModelType=="PRSCC"))
table(pcta$Key)
new=setdiff(pcta$geo_accession,original$ID)
#now manual inspection
add=pcta[pcta$geo_accession %in% new,]
write.csv(add,"add_samples.csv")
# now load the confirmed 
add=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/add_samples.csv")
