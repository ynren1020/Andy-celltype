############2020-06-04###################
# create gene sets file for g:Profiler###
# gene list are snf 3 groups gene list ##
# logFC >1 and adj.p <0.05 ##############

library(readxl)
library(dplyr)
library(tidyverse)

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
    
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    humanx <- unique(genesV2[, 2])
    
    # Print the first 6 genes found to the screen
    print(head(humanx))
    return(humanx)
}


# data ---
path <- "Pijuan.xlsx"
path <- "Brink.xlsx"

list_all <- path %>% 
    excel_sheets() %>% 
    set_names() %>% 
    map(read_excel, path = path)


genes <- NULL

for (i in 1:length(list_all)) {

#genes[i] <- paste0(convertMouseGeneList(list_all[[i]]$...2), collapse = ";")
genes[i] <- paste0(convertMouseGeneList(list_all[[i]]$names), collapse = ";")


}

df4geneset <- data.frame(genesets = names(list_all), description = "na", genes = genes)

write.table(df4geneset, "Pijuan_geneset.txt", sep = "\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE)

write.table(df4geneset, "Brink_geneset.txt", sep = "\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE)
# Brink study ---
path <- "Brink.xlsx"

list_all2 <- path %>% 
    excel_sheets() %>% 
    set_names() %>% 
    map(read_excel, path = path)


nmp <- intersect(convertMouseGeneList(list_all$`27-NMP`$...2), convertMouseGeneList(list_all2$cl_07$names))
#> nmp
#[1] "LHPP"  "CHD3"  "HOXA3" "PSIP1" "ACOT7" "SP8"   "SMIM3" "EPHA5" "HOXB9"
#[10] "HOXC8" "HOXC4" "HOXB8"


# venn diagram ---
library(VennDiagram)

# Generate 3 sets of 200 words
pijuan_27 <- convertMouseGeneList(list_all$`27-NMP`$...2)
brink_7 <- convertMouseGeneList(list_all2$cl_07$names)

# Chart
venn.diagram(
    x = list(pijuan_27, brink_7),
    category.names = c("pijuan_27" , "brink_7 "),
    filename = 'pijuan27&brink7_venn_diagramm.png',
    output=TRUE
)





