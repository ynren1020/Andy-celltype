###############2020-06-20###################
## gProfiler_snf*_Brinks_human_allres.txt ##
## Brinks cell type ~ -log10(fdr) ##########
############################################


pre_plot <- function(data, group){
    # read data ---
    snf1 <- read.delim(data, header = TRUE,
                      stringsAsFactors = FALSE)
    # log transform FDR ---
    snf1[[group]] <- -log10(snf1$FDR)
    # choose cell type and logFDR
    snf1.sub <- snf1[,c(1,7)]
    # cell type as rownames ---
    rownames(snf1.sub)<- snf1.sub$GO.ID
    snf1.sub$GO.ID <- NULL
    
    # find missing cell type ---
    cell <- c(paste0("cl_0",1:9), paste0("cl_",10:13))
    miss <- setdiff(cell, rownames(snf1.sub))
    # create df for missing cell type, logFDR as -1 ---
    snf1.miss <- data.frame(C1 = rep(-1, length(miss)))
    colnames(snf1.miss) <- group
    rownames(snf1.miss) <- miss
    # join dfs ---
    snf1.join <- rbind(snf1.sub,snf1.miss)
    
    snf1.join.T <- as.data.frame(t(snf1.join))
    return(snf1.join.T)
    
}


snf1 <- pre_plot("gProfiler_snf1_Brinks_human_allres.txt", group = "C1")
snf2 <- pre_plot("gProfiler_snf2_Brinks_human_allres.txt", group = "C2")
snf3 <- pre_plot("gProfiler_snf3_Brinks_human_allres.txt", group = "C3")

snf.all <- do.call("rbind", list(snf1.join.T,snf2,snf3))
snf.all <- snf.all[,order(names(snf.all))]
names(snf.all) <- c("Cardiac", "Paraxial MD", "Diff somite", "Somite", "Diff front",
                    "Presomitic MD", "NMPs", "Spinal cord", "Mesenchyme", "Endothelium",
                    "allantois", "ExE ectoderm", "endoderm")

# heatmap ---
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
col_fun = colorRamp2(c(0, 1.3, 2), c("#deebf7", "#9ecae1", "#3182bd"))
pdf("Figure19.SNFallfeature_Brinks_neglog10FDR.pdf",width=10,height=2)
Heatmap(as.matrix(snf.all), 
        name = "-log10(FDR)",
        col = col_fun,
        cluster_rows=FALSE,
        cluster_columns =FALSE,
        #right_annotation = c("C1","C2","C3"),
        #split = snf3.sub$SNF_allfeatures_sc3,
        row_names_gp = gpar(fontsize = 8,fontface="bold"),
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 10),
        column_names_side = "top",
        row_names_side = "left"
)
#)
dev.off()
