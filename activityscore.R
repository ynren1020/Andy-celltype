################2020-06-17###################
## Brinks cell type activity score ##########
#############################################

expr <- read.delim("../rnaseq/chordoma_logrpkm_normalized_108samples.txt", header = TRUE,
                   stringsAsFactors = FALSE, check.names = FALSE)

snf3 <- read.delim("../multiomics/SNF_spectralcluster_allfeature_3clusters.txt", header = TRUE, 
                   stringsAsFactors = FALSE, check.names = FALSE)
snf3.sub <- snf3[grepl("a|[0-9]$",snf3$Chordoma.ID),]
snf3.sub <- snf3.sub[order(snf3.sub$SNF_allfeatures_sc3),]
#write.table(snf3.sub, "../multiomics/SNF_spectralcluster_allfeature_3clusters_91patients_ordered.txt",
#            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

expr.sub <- expr[,names(expr)%in%snf3.sub$Chordoma.ID]

# mean and sd of expression ---
expr.avg <- mean(as.matrix(expr.sub[])) # 1.448046
expr.sd <- sd(as.matrix(expr.sub[])) # 2.168496

# zscore ---
expr.z <- matrix(data=NA, nrow=dim(expr.sub)[1],ncol=dim(expr.sub)[2])
for (i in 1:ncol(expr.z)){
    for (j in 1:nrow(expr.z)){
        expr.z[j,i] <- (expr.sub[j,i] - expr.avg) / expr.sd
    }
    
}
rownames(expr.z) <- rownames(expr.sub)
colnames(expr.z) <- colnames(expr.sub)

write.table(expr.z, "chordoma_logrpkm_zscore_SNF3groups_samples.txt",
                        quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
            
#cell type -------------
brinks <- read.delim("Brink_geneset.gmt", header = FALSE, stringsAsFactors = FALSE)
brinks <- brinks[,-2]

res <- list()
for (i in 1:nrow(brinks)){
    cl1<-unlist(brinks[i,-1])
    cl1.z <- expr.z[rownames(expr.z)%in%cl1,]
    res[[i]] <-data.frame(Chordoma.ID = colnames(cl1.z), AS=colMeans(cl1.z))
}

names(res) <- brinks$V1

res.df <- rlist::list.cbind(res)
discard <- seq(1,26,by=2)
res.df <- res.df[,-discard]

res.df1 <- res.df[rownames(res.df)%in%snf3.sub$Chordoma.ID[snf3.sub$SNF_allfeatures_sc3==1],]
res.df2 <- res.df[rownames(res.df)%in%snf3.sub$Chordoma.ID[snf3.sub$SNF_allfeatures_sc3==2],]
res.df3 <- res.df[rownames(res.df)%in%snf3.sub$Chordoma.ID[snf3.sub$SNF_allfeatures_sc3==3],]

res.df.order <- do.call("rbind", list(res.df1, res.df2, res.df3))
write.table(res.df.order, "SNF_3groups_brinks_activity.score.txt",sep="\t",
            quote = FALSE, col.names = TRUE, row.names = TRUE)

# heatmap ---
library(ComplexHeatmap)
df = data.frame(SNFtype = c(rep("C1", 33), rep("C2", 46), rep("C3", 12)))
# column ----
#ha = HeatmapAnnotation(df = df,col = list(SNFtype = c("C1" =  "red", "C2" = "blue",
#                                                      "C3" = "green")))
# row -----
ha = rowAnnotation(df = df,col = list(SNFtype = c("C1" =  "red", "C2" = "blue",
                                                  "C3" = "green")),
                   show_annotation_name = FALSE)

pdf("Figure18.SNFallfeature_Brinks_activityscore.pdf",width=10,height=10)
Heatmap(as.matrix(res.df.order), 
        name = "Activity Score", 
        cluster_rows=FALSE,
        cluster_columns =FALSE,
        right_annotation = ha,
        split = snf3.sub$SNF_allfeatures_sc3,
        #row_names_gp = gpar(fontsize = 8,fontface="bold"),
        show_row_names = FALSE,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 10),
        column_names_side = "top"
)
#)
dev.off()

