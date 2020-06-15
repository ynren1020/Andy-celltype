################ 2020-06-15 #################
# GSEA results for Pijuan and Brinks ########
# NES vs rank scatter plot for SNF 3 groups##
#############################################

pre_plot <- function(pos, neg, n = 13){
    pos <- read.delim(pos, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    if (nrow(pos) == 0) {
        neg <- read.delim(neg, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        neg <- neg[, c("NAME","NES","FDR.q.val")]
        neg$RANK <- order(neg$NES)
        neg$logP <- -log10(neg$FDR.q.val + 0.00001)
        return(neg)
    } else if (nrow(pos) == n){
        pos <- pos[, c("NAME","NES","FDR.q.val")]
        pos$RANK <- order(pos$NES)
        pos$logP <- -log10(pos$FDR.q.val + 0.00001)
        return(pos)
    } else {
        neg <- read.delim(neg, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        pos.neg <- rbind(neg, pos)
        pos.neg <- pos.neg[, c("NAME","NES","FDR.q.val")]
        pos.neg$RANK <- order(pos.neg$NES)
        pos.neg$logP <- -log10(pos.neg$FDR.q.val + 0.00001)
        return(pos.neg)
    } 
    
}


pijuan1 <- pre_plot("Pijuan_1vsothers_pos.txt", "Pijuan_1vsothers_neg.txt")
pijuan2 <- pre_plot("Pijuan_2vsothers_pos.txt", "Pijuan_2vsothers_neg.txt")
pijuan3 <- pre_plot("Pijuan_3vsothers_pos.txt", "Pijuan_3vsothers_neg.txt")

brinks1 <- pre_plot("Brinks_1vsothers_pos.txt", "Brinks_1vsothers_pos.txt")
brinks2 <- pre_plot("Brinks_2vsothers_pos.txt", "Brinks_2vsothers_neg.txt")
brinks3 <- pre_plot("Brinks_3vsothers_pos.txt", "Brinks_3vsothers_neg.txt")

# plot ----

library(ggplot2)
library(RColorBrewer)
library(ggrepel)


ptest<-ggplot2::ggplot(brinks1,aes(x=RANK,y=NES,label = ifelse(logP >4,as.character(NAME),'')))+
    ggplot2::geom_point(aes(color = logP))+
    labs(color='-log10(q)') +
    #geom_text(aes(label = ifelse(color=="dark",as.character(NAME),'')), vjust= -1)+
    #ggplot2::geom_point(aes(color=ifelse(FDR.q.val < 0.001, "FDR < = 0.001", ifelse(FDR.q.val < 0.01, "FDR < 0.01", ifelse(FDR.q.val ))))+ #mapping aes(color),here define different color based on isSuper!!!
    #ggplot2::scale_colour_manual(name='',values=c('FDR > 0.05'='grey','FDR < = 0.05'='red'))+ #add color manually!!!!
    scale_color_gradient(low = "blue", high = "red") +
    #geom_vline(xintercept = 25982,linetype="dashed",size=0.5,alpha=0.2)+
    #geom_hline(yintercept = 1.3,linetype="dashed",size=0.5,alpha=0.2)+
    #geom_hline(yintercept = 2,linetype="dashed",size=0.5,alpha=0.2)+
    ggplot2::labs(x="Rank",y="Normal Enrichement Score")+
    ggplot2::theme_bw()+
    ggplot2::theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())+
    ggplot2::theme(legend.position = c(0.9, 0.2))
   
set.seed(00)
ppd<-ptest+geom_text_repel(data=brinks1,nudge_x=-1,direction="y",force=2,max.iter=4000,hjust=1,segment.size=0.2)

ggsave("Brinks1_gsea.pdf", width = 6, height = 6)



#pijuan plot ---

ptest<-ggplot2::ggplot(pijuan3,aes(x=RANK,y=NES,label = ifelse((logP >4 | NAME == "27-NMP"),as.character(NAME),'')))+
    ggplot2::geom_point(aes(color = logP))+
    labs(color='-log10(q)') +
    #geom_text(aes(label = ifelse(color=="dark",as.character(NAME),'')), vjust= -1)+
    #ggplot2::geom_point(aes(color=ifelse(FDR.q.val < 0.001, "FDR < = 0.001", ifelse(FDR.q.val < 0.01, "FDR < 0.01", ifelse(FDR.q.val ))))+ #mapping aes(color),here define different color based on isSuper!!!
    #ggplot2::scale_colour_manual(name='',values=c('FDR > 0.05'='grey','FDR < = 0.05'='red'))+ #add color manually!!!!
    scale_color_gradient(low = "blue", high = "red") +
    #geom_vline(xintercept = 25982,linetype="dashed",size=0.5,alpha=0.2)+
    #geom_hline(yintercept = 1.3,linetype="dashed",size=0.5,alpha=0.2)+
    #geom_hline(yintercept = 2,linetype="dashed",size=0.5,alpha=0.2)+
    ggplot2::labs(x="Rank",y="Normal Enrichement Score")+
    ggplot2::theme_bw()+
    ggplot2::theme(panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank())+
    ggplot2::theme(legend.position = c(0.9, 0.2))

set.seed(00)
ppd<-ptest+geom_text_repel(data=pijuan3,nudge_x=-10,direction="y",force=2,max.iter=4000,hjust=1,segment.size=0.2)

ggsave("Pijuan3_gsea.pdf", width = 6, height = 6)




