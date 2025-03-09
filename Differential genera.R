library(ggplot2)
library(openxlsx)
library(ggrepel)
# Fig.3C
data<-read.xlsx('Table S3.xlsx',)
diff.plot <- data[,c('feature','coef','pval','qval.fdr')]
diff.plot$label = ifelse(diff.plot$pval<0.05&diff.plot$qval.fdr<0.1, ifelse(diff.plot$coef>0,paste0("Upregulated in ", 'BC'),paste0("Downregulated in ", 'BC')),"Not significant")
diff.plot$significant <- with(diff.plot, pval<0.05&qval.fdr<0.1)

ggplot(diff.plot, aes(x = coef, y = -log10(pval), colour=label)) +
  geom_point(alpha=0.8, size=3.5) +
  scale_color_manual(values=c("#0072B2", "#C0C6CC","#D55E00"))+
  geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="Coef",
       y="-log10 (p-value)")+
  geom_text_repel(
    data = subset(diff.plot, significant == TRUE),   # 只对显著的点添加标签
    aes(label = feature), 
    size = 5, 
    box.padding = 0.3,
    max.overlaps = 10)+
  theme_bw()+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),  
        legend.position="none", 
        legend.title = element_blank())

