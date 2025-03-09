library(reshape2)
library(ggplot2)
library(openxlsx)
library(cowplot)

data<-read.xlsx('F5B.xlsx',colNames=F)
d4p<-melt(data[c(1,3:12)],id.vars = c('features'))
importance_plot <- ggplot(d4p, aes(x = features, y = value)) +
  geom_boxplot(fill='#FFC073',alpha=1) +
  theme_minimal() +
  labs( x = "Genus", y = "Importance")+coord_flip() 

foldchange_results <- markers %>%
  group_by(Group) %>%
  summarise(across(colnames(markers)[1:20], mean, .names = "mean_{.col}")) %>%
  pivot_longer(
    cols = -Group,
    names_to = c("Statistic", "Feature"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = c("Group"),
    values_from = "Value"
  ) %>%
  mutate(
    FoldChange = log2(`BC` / `control`)  # 计算 fold change
  )%>%as.data.frame()
data$FC<-foldchange_results[match(data$features,foldchange_results$Feature),'FoldChange']
auc_plot<-ggplot(data, aes(x = features, y = 1, size = AUC, fill = FC)) +
  geom_point(shape = 21, color = "black", stroke = 0.5, alpha = 0.8) +
  scale_size(range = c(3, 10)) + 
  scale_fill_gradient2(low = "#4169E1", mid = "white", high = "#FF4500", midpoint = 0)+
  theme_minimal() +
  labs(x = "Feature",y = "",fill = "Fold Change",size = "AUC") +
  theme(axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    panel.grid.major.y = element_blank()) +
  coord_flip()

plot_grid(importance_plot, auc_plot,
          ncol=2,
          rel_widths = c(0.75, 0.25))
