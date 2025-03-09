library(ggplot2)
library(dplyr)
library(ggprism)
library(tidyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggprism)

# Fig. 4A
data<-read.xlsx('Table S4.xlsx',sheet = 1)
meta<-read.xlsx('Table S2.xlsx')
meta<-meta[meta$CohortType=='Discovery cohort',]
rownames(meta)<-meta$SampleID
meta$Group<-factor(meta$Group,levels = (c('BC','control')))

daa_results_df=data
Group=meta$Group
ko_to_kegg = FALSE
p_value_bar = TRUE
x_lab = "KOgene"
select = c(data[data$coef>0&data$PathwayL1=='Metabolism'&data$PathwayL2=='Nucleotide metabolism','feature'],
           data[data$coef<0&data$PathwayL1=='Metabolism'&data$PathwayL2=='Metabolism of cofactors and vitamins','feature'][2:13])
daa_results_df <- daa_results_df[!is.na(daa_results_df[,x_lab]), ]
colors <- c("#d93c3e", "#3685bc", "#6faa3e", "#e8a825", 
            "#c973e6", "#ee6b3d", "#2db0a7", "#f25292")[1:nlevels(as.factor(Group))]
rownames(daa_results_df)<-daa_results_df$feature
daa_results_filtered_df <- daa_results_df[select, ]
daa_results_filtered_sub_df <- daa_results_filtered_df[daa_results_filtered_df$feature, ]
relative_abundance_mat <- as.matrix(feat_abd)
sub_relative_abundance_mat <- relative_abundance_mat[rownames(relative_abundance_mat) %in% 
                                                       daa_results_filtered_sub_df$feature, ]
error_bar_matrix <- cbind(sample = colnames(sub_relative_abundance_mat), 
                          group = Group, t(sub_relative_abundance_mat))
error_bar_df <- as.data.frame(error_bar_matrix)
error_bar_df$group <- factor(Group, levels = levels(as.factor(Group)))
error_bar_pivot_longer_df <- tidyr::pivot_longer(error_bar_df,
                                                 -c(sample, group))
error_bar_pivot_longer_tibble <- mutate(error_bar_pivot_longer_df, 
                                        group = as.factor(group))
error_bar_pivot_longer_tibble$sample <- factor(error_bar_pivot_longer_tibble$sample)
error_bar_pivot_longer_tibble$name <- factor(error_bar_pivot_longer_tibble$name)
error_bar_pivot_longer_tibble$value <- as.numeric(error_bar_pivot_longer_tibble$value)
error_bar_pivot_longer_tibble_summarised <- error_bar_pivot_longer_tibble %>% 
  group_by(name, group) %>% summarise(mean = mean(value), 
                                      sd = stats::sd(value))
error_bar_pivot_longer_tibble_summarised <- error_bar_pivot_longer_tibble_summarised %>% 
  mutate(group2 = "nonsense")

daa_results_filtered_sub_df$pro<-ifelse(daa_results_filtered_sub_df$coef>0,'BC','control')
daa_results_filtered_sub_df <- daa_results_filtered_sub_df[order(daa_results_filtered_sub_df$pval), ]
daa_results_filtered_sub_df <- daa_results_filtered_sub_df[order(abs(daa_results_filtered_sub_df$coef),decreasing=T), ]
daa_results_filtered_sub_df <- daa_results_filtered_sub_df[order(daa_results_filtered_sub_df$pro), ]
error_bar_pivot_longer_tibble_summarised_ordered <- data.frame(name = NULL, 
                                                               group = NULL, mean = NULL, sd = NULL)
for (i in daa_results_filtered_sub_df$feature) {
  error_bar_pivot_longer_tibble_summarised_ordered <- rbind(error_bar_pivot_longer_tibble_summarised_ordered, 
                                                            error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name == 
                                                                                                       i, ])
}
error_bar_pivot_longer_tibble_summarised_ordered[, x_lab] <- rep(daa_results_filtered_sub_df[,x_lab], each = length(levels(factor(error_bar_pivot_longer_tibble_summarised_ordered$group))))


error_bar_pivot_longer_tibble_summarised_ordered$name <- factor(error_bar_pivot_longer_tibble_summarised_ordered$name, 
                                                                levels = rev(daa_results_filtered_sub_df$feature))
bar_errorbar <- ggplot(error_bar_pivot_longer_tibble_summarised_ordered, 
                       aes(mean, name, fill = group)) + 
  geom_errorbar(aes(xmax = mean +sd, xmin = 0), 
                position = position_dodge(width = 0.8),width = 0.5, size = 0.5, color = "black") + 
  geom_bar(stat = "identity",position = position_dodge(width = 0.8), width = 0.8) + 
  GGally::geom_stripped_cols(width = 10) + scale_fill_manual(values = colors) + 
  scale_color_manual(values = colors) + ggprism::theme_prism() + 
  scale_x_continuous(expand = c(0, 0), guide = "prism_offset_minor", 
  ) + scale_y_discrete(labels = rev(daa_results_filtered_sub_df[,x_lab])) + labs(x = "Relative Abundance", y = NULL) + 
  theme(axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), axis.line.x = element_line(size = 0.5), 
        axis.ticks.x = element_line(size = 0.5), 
        panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(), 
        axis.text = element_text(size = 10, color = "black"), 
        axis.text.x = element_text(margin = margin(r = 0)), 
        axis.text.y = element_text(size = 10, color = "black", 
                                   margin = margin(b = 6)), 
        axis.title.x = element_text(size = 10,color = "black", hjust = 0.5), legend.position = "top", 
        legend.key.size = unit(0.1, "cm"), legend.direction = "vertical", 
        legend.justification = "left", legend.text = element_text(size = 8,face = "bold"), 
        legend.box.just = "right", plot.margin = margin(0,0.5, 0.5, 0, unit = "cm")) + coord_cartesian(clip = "off")

daa_results_filtered_sub_df <- cbind(daa_results_filtered_sub_df, 
                                     negative_log10_p = -log10(daa_results_filtered_sub_df$pval), 
                                     group_nonsense = "nonsense", log_2_fold_change = NA)
for (i in daa_results_filtered_sub_df$feature) {
  mean <- error_bar_pivot_longer_tibble_summarised_ordered[error_bar_pivot_longer_tibble_summarised_ordered$name %in% 
                                                             i, ]$mean
  daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature == 
                                i, ]$log_2_fold_change <- log2(mean[1]/mean[2])
}

daa_results_filtered_sub_df$feature <- factor(daa_results_filtered_sub_df$feature, 
                                              levels = rev(daa_results_filtered_sub_df$feature))
p_values_bar <- daa_results_filtered_sub_df %>% ggplot(aes(feature, coef, fill = group_nonsense)) + 
  geom_bar(stat = "identity",position = position_dodge(width = 0.8), width = 0.8) + 
  labs(y = "coef", x = NULL) + GGally::geom_stripped_cols() + 
  scale_fill_manual(values = "#87ceeb") + scale_color_manual(values = "#87ceeb") + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", 
             color = "black") + ggprism::theme_prism() + 
  scale_y_continuous(expand = c(0,0), guide = "prism_offset_minor") + 
  theme(axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.5),
        axis.ticks.x = element_line(size = 0.5), 
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(), 
        axis.text = element_text(size = 10,color = "black"), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black",margin = margin(b = 6)), 
        axis.title.x = element_text(size = 11, color = "black", hjust = 0.5), legend.position = "non") + 
  coord_flip()

daa_results_filtered_sub_df$pval <- as.character(daa_results_filtered_sub_df$pval)
daa_results_filtered_sub_df$unique <- nrow(daa_results_filtered_sub_df) - 
  seq_len(nrow(daa_results_filtered_sub_df)) + 1
daa_results_filtered_sub_df$ann<-NA
for (i in 1:nrow(daa_results_filtered_sub_df)){
  if (daa_results_filtered_sub_df$pval[i]<0.001){
    daa_results_filtered_sub_df[i,'ann']<-'***'
  }else if(daa_results_filtered_sub_df$pval[i]<0.01){
    daa_results_filtered_sub_df[i,'ann']<-'**'
  }else if(daa_results_filtered_sub_df$pval[i]<0.05){
    daa_results_filtered_sub_df[i,'ann']<-'*'
  }
}
p_annotation <- daa_results_filtered_sub_df %>% ggplot(aes(group_nonsense, pval)) +
  geom_text(aes(group_nonsense, unique, label = ann), 
            size = 3.5, color = "black",fontface = "bold", family = "sans") +
  labs(y = "p-value") + 
  scale_y_discrete(position = "right") + ggprism::theme_prism() + 
  theme(axis.ticks = element_blank(), 
        axis.line = element_blank(), panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), panel.background = element_blank(), 
        axis.text = element_blank(), plot.margin = unit(c(0,0.2, 0, 0), "cm"),
        axis.title.y = element_text(size = 11,color = "black", vjust = 0), 
        axis.title.x = element_blank(), 
        legend.position = "non")

combination_bar_plot <- bar_errorbar + p_values_bar + 
  p_annotation + patchwork::plot_layout(ncol = 3, 
                                        widths = c(2.3, 0.7, 0.3))

# Fig. 4B
format_pvalue <- function(p) {
  if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else if (p < 0.001) {
    return("***")
  }else {
    return("")
  }
}
diff_sig<-read.xlsx( "Table S4.xlsx",sheet = 2,rowNames = T)
rdata<-read.xlsx('./F4E.xlsx',sheet = 1)
pdata<-read.xlsx('./F4E.xlsx',sheet = 2)
rdata<-rdata[,select_pwy]
pdata<-pdata[,select_pwy]
sig_matrix <- (apply(pdata, c(1, 2), format_pvalue))
sig_matrix<-as.data.frame(t(sig_matrix))
annotable<-diff_sig[rownames(rdata),]
annotable$diff =ifelse(annotable$coef>0,'up','down')
annotable2<-read.xlsx('Table S4.xlsx',sheet = 2)
rownames(annotable2)<-annotable2$description
annotable2<-annotable2[colnames(rdata),]
column_ha = HeatmapAnnotation(sig_diff=annotable$diff,
                              col = list(sig_diff=c(up = "#F75000", down = "#0280FB",ns='grey')))
summary(annotable2$coef)
row_ha=rowAnnotation(coef=annotable2$coef,class=annotable2$class,
                     col = list(coef=colorRamp2(c(-0.013, 0, 0.013), c("#0B8B42", "white", "#EB8B2D"))))
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
Heatmap(t(rdata),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sig_matrix[i, j], x, y, gp = gpar(fontsize = 10))},
            bottom_annotation = column_ha,
            right_annotation = row_ha,
            width = ncol(t(rdata))*unit(15, "mm"), height = nrow(t(rdata))*unit(7.5, "mm"),
            column_order = c('Finegoldia','Bacteroides','Anaerococcus','Muribaculaceae','Fenollaria','Escherichia-Shigella'),
            col=col_fun,cluster_columns = F,cluster_rows = F
)
