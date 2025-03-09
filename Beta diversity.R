library(openxlsx)
library(vegan)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)


#Data
pcoa<-read.xlsx('F2B.xlsx')
pc<-c(16.53,8.46 )

#Fig.1B
# main plot
g.main <- pcoa %>% 
  ggplot(aes(x=PCoA1, y=PCoA2, shape= Group ,col= Study)) +
  stat_ellipse(aes(linetype = Group), level = 0.95) +
  geom_point(size = 1.5 )  +
  xlab(paste0("PCoA1 [",pc[1],"%]")) + ylab(paste0("PCoA2 [",pc[2],"%]")) +
  scale_colour_manual(values=c('CHI1' = '#FF5959', 'CHI2' = '#FACF5A',
                               'POL' = '#49BEB7','CRO'='#085F63'),guide = FALSE) +
  scale_shape_manual(values=c(16,  17)) +
  scale_x_continuous(position='top') +
  theme(panel.background = element_rect(fill='white', color = 'black'),
        legend.position = 'none',
        axis.ticks=element_blank(), axis.text = element_blank(),axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        panel.grid = element_blank())
# study boxplot axis 1
g.s.1 <- pcoa %>% 
  mutate(Study=factor(Study, levels=names(c('CHI1' = '#FF5959', 'CHI2' = '#FACF5A',
                                            'POL' = '#49BEB7')))) %>% 
  ggplot(aes(y=PCoA1, x=Study, fill=Study)) +
  xlab(paste0('Study','')) +
  geom_boxplot() +
  scale_fill_manual(values=c('CHI1' = '#FF5959', 'CHI2' = '#FACF5A',
                             'POL' = '#49BEB7'), guide=FALSE) +
  theme(axis.ticks = element_blank(),
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y  = element_text(size = 10),
        panel.grid = element_blank()) + 
  coord_flip()
# study boxplot axis 2
g.s.2 <- pcoa %>% 
  mutate(Study=factor(Study, levels=names(c('CHI1' = '#FF5959', 'CHI2' = '#FACF5A',
                                            'POL' = '#49BEB7')))) %>% 
  ggplot(aes(y=PCoA1, x=Study, fill=Study)) + 
  xlab(paste0('Study','')) +
  geom_boxplot() + 
  scale_fill_manual(values=c('CHI1' = '#FF5959', 'CHI2' = '#FACF5A',
                             'POL' = '#49BEB7'), guide = FALSE) +
  
  scale_x_discrete(position='top') +
  theme(axis.ticks=element_blank(), 
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x  = element_text(size = 10),
        panel.grid = element_blank())
# group plot axis1
g.g.1 <- pcoa %>% 
  ggplot(aes(x=Group, y=PCoA1, fill=Group)) +
  xlab(paste0('Group','')) +
  geom_boxplot() +
  scale_fill_manual(values= c('control' = '#3C84B8', 'BC' = '#FFC073'),guide= FALSE) + 
  ylab(paste0("PCoA1 [",pc[1],"%]")) +
  theme(axis.ticks.y=element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y=element_blank(),
        axis.title.y = element_text(size = 10),
        axis.title.x  = element_blank(),
        legend.title=element_text(size =15),legend.text=element_text(size = 15),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank()) + 
  coord_flip()
# group plot axis2
g.g.2 <- pcoa %>% 
  ggplot(aes(x=Group, y=PCoA2, fill=Group)) +
  xlab(paste0('Group','')) +
  geom_boxplot() +
  scale_fill_manual(values=c('control' = '#3C84B8', 'BC' = '#FFC073'), guide=FALSE) + 
  scale_x_discrete(position='top') + 
  scale_y_continuous(position = 'right') +
  
  ylab(paste0("PCoA2 [",pc[2],"%]")) + 
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank())

plot_grid(g.main, g.s.2, g.g.2, g.s.1, NULL,NULL,g.g.1, NULL, NULL,
          nrow=3,
          rel_widths = c(0.8, 0.2, 0.2), rel_heights = c(0.8, 0.2, 0.2))


# Confound analysis
phenotype_microbiota_corr<- function(dist,microbiota){
  df_result <- matrix(nrow = 2, ncol = ncol(meta),0)
  row.names(df_result) <- c('R2','Pvalue')
  colnames(df_result) <- colnames(meta)
  for (i in 1:ncol(meta)){
    tmp = meta[,i]
    if (distance == "bc"){
      permanova = adonis2(dist~tmp, permutations=count, method = "bray")
    } else if (distance == "euclidean"){
      permanova = adonis2(dist~tmp, permutations=count, method = "euclidean")
    }
    df_result['R2',colnames(meta)[i]] <- as.data.frame(permanova)['Model','R2']
    df_result['Pvalue',colnames(meta)[i]] <- as.data.frame(permanova)['Model','Pr(>F)']
  }
  return(df_result)    
}
meta<-read.xlsx('Table S1.xlsx')
meta<-meta[meta$CohortType=='Discovery cohort',]
meta<-meta[,c('Group','Study')]
group <-'Group'
count <- 999
distance<-'bc'
dist = vegdist(feat_abd, method = 'bray', correction = 'lingoes', na.rm = TRUE)
result <- phenotype_microbiota_corr(dist,feat_abd)
print(result)







