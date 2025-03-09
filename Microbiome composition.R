library(openxlsx)
library(microeco)
library(ggplot2)
library(magrittr)
library(dplyr)

# Fig.2C
col.hm <- rev(c('#9D9B9C','#550A46','#902044','#CB414B',
                '#D16F7C','#E56B46','#F4A862','#F6DB86',
                '#DFE899','#A4D5A0','#62BA9F','#3681AD','#5A4C97'))

data<-read.xlsx('./F2C.xlsx')
top12<-unique(data$Taxonomy)[1:12]
data$Phylum<-ifelse(data$Taxonomy%in%top12,data$Taxonomy,'others')
data$Phylum<-factor(data$Phylum,levels = c(top12,'others'))
data <- data %>%
  group_by(Sample,Phylum) %>%
  summarise(value = sum(Abundance))
head(data)
data$Study<-sapply(strsplit(as.character(data$Sample),'-'), "[", 1)
data$Study<-factor(data$Study,levels = c('CHI1','CHI2','POL'))
data$Group<-sapply(strsplit(as.character(data$Sample),'-'), "[", 2)
data$Group<-factor(data$Group,levels = c('BC','control'))
ggplot(data, aes(x=Group, y=value, fill=Phylum))  +
  geom_bar(stat = "identity", width=0.8, col='black')  +
  #  theme_pander() +
  ylab('Relative abundance')+
  scale_fill_manual(values = col.hm) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(hjust=1, angle=90),
        panel.background = element_rect(fill=NULL, colour = 'white'))+
  facet_grid(. ~ Study,scales="free_x", space="free_x")


# Fig.2D-E
meta<-read.xlsx('Table S2.xlsx')
meta<-meta[meta$CohortType=='Discovery cohort',]
rownames(meta)<-meta$SampleID
meta$Group<-factor(meta$Group,levels = (c('BC','control')))
meta.BC<-meta[meta$Group=='BC',]
meta.C<-meta[meta$Group=='control',]
data_BC<-feat_abd[rownames(meta.BC),];data_BC<-data_BC[,colSums(data_BC)>0]
data_C<-feat_abd[rownames(meta.C),];data_C<-data_C[,colSums(data_C)>0]
taxon_BC<-taxon[colnames(data_BC),];taxon_BC %<>% tidy_taxonomy
dataset_BC <- microtable$new(otu_table = as.data.frame(t(data_BC)),
                                 sample_table = meta.BC,
                                 tax_table = taxon_BC)
dataset_BC$cal_abund()
taxon_C<-taxon[colnames(data_C),];taxon_C %<>% tidy_taxonomy
dataset_C <- microtable$new(otu_table = as.data.frame(t(data_C)),
                             sample_table = meta.C,
                             tax_table = taxon_C)
dataset_C$cal_abund()
t1 <- trans_abund$new(dataset = dataset_BC, taxrank = "Genus", ntaxa = 10, groupmean = "Study")
t1$plot_tern()

t2 <- trans_abund$new(dataset = dataset_C, taxrank = "Genus", ntaxa = 10, groupmean = "Study")
t2$plot_tern()







