library(openxlsx)
library(ggplot2)
library(ggpubr)
library(vegan)
library(Rmisc)

#Function
GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)
geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}
#Data
data<-read.xlsx('F2A.xlsx')
all<-data;all$Study<-'ALL'
data<-rbind(data,all)
data$Study<-factor(data$Study,levels = c('CHI1','CHI2','POL','ALL'))

#Plot
Data_summary <- summarySE(data, measurevar='chao1', groupvars=c("Study","Group"))
ggplot(data=data, aes(x=Study, y=chao1,fill=Group)) + 
  geom_split_violin(trim=FALSE,color="white") + 
  geom_point(data = Data_summary,aes(x=Study, y=chao1),pch=19,position=position_dodge(0.1),size=0.5)+ 
  geom_errorbar(data = Data_summary,aes(ymin = chao1-se, ymax=chao1+se), 
                width=0.12,  position=position_dodge(0.1), 
                color="black",alpha = 0.7,size=0.5) +
  scale_fill_manual(values= c('control' = '#3C84B8', 'BC' = '#FFC073')) +
  theme_bw()+
  ylim(c(0,150))+
  theme(axis.text.x=element_text(angle=0,hjust = 0.5,colour="black",size=15), 
        axis.text.y=element_text(size=16,face="plain",colour="black"), 
        axis.line.x.top = element_line(colour = "black",size=1),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text=element_text(face="plain",  colour="black", size=16),
        legend.title=element_text(face="plain", colour="black",size=18),
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5, vjust = 1),
        panel.grid.minor = element_blank())+ 
  ylab("")+xlab("")+ggtitle('Chao1')

Data_summary <- summarySE(data, measurevar='shannon', groupvars=c("Study","Group"))
ggplot(data=data, aes(x=Study, y=shannon,fill=Group)) + 
  geom_split_violin(trim=FALSE,color="white") + 
  geom_point(data = Data_summary,aes(x=Study, y=shannon),pch=19,position=position_dodge(0.1),size=0.5)+ 
  geom_errorbar(data = Data_summary,aes(ymin = shannon-ci, ymax=shannon+ci), 
                width=0.12,  position=position_dodge(0.1), 
                color="black",alpha = 0.7,size=0.5) +
  scale_fill_manual(values= c('control' = '#3C84B8', 'BC' = '#FFC073')) +
  theme_bw()+
  theme(axis.text.x=element_text(angle=0,hjust = 0.5,colour="black",size=15), 
        axis.text.y=element_text(size=16,face="plain",colour="black"), 
        axis.line.x.top = element_line(colour = "black",size=1),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text=element_text(face="plain",  colour="black", size=16),
        legend.title=element_text(face="plain", colour="black",size=18),
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5, vjust = 1),
        panel.grid.minor = element_blank())+ 
  ylab("")+xlab("")+ggtitle('Shannon')




