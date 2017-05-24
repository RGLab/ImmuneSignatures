rm(list=ls())
library("ggplot2")
library(plyr)

setwd('C:/Hailong/Projects/HIPC/manuscript/ScienceImunology/figures/')

data=read.table("moduleActivity.txt", header=TRUE, sep="\t", as.is= TRUE)
young.sig = c("platelet activation (III) (M42)",
              "enriched in T cells (II) (M223)",
              "E2F1 targets (Q4) (M10.1)",
              "TBA (M72.0)",
              "inflammatory response (M33)",
              "E2F1 targets (Q3) (M10.0)",
              "TBA (M72.1)",
              "enriched in activated dendritic cells/monocytes (M64)",
              "TBA (M198)",
              "transmembrane transport (II) (M191)",
              "BCR signaling (M54)"
)
old.sig = c("lysosomal/endosomal proteins (M139)",
            "enriched in antigen presentation (II) (M95.0)")
all.genemodule = data$Module

signatures = rep("Other", length(all.genemodule))
signatures[match(young.sig, all.genemodule)]="Young"
signatures[match(old.sig, all.genemodule)]="Older"
signatures = factor(signatures, levels = c("Other","Young", "Older") )

cor(data$Young, data$Older)
# calculate p value of this correlation
#cor.test(data$Young, data$Older)

data.t = data.frame(Young=data$Young, Older=data$Older, signtures=signatures)
data.g = data.t %>%
         arrange(signatures)
       
ggplot()  +
  geom_vline(xintercept=0,size=0.1, colour="black")+
  geom_hline(yintercept=0,size=0.1, colour="black")+
  geom_point(data=filter(data.g, signatures=="Other"),aes(x=Young, y=Older), shape = 19, color ="darkgray", size=1.4)+
  geom_point(data=filter(data.g, signatures=="Young"),aes(x=Young, y=Older), shape = 15, color ="Black", size=2.5)+
  geom_point(data=filter(data.g, signatures=="Older"),aes(x=Young, y=Older), shape = 17, color ="Black", size=2.5)+
  theme_bw()+
  theme(axis.title=element_text(size=16,face="bold"),
        #legend.title = element_blank(),legend.key = element_rect(colour = "white"),
        #legend.position = c(0.25, 0.9), legend.text = element_text(size = 15,face="bold"),legend.background=element_rect(colour = "lightgray")
        legend.position="none"
        ) +
  xlab("Gene module activity (young)")+
  ylab("Gene module activity (older)")+
  scale_x_continuous(limits = c(-0.26, 0.26))   +
  scale_y_continuous(limits = c(-0.24, 0.24)) +
  annotate("text", label = "r = -0.65", x=-0.15, y=-0.2, size = 7, colour = "black",fontface = "bold")
                

      

  
