# For NatComms rev1

# Housekeeping

setwd("~/code/PolySyn_MSNs/")

library(ggplot2)
library(reshape2)

# Load in single feature contrast maps (not necessary but good to store with these other data)

SCAx_morph_nodes_effects_lh <- data.frame(CT=read.table("Data/CNVs/SCAx_CT_DXvNV.txt")[1:152,],
                                           GM=read.table("Data/CNVs/SCAx_CT_DXvNV.txt")[1:152,],
                                           SA=read.table("Data/CNVs/SCAx_CT_DXvNV.txt")[1:152,],
                                           MC=read.table("Data/CNVs/SCAx_CT_DXvNV.txt")[1:152,],
                                           IC=read.table("Data/CNVs/SCAx_CT_DXvNV.txt")[1:152,])
SCAy_morph_nodes_effects_lh <- data.frame(CT=read.table("Data/CNVs/SCAy_CT_DXvNV.txt")[1:152,],
                                           GM=read.table("Data/CNVs/SCAy_GM_DXvNV.txt")[1:152,],
                                           SA=read.table("Data/CNVs/SCAy_SA_DXvNV.txt")[1:152,],
                                           MC=read.table("Data/CNVs/SCAy_MC_DXvNV.txt")[1:152,],
                                           IC=read.table("Data/CNVs/SCAy_IC_DXvNV.txt")[1:152,])
Downs_morph_nodes_effects_lh <- data.frame(CT=read.table("Data/CNVs/Downs_CT_DXvNV.txt")[1:152,],
                                           GM=read.table("Data/CNVs/Downs_GM_DXvNV.txt")[1:152,],
                                           SA=read.table("Data/CNVs/Downs_SA_DXvNV.txt")[1:152,],
                                           MC=read.table("Data/CNVs/Downs_MC_DXvNV.txt")[1:152,],
                                           IC=read.table("Data/CNVs/Downs_IC_DXvNV.txt")[1:152,])
Turner_morph_nodes_effects_lh <- data.frame(CT=read.table("Data/CNVs/Turner_CT_DXvNV.txt")[1:152,],
                                           GM=read.table("Data/CNVs/Turner_GM_DXvNV.txt")[1:152,],
                                           SA=read.table("Data/CNVs/Turner_SA_DXvNV.txt")[1:152,],
                                           MC=read.table("Data/CNVs/Turner_MC_DXvNV.txt")[1:152,],
                                           IC=read.table("Data/CNVs/Turner_IC_DXvNV.txt")[1:152,])
VELOC_morph_nodes_effects_lh <- data.frame(CT=read.table("Data/CNVs/VELOC_CT_DXvNV.txt")[1:152,],
                                           GM=read.table("Data/CNVs/VELOC_GM_DXvNV.txt")[1:152,],
                                           SA=read.table("Data/CNVs/VELOC_SA_DXvNV.txt")[1:152,],
                                           MC=read.table("Data/CNVs/VELOC_MC_DXvNV.txt")[1:152,],
                                           IC=read.table("Data/CNVs/VELOC_IC_DXvNV.txt")[1:152,])
WAGR_morph_nodes_effects_lh <- data.frame(CT=read.table("Data/CNVs/WAGR_CT_DXvNV.txt")[1:152,],
                                           GM=read.table("Data/CNVs/WAGR_GM_DXvNV.txt")[1:152,],
                                           SA=read.table("Data/CNVs/WAGR_SA_DXvNV.txt")[1:152,],
                                           MC=read.table("Data/CNVs/WAGR_MC_DXvNV.txt")[1:152,],
                                           IC=read.table("Data/CNVs/WAGR_IC_DXvNV.txt")[1:152,])

# read in tables if leave-one-feature-out (LOFO) MS change maps

SCAx_LOFO <- read.table("Data/CNVs/SCAx_MSNs_LOFO.txt")
SCAy_LOFO <- read.table("Data/CNVs/SCAy_MSNs_LOFO.txt")
Downs_LOFO <- read.table("Data/CNVs/Downs_MSNs_LOFO.txt")
Turner_LOFO <- read.table("Data/CNVs/Turner_MSNs_LOFO.txt")
VELOC_LOFO <- read.table("Data/CNVs/VELOC_MSNs_LOFO.txt")
WAGR_LOFO <- read.table("Data/CNVs/WAGR_MSNs_LOFO.txt")

# make heatmap of correlations between LOFO maps and empirical MS change map

hm <- data.frame("+X"=cor(SCAx_LOFO)[1,],"-X"=cor(Turner_LOFO)[1,],
                 "+Y"=cor(SCAy_LOFO)[1,],"+21"=cor(Downs_LOFO)[1,],
           "-22q11"=cor(VELOC_LOFO)[1,],"-11p13"=cor(WAGR_LOFO)[1,])
colnames(hm) <- c("+X","-X","+Y","+21","-22q11","-11p13")
hm <- melt(t(hm))
colnames(hm) <- c("Group","Feature","r")
hm_mins <- rep("",36)
hm_mins[hm$r %in% sapply(levels(hm$Group),function(x) min(hm$r[which(hm$Group==x)]))] <- "#"
ggplot(hm,aes(Group,Feature,fill=r))+geom_tile()+geom_text(label=hm_mins,size=7)+theme_minimal()+
  theme(text=element_text(size=24))+scale_fill_gradient2(low="blue",mid="white",high="red")+
  ylab("Feature (left out)")

