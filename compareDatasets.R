library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(eulerr)

outPath="."
scriptPath="."
source(paste0(scriptPath,"/functions.R"))

theme_set(
  theme_bw()+
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
)

scriptName="compareDatasets"
fileNamePrefix=paste0(scriptName,"/")
makeDirs(outPath,dirNameList=c(paste0(c("plots/"),scriptName)))

padjVal=0.05
lfcVal=0.5
# DEseq2 results files to use
RNAseqDir="/Users/semple/Documents/MeisterLab/otherPeopleProjects/Bolaji/BolajiRNAseq_20211216"

fileList<-data.frame(
  sampleName=c("COH1cs"),
  filePath=paste0(RNAseqDir,
                  c("/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")))








###################-
## Compare coh-1cs analysis
###################-

mdas<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/2021_RNAseq_MDas/rds/p0.05_lfc0.5_filtChrAX/filtChrAX_X.wt.wt.0mM_coh1cs_vs_wt_DESeq2_fullResults_p0.05.rds")

salmon<-readRDS(file=fileList[fileList$sampleName=="COH1cs","filePath"])

sigGenesUp<-list()
sigGenesUp[["mdas"]]<-mdas$wormbaseID[!is.na(mdas$padj) & mdas$padj<0.05 & mdas$log2FoldChange>0]
sigGenesUp[["fount"]]<-salmon$wormbaseID[!is.na(salmon$padj) & salmon$padj<0.05 & salmon$log2FoldChange>0]


fit<-euler(sigGenesUp)
p1<-plot(fit, quantities=list(type="counts"),
         main=list(label=paste0("Up regulated genes: ", "padj<",padjVal,"\n",
                                paste(lapply(row.names(fit$ellipses), function(x){
                                  paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                                }), collapse="  ")), fontsize=8))
p1


sigGenesDown<-list()
sigGenesDown[["mdas"]]<-mdas$wormbaseID[!is.na(mdas$padj) & mdas$padj<0.05 & mdas$log2FoldChange<0]
sigGenesDown[["fount"]]<-salmon$wormbaseID[!is.na(salmon$padj) & salmon$padj<0.05 & salmon$log2FoldChange<0]

fit<-euler(sigGenesDown)
p2<-plot(fit, quantities=list(type="counts"),
         main=list(label=paste0("Down regulated genes: ", "padj<",padjVal,"\n",
                                paste(lapply(row.names(fit$ellipses), function(x){
                                  paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                                }), collapse="  ")), fontsize=8))
p2

# p<-ggarrange(p1,p2,nrow=2)
# ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_coh1_mdasVfount.pdf"),plot=p, device="pdf",width=19,height=29,units="cm")

c1<-inner_join(data.frame(mdas),data.frame(salmon),by="wormbaseID")
dim(c1)
minScale=-3
maxScale=3
p3<-ggplot(c1,aes(x=log2FoldChange.x,y=log2FoldChange.y)) +
  geom_point(size=1,alpha=0.4) +
  coord_cartesian(xlim=c(minScale,maxScale), ylim=c(minScale,maxScale)) +
  geom_smooth(method=lm,se=F,fullrange=T, size=0.7,show.legend = T) +
  geom_hline(yintercept=0,lty=3,col="grey70",) +
  geom_vline(xintercept=0,lty=3,col="grey70") +
  ggpubr::stat_cor(aes(label = ..r.label..), method="spearman",
                   cor.coef.name = c("Rs"), output.type = "text",
                   show.legend=F,size=3,
                   label.x=minScale+0.1, label.y=c(maxScale-0.1,maxScale-0.5)) +
  ggpubr::stat_cor(aes(label = ..r.label..), method="pearson",
                   cor.coef.name = c("Rp"), output.type = "text",
                   show.legend=F,size=3,
                   label.x=minScale+0.1, label.y=c(maxScale-0.6,maxScale-1.1),colour="red") +
  xlab(label="mdas LFC") + ylab("fount LFC")

p<-ggarrange(p1,p2,p3,nrow=2,ncol=2)
p
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_coh1_mdasVfount.pdf"),plot=p, device="pdf",width=19,height=19,units="cm")
