#https://hbctraining.github.io/DGE_workshop_salmon/lessons/09_sleuth.html
#https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/
#https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html

library(wasabi)
library(sleuth)
library(annotables)
library(tidyverse)
library(AnnotationDbi)
library(ggplot2)
library(ggpubr)
library(BSgenome.Celegans.UCSC.ce11)
library(dplyr)
library(ggbio)
library(GENOVA)
library(patchwork)

theme_set(
  theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(),
          axis.title.x=ggtext::element_markdown())
)

outPath="."
genomeVer="WS285"
genomeDir <- paste0("/Users/semple/Documents/MeisterLab/GenomeVer/",genomeVer)
fastqList_file <- "./fastqList_coh1.csv"
fastqList<-read.csv(fastqList_file)

grp="coh1"

source(paste0(outPath,"/functions.R"))
source(paste0(outPath,"/DESeq2_functions.R"))
source(paste0(outPath,"/Sleuth_functions.R"))

makeDirs(outPath,c("sleuth"))
clrs<-getColours()
fileNamePrefix=paste0("sleuth/")

print(paste("fastqList_file is: ",fastqList_file))
print(paste("outPath is: ",outPath))
print(paste("genomeDir is: ",genomeDir))





### genomic plot ------

# get hic data
tev<-load_contacts(signal_path="/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/2021_HiCworked/366_cis1000bins_1000.cool", sample_name= "TEVonly", centromeres=F, balancing = T, verbose = T)

coh1<-load_contacts(signal_path="/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/2021_HiCworked/828_cis1000bins_1000.cool", sample_name= "coh1", centromeres=F, balancing = T, verbose = T)
# get fountains
fountains<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/detected_fountains_equalQ.RDS")
fountains$fountainName<-paste0("fount",1:length(fountains))
# get daugherty enhancers
daugherty<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/publicData/daugherty2017_L3enhancers_ce11.rds")
activedaugherty<-daugherty[daugherty$L3_chromHMMState=="L3_activeEnhancer"]
repdaugherty<-daugherty[daugherty$L3_chromHMMState=="L3_repressedEnhancer" | daugherty$L3_chromHMMState=="L3_H3K27me3Repressed"]
# get jaenes enhancers
jaenes<-import("/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/publicData/Jaenes2018_enhancers_ce11_stages_L3chromHMM.bed")
activejaenes<-jaenes[jaenes$name=="Active enhancer"]
repjaenes<-jaenes[jaenes$name=="Repressed enhancer" | jaenes$name=="H3K27me3 repressed"]

# get transcripts
if(!file.exists(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                       "_ce11.annotations.sqlite"))){
  makeTxDbsqlite_ce11(genomeDir,genomeVer)
}
txdb<-loadDb(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                    "_ce11.annotations.sqlite"))

dir.create(paste0(outPath,"/sleuth/sleuthDTE"))

theme_bed = function(){
  theme(title=element_text(size=8),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank())
}

#' Try bed track plot but if no overlaps in this region - do empty plot keeping title
tryPlot<-function(gr1,gr2,title,fill){
  p1<-autoplot(subsetByOverlaps(gr1,gr2,ignore.strand=T), geom="rect",fill=fill,legend=F) +
    xlim(gr2) + ggtitle(title)
  p2<-autoplot(Celegans, which = gr2) + xlim(gr2)+ ggtitle(title)
  if(is(try(print(p1)),"try-error")) p2 else p1
}

#wbid<-so$wormbaseID[so$publicID=="cat-1"][1]

for (wbid in unique(so$wormbaseID[so$diffuptx==T])) {
  #wbid="WBGene00020131"
  gr<-GenomicRanges::reduce(so[so$wormbaseID==wbid],drop.empty.ranges=T)
  #f<-nearest(gr,fountains,ignore.strand=T)
  #gr<-reduce(fountains[f],gr)
  gr<-trim(resize(gr,width=100000,fix="center"))
  strand(gr)<-"*"
  p.bg <- autoplot(Celegans, which = gr) + xlim(gr)
  p.txdb<-autoplot(txdb,aes(fill=strand),which=gr,label.size=1.5) + xlim(gr) +
    scale_fill_manual(values=c("#A58AFF","#cccccccc","#00C094"))
  p.fount <- tryPlot(fountains,gr,title="Detected fountains",fill="darkblue") + theme_bed()
  p.activeD<-tryPlot(activedaugherty,gr,title="Active enhancers (Daugherty et al.)",fill="darkgreen") + theme_bed()
  p.activeJ<-tryPlot(activejaenes,gr,fill="darkgreen",title="Active enhancers (Jaenes et al.)") + theme_bed()
  p.repD<-tryPlot(repdaugherty,gr,title="Repressed enhancers (Daugherty et al.)",fill="darkgray")  + theme_bed()
  p.repJ<-tryPlot(repjaenes,gr,title="Repressed enhancers (Jaenes et al.)",fill="darkgray")  + theme_bed()

  p.hic<-pyramid(exp=tev,chrom=as.character(seqnames(gr)),start=start(gr),end=end(gr),crop_y=100000) +
    ggtitle("TEVonly") + theme(legend.key.size = unit(0.2,'cm'))
  p.hic1<-pyramid(exp=coh1,chrom=as.character(seqnames(gr)),start=start(gr),end=end(gr),crop_y=100000) +
    ggtitle("COH-1cs") + theme(legend.key.size = unit(0.2,'cm'))
  p.hicdiff<-pyramid_difference(coh1,tev,chrom=as.character(seqnames(gr)),start=start(gr),end=end(gr),crop_y=100000) +
    ggtitle("COH-1cs - TEV") + theme(legend.position=c(0.9,0.8))
  # look at transcripts
  pptxdb<-autoplot(subsetByOverlaps(so,gr,ignore.strand=T),aes(fill=b),names.expr="target_id") +xlim(gr)+
    scale_fill_gradient2(na.value="grey90") + theme(legend.key.size = unit(0.2,'cm'))+
    labs(fill="LFC")

  p1<-(p.hic+p.hic1)/p.hicdiff/p.txdb@ggplot/pptxdb@ggplot/p.fount@ggplot/plot_spacer()/p.activeD@ggplot/plot_spacer()/p.activeJ@ggplot/plot_spacer()/p.repD@ggplot/plot_spacer()/p.repJ@ggplot/plot_spacer()/p.bg@ggplot +
    plot_layout(heights=c(10,25,10,5,0.5,-2,0.5,-2,0.5,-2,0.5,-2,0.5,-2,0.1)) +
    plot_annotation(title=paste0(wbid,": ",so$publicID[so$wormbaseID==wbid]))
  ggbio::ggsave(paste0(outPath,"/sleuth/sleuthDTE/hic_",wbid,".pdf"),p1,device="pdf",height=29,width=19,units="cm")
}

plotList=list()
for (wbid in unique(so$wormbaseID[so$diffuptx==T])) {
  tmp<-data.frame(so[so$wormbaseID==wbid])
  tmp$txt<-paste0("p=",round(tmp$qval,2))
  tmp$txt[tmp$qval>0.05]<-"n.s."

  plotList[[wbid]]<-ggplot(tmp,aes(x=target_id,y=b,fill=b)) +
    geom_bar(stat="identity",color="black") +
    ggtitle(paste0(tmp$wormbaseID[1],": ",tmp$publicID)) + ylab(label="Log<sub>2</sub>FC") +
    geom_errorbar(aes(ymin=b-se_b,ymax=b+se_b),width=0.2) +
    geom_hline(yintercept=0) + scale_fill_gradient2() +
    geom_text(aes(y=0,label=txt)) +
    theme(legend.key.size = unit(0.8,'cm'),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(fill="LFC")
}

library(gridExtra)
p<-gridExtra::marrangeGrob(grobs=plotList,ncol=2,nrow=2)
p
pdf(paste0(outPath,"/sleuth/sleuthDTE/sleuthDTE_barplot.pdf"),paper="a4",height=10,width=10)
print(p)
dev.off()
