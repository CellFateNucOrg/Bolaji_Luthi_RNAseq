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
library(eulerr)
library(ggvenn)

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


################-
## sleuth analysis ------
################-
sf_dirs<-file.path("salmon/mRNA",fastqList$sampleName)


## convert file from salmon to kallisto format
prepare_fish_for_sleuth(sf_dirs)

# prepare sampleTable and design formula
sampleTable<-data.frame(sample=fastqList$sampleName,
                        condition=factor(fastqList$strain,levels=c("366","828")),
                        path=sf_dirs,
                        lane=fastqList$lane,
                        replicate=fastqList$repeatNum)

design<- formula("~ lane + replicate + condition")
reduced<- formula("~ lane  + replicate")

# prepare metadata
tx2gene<-getTx2Gene(genomeDir,genomeVer)
tx2gene$GENEID<-gsub("Gene:","",tx2gene$GENEID)
colnames(tx2gene)<-c("target_id","gene_id")

txMeta<-getTxMetadataGR(genomeDir,genomeVer)
#txMeta<-txMeta[txMeta$class=="protein_coding_gene"]
txMeta<-tagOscillating(txMeta)


############-
# Sleuth
############-

# read in data to sleuth
so<-sleuth_prep(sampleTable,normalize=F)

# filter input genes by counts and using metadata
counts<-sleuth_to_matrix(so,"obs_raw","est_counts")
filter<-apply(counts,1,basic_filter)

bioFilter<-txMeta$txptSeqID[txMeta$class=="protein_coding_gene" &
                              txMeta$oscillating=="no"]
txFilter<-bioFilter[which(bioFilter %in% names(filter)[filter])]

so<-sleuth_prep(sampleTable,
                full_model=design,
                target_mapping=tx2gene,
                num_cores=1,
                read_bootstrap_tpm=T,
                extra_bootstrap_summary=T,
                transformation_function=function(x) log2(x + 0.5),
                filter_target_id=txFilter)

# 18363 targets passed the filter without gene list
# 13198 targets passed the filter when removing oscillating

# fit model
so <- sleuth_fit(so)

# get names of coefficients
models(so)

# Wald test for differential expression of isoforms
oe <- sleuth_wt(so,
                which_beta = 'condition828')

saveRDS(oe,file=paste0(outPath,"/sleuth/sleuthObject.rds"))

p1<-plot_pc_variance(oe)
p2<-plot_loadings(oe,pc_input=1L)
p2a<-plot_loadings(oe,pc_input=2L)

p3<-plot_pca(oe,
         color_by = 'replicate',
         text_labels = TRUE)

p4<-plot_pca(oe,
         color_by = 'condition',
         text_labels = TRUE)

p4a<-plot_pca(oe, pc_x=3L, pc_y=4L,
             color_by = 'condition',
             text_labels = TRUE)

#plot_sample_heatmap(oe)

p5<-plot_group_density(oe,
                   use_filtered = FALSE,
                   units = "est_counts",
                   trans = "log",
                   grouping = "condition")

p6<-plot_group_density(oe,
                   use_filtered = TRUE,
                   units = "est_counts",
                   trans = "log",
                   grouping = "condition")


p<-ggpubr::ggarrange(p1,p2,p2a,p3,p4,p4a,p5,p6,nrow=3,ncol=3)
p
ggplot2::ggsave(paste0(outPath,"/sleuth/QC_sleuth.pdf"),p,dev="pdf",
       height=13,width=15)

# output results
res <- sleuth_results(oe,
                      test = 'condition828',
                      show_all = TRUE)

res<-left_join(res,data.frame(txMeta),by=join_by(target_id==txptSeqID))

plotTxPvalQC(res,contrastName=grp,outPath=outPath, fileNamePrefix=fileNamePrefix)


# look at skn-1
skn1<-res$target_id[res$publicID=="skn-1"]
tmp<-res[res$publicID=="skn-1",]
tmp<-tmp[complete.cases(tmp),]
tmp$txt<-round(tmp$qval,2)
tmp$txt[tmp$txt>0.05]<-"n.s."
tmp$txtPos<-0.2
tmp$label<-c("a","b")

p<-ggplot(tmp,aes(x=label,y=b)) +
  geom_bar(stat="identity",fill="lightgrey",color="black") +
  xlab("*skn-1* isoform") + ylab(label="Log<sub>2</sub>FC") +
  geom_errorbar(aes(ymin=b-se_b,ymax=b+se_b),width=0.2) +
  geom_hline(yintercept=0) + coord_cartesian(ylim=c(-0.23,0.21)) +
  geom_text(aes(y=txtPos,label=txt))
p
ggplot2::ggsave(paste0(outPath,"/sleuth/skn1_bar.pdf"),p,dev="pdf",
       height=3,width=1.5)

# counts1<-sleuth_to_matrix(so,"obs_norm","est_counts")
# tmp<-counts1[grep("T19E7.2",rownames(counts1))[1:2],]
# tmp
# ratios<-data.frame(tmp[1,]/tmp[2,])
# ratios$name<-gsub("_.*","",rownames(ratios))
# colnames(ratios)<-c("ratio","sample")
# ggplot(ratios,aes(x=sample,y=ratio)) + geom_jitter(width=0.1)

# plotList<-list()
# for(tid in skn1){
#   if(tid %in% res$target_id & !is.na(res[res$target_id==tid,"pval"])){
#     plotList[[tid]]<-
#       plot_bootstrap(oe,
#                target_id = tid,
#                units = "tpm",
#                color_by = "condition")
#   }
# }
#
# p<-ggpubr::ggarrange(plotlist=plotList,nrow=2,ncol=2)
# ggsave(paste0(outPath,"/sleuth/skn1.pdf"),p,dev="pdf",
#        height=13,width=15)



sigtxpts <- res %>%
  filter(qval < 0.05)


#plot_transcript_heatmap(oe,
#                        transcripts = sigtxpts$target_id[1:20])

plot_bootstrap(oe,
               target_id = sigtxpts$target_id[3],
               units = "est_counts",
               color_by = "condition")
#so$sample_to_covariates

plot_mean_var(oe)
plot_qq(oe,test="condition828",test_type="wt")
plot_ma(oe,test="condition828",test_type="wt")
#sleuth_live(oe)

write.csv(res,paste0(outPath,"/sleuth/coh1cs_DTE.csv"))
gr<-makeGRangesFromDataFrame(res,keep.extra.column=T)
seqinfo(gr)<-seqinfo(Celegans)
gr<-sort(gr)
saveRDS(gr,paste0(outPath,"/sleuth/coh1cs_DTE.RDS"))

#######################-
# find DTE ------
#######################-
so<-readRDS(paste0(outPath,"/sleuth/coh1cs_DTE.RDS"))

data.frame(so) %>% dplyr::filter(!is.na(qval),qval<0.05,b>0) %>% nrow()
upGenes<-data.frame(so) %>% dplyr::filter(!is.na(qval),qval<0.05,b>0) %>% dplyr::select(wormbaseID) %>% unique()
dim(upGenes)
write.table(upGenes,"sigGenesUp_DTEsleuth.txt")

data.frame(so) %>% dplyr::filter(!is.na(qval),qval<0.05,b<0) %>% nrow()
downGenes<-data.frame(so) %>% dplyr::filter(!is.na(qval),qval<0.05,b<0) %>% dplyr::select(wormbaseID) %>% unique()
dim(downGenes)
write.table(downGenes,"sigGenesDown_DTEsleuth.txt")



plot(euler(c(up=upGenes,down=downGenes)))
ggvenn(list(up=upGenes$wormbaseID,down=downGenes$wormbaseID))
sum(upGenes$wormbaseID %in% downGenes$wormbaseID)


uptx<-data.frame(so) %>% dplyr::filter(wormbaseID %in% upGenes$wormbaseID) %>%
  dplyr::group_by(wormbaseID) %>%
  dplyr::mutate(count=n(),notNA=which(!is.na(qval))[1],sig1=any(qval<0.05)) %>%
  dplyr::filter(count>1) %>% mutate(isoRatio=b-b[notNA]) %>%
  dplyr::filter(!is.na(isoRatio),isoRatio>0.1,sig1==T)

print(uptx,width=Inf)
dim(uptx)
length(unique(uptx$wormbaseID))
diffuptx<-unique(uptx$wormbaseID)
so$diffuptx<-F
so$diffuptx[so$wormbaseID %in% diffuptx]<-T
length(so[so$diffuptx==T])
length(unique(so[so$diffuptx==T]$wormbaseID))
saveRDS(so[so$diffuptx==T],file=paste0(outPath,"/sleuth/coh1cs_DTE_isoformDiffratio.RDS"))

write.table(unique(so[so$diffuptx==T]$wormbaseID),"sigGenesUp_multiTxptDiffReg.txt")
length(unique(so[so$diffuptx==T]$wormbaseID))

so<-readRDS(paste0(outPath,"/sleuth/coh1cs_DTE_isoformDiffratio.RDS"))
length(so)


fountains<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/detected_fountains_equalQ.RDS")
fountains$fountainName<-paste0("fount",1:length(fountains))
fountains<-resize(fountains, width=6000,fix="center")

nonFount<-gaps(fountains)
nonFount<-resize(nonFount,width=6000,fix="center")

so$olFount<-F
ol<-findOverlaps(so[so$diffuptx],fountains,ignore.strand=T)
so$olFount[queryHits(ol)]<-T
table(so$olFount,so$diffuptx)

so$olGap<-F
ol<-findOverlaps(so[so$diffuptx],nonFount,ignore.strand=T)
so$olGap[queryHits(ol)]<-T
table(so$olGap,so$diffuptx)



write.csv(so[so$publicID %in% c("daf-16","par-1","hlh-30","casy-1","egl-44","src-1","daf-3"),],
          paste0(outPath,"/sleuthDTE_genesOfInterest.csv"),quote=F,row.names=F)
