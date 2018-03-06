# 
# Host gene expression analysis
# Clustering analysis of tissue gene expression.
# Xiao et al. (2017) 
# Patrick T. Dolan, Ph.D. (ptdolan@stanford.edu)
#

# # # # # # #
# Requirements
# # # # # # #

library(ggplot2)
library(tidyr)
library(limma)
library(reshape2)

#Fix path do data...
path<-"~/Path/To/InVivoPVRNAseq/"

#limits
nucLims<-1:7404

Names<-c("UTR5","VP4","VP2","VP3","VP1",  "2A", "2B", "2C", "3A", "Vpg", "3C","3D","UTR3")
LimitsRes<-c(1,   70,   342,  580,   882, 1031,  1129, 1457, 1544, 1565, 1748,  2209)
LimitsNuc<-c(1,743+(3*(LimitsRes)),7404)

elem<-1
annotV<-NULL
for(Reg in Names){
  elem<-elem+1
  annotV<-c(annotV,rep(Reg, LimitsNuc[elem]-LimitsNuc[elem-1]))
}

### Functions
coverage<-function(X){ #coverage based on wide data
  apply(X,1,FUN=function(S){
    sums<-sum(as.integer(S[4:7]))
    return(sums)})}

entropies<-function(X){ #coverage based on wide data
  apply(X,1,FUN=function(S){
    p<-as.integer(S[4:7])/sum(as.integer(S[4:7]))
    h<-sum(-p*log(p,base=4),na.rm = T)
    return(h)})}

generateCounts<-function(all){
  wtVec<-NULL
  sumVec<-NULL
  for(element in 1:nrow(all)){
    wtVec<-c(wtVec,all[element,as.character(all$Nucl[element])])
    sumVec<-c(sumVec,sum(all[element,4:7]))
  }
  all$wtC<-wtVec
  all$muC<-sumVec-wtVec
  return(all)
}

### Collate files

fileList<-list.files(full.names = T,pattern = "BINOM95.txt",path)
all=NULL

for (Fl in fileList){
  nameVec<-strsplit2(Fl,split = "_")
  sample=nameVec[5]
  genotype=nameVec[6]
  time=nameVec[7]
  tissue=nameVec[8]
  currFile<-read.delim(Fl,header = T,sep = "\t")
  currFileClean<-currFile[,1:7]
  currFile<-currFileClean[currFileClean$Pos<7441,]
  currFile$genotype<-as.factor(genotype)
  currFile$sample<-as.factor(sample)
  currFile$time<-as.factor(time)
  currFile$tissue<-as.factor(tissue)
  all<-rbind(all,currFile)
}

cover<-coverage(all)

all$coverage<-as.vector(cover)

entropy<-entropies(all)
all$entropy<-as.numeric(entropy)

counts<-gather(data = all,mut,count,A:T)

# # # Remove reference differences and Z
clean<-counts[counts$Nucl!=counts$mut&(!counts$Pos%in%c(2133,6176,6177,6178,6221)),]
clean$freq<-clean$count/clean$coverage

write.csv(file = "topmuts_RNAseqTissues.csv",clean[ clean$freq>.05, ])
# # # Filter
#clean_filt<-clean[clean$freq>0.01&clean$freq<0.05,]
clean_filt<-clean[clean$count>10,]
clean_filt$exp<-as.factor(paste(clean_filt$genotype,clean_filt$tissuetime))
# # # ggplot block
ggplot(clean)+geom_density(aes(count))+scale_x_log10()+facet_wrap(genotype+tissue~sample)

Colors<- brewer.pal(4,"Dark2")[c(3,1,2,4)]

C<-ggplot(clean_filt)+theme_classic()
C+geom_vline(xintercept = LimitsNuc, lwd = 0.25)+stat_summary(aes(Pos,entropy,group=tissue,color=tissue),cex=0.1,fun.data = mean_se )+facet_wrap(genotype~tissue)
ggsave("ComparisonOfEntropyValues.pdf",C+geom_vline(xintercept = LimitsNuc, lwd = 0.25)+stat_summary(aes(Pos,entropy,color=tissue),lwd=0.5,geom = "bar",fun.y = mean )+facet_grid(genotype~time)+scale_color_manual(values = Colors[4:2]))
C+geom_vline(xintercept = LimitsNuc)+stat_summary(aes(Pos,freq,group=tissue,color=tissue),cex=0.2,fun.data = mean_se )+scale_y_log10()+facet_wrap(tissue~genotype)C+geom_point(aes(Pos,entropy,color=sample))

melted<-melt(clean_filt[clean_filt$coverage>50000,], measure.vars = c(8,9,11,12))

thing<-dcast(melted[melted$variable=="entropy",],formula = Pos+tissue~genotype,fun.aggregate = mean,value.var = "value")

logENT<-ggplot(thing)+scale_y_log10()+scale_x_log10()+scale_color_manual(values = Colors[4:2])+theme_bw()
ENT<-ggplot(thing)+ylim(0,0.5)+xlim(0,0.5)+scale_color_manual(values = Colors[4:2])+theme_bw()
lmCG<-with(thing,lm(WTC~GD))
with(thing,cor.test(WTC,WTB,method = 'kendall'))
summary(lmCG)
ENT+geom_smooth(aes(WTC,GD,color=tissue),method = 'lm')+geom_point(aes(WTC,GD,color=tissue),alpha=0.6)
lmBG<-with(thing,lm(WTB~GD))
summary(lmBG)
ENT+geom_smooth(aes(WTB,GD,color=tissue),method = 'lm')+geom_point(aes(WTB,GD,color=tissue),alpha=0.6)
ENT+geom_smooth(aes(WTC,WTB,color=tissue),method = 'lm')+geom_point(aes(WTC,WTB,color=tissue),alpha=0.6)

logENT+geom_point(aes(WTC,GD,color=tissue))
logENT+geom_point(aes(WTB,GD,color=tissue))
logENT+geom_point(aes(WTC,WTB,color=tissue))

C+geom_tile(aes(Pos,genotype,fill=entropy))+facet_grid(tissue~.)

Bh<-C+geom_boxplot(aes(genotype,entropy,fill=genotype))+scale_y_log10()+facet_wrap(tissue~.)
Vh<-C+geom_violin(draw_quantiles = c(.025,.25,.5,.75,.975),aes(tissue,entropy,fill=genotype))+scale_y_log10()+facet_wrap(~genotype)
Vf<-C+geom_violin(draw_quantiles = c(.025,.25,.5,.75,.975),aes(sample,freq,fill=genotype))+scale_y_log10()
Pf<-C+geom_point(aes(Pos,freq,color=sample),alpha=0.7)+scale_y_log10()+facet_grid(genotype~tissue,drop = T)+theme(legend.position="none")
Pc<-C+geom_point(aes(Pos,coverage,color=sample),alpha=0.7)+scale_y_log10()+facet_grid(genotype~tissue,drop = T)+theme(legend.position="none")

ggsave(Bh,path=path,width = 6.3,height=5,filename = "EntropyBox_filt.pdf")
ggsave(Vh,path=path,width = 6.3,height=5,filename = "EntropyViolin_filt.pdf")
ggsave(Vf,path=path,width = 6.3,height=5,filename = "FrequencyViolin_filt.pdf")
ggsave(Pf,path=path,width = 11,height=5,filename = "FrequencyManhattanPlot_filt.pdf")
ggsave(Pc,path=path,width = 11,height=5,filename = "CoverageManhattanPlot_filt.pdf")
######

#Cast freq. matrix for dimension reduction
freqmat<-dcast(fill = 0.0,clean_filt,formula = Pos+Nucl+mut~genotype+tissue+sample,value.var = "freq")

# # # MDS
popDist<-dist(t(cbind(freqmat[,4:ncol(freqmat)],Ref_Ref_Ctrl=rep(0,nrow(freqmat)))))
mdsPops<-as.data.frame(cmdscale(popDist))
mdsMuts<-as.data.frame(cmdscale(dist(cbind(freqmat[1:1000,4:ncol(freqmat)],Reference=rep(0,100)))))
plot(mdsMuts)
colnames(mdsPops)<-c("Genetic.Dist..1", "Genetic.Dist..2")
mdsDF<-data.frame(mdsPops)
nameInfo<-matrix(unlist(lapply(rownames(mdsPops),FUN = strsplit2,split="_")),byrow=T,ncol=3,dimnames = list(NULL,c("genotype","tissue","sample")))
nameDF<-as.data.frame(nameInfo,stringsAsFactors = T)
mdsDF<-cbind(nameDF,mdsDF)

#Plot MDS
mdsGenotype<-ggplot(mdsDF)+geom_point(aes(`Genetic.Dist..1`,`Genetic.Dist..2`,pch=genotype),alpha=.9,cex=3)
mdsTissue<-ggplot(mdsDF)+geom_point(aes(`Genetic.Dist..1`,`Genetic.Dist..2`,pch=genotype,color=tissue),alpha=.9,cex=3)+
  ylab(label = "Viral Genetic Dist. 2")+xlab("Viral Genetic Dist. 1")+
  scale_color_brewer(palette = "Dark2")
ggplot(mdsDF)+geom_point(aes(`Genetic.Dist..1`,`Genetic.Dist..2`,color=tissue, pch=genotype))

ggsave(mdsGenotype,path=path,width = 7,height=4,filename = "geneticDistGenotype_filt.pdf")
ggsave(mdsTissue,path=path,width = 7,height=4,filename = "geneticDistTissue_filt.pdf")