# 
# Host gene expression analysis
# Clustering analysis of tissue gene expression.
# Xiao et al. (2017) 
# Patrick T. Dolan, Ph.D. (ptdolan@stanford.edu)
#

# # # # # # #
# Requirements
# # # # # # #

# source("https://bioconductor.org/biocLite.R")
# biocLite("cummeRbund")
# biocLite("limma")

library(ggplot2)
library(RColorBrewer) 
library(limma)
library(reshape2)
library(IDPmisc)
library(viridis)
library(boot)
library(gplots)

# # # # # Mus musculus GO table (http://www.informatics.jax.org/)
mmusAnnots<-read.delim("~/Downloads/gene_association.mgi",sep = "\t",comment.char = "!",header=F)
colnames(mmusAnnots)<-c('db',"db.geneid","gene","none","GO.ID","ID","Evidence","GO Type","gene.info","none2","elementtype","taxon","gene.num","db2")

# # # # # # # # # # # # # # # #
# # # # # Read Data # # # # # #
# # # # # # # # # # # # # # # #
root="~/GitHub/Xiaoetal_2017_analysis/RawCuffDiffData/" #path to cuffdiff files. 

#Anti-PV genes from Schoggins, et al., 
SchogginsGenes<-c("IRF1",
                  "APOL6",
                  "TRIM25",
                  "EHD4",
                  "TRAFD1",
                  "ELF1",
                  "THBD",
                  "RARRES3",
                  "EPSTI1",
                  "TLR3",
                  "ZBP1",
                  "TDRD7",
                  "MX1",
                  "DDX58")

#save for quick loading later

#save(mmusAnnots,file = "~/GitHub/GeneExpression/mmusAnnots.Rdata")
#load("~/GitHub/GeneExpression/mmusAnnots.Rdata")

###########
#FUNCTIONS#
###########

#Takes GO ID (e.g. "GO:174645") and returns member genes. 
GOselect<-function(terms){ 
  select<-allgenes$gene %in% mmusAnnots$gene[mmusAnnots$GO.ID%in%terms]
  return(select)
}

#Computes genewise distance matrix, Dmat, with selected method for dist() (e.g "euclidian","manhattan").
genewiseDist<-function(X,methodSel){ 
  print(rownames(X))
  print(paste("Pairwise gene distances by",methodSel,"..."))
  print("parallel distance function")
  Dmat<-dist(X,methodSel)
  heatmap.2(as.matrix(Dmat), scale ='none',col=viridis(100),trace = "none",labCol = rownames(Dmat), rowCol = rownames(Dmat), symm = T)
  return(Dmat)
}

#Computes sample-wise distance with selected method (e.g "euclidian","manhattan").
sampleDist<-function(X,methodSel){ 
  print(paste("Samplewise distances by",methodSel,"..."))
  Dmat<-dist(t(X),methodSel)
  heatmap.(as.matrix(Dmat), scale ='none')
  return(Dmat)
}

# Ward clustering and cutree to choose k clusters.
ward2C<-function(Dmat,kSel){ 
  print(paste("...Clustering with",kSel))
  HD<-cluster::agnes(Dmat,method = 'ward')
  CT<-cutree(HD,k = kSel)
  return(CT)
}

plotClusters<-function(FCTable, plotpath="/"){ 
  MDSoutput<-FCTable[FCTable$variable %in% c("X1","X2"),]
  meltedMDS<-melt(MDSoutput)
  
  MDSTable<-dcast(MDSoutput,formula = gene+cluster+Annot~variable,fun.aggregate = mean)
  
  #generate ggplot object with MDS table data 
  MDS<-ggplot(MDSTable)+scale_color_brewer(palette = "Spectral")
  
  ggsave(MDS+geom_point(aes(X1,X2,col=cluster))+
           geom_text(aes(X1,X2,label=toupper(gene)),cex = 1,nudge_x=.2)+  # can add to filter: FCTable<-clustAnalysis(selectSig,7,"manhattan")
           ylab("Expression Pattern Dim. 2")+
           xlab("Expression Pattern Dim. 1"),filename = paste(plotpath,"MDSgeneplot_withClustersall.pdf",sep="/"),width=5.5,height=5)
  
  MDS<-ggplot(MDSTable)+scale_color_brewer(palette = "Set2")
  
  MDS+geom_point(aes(X1,X2,col=Annot))+
    geom_text(aes(X1,X2,label=toupper(gene)),cex = 1,nudge_x=.2)+ 
    ylab("Expression Pattern Dim. 2")+
    xlab("Expression Pattern Dim. 1")
  ggsave(paste(plotpath,"MDSgeneplot_annot.pdf",sep="/"),width=5.5,height=5)
  
  FCTable<-FCTable[!FCTable$variable %in% c("X1","X2")&!FCTable$genotype %in% c("Ref","logq.GD","logq.WTB","logq.WTC"),]
  
  GG<-ggplot(FCTable)+scale_color_brewer(palette = "Spectral")
  GG+geom_hline(yintercept=0)+geom_line(aes(numLab,value,group=gene,color=cluster),alpha=0.3)+
    coord_cartesian(ylim = c(-5,5))+
    stat_summary(aes(numLab,value,group=cluster,color=cluster), geom = "line", fun.y = median, size = 1)+facet_grid(tissue~genotype)+theme_minimal()+
    ylab(expression('log'[2]*' Fold Change'))+
    xlab("Days Post-infection")
  ggsave(paste(plotpath,"TrajectoryPlot_withClustersall.pdf",sep=""),width=7,height=4)
  
  IFNmat<-acast(FCTable, variable~gene ,value.var = "value" )
  IFN_D<-dist(IFNmat,method = "manhattan")
  MDS_IFN<-cmdscale(IFN_D)
  MDS.DF<-data.frame(MDS_IFN,strsplit2(rownames(MDS_IFN),"_"))
  colnames(MDS.DF)<-c("D1","D2","genotype","time","tissue")
  plot(ggplot(MDS.DF)+geom_point(aes(D1,D2,color=tissue,pch=genotype),size=4)+
         geom_text(aes(D1,D2,label=time),nudge_x = 6)+
         scale_color_brewer(palette = "Dark2"))
  ggsave("ExpressionDistanceall.pdf",path=plotpath,width = 5.5,height=5)
}

#perform clustering on the 
clustAnalysis<-function(ingenes, kSel, distMeth,suffix=""){
  Qmat<-dcast(ingenes,gene_id+gene+Annot~genotype+time+tissue,value.var = "q_value")
  FCmat<-dcast(ingenes,gene_id+gene+Annot~genotype+time+tissue,value.var = "log2.fold_change.")
  FCQmat<-data.frame(FCmat,Ref_d0_Ref=0,WTB_d0_LV=0,WTC_d0_LV=0,GD_d0_LV=0,WTB_d0_KD=0,WTC_d0_KD=0,GD_d0_KD=0,logq=log10(Qmat[,-(1:4)]))
  FCQmat<-FCmat[complete.cases(FCmat[,4:ncol(FCmat)]),]
  
  print(head(FCQmat))
  
  Dmat<-genewiseDist(FCQmat[,4:ncol(FCQmat)],distMeth)
  DC<-ward2C(Dmat,kSel)
  Dscale<-cmdscale(Dmat)
  
  FCmat<-data.frame(FCQmat,cluster=as.factor(DC),Dscale) #add both clusterings
  
  FCs<-melt(FCQmat)
  labels<-data.frame(strsplit2(FCs$variable,"_"),stringsAsFactors = T)
  colnames(labels)<-c("genotype","time","tissue")
  labels<-data.frame(labels,numLab=as.numeric(strsplit2(labels$time,"d")[,2]))
  output<-data.frame(FCs,labels)
  return(output)
}
# # generate input datasets # #
allgenes<-NULL
for(direc in list.dirs(full.names = F,path = root)[2:(length(list.dirs(full.names = T,path = root)))]){
  path=paste(root,direc,sep="/")
  label<-strsplit(direc,"_")[[1]][1]
  print(label)
  labelvector<-strsplit(label,split = "")[[1]]
  veclength<-length(labelvector)
  tissue<-as.factor(paste(labelvector[(veclength-1):veclength],collapse = ""))
  time<-as.factor(paste(labelvector[(veclength-3):(veclength-2)],collapse = ""))
  genotype<-as.factor(paste(labelvector[1:(veclength-4)],collapse=""))
  diffFile=paste(path,"gene_exp.diff",sep="/")
  diffs<-read.delim(diffFile)
  diffs<-data.frame(diffs,label,genotype,tissue,time)
  print(dim(diffs))
  allgenes<-rbind(allgenes,diffs)
}
input<-allgenes

save(list = "input",file = "~/GitHub/Xiaoetal_2017_analysis/allGene_cuffdiff_3-27.Rdata")

#Save time and reload dataset as Rdata
load("~/GitHub/Xiaoetal_2017_analysis/allGene_cuffdiff_3-27.Rdata")

allgenes<-input

#Annotate antiviral genes
allgenes$Annot<-NULL

#GOterm classes
NegReg<-GOselect("GO:0048525")
Immune<-GOselect("GO:0002376")
Innate<-GOselect(c("GO:0045087","GO:0045088"))

#Immune and Innate immune classes all included in Fig. 
allgenes$Annot[Immune]<-"Immune"
allgenes$Annot[NegReg]<-"Virus Restriction"
allgenes$Annot[Innate]<-"Innate Immunity"
allgenes$Annot[toupper(allgenes$gene) %in% SchogginsGenes] <- "AntiPV"
allgenes$Annot[is.na(allgenes$Annot)]<-"None"
allgenes$Annot<-factor(allgenes$Annot)

#Filter Criteria
allgenes<-allgenes[allgenes$status=="OK", ]
allgenes$log2.fold_change.[!is.finite(allgenes$log2.fold_change.)&allgenes$log2.fold_change.>0] <- NA
allgenes$log2.fold_change.[!is.finite(allgenes$log2.fold_change.)&allgenes$log2.fold_change.<0] <- NA
write.csv(file = "hostgeneexpressionvalues_PVinfection.csv", allgenes)

sig<-allgenes[allgenes$q_value<0.01&abs(allgenes$log2.fold_change.)>2,]
sigLV<-allgenes[allgenes$tissue=="LV"&allgenes$q_value<0.01&abs(allgenes$log2.fold_change.)>2,]
sigKD<-allgenes[allgenes$tissue=="KD"&allgenes$q_value<0.01&abs(allgenes$log2.fold_change.)>2,]
dim(sig)

#Join Fold change, q-value and ctrl measurements
FCcast<-dcast(allgenes,formula = gene+gene_id~tissue+genotype,value.var = "log2.fold_change.",fun.aggregate = function(X){mean(X,na.rm=TRUE)})
Qcast<-dcast(allgenes,formula = gene+gene_id~tissue+genotype,value.var = "q_value",fun.aggregate = function(X){mean(X,na.rm=TRUE)})
ctrl<-dcast(allgenes,formula = gene+gene_id~tissue+genotype,value.var = "value_1",fun.aggregate = function(X){mean(X,na.rm=TRUE)})

jointcast<-data.frame(FC=FCcast,q=Qcast[,-c(1:2)],ctrl=ctrl[,-c(1:2)])

PC<-princomp(NaRV.omit(jointcast[,-c(1:2)]),cor = T)

data<-jointcast[complete.cases(jointcast[,-c(1:2)]),]
combined<-data.frame(data,PC$scores)
#Select significant genes
siggenes<-as.character(unique(sig$gene)) #565

#ggplot(allgenes)+geom_boxplot(aes(Annot,log2.fold_change.))+ylim(-3,3)+facet_wrap(tissue+genotype~time)
Init=1

# # # # # # # # #
###### MAIN #####
# # # # # # # # #

FCTable<-clustAnalysis(selectSig,5,"euclidian")
save(file = "~/GitHub/GeneExpression/ClusteredTable.RData",FCTable)
load("~/GitHub/GeneExpression/ClusteredTable.RData")
suffix=""
plotpath = paste(root,suffix,sep="")
plotClusters(FCTable,plotpath)
str(FCTable)
#FCall<-clustAnalysis(allgenes,7,"euclidian")

IFNmat<-acast(FCTable[!FCTable$variable %in% c("X1","X2"),], variable~gene ,value.var = "value",fun.aggregate = mean )

geneMat<-t(IFNmat)
geneDF<-data.frame(gene=rownames(geneMat),geneMat[,-c(13:19)])
melted<-melt(geneDF)
geneDF<-data.frame(melted,stringsAsFactors = T)
casted<-dcast(geneDF,gene~variable)

# # # # # # # # # # #
# # Generate and save all plots and tables
# # # # # # # # # # #

melty<-melt(FCTable,measure.vars = "value")

castedTiss<-dcast(melty,gene+genotype+time+cluster+Annot~tissue,value.var = "value",fun.aggregate = mean)

ClusterTissueExpression<-ggplot(castedTiss[castedTiss$genotype%in%c("GD","WTB","WTC")&castedTiss$time!="d0",],aes(KD,LV,color=cluster))+
  coord_cartesian(ylim = c(-5,5))+geom_abline(slope=1)+geom_point(cex=0.5)+
  #geom_text(data=aes(KD,LV,label=ifelse(abs(KD-LV)>1,toupper(as.character(gene)),""),cex=1),color='black')+
  facet_grid(genotype~time)+
  theme_bw()+ylab("expression('Liver Expression log'[2]*' Fold Change')")+xlab(expression('Kidney Expression log'[2]*' Fold Change'))

AnnotTissueExpression<-ggplot(castedTiss[castedTiss$genotype%in%c("GD","WTB","WTC")&castedTiss$time!="d0",],aes(KD,LV,color=Annot))+
  coord_cartesian(ylim = c(-5,5))+geom_abline(slope=1)+geom_point(cex=0.5,alpha=0.5)+
  #geom_text(data=aes(KD,LV,label=ifelse(abs(KD-LV)>1,toupper(as.character(gene)),""),cex=1),color='black')+
  facet_grid(genotype~time)+
  theme_bw()+ylab("expression('Liver Expression log'[2]*' Fold Change')")+xlab(expression('Kidney Expression log'[2]*' Fold Change'))

castedTiss<-dcast(melty,gene+genotype+time+cluster+tissue+Annot+numLab~.,value.var = "value",fun.aggregate = mean)  

outputformat<-dcast(castedTiss,formula = toupper(gene)+cluster~time+tissue+genotype,value.var = ".")

write.table(file = "~/GitHub/GeneExpression/allgenetable.txt",outputformat,sep = "\t",quote = F,row.names = F)

ClusterExpressionBoxplot<-
  ggplot(castedTiss[castedTiss$genotype%in%c("WTB","WTC","GD")&castedTiss$time!="d0",],aes(cluster,.))+
  geom_hline(aes(yintercept=0))+
  geom_point(aes(color=cluster),position="jitter",cex=.5,alpha=1)+geom_boxplot(outlier.shape = "")+
  facet_grid(genotype~time+tissue)+
  theme_bw()+ylab(expression('log'[2]*' Fold Change'))+xlab("Cluster")

ClusterExpressionLines<-
  ggplot(castedTiss[castedTiss$genotype%in%c("WTB","WTC","GD")&castedTiss$time!="d0",],aes(time,.))+
  coord_cartesian(ylim=c(-5,5))+
  geom_hline(aes(yintercept=0))+
  geom_line(aes(group=gene,color=cluster),cex=.25,alpha=0.25)+
  geom_smooth(aes(group=cluster,color=cluster),fill='white',cex=1,method = 'gam')+
  facet_grid(genotype~tissue)+
  theme_bw()+ylab("Time")+xlab("Cluster")

ggsave(paste(root,"TissueExpressionClustersall.pdf",sep=""),ClusterTissueExpression+scale_color_brewer(palette = "Spectral"),width = 5,height = 3.5)
ggsave(paste(root,"TissueExpressionClustersSmooth.pdf",sep=""),ClusterExpressionLines+scale_color_brewer(palette = "Spectral"),width = 5,height = 3.5)
ggsave(paste(root,"TissueExpressionBoxplotall.pdf",sep=""),ClusterExpressionBoxplot+scale_color_brewer(palette = "Spectral"),width = 6,height = 3.5)
ggsave(filename = 'pairplot.pdf',pairPlot(sig,ggstuff = 'scale_y_log10()+scale_x_log10()'))

clusters<-unique.data.frame(data.frame(castedTiss$gene,castedTiss$cluster))
write.table(file = paste(root,"ClusterAssignments.txt"),clusters,quote = F,sep = "\t",row.names = F)

allresting=dcast(allgenes,formula = gene~genotype+tissue+time,value.var="value_1",mean)   
allcounts=dcast(allgenes,formula = gene~genotype+tissue+time,value.var="value_2",mean)
allqval=dcast(allgenes,formula = gene~genotype+tissue+time,value.var="q_value",mean)
allFC=dcast(allgenes,formula = gene~genotype+tissue+time,value.var="log2.fold_change.",mean)
allcounts$gene<-toupper(allcounts$gene)

alldata<-data.frame(allcounts,fc= allFC[-1], resting=allresting[,-1],qval=allqval[,-1])
write.table(x = alldata,quote = F,row.names = F,file = "cuffDiffExpVals.txt",sep = "\t")

