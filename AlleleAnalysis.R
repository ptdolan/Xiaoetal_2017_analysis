IFN<-read.delim("~/Research/CirSeq/Hong/Download032716/master_Q20threshold_PV_IFNP6_annot.txt",sep="\t")
IFN<-IFN[IFN$wtNT!=IFN$mutNT&IFN$ntpos%in%1:7439,]
WT1<-read.delim("~/Research/CirSeq/Hong/newcopiedset/PV_replicates_1/p6/WTrep1_6-Q20threshold_annot.txt",sep="\t")
WT1<-WT1[WT1$wtNT!=WT1$mutNT&WT1$ntpos%in%1:7439,]
WT2<-read.delim("~/Research/CirSeq/Hong/newcopiedset/PV_replicates_2/p6/WTrep2_6-Q20threshold_annot.txt",sep="\t")
WT2<-WT2[WT2$wtNT!=WT2$mutNT&WT2$ntpos%in%1:7439,]
GD<-read.delim("~/Research/CirSeq/Hong/Processed_Data/Renamed/G64SD79H/Rep1/P6/Q20threshold_annot.txt",sep="\t")
GD<-GD[GD$wtNT!=GD$mutNT&GD$ntpos%in%1:7439,]

cIFN<-dcast(IFN,ntpos+wtNT+resPos+wtRes+muRes+mutNT+count+coverage+synNonsyn~.,value.var="freq")
cIFN<-data.frame(cIFN,sample="IFN",order=rank(cIFN$.))
cWT1<-dcast(WT1,ntpos+wtNT+resPos+wtRes+muRes+mutNT+count+coverage+synNonsyn~.,value.var="freq")
cWT1<-data.frame(cWT1,sample="WT1",order=rank(cWT1$.))
cWT2<-dcast(WT2,ntpos+wtNT+resPos+wtRes+muRes+mutNT+count+coverage+synNonsyn~.,value.var="freq")
cWT2<-data.frame(cWT2,sample="WT2",order=rank(cWT2$.))
cGD<-dcast(GD,ntpos+wtNT+resPos+wtRes+muRes+mutNT+count+coverage+synNonsyn~.,value.var="freq")
cGD<-data.frame(cGD,sample="GD",order=rank(cGD$.))

stack<-rbind(cIFN,cWT1,cWT2,cGD)
stack$wtNT<-factor(tolower(as.character(stack$wtNT)))
stack$mutNT<-factor(tolower(as.character(stack$mutNT)))

DF<-data.frame(IFN=cIFN,WT1=cWT1[,c(7,8,10)],WT2=cWT2[,c(7,8,10)])

#alleles
ggplot(stack)+geom_density(aes(.,color=sample))+scale_x_log10()+coord_cartesian(xlim=c(1e-6,1))

#AlleleFrequencySpectrum
stack$order<-22317-stack$order

ggplot(stack)+scale_color_manual(values = c("goldenrod","black",'black','red3'))+
  geom_line(aes(group=sample,order,.,color=sample),lwd = 2,alpha=0.8)+scale_y_log10()+
  xlab("Allele Rank")+ylab("Allele Frequency")+coord_cartesian(xlim=c(1,10000),ylim=c(1E-5,1E-1))


G<-ggplot(DF)+theme_bw()+scale_y_log10()+scale_x_log10()+geom_abline()+coord_fixed()
ggsave(width = 7,height=7,filename = "WT1_2_plots.pdf",G+geom_point(aes(WT2..,WT1..,color=IFN.synNonsyn))+geom_text(aes(WT2..,WT1..,label=paste(IFN.wtRes,IFN.resPos,IFN.muRes,sep="")),check_overlap = T))
ggsave(width = 7,height=7,filename = "WT1_IFN_plot.pdf",G+geom_point(aes(IFN..,WT1..,color=IFN.synNonsyn))+geom_text(aes(IFN..,WT1..,label=paste(IFN.wtRes,IFN.resPos,IFN.muRes,sep="")),check_overlap = T))
ggsave(width = 7,height=7,filename = "WT2_IFN_plot.pdf",G+geom_point(aes(IFN..,WT2..,color=IFN.synNonsyn))+geom_text(aes(IFN..,WT2..,label=paste(IFN.wtRes,IFN.resPos,IFN.muRes,sep="")),check_overlap = T))


###### plot entropy or frequency of alleles. 

