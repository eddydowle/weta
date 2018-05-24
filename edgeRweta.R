#edgeR for preliminary weta data 

library(edgeR)
library(stringr)
library(locfit)
library(statmod)
library(plyr)
library(ggplot2)
library(reshape)
library(heatmap3)
setwd("/Users/edwinadowle/Documents/Earwigs/Weta/FirstRun/weta/salmon/")
#setwd("/Users/edwinadowle/Documents/Earwigs/FirstRuns/March2018/Earwigs_salmon/")

#all cerosi work is now done with file that have had funny characters that end the sequence name removed, causing too many issues

warnings()
file_list<-list.files(pattern='*.sf') #gernam
#file_list<-list.files(pattern='.trinity.quant.sf') #trinity
file_list
table<-c()

for (file in file_list){
  sampleID<-str_match(file,"(^[A-Za-z0-9]*)")[,1]
  #  sampleID<-str_match(file,"(.*)?_.trinity.ours.ncbi.quant.sf")[,2]
  #  if (is.na(sampleID)){
  #    sampleID<-str_match(file,"(^[a-z]*_[A-Za-z]*_[A-Z0-9]*_[0-9]*)")[,1]
  #  }
  print(sampleID)
  df<-read.table(file,header=T,sep="\t",row.names=NULL)
  df2<-df[,c(1,5)]
  colnames(df2)[colnames(df2)=="NumReads"] <- sampleID
  if (length(table) ==0){table<-df2
  } else {
    table<-merge(table,df2,by="Name") }
}  

table2<-table

boxplot(table[-1],las=2)
colnames(table)
?boxplot
filterMinCount<- function(x) {
  pres<-x >=1
  out=F
  if ((sum(pres)/length(pres)) >= 0.5) {out=T}
  return(out)
}

filterInd2<-apply(table2[,(-1)],1,filterMinCount)
table2 <- table2[filterInd2,] #9267 #7910

table.dge<-DGEList(counts=table2[,2:4],genes=table2[,(1)])
table.dge<- calcNormFactors(table.dge)
table.dge$samples


plotMDS(table.dge)
?plotMDS
colnames(table.dge)

#plotMDS(table2.dge, top=500,pch = c(1,1,rep(2,6),rep(3,5),rep(4,5)),col=c("green4","green4",rep("purple4",6),rep("red",5),rep("blue",5)))
plotMDS(table.dge, top=500,pch = 1,col='blue4')

colnames(table.dge)
table.ordered<-table %>% arrange(desc(rowSums(.[-1]))) 

#?read.table
anottate<-read.table('../../../FirstRun/Trinotate.weta.txt',sep="\t",row.names=NULL,header=T,quote="",na.strings='NA')
colnames(anottate)[colnames(anottate)=="transcript_id"] <- "Name"
#?merge
table.ordered.anot<-merge(table.ordered,anottate,by="Name")
colnames(table.ordered.anot)
table.ordered.anot.ordered<-table.ordered.anot %>% arrange(desc(rowSums(.[2:4]))) 
write.table(table.ordered.anot.ordered,"salmon_results.annot.txt",sep="\t",row.names=F,quote=F)

#stage<-factor(c(rep("After_Emergence",2),rep("Before_Emergence",6),rep("Healthy",5),rep("Resistant",5)))
stage<-factor(c(rep("infected",1),rep("uninfected",1),rep("infected",1)))
#altitute<-factor(c(rep("High",20),rep("Low",19)))
stage
table.sampleData<-data.frame(Sample=colnames(table.dge),stage)
table.sampleData
#model
levels(stage)
#month_<-factor(month, levels = c("2M", "2_5M","3_5M","4M","4_5M"))
stage<-factor(stage,levels = c ("uninfected","infected"))
#levels(month_)
#levels(altitute)
design <- model.matrix(~stage)
rownames(design) <- colnames(table.dge)
design

#> design
#(Intercept) stageinfected
#W3infected             1             1
#W4uninfected           1             0
#W8infected             1             1


stage<-factor(c(rep("dualinfected",1),rep("uninfected",1),rep("singleinfected",1)))
levels(stage)
stage<-factor(stage,levels = c ("uninfected","singleinfected","dualinfected"))
design <- model.matrix(~stage)
rownames(design) <- colnames(table.dge)
design

#(Intercept) stagesingleinfected stagedualinfected
#W3infected             1                   0                 1
#W4uninfected           1                   0                 0
#W8infected             1                   1                 0

table.dge<- estimateDisp(table.dge, design, robust=TRUE)
table.dge$common.dispersion
#
sqrt(1.350319)
#0.7 fucken high coefficent of variation


#plot
plotBCV(table.dge)

#can look at trended dispersion that is used for QL (Quasi-Likelihood pipeline)
fit <- glmFit(table.dge, design)
colnames(design)

#fit <- glmQLFit(table.dge, design)
#plotQLDisp(fit)


##############contrasts##################

fit <- glmFit(table.dge, design)

colnames(design)

infected.uninfected<- glmLRT(fit,contrast=c(0,1))
infected.uninfected.toptags <- data.frame(topTags(infected.uninfected,n=nrow(table),sort="none"))
infected.uninfected.toptags

infected.uninfected.toptags.sig<-infected.uninfected.toptags[infected.uninfected.toptags$FDR < 0.05,]
infected.uninfected.toptags.sig

infected.uninfected.toptags.sig %>% arrange(desc(rowSums(.[-1]))) 
library(data.table)
infected.uninfected.toptags.sig.ordered<-setDT(infected.uninfected.toptags.sig)[order(-abs(logFC)), .SD]


write.table(infected.uninfected.toptags.sig.ordered,"InfectedUninfectedContrastedgeR.toptags.sig.ordered.txt",sep="\t",row.names=F,quote=F)



infected.uninfected.toptags.ordered<-setDT(infected.uninfected.toptags)[order(-abs(logFC)), .SD]


write.table(infected.uninfected.toptags.ordered,"InfectedUninfectedContrastedgeR.toptags.ordered.txt",sep="\t",row.names=F,quote=F)




