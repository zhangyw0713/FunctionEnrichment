# source("http://bioconductor.org/biocLite.R")
BiocManager::install("KEGG.db")
#setwd("D:\\GO")
library(KEGG.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)

x<-org.Hs.egPATH
#GettheentrezgeneidentifiersthataremappedtoaKEGGpathwayID
mapped_genes<-mappedkeys(x)
#Converttoalist
xx<-as.list(x[mapped_genes])
name<-matrix(ncol=1,nrow=1)
s<-unlist(xx,use.names = FALSE)
for(i in 1:length(xx))
{

	n<-length(xx[i][[1]])

	name<-rbind(name,as.matrix(replicate(n,names(xx)[i]),ncol=1))
}
name<-na.omit(name)
kegg<-cbind(name,as.matrix(s,ncol=1))
write.table(kegg,file="Hs.gene2kegg",sep="\t",row.names=F,col.names=F)


x<-org.Mm.egPATH
#GettheentrezgeneidentifiersthataremappedtoaKEGGpathwayID
mapped_genes<-mappedkeys(x)
#Converttoalist
xx<-as.list(x[mapped_genes])
name<-matrix(ncol=1,nrow=1)
s<-unlist(xx,use.names = FALSE)
for(i in 1:length(xx))
{

	n<-length(xx[i][[1]])

	name<-rbind(name,as.matrix(replicate(n,names(xx)[i]),ncol=1))
}
name<-na.omit(name)
kegg<-cbind(name,as.matrix(s,ncol=1))
write.table(kegg,file="Mm.gene2kegg",sep="\t",row.names=F,col.names=F)


x<-org.Rn.egPATH
#GettheentrezgeneidentifiersthataremappedtoaKEGGpathwayID
mapped_genes<-mappedkeys(x)
#Converttoalist
xx<-as.list(x[mapped_genes])
name<-matrix(ncol=1,nrow=1)
s<-unlist(xx,use.names = FALSE)
for(i in 1:length(xx))
{

	n<-length(xx[i][[1]])

	name<-rbind(name,as.matrix(replicate(n,names(xx)[i]),ncol=1))
}
name<-na.omit(name)
kegg<-cbind(name,as.matrix(s,ncol=1))
write.table(kegg,file="Rn.gene2kegg",sep="\t",row.names=F,col.names=F)




xx <- as.list(KEGGEXTID2PATHID)
grep("mm",xx)->mmid
mmxx<-xx[mmid]
lapply(mmxx,length)->len
names(mmxx)->n
alln<-""
for(i in 1:length(n))
{
	alln<-c(alln,rep(n[i],len[i]))
}
unlist(mmxx)->mmkegg
cbind(alln[2:length(alln)],mmkegg)->mmg2kegg
write.table(mmg2kegg,file="2Mm.gene2kegg",sep="\t",row.names=F,col.names=F)


xx<-as.list(KEGGPATHID2NAME)
s<-unlist(xx)
s<-as.matrix(s)
write.table(s,file="kegg2name",sep="\t",row.names=T,col.names=F)


