  
BiocManager::install("org.Hs.eg.db")       
 
BiocManager::install("org.Mm.eg.db")                                 

BiocManager::install("AnnotationDbi")
#source("http://www.bioconductor.org/BiocManager::install.R") 
#BiocManager::install("GenomeInfoDb")

BiocManager::install("org.Rn.eg.db")                                 
  
BiocManager::install("GO.db")       

#setwd("D:\\GO")

library(org.Hs.eg.db)

x<-org.Hs.egGO
#GettheentrezgeneidentifiersthataremappedtoaGOID
mapped_genes<-mappedkeys(x)
#Converttoalist
xx<-as.list(x[mapped_genes])






s<-xx
id<-t(matrix(unlist(strsplit(names(unlist(s)),".GO:")),nrow=6))
id<-id[,1]
goinfo<-t(matrix(unlist(s), nrow=3))
allgo<-cbind(id,goinfo)
go_bp<-allgo[which(allgo[,4]=="BP"),]
go_cc<-allgo[which(allgo[,4]=="CC"),]
go_mf<-allgo[which(allgo[,4]=="MF"),]

library("GO.db")
xx<-as.list(GOBPANCESTOR)
s<-xx
p_go<-t(matrix(unlist(s),nrow=1))
o_go<-names(unlist(s))
o_go<-as.matrix(substr(o_go,1,10),ncol=1)
go<-cbind(o_go,p_go)
go<-go[which(go[,2]!="all"),]
colnames(go)<-c("go","parent_go")
colnames(go_bp)<-c("id","go","evidence","ontology")
all_parent_gobp<-merge(go_bp,go)
tmp<-cbind(go_bp[,c(2,1,3,4)],go_bp[,2])
colnames(tmp)<-c("go","id","evidence","ontology","parent_go")
all_parent_gobp<-rbind(all_parent_gobp,tmp)



xx<-as.list(GOCCANCESTOR)
s<-xx
p_go<-t(matrix(unlist(s),nrow=1))
o_go<-names(unlist(s))
o_go<-as.matrix(substr(o_go,1,10),ncol=1)
go<-cbind(o_go,p_go)
go<-go[which(go[,2]!="all"),]
colnames(go)<-c("go","parent_go")
colnames(go_cc)<-c("id","go","evidence","ontology")
all_parent_gocc<-merge(go_cc,go)
tmp<-cbind(go_cc[,c(2,1,3,4)],go_cc[,2])
colnames(tmp)<-c("go","id","evidence","ontology","parent_go")
all_parent_gocc<-rbind(all_parent_gocc,tmp)


xx<-as.list(GOMFANCESTOR)
s<-xx
p_go<-t(matrix(unlist(s),nrow=1))
o_go<-names(unlist(s))
o_go<-as.matrix(substr(o_go,1,10),ncol=1)
go<-cbind(o_go,p_go)
go<-go[which(go[,2]!="all"),]
colnames(go)<-c("go","parent_go")
colnames(go_mf)<-c("id","go","evidence","ontology")
all_parent_gomf<-merge(go_mf,go)
tmp<-cbind(go_mf[,c(2,1,3,4)],go_mf[,2])
colnames(tmp)<-c("go","id","evidence","ontology","parent_go")
all_parent_gomf<-rbind(all_parent_gomf,tmp)

write.table(all_parent_gobp,file="hs.gene2go_bp_2.txt",sep="\t",row.names=FALSE,col.names=TRUE)
write.table(all_parent_gocc,file="hs.gene2go_cc_2.txt",sep="\t",row.names=FALSE,col.names=TRUE)
write.table(all_parent_gomf,file="hs.gene2go_mf_2.txt",sep="\t",row.names=FALSE,col.names=TRUE)



#####mouse





library(org.Mm.eg.db)
x<-org.Mm.egGO
#GettheentrezgeneidentifiersthataremappedtoaGOID
mapped_genes<-mappedkeys(x)
#Converttoalist
xx<-as.list(x[mapped_genes])






s<-xx
id<-t(matrix(unlist(strsplit(names(unlist(s)),".GO:")),nrow=6))
id<-id[,1]
goinfo<-t(matrix(unlist(s), nrow=3))
allgo<-cbind(id,goinfo)
go_bp<-allgo[which(allgo[,4]=="BP"),]
go_cc<-allgo[which(allgo[,4]=="CC"),]
go_mf<-allgo[which(allgo[,4]=="MF"),]

library("GO.db")
xx<-as.list(GOBPANCESTOR)
s<-xx
p_go<-t(matrix(unlist(s),nrow=1))
o_go<-names(unlist(s))
o_go<-as.matrix(substr(o_go,1,10),ncol=1)
go<-cbind(o_go,p_go)
go<-go[which(go[,2]!="all"),]
colnames(go)<-c("go","parent_go")
colnames(go_bp)<-c("id","go","evidence","ontology")
all_parent_gobp<-merge(go_bp,go)
tmp<-cbind(go_bp[,c(2,1,3,4)],go_bp[,2])
colnames(tmp)<-c("go","id","evidence","ontology","parent_go")
all_parent_gobp<-rbind(all_parent_gobp,tmp)



xx<-as.list(GOCCANCESTOR)
s<-xx
p_go<-t(matrix(unlist(s),nrow=1))
o_go<-names(unlist(s))
o_go<-as.matrix(substr(o_go,1,10),ncol=1)
go<-cbind(o_go,p_go)
go<-go[which(go[,2]!="all"),]
colnames(go)<-c("go","parent_go")
colnames(go_cc)<-c("id","go","evidence","ontology")
all_parent_gocc<-merge(go_cc,go)
tmp<-cbind(go_cc[,c(2,1,3,4)],go_cc[,2])
colnames(tmp)<-c("go","id","evidence","ontology","parent_go")
all_parent_gocc<-rbind(all_parent_gocc,tmp)


xx<-as.list(GOMFANCESTOR)
s<-xx
p_go<-t(matrix(unlist(s),nrow=1))
o_go<-names(unlist(s))
o_go<-as.matrix(substr(o_go,1,10),ncol=1)
go<-cbind(o_go,p_go)
go<-go[which(go[,2]!="all"),]
colnames(go)<-c("go","parent_go")
colnames(go_mf)<-c("id","go","evidence","ontology")
all_parent_gomf<-merge(go_mf,go)
tmp<-cbind(go_mf[,c(2,1,3,4)],go_mf[,2])
colnames(tmp)<-c("go","id","evidence","ontology","parent_go")
all_parent_gomf<-rbind(all_parent_gomf,tmp)

write.table(all_parent_gobp,file="Mm.gene2go_bp.txt",sep="\t",row.names=FALSE,col.names=TRUE)
write.table(all_parent_gocc,file="Mm.gene2go_cc.txt",sep="\t",row.names=FALSE,col.names=TRUE)
write.table(all_parent_gomf,file="Mm.gene2go_mf.txt",sep="\t",row.names=FALSE,col.names=TRUE)






#####rat





library(org.Rn.eg.db)
x<-org.Rn.egGO
#GettheentrezgeneidentifiersthataremappedtoaGOID
mapped_genes<-mappedkeys(x)
#Converttoalist
xx<-as.list(x[mapped_genes])






s<-xx
id<-t(matrix(unlist(strsplit(names(unlist(s)),".GO:")),nrow=6))
id<-id[,1]
goinfo<-t(matrix(unlist(s), nrow=3))
allgo<-cbind(id,goinfo)
go_bp<-allgo[which(allgo[,4]=="BP"),]
go_cc<-allgo[which(allgo[,4]=="CC"),]
go_mf<-allgo[which(allgo[,4]=="MF"),]

library("GO.db")
xx<-as.list(GOBPANCESTOR)
s<-xx
p_go<-t(matrix(unlist(s),nrow=1))
o_go<-names(unlist(s))
o_go<-as.matrix(substr(o_go,1,10),ncol=1)
go<-cbind(o_go,p_go)
go<-go[which(go[,2]!="all"),]
colnames(go)<-c("go","parent_go")
colnames(go_bp)<-c("id","go","evidence","ontology")
all_parent_gobp<-merge(go_bp,go)
tmp<-cbind(go_bp[,c(2,1,3,4)],go_bp[,2])
colnames(tmp)<-c("go","id","evidence","ontology","parent_go")
all_parent_gobp<-rbind(all_parent_gobp,tmp)



xx<-as.list(GOCCANCESTOR)
s<-xx
p_go<-t(matrix(unlist(s),nrow=1))
o_go<-names(unlist(s))
o_go<-as.matrix(substr(o_go,1,10),ncol=1)
go<-cbind(o_go,p_go)
go<-go[which(go[,2]!="all"),]
colnames(go)<-c("go","parent_go")
colnames(go_cc)<-c("id","go","evidence","ontology")
all_parent_gocc<-merge(go_cc,go)
tmp<-cbind(go_cc[,c(2,1,3,4)],go_cc[,2])
colnames(tmp)<-c("go","id","evidence","ontology","parent_go")
all_parent_gocc<-rbind(all_parent_gocc,tmp)


xx<-as.list(GOMFANCESTOR)
s<-xx
p_go<-t(matrix(unlist(s),nrow=1))
o_go<-names(unlist(s))
o_go<-as.matrix(substr(o_go,1,10),ncol=1)
go<-cbind(o_go,p_go)
go<-go[which(go[,2]!="all"),]
colnames(go)<-c("go","parent_go")
colnames(go_mf)<-c("id","go","evidence","ontology")
all_parent_gomf<-merge(go_mf,go)
tmp<-cbind(go_mf[,c(2,1,3,4)],go_mf[,2])
colnames(tmp)<-c("go","id","evidence","ontology","parent_go")
all_parent_gomf<-rbind(all_parent_gomf,tmp)



write.table(all_parent_gobp,file="Rn.gene2go_bp.txt",sep="\t",row.names=FALSE,col.names=TRUE)
write.table(all_parent_gocc,file="Rn.gene2go_cc.txt",sep="\t",row.names=FALSE,col.names=TRUE)
write.table(all_parent_gomf,file="Rn.gene2go_mf.txt",sep="\t",row.names=FALSE,col.names=TRUE)

