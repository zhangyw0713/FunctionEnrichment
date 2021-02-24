args=commandArgs(T)
datapath=args[1]
inputpath=args[2]
outpath=args[3]

setwd(datapath)
files<-dir(datapath)
grep("HALLMARK",files,value=T)->hallfiles
if(file.info(inputpath)$size>0)
{
	read.table(inputpath,sep="\t",header=F)->qrna
	if(nrow(qrna)>0)
	{
		qrna[2:nrow(qrna),1]->query
		length(query)->a
		matrix(nrow=length(hallfiles),ncol=3)->res
		dd=1
		for(i in hallfiles)
		{
			read.table(i,sep="\t",header=F)->glist
			glist[3:nrow(glist),1]->glist
			length(glist)->b
			length(intersect(query,glist))->c
			N=20377
			phyper(c,a,N-a,b,lower.tail=FALSE)->pvalue
			i->res[dd,1]
			pvalue->res[dd,2]
			dd=dd+1
		}

		p.adjust(as.numeric(as.character(res[,2])),method="fdr",n=nrow(res))->res[,3]
		# outpath=paste(outpath,"_MsigDB.tsv",sep="")
		colnames(res)<-c("Hallmarks","p-value","FDR")
		write.table(res,outpath,sep="\t",row.names=F,col.names=T,quote=F)
	}
}
