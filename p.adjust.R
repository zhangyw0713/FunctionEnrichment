#Rscript  --slave  --no-restore --no-save --no-init-file p.adjust.R indir outdir
#Rscript  --slave  --no-restore --no-save --no-init-file p.adjust.R tmp_enrich tmp_adjust
#perl get.goterm.pl tmp_adjust tmp_go
#Rscript  --slave  --no-restore --no-save --no-init-file p.adjust.R tmp_enrich tmp_adjust
#Rscript  --slave  --no-restore --no-save --no-init-file p.adjust.R hub/guo_co/hub_enrich hub/guo_co/hub_adjust 
#perl get.goterm.pl hub/guo_co/hub_adjust  hub/guo_co/hub_go


argv <- commandArgs(TRUE)
indir<-argv[1]
outdir<-argv[2]
files<-list.files(indir)

library(multtest)
for(i in 1:length(files))
{
	g<-read.table(paste(indir,"/",files[i],sep=""),sep="\t")
	p.adjust(as.numeric(g[,ncol(g)]),method="fdr")->q
	cbind(g,q)->gg
	write.table(gg,file=paste(outdir,"/",files[i],sep=""),sep="\t",row.names=FALSE,col.names=FALSE)
}



