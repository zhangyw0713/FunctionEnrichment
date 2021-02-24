open(go,"GO_level.txt");
while(<go>)
{
	chomp;
	chop;
	@token=split(/\t/,$_);
	$go{$token[0]}=$token[1];
	$level{$token[0]}=$token[3];
}
close(go);
open(go,"kegg2name");
while(<go>)
{
	chomp;
	$_=~s/"//g;

	@token=split(/\t/,$_);
	$go{$token[0]}=$token[1];
	$token[0]=~s/^0//g;
	$token[0]=~s/^0//g;
	$token[0]=~s/^0//g;
	$token[0]=~s/^0//g;
	$token[0]=~s/^0//g;
	$go{$token[0]}=$token[1];
	

}
close(go);
opendir(dir,"$ARGV[0]");
system("mkdir $ARGV[1]");
while(my $f=readdir(dir))
{
	open(one,"$ARGV[0]/$f");
	open(res,">$ARGV[1]/$f");

	print res "GO Terms	GO level	GO ID	Total Number of genes with this GO 	Number of genes with this GO in this dataset 	Number of genes in this dataset 	Total number of genes	Gene Names with this GO in the dataset	P-value	Adjust p-value\n" if(!grep(/kegg/,$f));
	print res "KEGG Pathway	No information	KEGG ID	Total Number of genes with this pathway 	Number of genes with this GO in this dataset 	Number of genes in this dataset 	Total number of genes	Gene Names with this pathway in the dataset	P-value	Adjust p-value\n" if(grep(/kegg/,$f));
	while(<one>)
	{
		chomp;
		$_=~s/"//g;
		@token=split(/\t/,$_);
	
		print res $go{$token[0]}."\t$level{$token[0]}\t$_\n";
	}
	close(one);
	close(res);
}
closedir(dir);
