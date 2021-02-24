use Math::GSL::CDF qw /:all/;
#/home/zhengqi/coexpression/tmp 

$gofile="/home/data/GO/Rn.gene2go_bp.txt" if ($ARGV[3]==3);
$gofile="/home/data/GO/Mm.gene2go_bp.txt" if ($ARGV[3]==2);
$gofile="/home/data/GO/hs.gene2go_bp.txt" if ($ARGV[3]==1);

open(gos,"$gofile"); #Hs Rn Mm


while(<gos>)
{
	chomp;
#	chop;
	$_=~s/"//g;
	@token=split(/\t/,$_);
	
	
	$upgo=$token[1]."_up";
	$allgenes{$token[1]}=0;
	${$upgo}{$token[0]}=0;
	${$upgo}{$token[4]}=0;
	$gonum=$token[0]."_go";
	${$gonum}{$token[1]}=0;
	$gonum=$token[4]."_go";
	${$gonum}{$token[1]}=0;

}
close(gos);
$orgfile="/home/data/Rattus_norvegicus.gene_info" if ($ARGV[3]==3);
$orgfile="/home/data/Mus_musculus.gene_info" if ($ARGV[3]==2);
$orgfile="/home/data/Homo_sapiens.gene_info" if ($ARGV[3]==1);
open(gee,"$orgfile"); #Rattus_norvegicus  Homo_sapiens Mus_musculus.gene_info.
while(<gee>)
{
	$_="\U$_";
	@token=split(/\t/,$_);
	$allg{$token[1]}=0;
	$ginfo{$token[1]}=$token[2];
	$token[2]="\U$token[2]";
	$n2gid{$token[2]}=$token[1];
	$gtype{$token[2]}=$token[9]; 

#	print $token[9]."\n";
}
close(gee);



#$target_dir="target_filter_by_dhs";
#$enrich_dir="enrich_filter_by_dhs";
$target_dir="$ARGV[0]"; #INPUT directory
$enrich_dir="$ARGV[1]"; # output diectory
#system("mkdir $enrich_dir");

opendir(tar,"$target_dir");
while(my $file=readdir(tar))
{
	if(!grep(/^\./,$file))
	{
%mygene=();
%mygonum=();
%info=();
%gog=();
	open(gene,"$target_dir/$file");
while(<gene>)
{
	chomp;
	print $_."\n";
	chop if($ARGV[2]==1);
	print $_."\n";
	$_=~s/"//g;
	$_="\U$_";
	@token=split(/\t/,$_);
#	if(grep(/_coding/,$token[2]))
#	{
		
		$name=$token[0];
		$name=~s/_coding//g;
			
		if(exists($n2gid{$name}))
		{#	print "aa\n$name\n";
		
			$t=$n2gid{$name};
			if(exists($allgenes{$t}))
			{
			
			$mygene{$t}=0  ;
			$upgo=$t."_up";
			@gos=keys(%{$upgo});
			
			foreach my $g(@gos)
			{
				#	print $g;
				$mygonum{$g}++;
				$gog{$g}=$gog{$g}.";".$name;
			}
		}
		}
		
		
#	}
}
close(gene);

	open(res,">$enrich_dir/$file");


$n=keys(%mygene);
$total_gene=keys(%allgenes);
while(my($k,$v)=each(%mygonum))
{
	$gonum=$k."_go";
	$has_go_gene_num=keys(%{$gonum});
	$in_num=$v;
	$gogs=substr($gog{$k},1);
	$p= gsl_cdf_hypergeometric_Q($v, $has_go_gene_num,$total_gene-$has_go_gene_num, $n); 
	$info{$k."\t".$has_go_gene_num."\t$in_num\t$n\t$total_gene\t$gogs\t$p\n"}=$p if($p!=0 and $in_num>=3);


}

@sorted = map { { ($_ => $info{$_}) } }
        sort { $info{$a} <=> $info{$b}
              or $a <=> $b
            } keys %info;
foreach $hashref (@sorted) {
  ($key, $value) = each %$hashref;
  print res  "$key";
}
close(res);

}
}


