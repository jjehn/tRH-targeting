$spec="dmel";
$geneID_format="FBgn";
$transcriptID_format="FBtr";
$genome="BDGP6.genome";
$gtf="BDGP6.94.gtf";
$|=1;

$skip=0;
open(GENOME,$genome)||die print$!;
while(<GENOME>)
	{
	$_=~s/\s*$//;
	if($_=~s/^>//)
		{
		$chr=$_;
		print"\nLoad $chr...";
		}
	else
		{
		$genome{$chr}.=uc$_;
		}
	}
close GENOME;

open(GTF,$gtf)||die print$!;
open(UTR_3,">$spec-3pUTR.transcripts.fas")||die print$!;
open(UTR_5,">$spec-5pUTR.transcripts.fas")||die print$!;
open(CDS,">$spec-CDS.transcripts.fas")||die print$!;
$prev3ID="";
$prev5ID="";
$prevCDSID="";
while(<GTF>)
	{
	@d=split("\t",$_);
	if($d[2]eq'three_prime_utr')
		{
		get();
		if($prev3ID ne $trID)
			{
			print UTR_3"\n>$trID\n$seq";
			}
		else
			{
			print UTR_3"$seq";
			}
		$prev3ID=$trID;
		}
	elsif($d[2]eq'five_prime_utr')
		{
		get();
		if($prev5ID ne $trID)
			{
			print UTR_5"\n>$trID\n$seq";
			}
		else
			{
			print UTR_5"$seq";
			}
		$prev5ID=$trID;
		}
	elsif($d[2]eq'CDS')
		{
		get();
		if($prevCDSID ne $trID)
			{
			print CDS"\n>$trID\n$seq";
			}
		else
			{
			print CDS"$seq";
			}
		$prevCDSID=$trID;
		}
	}
close GTF;
close UTR_3;
close UTR_5;
close CDS;

# write transcripts with maximum length for each gene to a new file
%max=();
foreach$gene_ID(keys%tr_length)
	{
	$max=0;
	foreach$trID(keys%{$tr_length{$gene_ID}})
		{
		if(!$max{$gene_ID}||$tr_length{$gene_ID}{$trID}>$tr_length{$gene_ID}{$max{$gene_ID}})
			{
			$max{$gene_ID}=$trID;
			}
		}
	}
%tr_2_gene=();
foreach$gene_ID(keys%max)
	{
	$tr_2_gene{$max{$gene_ID}}=$gene_ID;
	}

$max=0;
foreach$file("3pUTR","5pUTR","CDS")
	{
	open(IN,"$spec-$file.transcripts.fas")||die print$!;
	open(OUT,">$spec-$file.genes.fas")||die print$!;
	while(<IN>)
		{
		$_=~s/\s*$//;
		if($_=~s/^>//)
			{
			if($tr_2_gene{$_})
				{
				$max=1;
				$gene_id=$tr_2_gene{$_};
				}
			else
				{
				$max=0;
				}
			}
		elsif($max==1)
			{
			print OUT">$gene_id\n$_\n";
			}
		}
	close IN;
	close OUT;
	}

sub get
	{
	$d[8]=~s/$transcriptID_format\d+//;
	$trID=$&;
	$d[8]=~s/$geneID_format\d+//;
	$geneID=$&;
	$seq=substr($genome{$d[0]},$d[3]-1,($d[4]-$d[3])+1);
	if($d[6]eq'-')
		{
		$seq=reverse$seq;
		$seq=~tr/ATGC/TACG/;
		}
	$tr_length{$geneID}{$trID}+=length$seq;
	}
