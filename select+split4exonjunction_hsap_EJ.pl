$spec="hsap";
$geneID_format="ENSG";
$transcriptID_format="ENST";
$genome="GRCh38.genome";
$gtf="GRCh38.94.gtf";
$mask_size=8; # needs to be an even number
$|=1;

$skip=0;
open(GENOME,$genome)||die print$!;
while(<GENOME>)
	{
	$_=~s/\s*$//;
	if($_=~s/^.+chromosome //)
		{
		if($_=~/unlocalized/||$_=~/genomic patch/||$_=~/alternate locus/)
			{
			$skip=1;
			print"\nSkip $_...";
			}
		else
			{
			$skip=0;
			$_=~s/\,.+$//;
			$chr=$_;
			print"\nLoad $chr...";
			}
		}
	elsif($skip==0)
		{
		$genome{$chr}.=uc$_;
		}
	}
close GENOME;

open(GTF,$gtf)||die print$!;

# REGULAR OUTPUT
open(UTR_3,">$spec-3pUTR.transcripts.fas")||die print$!;
open(UTR_5,">$spec-5pUTR.transcripts.fas")||die print$!;
open(CDS,">$spec-CDS.transcripts.fas")||die print$!;

# ONLY EXON JUNCTIONS
open(UTR_3_EJ,">$spec-3pUTR.transcripts.EJ.fas")||die print$!;
open(UTR_5_EJ,">$spec-5pUTR.transcripts.EJ.fas")||die print$!;
open(CDS_EJ,">$spec-CDS.transcripts.EJ.fas")||die print$!;

# EXON STARTS ONLY
open(UTR_3_ES,">$spec-3pUTR.transcripts.ES.fas")||die print$!;
open(UTR_5_ES,">$spec-5pUTR.transcripts.ES.fas")||die print$!;
open(CDS_ES,">$spec-CDS.transcripts.ES.fas")||die print$!;

# EXON ENDS ONLY
open(UTR_3_EE,">$spec-3pUTR.transcripts.EE.fas")||die print$!;
open(UTR_5_EE,">$spec-5pUTR.transcripts.EE.fas")||die print$!;
open(CDS_EE,">$spec-CDS.transcripts.EE.fas")||die print$!;


$prev3ID="";
$prev5ID="";
$prevCDSID="";
while(<GTF>)
	{
	@d=split("\t",$_);
	if($d[2]eq'three_prime_utr')
		{
		get();
		$seq_EJ=substr($seq,0,$mask_size/2). "N" x ((length$seq)-$mask_size) . substr($seq,-($mask_size/2));
		$seq_ES=substr($seq,0,$mask_size). "N" x ((length$seq)-$mask_size);
		$seq_EE="N" x ((length$seq)-$mask_size) . substr($seq,-$mask_size);
		if($prev3ID ne $trID)
			{
			print UTR_3"\n>$trID\n$seq";
			print UTR_3_EJ"\n>$trID\n$seq_EJ";
			print UTR_3_ES"\n>$trID\n$seq_ES";
			print UTR_3_EE"\n>$trID\n$seq_EE";
			}
		else
			{
			print UTR_3"$seq";
			print UTR_3_EJ"$seq_EJ";
			print UTR_3_ES"$seq_ES";
			print UTR_3_EE"$seq_EE";
			}
		$prev3ID=$trID;
		}
	elsif($d[2]eq'five_prime_utr')
		{
		get();
		if($prev5ID ne $trID)
			{
			print UTR_5"\n>$trID\n$seq";
			print UTR_5_EJ"\n>$trID\n$seq_EJ";
			print UTR_5_ES"\n>$trID\n$seq_ES";
			print UTR_5_EE"\n>$trID\n$seq_EE";
			}
		else
			{
			print UTR_5"$seq";
			print UTR_5_EJ"$seq_EJ";
			print UTR_5_ES"$seq_ES";
			print UTR_5_EE"$seq_EE";
			}
		$prev5ID=$trID;
		}
	elsif($d[2]eq'CDS')
		{
		get();
		if($prevCDSID ne $trID)
			{
			print CDS"\n>$trID\n$seq";
			print CDS_EJ"\n>$trID\n$seq_EJ";
			print CDS_ES"\n>$trID\n$seq_ES";
			print CDS_EE"\n>$trID\n$seq_EE";
			}
		else
			{
			print CDS"$seq";
			print CDS_EJ"$seq_EJ";
			print CDS_ES"$seq_ES";
			print CDS_EE"$seq_EE";
			}
		$prevCDSID=$trID;
		}
	}
close GTF;
close UTR_3;close UTR_3_EJ;close UTR_3_ES;close UTR_3_EE;
close UTR_5;close UTR_5_EJ;close UTR_5_ES;close UTR_5_EE;
close CDS;close CDS_EJ;close CDS_ES;close CDS_EE;

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
	
	open(IN_EJ,"$spec-$file.transcripts.EJ.fas")||die print$!;
	open(OUT_EJ,">$spec-$file.genes.EJ.fas")||die print$!;
	
	open(IN_ES,"$spec-$file.transcripts.ES.fas")||die print$!;
	open(OUT_ES,">$spec-$file.genes.ES.fas")||die print$!;
	
	open(IN_EE,"$spec-$file.transcripts.EE.fas")||die print$!;
	open(OUT_EE,">$spec-$file.genes.EE.fas")||die print$!;
	
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
	
	while(<IN_EJ>)
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
			print OUT_EJ">$gene_id\n$_\n";
			}
		}
	close IN_EJ;
	close OUT_EJ;
	
	while(<IN_ES>)
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
			print OUT_ES">$gene_id\n$_\n";
			}
		}
	close IN_ES;
	close OUT_ES;
	
	while(<IN_EE>)
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
			print OUT_EE">$gene_id\n$_\n";
			}
		}
	close IN_EE;
	close OUT_EE;
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
