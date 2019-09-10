# tRF-5p-Glu-CTC	TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT	l=33
# tRF-5p-Gly-GCC	GCATTGGTGGTTCAGTGGTAGAATTCTCGCCT	l=32
# T3-AspGUC-4-23	TCGATAGTATAGTGGTTAGT				l=20
# T6-GluCUC-6-25	TATTGTCTAGTGGTTAGGAT				l=20
# T10-LysUUU-4-23	CGGATAGCTCAGTCGGTAGA				l=20
$|=1;

$DESeq2result="res_expressed_Glu.txt"; # Result file of DESeq2 (depleted by NA values; filtered for fpkm >1). Has a header. GeneID: column 0 | log2FC: column 2 | p_adj: column 6


@kmer_maps=
	(
	# human HEK293 Glu-CTC
	"map_kmers_strict/hsap-3pUTR.genes.fas.tRF-5p-Glu-CTC.seqmap.merge.clean",
	#~ "map_kmers_strict/hsap-5pUTR.genes.fas.tRF-5p-Glu-CTC.seqmap.merge.clean",
	#~ "map_kmers_strict/hsap-CDS.genes.fas.tRF-5p-Glu-CTC.seqmap.merge.clean",
	
	# human HEK293 Gly-GCC
	#~ "map_kmers_strict/hsap-3pUTR.genes.fas.tRF-5p-Gly-GCC.seqmap.merge.clean",
	#~ "map_kmers_strict/hsap-5pUTR.genes.fas.tRF-5p-Gly-GCC.seqmap.merge.clean",
	#~ "map_kmers_strict/hsap-CDS.genes.fas.tRF-5p-Gly-GCC.seqmap.merge.clean",
	
	# fly S2 T3
	#~ "8-14-20-26-32.14-20-24.strict/dmel-3pUTR.genes.fas.T3-AspGUC-4-23.seqmap.merge.clean",
	#~ "8-14-20-26-32.14-20-24.strict/dmel-5pUTR.genes.fas.T3-AspGUC-4-23.seqmap.merge.clean",
	#~ "8-14-20-26-32.14-20-24.strict/dmel-CDS.genes.fas.T3-AspGUC-4-23.seqmap.merge.clean",
	
	# fly S2 T6
	#~ "8-14-20-26-32.14-20-24.strict/dmel-3pUTR.genes.fas.T6-GluCUC-6-25.seqmap.merge.clean",
	#~ "8-14-20-26-32.14-20-24.strict/dmel-5pUTR.genes.fas.T6-GluCUC-6-25.seqmap.merge.clean",
	#~ "8-14-20-26-32.14-20-24.strict/dmel-CDS.genes.fas.T6-GluCUC-6-25.seqmap.merge.clean",
	
	# fly S2 T10
	#~ "8-14-20-26-32.14-20-24.strict/dmel-3pUTR.genes.fas.T10-LysUUU-4-23.seqmap.merge.clean",
	#~ "8-14-20-26-32.14-20-24.strict/dmel-5pUTR.genes.fas.T10-LysUUU-4-23.seqmap.merge.clean",
	#~ "8-14-20-26-32.14-20-24.strict/dmel-CDS.genes.fas.T10-LysUUU-4-23.seqmap.merge.clean",
	);

$results_file="DESeq2result_tRF-5p-Glu-CTC_3pUTR.strict.txt";

$min_target_site=5;
$seq_length=33;


$up=0;
$down=0;
$unchanged=0;
%up=();
%down=();
%unchanged=();
%expressed=();
print"\nLoad DESeq2 result table...";
open(DE,$DESeq2result)||die print$!;
$head=<DE>;
while(<DE>)
	{
	$_=~s/\s*$//;
	@de=split("\t",$_);
	
	$log2FC = $de[2];
	$p_adj = $de[6];
	$totalgenes{$de[0]}=1;
	
	if($p_adj>0.01) # same TH like volcano plot
		{
		$unchanged++; # meaning in this context not significant
		$unchanged{$de[0]}=1;
		}
	elsif($log2FC>0) #up
		{
		$up++;
		$up{$de[0]}=1;
		}
	elsif($log2FC<0) # down
		{
		$down++;
		$down{$de[0]}=1;
		}
	}
close DE;
$totalgenes=keys%totalgenes;
print" done.
Genes in DESeq2 result file: $totalgenes
With threshold for significant differential expression set to adjusted p-value < 0.01:
up-regulated:   $up
down-regulated: $down
unchanged:      $unchanged
";

if($up==0||$down==0)
	{
	print"\nFurther analysis makes no sense with $up up- and $down down-regulated genes.\n\n";exit;
	}

%targets=(); # key1:length key2:start key3:up/down value:number_of_genes
foreach$kmer_map(@kmer_maps)
	{
	print"\nRead kmer map file $kmer_map...";
	open(MAP,$kmer_map)||die print$!;
	while(<MAP>)
		{
		@map=split("\t",$_);
		@d=split(':',$map[3]);
		$start=$seq_length-$d[1]-$d[2]+1;
		
		# check for canonical miRNA target site (2-8, max 1 G:U whobble)
		if($start==2&&$d[2]>=7)
			{
			$seed_ok=1;
			$GU=0;
			foreach$p(1..7)
				{
				if(substr($map[2],$p,1) ne substr($map[4],$p,1))
					{
					if(substr($map[2],$p,1).substr($map[4],$p,1)eq'TC'||substr($map[2],$p,1).substr($map[4],$p,1)eq'GA')
						{
						$GU++;
						if($GU>1)
							{
							$seed_ok=0;
							last;
							}
						}
					else
						{
						$seed_ok=0;
						last;
						}
					}
				}
			if($seed_ok==1)
				{
				# this is a canonical miRNA target site
				if($up{$map[0]}==1)
					{
					$miR_up{$map[0]}++;
					}
				elsif($down{$map[0]}==1)
					{
					$miR_down{$map[0]}++;
					}
				elsif($unchanged{$map[0]}==1)
					{
					$miR_unchanged{$map[0]}++;
					}
				}
			}
		
		if($up{$map[0]}==1)
			{
			$targets{$d[2]}{$start}{1}{$map[0]}++;
			}
		elsif($down{$map[0]}==1)
			{
			$targets{$d[2]}{$start}{2}{$map[0]}++;
			}
		elsif($unchanged{$map[0]}==1)
			{
			$targets{$d[2]}{$start}{0}{$map[0]}++;
			}
		}
	close MAP;
	print" done.";
	}

print"\nSave results...";
open(RESULTS,">$results_file")||die print$!;
print RESULTS"length\tstart\ttargeted\tup\tdown\tunchanged\tup/down\tup/total\tdown/total\tnot-targeted\tup\tdown\tunchanged\tup/down\tup/total\tdown/total\tfraction-targeted\n";
foreach$l($min_target_site..$seq_length)
	{
	foreach$s(1..($seq_length-$min_target_site)+1)
		{
		$targeted_up=keys%{$targets{$l}{$s}{1}};
		$targeted_down=keys%{$targets{$l}{$s}{2}};
		$targeted_unchanged=keys%{$targets{$l}{$s}{0}};
		$targeted=$targeted_up+$targeted_down+$targeted_unchanged;
		$up_vs_down="";
		$up_vs_total="";
		$down_vs_total="";
		if($targeted_down>0){$up_vs_down=$targeted_up/$targeted_down;}
		if($targeted>0){$up_vs_total=$targeted_up/$targeted;$down_vs_total=$targeted_down/$targeted;}
		
		$nt_up=$up-$targeted_up;
		$nt_down=$down-$targeted_down;
		$nt_unchanged=($totalgenes-$up-$down)-$targeted_unchanged;
		$nt=$totalgenes-$targeted;
		$nt_up_vs_down="";
		$nt_up_vs_total="";
		$nt_down_vs_total="";
		if($nt_down>0){$nt_up_vs_down=$nt_up/$nt_down;}
		if($nt>0){$nt_up_vs_total=$nt_up/$nt;$nt_down_vs_total=$nt_down/$nt;}
		
		$fraction_targeted=$targeted/$totalgenes;
		
		print RESULTS"$l\t$s\t$targeted\t$targeted_up\t$targeted_down\t$targeted_unchanged\t$up_vs_down\t$up_vs_total\t$down_vs_total\t$nt\t$nt_up\t$nt_down\t$nt_unchanged\t$nt_up_vs_down\t$nt_up_vs_total\t$nt_down_vs_total\t$fraction_targeted\n";
		
		last if($s+$l>$seq_length);
		}
	}

# output for canonical miRNA target sites
$targeted_up=keys%miR_up;
$targeted_down=keys%miR_down;
$targeted_unchanged=keys%miR_unchanged;
$targeted=$targeted_up+$targeted_down+$targeted_unchanged;
$up_vs_down="";
$up_vs_total="";
$down_vs_total="";
if($targeted_down>0){$up_vs_down=$targeted_up/$targeted_down;}
if($targeted>0){$up_vs_total=$targeted_up/$targeted;$down_vs_total=$targeted_down/$targeted;}

$nt_up=$up-$targeted_up;
$nt_down=$down-$targeted_down;
$nt_unchanged=($totalgenes-$up-$down)-$targeted_unchanged;
$nt=$totalgenes-$targeted;
$nt_up_vs_down="";
$nt_up_vs_total="";
$nt_down_vs_total="";
if($nt_down>0){$nt_up_vs_down=$nt_up/$nt_down;}
if($nt>0){$nt_up_vs_total=$nt_up/$nt;$nt_down_vs_total=$nt_down/$nt;}

$fraction_targeted=$targeted/$totalgenes;
print RESULTS"\n8-miR\t2-miR\t$targeted\t$targeted_up\t$targeted_down\t$targeted_unchanged\t$up_vs_down\t$up_vs_total\t$down_vs_total\t$nt\t$nt_up\t$nt_down\t$nt_unchanged\t$nt_up_vs_down\t$nt_up_vs_total\t$nt_down_vs_total\t$fraction_targeted";


close RESULTS;
print" done.\n\n";
