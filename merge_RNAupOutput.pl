#!/usr/bin/perl
$|=1;

############################################################################################################################
# extract alignment position of transcript and tRF from RNAup output
############################################################################################################################

@regions=();
@targets=();

# files with target sequences
@regions=
	(
	"3pUTR",
	"5pUTR",
	"CDS",
	);

# ENSG ID of 34 potential targets (significantly HEK down & HepG2 up)	
@targets=
	(
	"ENSG00000007202",
	"ENSG00000029363",
	"ENSG00000070061",
	"ENSG00000082438",
	"ENSG00000086598",
	"ENSG00000093000",
	"ENSG00000100196",
	"ENSG00000101310",
	"ENSG00000101311",
	"ENSG00000102226",
	"ENSG00000103342",
	"ENSG00000111846",
	"ENSG00000113448",
	"ENSG00000116044",
	"ENSG00000117394",
	"ENSG00000119927",
	"ENSG00000120694",
	"ENSG00000120705",
	"ENSG00000124383",
	"ENSG00000125148",
	"ENSG00000126709",
	"ENSG00000126878",
	"ENSG00000145703",
	"ENSG00000155097",
	"ENSG00000155629",
	"ENSG00000158435",
	"ENSG00000163399",
	"ENSG00000167861",
	"ENSG00000174684",
	"ENSG00000179918",
	"ENSG00000183741",
	"ENSG00000184307",
	"ENSG00000204569",
	"ENSG00000213341",
	);


$i=0;
@d=();
foreach$region(@regions)
	{
	open(OUT, ">RNAup_in.hsap-$region.genes.fas.tRF-5p-Glu-CTC.merge.out");
	print OUT"Gene_ID\ttranscript_start\ttranscript_end\ttRF_start\ttRF_end\talignment_length\n";
	
	foreach$id(@targets)
		{		
		# get output for RNAup alignment
		open(IN,"RNAup_in.hsap-$region.genes.fas.tRF-5p-Glu-CTC.$id.fa.out");
			while(<IN>) 
				{
				@d=split("\n",$_);
				foreach $line (@d) 
					{
					if ($line =~ /^\(\(/) 
						{
						@info=split ' ', $line;
						@pos_transcript=split ",", $info[1];
						$transcript_start=$pos_transcript[0];
						$transcript_end=$pos_transcript[1];
						@pos_tRF=split ",", $info[3];
						$tRF_start=$pos_tRF[0];
						$tRF_end=$pos_tRF[1];
						$alignment_length=$tRF_end-$tRF_start+1;
						print OUT"$id\t$transcript_start\t$transcript_end\t$tRF_start\t$tRF_end\t$alignment_length\n";
						@info=();
						@pos_transcript==();
						$transcript_start=0;
						$transcript_end=0;
						@pos_tRF=0;
						$tRF_start=0;
						$tRF_end=0;
						$alignment_length=0;
						}
					}
				@d=();
				}
		close IN;
		}
	close OUT;
	}