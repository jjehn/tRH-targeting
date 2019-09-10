#!/usr/bin/perl
$|=1;

############################################################################################################################
# extract for each transcript region fasta-sequence of each potential target (34 genes; significantly HEK down & HepG2 up)
############################################################################################################################

@reference_files=();
@targets=();

# files with target sequences
@reference_files=
	(
	"hsap-3pUTR.genes.fas",
	"hsap-5pUTR.genes.fas",
	"hsap-CDS.genes.fas",
	);
	
# sequence of tRF-5p-Glu-CTC from 5' to 3'	
$tRF="TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT";
	
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
@seq_data=();
foreach$reference(@reference_files)
	{
	foreach$id(@targets)
		{
		# get sequence of potential target gene
		open(IN,"$reference");
				while(<IN>)
					{
					$_=~s/\s*$//;
					push(@seq_data,$_);
					$i+=0.5;
					if($i==1)
						{
						$i=0;
						for($seq_data[0] =~ /$id/)
							{
							# generate input file for RNAup (> tRF name \n tRF sequence \n > ENSG ID \n gene region sequence)
							open(TEMPIN,">RNAup_in.$reference.tRF-5p-Glu-CTC.$id.fa");
							print TEMPIN">tRF-5p-Glu-CTC\n$tRF\n$seq_data[0]\n$seq_data[1]\n";
							close TEMPIN;
							
							## run RNAup and remove temporary input file
							#system("RNAup -b -d2 --noLP -c 'S' < RNAup_in.$reference.tRF-5p-Glu-CTC.$id.fa > RNAup_$reference.tRF-5p-Glu-CTC.$id.out ");
							#unlink"RNAup_in.$reference.tRF-5p-Glu-CTC.$id.fa";
							}
						@seq_data=();
						}
					}
				close IN;	
		}
	}