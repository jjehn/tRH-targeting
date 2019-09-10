# map_kmers.pl maps all possible kmers of an input sequence
# to a collection of reference sequences. Nested alignments
# will be removed while retaining only the largest one.
# 
# DEPENDENCIES
# - seqmap executable in path. Download it here: http://www-personal.umich.edu/~jianghui/seqmap/#download
# - Perl modules that are not part of the core distribution: Parallel::ForkManager
#
# USAGE
# perl map_kmers.pl -ref reference.fasta -tRF name:sequence [-option value {...}]
#
use Parallel::ForkManager;
use Getopt::Long;

# default settings
$min_length=5; # minimum kmer length
$max_length=5; # maximum kmer length (restricted here for EJ-analysis)
$match_ends=2; # number of bp at both ends of the alignment without mismatch
$threads=1; # number of paralell threads

$mm1=8; # minimum length of a sequence to allow one mismatch
$mm2=14; # minimum length of a sequence to allow two mismatches
$mm3=20; # minimum length of a sequence to allow three mismatches
$mm4=26; # minimum length of a sequence to allow four mismatches
$mm5=32; # minimum length of a sequence to allow five mismatches
$insdel1=14; # minimum length of a sequence to allow one insertion/deletion
$insdel2=20; # minimum length of a sequence to allow two insertions/deletions
$insdel3=24; # minimum length of a sequence to allow three insertions/deletions

$|=1;

@reference_files=();
@tRFs=();
GetOptions
	(
	"min_length=i"=>\$min_length,
	"match_ends=i"=>\$match_ends,
	"threads=i"=>\$threads,
	"ref=s"=>\@reference_files,
	"trf=s"=>\@tRFs, # input format example: tRF-5p-Glu-CTC:TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT
	"mm1=i"=>\$mm1,
	"mm2=i"=>\$mm2,
	"mm3=i"=>\$mm3,
	"mm4=i"=>\$mm4,
	"mm5=i"=>\$mm5,
	"insdel1=i"=>\$insdel1,
	"insdel2=i"=>\$insdel2,
	"insdel3=i"=>\$insdel3,
	)||die print("\n\nError in command line arguments\n\n");

#~ # check if reference files exist
#~ if(@reference_files==0)
	#~ {
	#~ print"\n\nError. Provide at least one reference file with option -ref [filename].\n\n";exit;
	#~ }
#~ foreach$file(@reference_files)
	#~ {
	#~ if(!-e$file)
		#~ {
		#~ print"\n\nError. Cannot find file $file.\n\n";exit;
		#~ }
	#~ else
		#~ {
		#~ print"\nINPUT REFERENCE: $file";
		#~ }
	#~ }

#~ # check input format for tRFs
#~ if(@tRFs==0)
	#~ {
	#~ print"\n\nError. Provide at least one tRF name and sequence with option -trf [name:sequence].\n\n";exit;
	#~ }
#~ foreach$tRF(@tRFs)
	#~ {
	#~ @d=split(':',$$tRF);
	#~ if(@d!=2)
		#~ {
		#~ print"\n\nError. Format must be 'name:sequence'. E.g. tRF-1:TGCTTTGCTACGTCGAT\n\n";exit;
		#~ }
	#~ elsif($d[1]!~/^[ATGCatgc]$/)
		#~ {
		#~ print"\n\nError. Use [ATGCatgc] characters for tRF sequence.\n\n";exit;
		#~ }
	#~ elsif(length$d[0]==0||length$d[1]==0)
		#~ {
		#~ print"\n\nError. No sequence ('$d[1]') or sequence name ('$d[0']) provided.\n\n";exit;
		#~ }
	#~ $tRF{$d[0]}=$d[1];
	#~ print"\n$d[0] -> $d[1]";
	#~ }

# check mismatch settings
if($mm1>$mm2){print"\n\nError. value for mm1 must not be larger than value for mm2.\n\n";exit;}
elsif($mm2>$mm3){print"\n\nError. value for mm2 must not be larger than value for mm3.\n\n";exit;}
elsif($mm3>$mm4){print"\n\nError. value for mm3 must not be larger than value for mm4.\n\n";exit;}
elsif($mm4>$mm5){print"\n\nError. value for mm4 must not be larger than value for mm5.\n\n";exit;}
elsif($insdel1>$insdel2){print"\n\nError. value for insdel1 must not be larger than value for insdel2.\n\n";exit;}
elsif($insdel2>$insdel3){print"\n\nError. value for insdel2 must not be larger than value for insdel3.\n\n";exit;}


# files with target sequences
@reference_files=
	(
	# human
	"hsap-3pUTR.genes.fas",
	"hsap-3pUTR.genes.EE.fas",
	"hsap-3pUTR.genes.EJ.fas",
	"hsap-3pUTR.genes.ES.fas",
	
	"hsap-5pUTR.genes.fas",
	"hsap-5pUTR.genes.EE.fas",
	"hsap-5pUTR.genes.EJ.fas",
	"hsap-5pUTR.genes.ES.fas",
	
	"hsap-CDS.genes.fas",
	"hsap-CDS.genes.EE.fas",
	"hsap-CDS.genes.EJ.fas",
	"hsap-CDS.genes.ES.fas",
	);
	
# small RNA name and sequence from 5' to 3'	
%tRFs=
	(
	#human
	"tRF-5p-Glu-CTC"=>"TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT",
	#~ "tRF-5p-Gly-GCC"=>"GCATTGGTGGTTCAGTGGTAGAATTCTCGCCT",
	);

foreach$reference(@reference_files)
	{
	foreach$id(keys%tRFs)
		{
		# search for reverse complement (target sites)
		$tRF=$tRFs{$id};
		$tRF=reverse$tRF;
		$tRF=~tr/ATGC/TACG/;
		
		$pm=Parallel::ForkManager->new($threads);
		foreach$l($min_length..$max_length) # restricted due to EJ-analysis (instead of length$tRFs{$id})
			{
			$pm->start and next;
			$p=-1;
			while(1)
				{
				$p++;
				$kmer=substr($tRF,$p,$l);
				last if(length$kmer<$l);
				
				$mm=0;
				$insdel="";
				
				if($l>=$mm1){$mm=1;}
				if($l>=$mm2){$mm=2;}
				if($l>=$mm3){$mm=3;}
				if($l>=$mm4){$mm=4;}
				if($l>=$mm5){$mm=5;}
					
				if($l>=$insdel1){$insdel="/allow_insdel:1";}
				if($l>=$insdel2){$insdel="/allow_insdel:2";}
				if($l>=$insdel3){$insdel="/allow_insdel:3";}
				
				# create temporary input file
				$p_out=(length$tRF)-$p-$l;
				
				open(TEMPIN,">seqmap_in.$reference.$id.$p_out.$l.temp");
				print TEMPIN">$id:$p_out:$l\n$kmer\n";
				close TEMPIN;
				
				# run SeqMap and remove temporary input file
				$out_id++;
				system("seqmap $mm seqmap_in.$reference.$id.$p_out.$l.temp $reference $reference.$id.$p_out.$l.seqmap /forward_strand /available_memory:16000 /output_all_matches $insdel ");
				unlink"seqmap_in.$reference.$id.$p_out.$l.temp";
				}
			$pm->finish;
			}
		$pm->wait_all_children;
		
		# merge SeqMap output files
		$i=0;
		open(MERGE,">$reference.$id.seqmap.merge")||die print$!;
		foreach$l($min_length..length$tRF)
			{
			foreach$p_out(0..length$tRF)
				{
				last if(!-e"$reference.$id.$p_out.$l.seqmap");
				open(IN,"$reference.$id.$p_out.$l.seqmap");
				$head=<IN>;
				while(<IN>)
					{
					@d=split("\t",$_);
					
					# do not allow mismatch in first or last n positions of the alignment
					next if(substr($d[2],0,$match_ends) ne substr($d[4],0,$match_ends));
					next if(substr($d[2],-$match_ends) ne substr($d[4],-$match_ends));
					print MERGE$_;
					}
				close IN;
				unlink"$reference.$id.$p_out.$l.seqmap";
				}
			}
		close MERGE;
		
		# remove nested target sites, take only largest one
		@content=();
		open(INPUT,"$reference.$id.seqmap.merge")||die print$!;
		while(<INPUT>)
			{
			@d=split("\t",$_);
			$d[3]=~s/^[^:]+/x/;
			$l1=length$d[2];
			$l2=length$d[4];
			push(@content,"$d[0]\t$d[1]\t$l1");
			}
		close INPUT;
		
		$i=@content;
		%ok=();
		%sites=();
		while(@content>0)
			{
			$element=pop@content;
			$i=$i-1;
			@d=split("\t",$element);
			
			$ok=0;
			foreach$p($d[1]..$d[1]+$d[2])
				{
				if(!$sites{$d[0]}{$p})
					{
					foreach$p($d[1]..$d[1]+$d[2])
						{
						$sites{$d[0]}{$p}=1;
						}
					$ok{$i}=1;
					last;
					}
				}
			}
		
		open(INPUT,"$reference.$id.seqmap.merge")||die print$!;
		open(OUTPUT,">$reference.$id.seqmap.merge.clean")||die print$!;
		$i=-1;
		while(<INPUT>)
			{
			$i++;
			if($ok{$i}==1)
				{
				print OUTPUT$_;
				}
			}
		close INPUT;
		close OUTPUT;
		unlink"$reference.$id.seqmap.merge";
		undef@content;
		undef%ok;
		undef%sites;
		}
	}
