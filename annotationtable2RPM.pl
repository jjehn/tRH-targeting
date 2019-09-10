#!/usr/bin/perl -w

use strict;
use warnings;
#Initiate all variables, hashes and co


############### RPM ##########################
#Open folders in working directy

my @folders = glob("*"); #to get all folders in directory; extension ("*") as wildcard to get all names
foreach my $folder(@folders) #to speak to each element in directory
	{
	next if ($folder!~/^UNITAS_/); #skip elements which do not start with "UNITAS"
	opendir(DIR,$folder)||die print$!; #open folder, end script when opening is not possible (DIR is the "filehandle" for the directory)
	print"\n$folder";
	while( my $file=readdir(DIR)) #returns content of folder
		{
		next if($file!~/\.mapped_sequences$/); #get the mapped_sequences file we need to read out the reads
		print"\n$file"; #print out file names to make sure we get the right files

		my $reads = 0; #set the number of reads to 0 for each run

		open my $fileone, '<', "$folder/$file" or die "$folder/$file: $!";
		
		while(my $tocount=<$fileone>)#read file
			{
			chomp $tocount; #chomp to remove newlines at the end 
			$tocount =~ s/>//g; #remove all ">"
			next if ($tocount =~ /[A-Za-z]/); #skip lines which contain the sequence
					
			if ($tocount =~ /[0-9]/) #get the read-number
				{
				print"\n$tocount";
				$reads = ($reads + $tocount); # add up all reads
				}
			print"\nReads = $reads";

			}
		close $fileone;

		my %hash; #initiate empty hash

		my $trftable = 'unitas.tRF-table.txt'; #save file name in variable

		open my $trf, '<', "$folder/$trftable" or die "$folder/$trftable: $!";

		<$trf> for 1 .. 4;

		while( my $line=<$trf>)
			{
			chomp $line;
			my @line = split("\t",$line);
		
			my $tRNAname;

			if($line[0] =~ s/tRNA-[^-]+-...//) # "tRNA-"(matched tRNA und -) "[^-]+" beginning bis Ende, egal was "-..."(weiterer Strich bis Ende)
				{
				$tRNAname = $&; # "$&" = last pattern match
				print"\n$tRNAname";
				}
			else
				{
				$tRNAname = $line[0];
				$tRNAname =~ s/-ENS.+$//; # "-ENS.+$" ( matched allen die -ENS. bis Ende enthalten)
				print"\n$tRNAname";
				}
			my $hash = ($hash{$tRNAname} //= {} );
			$hash{$tRNAname}{"5p-tR-halves"} += $line[1]/$reads*1000000;
			$hash{$tRNAname}{"5p-tRFs"} += $line[3]/$reads*1000000;
			$hash{$tRNAname}{"3p-tR-halves"} += $line[5]/$reads*1000000;
			$hash{$tRNAname}{"3p-CCA-tRFs"} += $line[7]/$reads*1000000;
			$hash{$tRNAname}{"3p-tRFs"} += $line[9]/$reads*1000000;
			$hash{$tRNAname}{"tRF-1"} += $line[11]/$reads*1000000;
			$hash{$tRNAname}{"tRNA-leader"} += $line[13]/$reads*1000000;
			$hash{$tRNAname}{"misc-tRFs"} += $line[15]/$reads*1000000;
			}

		close $trf;

		open my $mergerpm,">>","$folder/mergerpm" or die "Could not open $folder/mergerpm : $!";
		
		my @tRF_types=("5p-tR-halves","5p-tRFs","3p-tR-halves","3p-CCA-tRFs","3p-tRFs","tRF-1","tRNA-leader","misc-tRFs");
		for my $tRNAname(sort{$a cmp $b}keys %hash) #sortiert die alphabetisch nach keys
			{
			print $mergerpm $tRNAname; # print tRNA name
			foreach my $tRF_type(@tRF_types) 
				{
				print $mergerpm"\t$hash{$tRNAname}{$tRF_type}"; # print counts for each tRF type separated by tab
				}
			print $mergerpm"\n";# print newline
			}
		close $mergerpm;
		}
	close DIR;
	}

############### RPT ##########################

my @folders = glob("*"); #to get all folders in directory; extension ("*") as wildcard to get all names
foreach my $folder(@folders) #to speak to each element in directory
	{
	next if ($folder!~/^UNITAS_/); #skip elements which do not start with "UNITAS"
	opendir(DIR,$folder)||die print$!; #open folder, end script when opening is not possible (DIR is the "filehandle" for the directory)
	print"\n$folder";
	while( my $file=readdir(DIR)) #returns content of folder
		{
		next if($file!~/\.mapped_sequences$/); #get the mapped_sequences file we need to read out the reads
		print"\n$file"; #print out file names to make sure we get the right files

		my $filename = 'unitas.annotation_summary.txt';

		my $sum = 0;

		open my $filetwo, '<', "$filename" or die "$filename: $!";
			while ( my $line=<$filetwo>)
			{
			if ($. == 8)
				{
				$sum = (split ' ', $line)[1];
				}
			print"\ntRFSum = $sum";
			}
		close $filetwo;

		my %hash; #initiate empty hash

		my $trftable = 'unitas.tRF-table.txt'; #save file name in variable

		open my $trf, '<', "$folder/$trftable" or die "$folder/$trftable: $!";

		<$trf> for 1 .. 4;

		while( my $line=<$trf>)
			{
			chomp $line;
			my @line = split("\t",$line);
		
			my $tRNAname;

			if($line[0] =~ s/tRNA-[^-]+-...//) # "tRNA-"(matched tRNA und -) "[^-]+" beginning bis Ende, egal was "-..."(weiterer Strich bis Ende)
				{
				$tRNAname = $&; # "$&" = last pattern match
				print"\n$tRNAname";
				}
			else
				{
				$tRNAname = $line[0];
				$tRNAname =~ s/-ENS.+$//; # "-ENS.+$" ( matched allen die -ENS. bis Ende enthalten)
				print"\n$tRNAname";
				}
			my $hash = ($hash{$tRNAname} //= {} );
			$hash{$tRNAname}{"5p-tR-halves"} += $line[1]/$sum*1000000;
			$hash{$tRNAname}{"5p-tRFs"} += $line[3]/$sum*1000000;
			$hash{$tRNAname}{"3p-tR-halves"} += $line[5]/$sum*1000000;
			$hash{$tRNAname}{"3p-CCA-tRFs"} += $line[7]/$sum*1000000;
			$hash{$tRNAname}{"3p-tRFs"} += $line[9]/$sum*1000000;
			$hash{$tRNAname}{"tRF-1"} += $line[11]/$sum*1000000;
			$hash{$tRNAname}{"tRNA-leader"} += $line[13]/$sum*1000000;
			$hash{$tRNAname}{"misc-tRFs"} += $line[15]/$sum*1000000;
			}

		close $trf;

		open my $mergerpt,">>","$folder/mergerpt" or die "Could not open $folder/mergerpt : $!";
		
		my @tRF_types=("5p-tR-halves","5p-tRFs","3p-tR-halves","3p-CCA-tRFs","3p-tRFs","tRF-1","tRNA-leader","misc-tRFs");
		for my $tRNAname(sort{$a cmp $b}keys %hash) #sortiert die alphabetisch nach keys
			{
			print $mergerpt $tRNAname; # print tRNA name
			foreach my $tRF_type(@tRF_types) 
				{
				print $mergerpt"\t$hash{$tRNAname}{$tRF_type}"; # print counts for each tRF type separated by tab
				}
			print $mergerpt"\n";# print newline
			}
		close $mergerpt;
		}
	close DIR;
	}
