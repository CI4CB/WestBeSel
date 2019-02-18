#! /usr/bin/perl -w

#####################################################################
## Check job status and read the sequence of protein, then do slide window and save to file
## Outfile: filename_prediction.txt
## Author: BO YAO   Jan 6 2012
## The correct command: perl SVMTriP.pl fasta_filename epilength
## Here epilength is only select from 10, 12, 14, 16, 18, 20
## fasta file only contains one sequence. More than one fasta sequence will only be dealt with the first one.

## Example: perl SVMTriP.pl data.fasta 20
#####################################################################

use strict;
use warnings;
use DBI;

my $filename = $ARGV[0];
my $epilength = $ARGV[1];



if (!(defined $filename) || !(defined $epilength))
{
	die ("Error: command line is wrong. The current commnad line is perl SVMTriP.pl fasta_filename epilength. \nHere epilength is only select from 10, 12, 14, 16, 18, 20.\n"
	. "The fasta file only contains one sequence. More than one fasta sequence will only be dealt with the first one. \nExample: perl SVMTriP.pl data.fasta 20. \n");
}

if($epilength != 10 && $epilength != 12 && $epilength != 14 && $epilength != 16 && $epilength != 18 && $epilength != 20)
{
  die ("Error: Here epilength is only select from 10, 12, 14, 16, 18, 20.");	
}

print "Input file name: " . $filename . "\n";
print "selected epitope length: " . $epilength . "aa\n";

# declare subroutines
sub sildewindow;
sub trim;
sub getunisequence;
sub getoptimalprediction;
sub updateuniseqindb;
sub checkunisequence;
sub copyprediction;


my $workdir="/code/westbesel/static/prediction/scripts/SVMTriP/SVMTriP/SVMTriP/"; # set up the work path where SVMTriP package exists


# Get the top protein sequence (delete symbol ">")
my $unisequence = getunisequence($filename);

# slide windwo protein sequence
slidewindow($filename, $unisequence, $epilength);


# transform protein sequence
system("perl ".$workdir."CreateTestData.pl " . $filename . " " . $epilength);

# sleep 5s
sleep(2);

# predict every slide window by "svm_classify"
system($workdir."svm_classify " . $workdir.$filename . "_testing.dat ".$workdir."model_" . $epilength . "aa.dat " . $workdir.$filename . "_" . $epilength . "aa_prediction.txt");

# remove the intermediate files "*_testing.dat" and "*_20aa.fas"
system("rm -rf " . $workdir.$filename . "_testing.dat");
system("rm -rf " . $workdir.$filename . "_" . $epilength . "aa.fas");

# pick up the most optimal prediction results
getoptimalprediction($filename, $unisequence, $epilength);


# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}


# delete illegal character, e.g., space, 1-9
sub short_sequence
{
	my $line;
	($line) = @_;

	my $short_str  = "";

	my %aminolist = ("A"=>1, "C"=>1, "D"=>1, "E"=>1, "F"=>1, "G"=>1, "H"=>1, "I"=>1, "K"=>1, "L"=>1, 
			"M"=>1, "N"=>1, "P"=>1, "Q"=>1, "R"=>1, "S"=>1, "T"=>1, "V"=>1, "W"=>1, "Y"=>1);

	my $trim_str = trim($line);

	for(my $index=0; $index<length($trim_str); $index++)
	{
		my $residue = substr($trim_str, $index, 1);

		if($aminolist{uc($residue)})
		{
			$short_str = $short_str . uc($residue);
		}

	}

	return $short_str;
}


sub getunisequence
{
	my($filename);
	($filename) = @_;
	
	open(IN_FILE, $workdir.$filename);

	#Only the top sequence will be dealt with (when more than one sequence are input)
	my $topsequence = "";

	my $head = <IN_FILE>;

	if($head ne ">")
	{
		$topsequence = $topsequence . short_sequence($head);
	}

	while(my $line = <IN_FILE>)
	{
		chomp($line);
		
		if(substr($line, 0, 1) eq ">")
		{
			last;
		}
	
		$topsequence = $topsequence . short_sequence($line);
		
	}
	
	close(IN_FILE);

	return $topsequence;
}



# Slide window on protein sequence
sub slidewindow
{
	my($filename, $sequence, $epilength);
	($filename, $sequence, $epilength) = @_;

	open(OUT, ">" . $workdir . $filename . "_" . $epilength . "aa.fas");

	if(length($sequence) < $epilength)
	{
		print OUT $sequence . "\n";
	}
	else
	{
		for(my $index = 0; $index <= length($sequence) - $epilength; $index++)
		{
			print OUT substr($sequence, $index, $epilength) . "\n";
		}
	}

	close(OUT);
}


sub round {
    my($number) = shift;
    return int($number + .5);
}


# Get the most optimal predictions from SVM result
sub getoptimalprediction
{
	my($filename, $sequence, $epilength);
	($filename, $sequence, $epilength) = @_;

	# The maximum number of optimal epitopes is 10
	my $top_number = 10;

	print "Start to pick up the most optimal predictions for input sequence: \n";

	open(IN, $workdir . $filename . "_" . $epilength . "aa_prediction.txt") or die("Not find the file " . $filename .  "_" . $epilength . "aa_prediction.txt\n");;

	my $count = 0;

	my @scorelist;
	my @loclist;
	my @iniscorelist;

	while(my $line = <IN>)
	{
		$count ++;

		chomp($line);

		push(@scorelist, $line);
		push(@loclist, $count);	
	}

	close(IN);

	# this array contains all the scores from AA1 - AAend
	@iniscorelist = @scorelist;

	# order scorelist, at the same time update the loclist using the same order
	for(my $i=0; $i<=$#scorelist; $i++)
	{
		for(my $j=$i+1; $j<=$#scorelist; $j++)
		{
			if($scorelist[$j] > $scorelist[$i])
			{
				my $tempscore = $scorelist[$i];
				$scorelist[$i] = $scorelist[$j];
				$scorelist[$j] = $tempscore;

				my $temploc = $loclist[$i];
				$loclist[$i] = $loclist[$j];
				$loclist[$j] = $temploc;
			}
		}
	}

	# this array contains all the top sequences we pick up
	my @outloc;

	for(my $index_loc_all = 0; $index_loc_all <= $#loclist; $index_loc_all++)
	{
		if($top_number <= ($#outloc + 1))
		{
			last;
		}

		my $new_loc = $loclist[$index_loc_all];
		my $is_new = 0;


		if(((($new_loc > 2) && $iniscorelist[$new_loc-1] > $iniscorelist[$new_loc-3]) || $new_loc <= 2)
        	     && ((($new_loc > 1) && $iniscorelist[$new_loc-1] > $iniscorelist[$new_loc-2]) || $new_loc <= 1)
		     && ((($new_loc <= $#iniscorelist) && $iniscorelist[$new_loc-1] > $iniscorelist[$new_loc]) || $new_loc > $#iniscorelist)
		     && ((($new_loc < $#iniscorelist) && $iniscorelist[$new_loc-1] > $iniscorelist[$new_loc+1]) || $new_loc >= $#iniscorelist)
		   )
		{
			$is_new = 1;

			for(my $index_outloc = 0; $index_outloc<=$#outloc; $index_outloc++)
			{
				my $old_loc = $outloc[$index_outloc];

				# if the two locations are too close, i.e., the segments are overlapped, forsake this one
				if($new_loc >= ($old_loc - $epilength - 1) && $new_loc <= ($old_loc + $epilength - 1))
				{
					$is_new = 0;
					last;
				}
			}
		}

		if($is_new)
		{
			push(@outloc, $new_loc);
		}
	}
	
	open(OUT_OPT, ">" . $filename . "_" . $epilength . "aa_optimal_predictions.txt");
	
	print "Location, Epitope, Score\n";
	print OUT_OPT "Location, Epitope, Score\n";

	for(my $index=0; $index<=$#outloc; $index++)
	{
		my $iniloc = $outloc[$index];
		my $location = ($iniloc - 2) . " - " . ($iniloc + $epilength - 3);
		my $epitope = substr($sequence, $iniloc - 1, $epilength);
		my $score = $iniscorelist[$iniloc - 1];
		
		print $location .  ", " . $epitope . ", " . $score . "\n";
		print OUT_OPT $location .  ", " . $epitope . ", " . $score . "\n";
		
	}
}
