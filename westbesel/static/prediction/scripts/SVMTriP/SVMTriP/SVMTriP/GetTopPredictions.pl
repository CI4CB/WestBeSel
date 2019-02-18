#! /usr/bin/perl -w

#Summary: get top segments for prediction results of Tri_Blosum62_tripropensity

#input file: the file containing the prediction score of slide window of 20aa
#output: the top $top_number 20aa segments with highest score

use strict;
use warnings;
use DBI;

# Read the argument from command line
my $filename = $ARGV[0];
my $top_number = 10;

my $workdir= "/code/westbesel/static/prediction/scripts/SVMTriP/SVMTriP/SVMTriP/";

print "Start to pick up the most optimal predictions for job id $filename\n";

open(IN, $workdir.$filename . "_prediction");

my $count = 0;

my @scorelist;
my @loclist;
my @iniscorelist;

while(<IN>)
{
	$count ++;

	chomp;
	my $line = $_;

	push(@scorelist, $line);
	push(@loclist, $count);
	
}

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

			# if the two locations are too close, i.e., the segments are overlapped, merge the two into one
			if($new_loc >= ($old_loc - 19) && $new_loc <= ($old_loc + 19))
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

# save the proteins sequence
for(my $index=0; $index<=$#outloc; $index++)
{
	print $outloc[$index] . ": ";
	my $index1 = $outloc[$index];
	print $iniscorelist[$index1-1] . "\n";
}


close(IN);


sub round {
    my($number) = shift;
    return int($number + .5);
}
