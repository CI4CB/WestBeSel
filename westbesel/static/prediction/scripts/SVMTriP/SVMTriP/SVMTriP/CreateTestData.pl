#! /usr/bin/perl

#Summary: Create Traint Data by SVM format of input
# kernel: tri-blosum62-tripropensity

#output file: Train_testing.dat
#output format: 1 feature:value for positive; -1 for negative

my $filename = $ARGV[0];
my $epilength = $ARGV[1];

my $workdir="./";  #set up the work path where the SVMTriP Package exists

print "Start to transform protein sequence into the data format for SVM classification\n";

@aminolist = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y");

#Get Hashtable containing blosum62 matrix scores
%ht_blosum62 = GetBlosum62();

#Get Hashtable containing single amino acid propensity
%ht_propensity = GetSingleAminoPropensity();

#Get Hashtable containing tri amino acides propensity
%ht_propensity_tri_residues = GetPropensityTriRes();


open(IN, $workdir.$filename . "_" . $epilength . "aa.fas");
open(OUT, ">" . $workdir.$filename . "_testing.dat");


#read lines from seq files
while(<IN>)
{
	chomp;
	my $line1 = $_;

	print OUT "0 ";

	$feature = 1;

	my $line = trim($line1);

	for(my $i=0; $i<=$#aminolist; $i++)
	{
		for(my $j=0; $j<=$#aminolist; $j++)
		{
			for(my $k=0; $k<=$#aminolist; $k++)
			{
				$sumscore = 0;
				$tempscore = 0;

				for(my $index=0; $index<length($line)-2; $index++)
				{
					$tempscore = GetTempScore_Tri_Resi(substr($line, $index, 3), $aminolist[$i], $aminolist[$j], $aminolist[$k]);
					$sumscore = $sumscore + $tempscore;
				}

				if($sumscore > 0)
				{
					my $temp = $sumscore / (length($line) - 2);

					if($temp > 0.01)
					{
						print OUT $feature . ":";
						printf OUT "%.2f", $temp ;
						print OUT " ";
					}
				}
				
				$feature++;
			}
		}
	}

	print OUT "\n";
}

print OUT "!";

close(IN);
close(OUT);

print "Transformation finished\n";


#Get Hashtable contains blosum62 matrix scores
sub GetBlosum62
{
	open(IN_BLOSUM62, $workdir."blosum62.csv");
	
	#The first line is seemed as the Head of blosum62 matrix
	$isHead = 1;

	while(<IN_BLOSUM62>)
	{
		chomp;
		$line1 = $_;
		$line = trim($line1);

		if($isHead)
		{
			@head = split(/,/, $line);
			
			$isHead = 0;
			next;
		}

		@collection = split(/,/, $line);

		for($i=0; $i<=$#head; $i++)
		{
			$str_key = GetKeyInHtBlosum62($collection[0], $head[$i]);

			if(!$temphash{$str_key})
			{
				$temphash{$str_key} = $collection[$i + 1];
			}
		}
	}

	close(IN_BLOSUM62);
	return %temphash;

}


# Create key for hashtable blosum62
# key is composed of two amino acid characters, always in alphabetical order
sub GetKeyInHtBlosum62
{
	my($residue1, $residue2) = @_;

	if($residue1 gt $residue2)
	{
		return $residue2 . $residue1;
	}

	return $residue1 . $residue2;
}


#Get Hashtable containing single amid acid propensity
sub GetSingleAminoPropensity()
{
	# Read the stat about epitope
	open(IN_Stat_Epitope, $workdir."stat_epitope.txt");
	while(<IN_Stat_Epitope>)
	{
		chomp;
		my $line = $_;
		my @collection = split(/\t/, trim($line));

		$ht_epitope{substr($collection[0], 0, 1)} = $collection[1];
	}

	close(IN_Stat_Epitope);

	# Read the stat about random protein
	open(IN_Stat_RandomProtein, $workdir."stat_randomprotein.txt");
	while(<IN_Stat_RandomProtein>)	
	{
		chomp;
		my $line = $_;
		my @collection = split(/\t/, trim($line));

		$ht_randomprotein{substr($collection[0], 0, 1)} = $collection[1];
	}

	close(IN_Stat_RandomProtein);

	foreach $key (keys(%ht_epitope))
	{
		$ht_propensity{$key} = $ht_epitope{$key} / $ht_randomprotein{$key};
	}

	return %ht_propensity;
}


# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}


#Get Hashtable containing tri amino acides propensity
sub GetPropensityTriRes
{
	open(IN_Stat_Tri_Res, $workdir."stat_tri_residues_" . $epilength . "aa.txt");

	while(<IN_Stat_Tri_Res>)
	{
		chomp;
		my $line = $_;
		my @collection = split(/\s/, trim($line));

		$ht_propensity_tri_resi{substr($collection[0], 0, 3)} = $collection[1];
	}

	close(IN_Stat_Tri_Res);

	return %ht_propensity_tri_resi;
}


sub GetTempScore_Tri_Resi
{
	my($seq, $amino1, $amino2, $amino3) = @_;

	$seq_trans = substr($seq, 2, 1) .substr($seq, 1, 1) . substr($seq, 0, 1);

	#Get the stat value of si and trans for seq1
	if($ht_propensity_tri_residues{$seq})
	{
		$propensity_tri_1_si = $ht_propensity_tri_residues{$seq};
	}
	else
	{
		$propensity_tri_1_si = 0;
	}

	if($ht_propensity_tri_residues{$seq_trans})
	{
		$propensity_tri_1_trans = $ht_propensity_tri_residues{$seq_trans};
	}
	else
	{
		$propensity_tri_1_trans = 0;
	}

	$seq_2_si = $amino1 . $amino2 . $amino3;
	$seq_2_trans = $amino3 . $amino2 . $amino1;
	
	#Get the stat value of si and trans for seq2
	if($ht_propensity_tri_residues{$seq_2_si})
	{
		$propensity_tri_2_si = $ht_propensity_tri_residues{$seq_2_si};
	}
	else
	{
		$propensity_tri_2_si = 0;
	}

	if($ht_propensity_tri_residues{$seq_2_trans})
	{
		$propensity_tri_2_trans = $ht_propensity_tri_residues{$seq_2_trans};
	}
	else
	{
		$propensity_tri_2_trans = 0;
	}
	
	my $tempscore1 = 0;
	my $tempscore2 = 0;
	my $tempscore3 = 0;

	my $mykey = GetKeyInHtBlosum62(substr($seq, 0, 1), $amino1);
	$tempscore1 = $ht_blosum62{$mykey};

	$mykey = GetKeyInHtBlosum62(substr($seq, 1, 1), $amino2);
	$tempscore2 = $ht_blosum62{$mykey};

	$mykey = GetKeyInHtBlosum62(substr($seq, 2, 1), $amino3);
 	$tempscore3 = $ht_blosum62{$mykey};

	$temp = ($tempscore1 + $tempscore2 + $tempscore3) / 3;

	if($temp > 0)
	{
		return $temp * ($propensity_tri_1_si + $propensity_tri_1_trans) * ($propensity_tri_2_si + $propensity_tri_2_trans) / 0.94; # 0.94 is the average value of propensity multiple.
	}

	return 0;
}


