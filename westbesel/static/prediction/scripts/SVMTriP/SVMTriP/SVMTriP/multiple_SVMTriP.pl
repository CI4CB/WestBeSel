#!/usr/bin/perl

use strict;
use warnings;

my $file= $ARGV[0];		#name of the file with many fasta sequences

open (my $fh, "<", "$file") or die "Can't open file : $!";
#`rm out.SVMT.pos`;


$/= ">";
my $dump= <$fh>;		#dump the first ">"

while(<$fh>)
{
    if(/^(.*?)\n(\w+)/)
    {
        my $fasid = $1;
        my $fasseq= $2;
        #print "$fasid\n$fasseq\n";
        open(OUT, ">", "out.SVMT$fasid");		#make a new file for each fasta
        print OUT ">$fasid\n$fasseq";
        system "perl SVMTriP.pl out.SVMT$fasid 20";
        system "rm out.SVMT$fasid";
        `(printf "$fasid\n";cat out.SVMT$fasid\_20aa_optimal_predictions.txt; echo) >> out.SVMT.pos`;
        system "rm out.SVMT$fasid\_20aa_optimal_predictions.txt";
        system "rm out.SVMT$fasid\_20aa_prediction.txt";
    }  
}    

