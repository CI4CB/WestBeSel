The command line is: perl SVMTriP.pl filename.fas epilenght. e.g.,
 
perl SVMTriP.pl example.fas 20
 
The current version of SVMTriP can only deal with a fasta sequence one time. If filename.fas contains more than one fasta sequenes, only the top one is dealt with.


Installation :

1.  extract SVMTriP and model files in the same directory
2.  cpan -i DBI
3.  change directory to working directory ($workdir) for SVMTrip.pl and GetTopPredictions.pl
4.  Correct default SVMTrip.pl line at the end " my $location = ($iniloc - 2) . " - " . ($iniloc + $epilength - 3); " 
