#!usr/bin/perl
use strict;
use warnings;

my $minThreshold = (shift @ARGV) + 1;
my $maxThreshold = 's';
my $inputFile = 'Matrices/1975_2015_window8';
my $outputFile = $inputFile . "_threshold" . ($minThreshold-1);
&thresholdMatrix($minThreshold, $maxThreshold, $inputFile, $outputFile);


#Applies a minimum and maximum number of co-occurrences threshold to a file by 
#copying the $inputFile to $outputFile, but ommitting lines that have less than
#$minThreshold number of co-occurrences, or greater than $maxThreshold
#If a value of 's' (for skip) is given for either threshold, that threshold 
#is not applied
sub thresholdMatrix {
    #grab the input
    my $minThreshold = shift;
    my $maxThreshold = shift;
    my $inputFile = shift;
    my $outputFile = shift;

    #open files
    open IN, $inputFile or die("ERROR: unable to open inputFile\n");
    open OUT, ">$outputFile" 
	or die ("ERROR: unable to open outputFile: $outputFile\n");

    print "Reading File\n";
    #threshold each line of the file
    my ($key, $copy, $cui1, $cui2, $val);
    while (my $line = <IN>) {
	#grab values and assume the line will be copied
	($cui1, $cui2, $val) = split(/\t/,$line);
	$copy = 1;

	#check minThreshold
	if ($minThreshold ne 's') {
	   if ($val < $minThreshold) {
		$copy = 0;
	   }
	}

	#check maxThreshold
	if ($maxThreshold ne 's') {
	   if ($val > $maxThreshold) {
		$copy = 0;
	   }
	}

	if ($copy == 1) {
	   print OUT $line;
	}
    }
    close IN;

    print "Done!\n";
}
