
=head1 NAME

cosineResults.pl

=head1 SYNOPSIS 

This program generates association scores and spearmans results for MiniMayo and UMNSRS, using cosine similarity with a given matrix of association scores.

=head1 USAGE

Usage: cosineResults.pl MATRIX

=head1 INPUT

=head2 [MATRIX] 

MATRIX: The file containing the matrix of co-occurrence data, NOT the matrix of association scores.  This matrix is needed to get a complete list of unique CUIs, and the association score matrix should be easy to find if named properly (for matrix "Matrix" it should be "Matrix_ImplicitMatrix").

=head1 OUTPUT

A directory named "CosineResults" which contains another directory, "Scores" as well as a file named "Data"
The "Scores" directory contains all the cosine similarity scores of all the CUI pairs within each dataset (MiniMayo coders & physicians, UMNSRS sim & rel).  There is one file per dataset.
The "Data" file contains the spearmans scores and N values for each dataset scorefile with its respective gold standard.

=head1 SYSTEM REQUIREMENTS

=over

=item * Perl (version 5.8.5 or better) - http://www.perl.org

=item * UMLS::Interface - http://search.cpan.org/dist/UMLS-Interface

=item * Text::NSP - http://search.cpan.org/dist/Text-NSP

=back

=head1 CONTACT US
    
    Alexander D. McQuilkin: alexmcq99@yahoo.com

=head1 AUTHOR

 Alexander D. McQuilkin, Virginia Commonwealth University

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to:

 The Free Software Foundation, Inc.,
 59 Temple Place - Suite 330,
 Boston, MA  02111-1307, USA.

=cut

###############################################################################

#                               THE CODE STARTS HERE
###############################################################################

#usr/bin/perl
use UMLS::Interface;
use UMLS::Association;
use strict;
use warnings;

my @datasets = qw(coders physicians sim rel);
my $matrix = "Matrices/1975_2015_window1";
my $implicitMatrix = $matrix . "_ImplicitMatrix";
my %matrix;
my %cuis;

#get list of unique cuis
open IN, $matrix or die "Cannot open $matrix for input: $!\n\n";
while (my $line = <IN>)
{
    chomp $line;
    my ($cui1,$cui2,$num) = split /\t/,$line;
    if(! exists $cuis{$cui1}){$cuis{$cui1} = 1;}
    if(! exists $cuis{$cui2}){$cuis{$cui2} = 1;}
}
close IN;

my @cuis = sort keys %cuis;
%cuis = ();

#create hash of association scores we need
open IN, $implicitMatrix or die "Cannot open $implicitMatrix for input $!\n\n";
while (my $line = <IN>)
{
    chomp $line;
    my ($cui1,$cui2,$num) = split /\t/,$line;
    $matrix{$cui1}{$cui2} = $num;
}
close IN;

#make directories if necessary
if (! -d "CosineResults"){`mkdir CosineResults`}
if (! -d "CosineResults/Scores"){`mkdir CosineResults/Scores`}

#generate scores and then data
my $datafile = "CosineResults/Data";
open DATA, ">$datafile" or die "Cannot open $datafile for output: $!\n\n";
foreach my $dataset (@datasets)
{
    #set up cuifile and goldfile depending on the dataset
    my $cuifile = "Data/";
    my $goldfile = "Data/";
    if($dataset eq "coders" || $dataset eq "physicians")
    {
	$cuifile .= "MiniMayoSRS.snomedct.cuis";
	$goldfile .= "MiniMayoSRS.snomedct." . $dataset;
    }
    else
    {
	$cuifile .= "UMNSRS_reduced_" . $dataset . ".cuis";
	$goldfile .= "UMNSRS_reduced_" . $dataset . ".gold";
    }

    #generate scores for each line of the cui file
    my $scorefile = "CosineResults/Scores/$dataset";
    open IN, $cuifile or die "Cannot open $cuifile for input: $!\n\n";
    open SCORES, ">$scorefile" or die "Canot open $scorefile for output: $!\n\n";
    while (my $line = <IN>)
    {
	chomp $line;
	my ($cui1,$cui2) = split /<>/,$line;
	my $score = _getScore($cui1,$cui2);
	print SCORES "$score<>$line\n";
    }
    close IN;
    close SCORES;
    
    #get spearmans data
    my @outputs = `perl spearmans.pl $goldfile $scorefile`;
			
    $outputs[0] =~ /(\d+)/g;
    my $n = $1;
    $outputs[1] =~ /(\d+\.\d+)/g;
    my $score = $1;

    print DATA _cutString($dataset) . "\t" . _cutString($score) . "\t$n\n";
}
close DATA;

sub _getScore
{
    my $concept1 = shift;
    my $concept2 = shift;

    my $dot = 0;
    my $magA = 0;
    my $magB = 0;

    foreach my $cui (keys %{$matrix{$concept1}})
    { 
        #check to see if a score exists for both $concept1 and $concept2 in a given column ($cui).  If so, add to the dot product.
	if (exists $matrix{$concept2}{$cui})
	{
	    $dot += $matrix{$concept1}{$cui} * $matrix{$concept2}{$cui};
	}
	
	#add the square this dimension's component of the A vector to the magnitude
	$magA += $matrix{$concept1}{$cui} ** 2;
    }

    #add the square this dimension's component of the A vector to the magnitude
    foreach my $cui (keys %{$matrix{$concept2}})
    {
	$magB += $matrix{$concept2}{$cui} ** 2;
    }
    
    #square root the magnitude variables to get the actual magnitudes
    $magA = sqrt($magA);
    $magB = sqrt($magB);

    if ($magA == 0 || $magB == 0)
    {
	return -1;
    }

    #get cosine between vectors
    return $dot / ($magA * $magB);
}

sub _cutString
{
    my $str = shift;
    my $max = 7; #maximum string length so the formatting won't mess up
    my $length = length $str;
    if(!defined $length)
    {
	$length = 0;
    }

    if($length > $max)
    {
	return substr($str,0,$max-$length);
    }
    else
    {
	return $str;
    }
}
