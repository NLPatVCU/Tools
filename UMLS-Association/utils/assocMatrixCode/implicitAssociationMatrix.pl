#!usr/bin/perl
use strict;
use warnings;

=head1 NAME

implicitAssociationMatrix.pl 

=head1 SYNOPSIS 

This program calculates an association matrix, calculating the association between all CUIs in the MiniMayo and UMNSRS datasets and all unique CUIs in a given matrix, functioning as the rows and columns, respectively.

=head1 USAGE

Usage: implicitAssociationMatrix.pl INPUTFILE MEASURE OUTPUTFILE

=head1 INPUT

=head2 [INPUTFILE] [MEASURE] [OUTPUTFILE]

INPUTFILE: The file containing the matrix to be worked with.  It is assumed that the implicit co-occurrence data files are already made, and stored in the same directory as the matrix, with its standard naming convention specified in getImplicitData.pl.  If not specified, the program will choose the matrix for 1975 onwards, window 1.

MEASURE: The association measure to be used in generating the matrix.  Use the same abbreviations as in UMLS::Interface.  If not specified, the program will choose Chi-Square (x2).

OUTPUTFILE: The file to which the implicit association matrix will be written.  If not specified, the output file will be the input file with "_ImplicitMatrix" concatenated to the end.

=head1 OUTPUT

A file, specified by the OUTPUTFILE option, that contains the aforementioned implicit association matrix. 

=head1 SYSTEM REQUIREMENTS

=over

=item * Perl (version 5.8.5 or better) - http://www.perl.org

=item * UMLS::Interface - http://search.cpan.org/dist/UMLS-Interface

=item * Text::NSP - http://search.cpan.org/dist/Text-NSP

=back

=head1 CONTACT US
    
    Bridget T. McInnes: btmcinnes at vcu.edu 
    Alexander D. McQuilkin: alexmcq99@yahoo.com

=head1 AUTHOR

 Bridget T. McInnes, Virginia Commonwealth University 
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

use UMLS::Interface;
use UMLS::Association;
use UMLS::Association::StatFinder;
use UMLS::Association::CuiFinder;

#change these for a different matrix and/or a different association measure
#Also specify the file the matrix is written to
#Change them from the command line
my $matrix = shift @ARGV;
if(! defined $matrix)
{
    $matrix = "Matrices/1975_2015_window1";
}

my $measure = shift @ARGV;
if(! defined $measure)
{
    $measure = "x2";
}

my $output = shift @ARGV;
if(! defined $output)
{
    $output = $matrix . "_ImplicitMatrix";
}

#don't change these
my %option_hash;
my %assoc_option_hash;
$option_hash{'t'} = 1;
$option_hash{'config'} = "UMLS-Association/utils/configuration";
my $umls = UMLS::Interface->new(\%option_hash);
$assoc_option_hash{"umls"} = $umls;
$assoc_option_hash{'database'} = "1809onward";
$assoc_option_hash{"matrix"} = $matrix;
$assoc_option_hash{"measure"} = $measure;

my $params = \%assoc_option_hash;
my $cuifinder = UMLS::Association::CuiFinder->new($params);
my $statfinder = UMLS::Association::StatFinder->new($params, $cuifinder);

my $npp = 0;
my $time;
my $start;

#get hash of all unique CUIs
my %cuis;
my %n1ps;
my %np1s;

my @cuiFiles = qw(UMLS-Association/utils/Data/MiniMayoSRS.snomedct.cuis UMLS-Association/utils/Data/UMNSRS_reduced_rel.cuis UMLS-Association/utils/Data/UMNSRS_reduced_sim.cuis);
my %importantCuis;
foreach my $cuiFile (@cuiFiles)
{
    open IN, $cuiFile or die "Cannot open $cuiFile for input: $!\n\n";
    while(my $line = <IN>)
    {
	chomp $line;
	my ($cui1, $cui2, $num) = split (/<>/,$line);

	$importantCuis{$cui1} = 1;
	$importantCuis{$cui2} = 1;
    }
    close IN;
}

open IN, $matrix or die "Cannot open $matrix for input: $!\n\n";
while(my $line = <IN>)
{
    chomp $line;
    my ($cui1, $cui2, $num) = split (/\t/,$line);

    if(!exists $cuis{$cui1})
    {
	$cuis{$cui1} = 1;
    }
    if(!exists $cuis{$cui2})
    {
	$cuis{$cui2} = 1;
    }
}
close IN;

$npp = scalar keys %cuis;

my $n1pFile = $matrix . "_ImplicitCUI1Data_sorted";
my $np1File = $matrix . "_ImplicitCUI2Data_sorted";

#If the implicit co-occurrence matrices are not built, build them
if(! -e $n1pFile || ! -e $np1File)
{
    `perl getImplicitData.pl $matrix`;
}

open IN, $n1pFile or die "Cannot open $n1pFile for input: $!\n\n";
while (my $line = <IN>)
{
    chomp $line;
    my ($cui, $data) = split /\t/, $line;
    $n1ps{$cui} = $data;
}
close IN;

open IN, $np1File or die "Cannot open $np1File for input: $!\n\n";
while (my $line = <IN>)
{
    chomp $line; 
    my ($cui, $data) = split /\t/, $line;
    $np1s{$cui} = $data;
}
close IN;

my @cuis = sort keys %cuis;
my @importantCuis = sort keys %importantCuis;
%cuis = ();
%importantCuis = ();
my $cuis = scalar @cuis;


open OUT, ">$output" or die "Cannot open $output for output: $!\n\n";
foreach my $i (0..scalar @importantCuis - 1)
{
    print STDERR "On Cui " . ($i+1) . " of " . (scalar @importantCuis) . "\n";
    my $row = _genRow($importantCuis[$i]);

    print OUT $row;
}
close OUT;

#  Method to get Implicit co-occurrence data for cui 1
#  input : $concept1 <- string containing cui 1
#  output: reference to %cui1Data, a hash with the keys as the cuis that appear second in a bigram with concept 1 being first and the values being 1
sub _getCui1Data
{
    my $concept1 = shift;
    my @cui1Data = split /,/,$n1ps{$concept1};
    my %cui1Data;
    foreach my $c (@cui1Data)
    {
	$cui1Data{$c} = 1;
    }

    return \%cui1Data;
}

#  Method to get Implicit co-occurrence data for cui 2
#  input : $concept2 <- string containing cui 2
#  output: reference to %cui2Data, a hash with the keys as the cuis that appear first in a bigram with concept 2 being second and the values being 1
sub _getCui2Data
{
    my $concept2 = shift;
    my @cui2Data = split /,/,$np1s{$concept2};
    my %cui2Data;
    foreach my $c (@cui2Data)
    {
	$cui2Data{$c} = 1;
    }

    return \%cui2Data; 
}

#  Method to generate a row of the implicit association matrix
#  input : $cui1 <- the cui for which this row is to be generated
#  output: $row <- a string contained each cui pair and implicit association score for this row, separated by a newline
sub _genRow
{
    my $cui1 = shift;
    my $row = "";
    my %cui1Data = %{_getCui1Data($cui1)};
    my $n1p = scalar keys %cui1Data;
    
    if(scalar keys %cui1Data == 0)
    {
	return $row;
    }

    foreach my $cui2 (@cuis)
    { 
	my %cui2Data = %{_getCui2Data($cui2)};

	if(scalar keys %cui2Data == 0)
	{
	    next;
	}
	my $np1 = scalar keys %cui2Data;
	
	my $n11 = 0;
	foreach my $c (keys %cui1Data)
	{
	    if (exists $cui2Data{$c})
	    {
		$n11++;
	    }
	}

	if($n11 == 0)
	{
	    next;
	}

	#print STDERR "CUI 1: $cui1 CUI 2: $cui2 N11: $n11 N1P: $n1p NP1: $np1 NPP: $npp\n";
	
	my $score = $statfinder->_calculateStatisticFromContingencyTable($n11, $n1p, $np1, $npp, $measure);
	if($score == 0)
	{
	    next;
	}
	
	$row .= "$cui1\t$cui2\t$score\n"; 
    }
    return $row;
}


