#UMLS::Association
#
# Perl module for scoring the semantic association of terms in the Unified
# Medical Language System (UMLS).
#
# Copyright (c) 2015
#
# Bridget T. McInnes, Virginia Commonwealth University
# btmcinnes at vcu.edu
#
# Keith Herbert, Virginia Commonwealth University
# herbertkb at vcu.edu
#
# Alexander D. McQuilkin, Virginia Commonwealth University 
# alexmcq99 at yahoo.com
#
# Sam Henry, Virginia Commonwealth University
# henryst at vcu.edu
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to
#
# The Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330,
# Boston, MA  02111-1307, USA.

package UMLS::Association::StatFinder;

use Fcntl;
use strict;
use warnings;
use DBI;
use bytes;
use File::Spec;

#  error handling variables
my $errorhandler = "";
my $cuifinder      = ""; 

my $pkg = "UMLS::Association::StatFinder";

#  debug variables
local(*DEBUG_FILE);

#  global variables
my $debug     = 0;
my $NPP       = 0; 
my $umls = undef;
my $precision = 4;
my $getdescendants = 0;
my $implicit = 0;
my $countMethod = 1;
my $noOrder = 0;
my $matrix = 0;

######################################################################
#  functions to initialize the package
######################################################################

#  method to create a new UMLS::Association::StatFinder object
#  input : $params <- reference to hash of database parameters
#          $handler <- reference to cuifinder object 
#  output: $self
sub new {
    my $self = {};
    my $className = shift;
    my $params = shift;
    my $handler = shift; 

    # bless the object.
    bless($self, $className);

    # initialize error handler
    $errorhandler = UMLS::Association::ErrorHandler->new();
    if(! defined $errorhandler) {
        print STDERR "The error handler did not get passed properly.\n";
        exit;
    }

    #  initialize the cuifinder
    $cuifinder = $handler; 

    #  initialize global variables
    $debug = 0; 

    # initialize the object.
    $self->_initialize($params);
    return $self;
}

#  method to initialize the UMLS::Association::StatFinder object.
#  input : $parameters <- reference to a hash of database parameters
#  output:
sub _initialize {

    my $self = shift;
    my $params = shift;
    my %params = %{$params};

    #set global variables using option hash
    $umls = $params{'umls'};
    $getdescendants = $params{'getdescendants'};
    $implicit = $params{'implicit'};
    $countMethod = !$params{'cooccurrence'};
    $noOrder = $params{'noorder'};
    $matrix = $params{'matrix'};
    my $function = "_initialize";
    &_debug($function);

    #  check self
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    $params = {} if(!defined $params);
    if(defined $params{'precision'})
    {
	$precision = $params{'precision'};
    }
}

sub _debug {
    my $function = shift;
    if($debug) { print STDERR "In UMLS::Association::StatFinder::$function\n"; }

}

######################################################################
#  functions to get statistical information about the cuis
######################################################################
#  Method to return the frequency of a concept pair
#  input : $concept1 <- string containing a cui 1
#          $concept2 <- string containing a cui 2
#  output: $frequency <- frequency of cui pair
sub _getFrequency {

    my $self = shift;
    my $concept1 = shift;
    my $concept2 = shift;

    my $function = "_getFrequency";

    #  check self
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #  check parameter exists
    if(!defined $concept1) {
        $errorhandler->_error($pkg, $function, "Error with input variable \$concept1.", 4);
    }
    if(!defined $concept2) {
        $errorhandler->_error($pkg, $function, "Error with input variable \$concept2.", 4);
    }

    #  check if valid concept
    if(! ($errorhandler->_validCui($concept1)) ) {
        $errorhandler->_error($pkg, $function, "Concept ($concept1) is not valid.", 6);
    }
    if(! ($errorhandler->_validCui($concept2)) ) {
        $errorhandler->_error($pkg, $function, "Concept ($concept2) is not valid.", 6);
    }
    
    #  check if concept exists
    if(! ($cuifinder->_exists($concept1)) ) {
        return -1; #$errorhandler->_error($pkg, $function, "Concept ($concept1) does not exist.", 6);
    }    
    #  check if concept exists
    if(! ($cuifinder->_exists($concept2)) ) {
        return -1; $errorhandler->_error($pkg, $function, "Concept ($concept2) does not exist.", 6);
    }    
    #  set up database
    my $db = $cuifinder->_getDB(); 
    
    my $freqRef;
    if($noOrder)
    {
	$freqRef = $db->selectcol_arrayref("select n_11 from N_11 where (cui_1='$concept1' and cui_2='$concept2') or (cui_1='$concept2' and cui_2='$concept1')"); 
    }
    else
    {
	$freqRef = $db->selectcol_arrayref("select n_11 from N_11 where cui_1='$concept1' and cui_2='$concept2'"); 
    }
    
    my $freq = shift @{$freqRef}; 
    
    if(defined $freq) 
    { 
	return $freq;
    } 
    else { return 0; }
}

#  Method to return the np1 of a concept 
#  input : $concept <- string containing a cui 1
#  output: $np1 <- number of times concept occurs in second bigram position
sub _getNp1 {

    my $self = shift;
    my $concept = shift; 

    my $function = "_getNp1"; 

    #  check self
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #  check parameter exists
    if(!defined $concept) {
        $errorhandler->_error($pkg, $function, "Error with input variable \$concept.", 4);
    }

    #  check if valid concept
    if(! ($errorhandler->_validCui($concept)) ) {
        $errorhandler->_error($pkg, $function, "Concept ($concept) is not valid.", 6);
    }   

    #  check if concept exists
    if(! ($cuifinder->_exists($concept)) ) {
        return -1; #$errorhandler->_error($pkg, $function, "Concept ($concept) does not exist.", 6);
    }    
    
    #  set up database
    my $db = $cuifinder->_getDB(); 

    my $np1Ref;
    if($noOrder)
    {
	$np1Ref = $db->selectcol_arrayref("select sum(n_11) from N_11 where cui_1='$concept' or cui_2 = '$concept'"); 
    }
    else
    {
	$np1Ref = $db->selectcol_arrayref("select sum(n_11) from N_11 where cui_2='$concept'"); 
    }
    
    my $np1 = shift @{$np1Ref}; 
    
    if(defined $np1) 
    {
	return $np1;
    } 
    else { return 0; }
}

#  Method to return the n1p of a concept 
#  input : $concept <- string containing a cui 1
#  output: $n1p <- number of times concept occurs in second bigram position
sub _getN1p {

    my $self = shift;
    my $concept = shift; 

    my $function = "_getN1p"; 

    #  check self
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #  check parameter exists
    if(!defined $concept) {
        $errorhandler->_error($pkg, $function, "Error with input variable \$concept.", 4);
    }

    #  check if valid concept
    if(! ($errorhandler->_validCui($concept)) ) {
        $errorhandler->_error($pkg, $function, "Concept ($concept) is not valid.", 6);
    }    
    
    #  check if concept exists
    if(! ($cuifinder->_exists($concept)) ) {
        return -1; #$errorhandler->_error($pkg, $function, "Concept ($concept) does not exist.", 6);
    } 
    
    #  set up database
    my $db = $cuifinder->_getDB(); 
    
    my $n1pRef;

    if($noOrder)
    {
	$n1pRef = $db->selectcol_arrayref("select sum(n_11) from N_11 where cui_1='$concept' or cui_2 = '$concept'"); 
    }
    else
    {
	$n1pRef = $db->selectcol_arrayref("select sum(n_11) from N_11 where cui_1='$concept'"); 
    }
    
    my $n1p = shift @{$n1pRef}; 
    
    if(defined $n1p) 
    { 
	return $n1p;
    } 
    else { return 0; }
}

#  Method to return npp
#  input : none
#  output: $npp <- number of total concept pairs
sub _getNpp {

    my $self = shift;
    
    my $function = "_getNpp"; 

    if($NPP > 0) { return $NPP; }
    
    #  check self
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    if($matrix)
    {
	my %cuis;

	if($matrix eq "Matrices/1975_2015_window1")
	{
	    #print STDERR "Got to fast NPP statement\n";
	    return 3644085387;
	}
	if($matrix eq "Matrices/1975_2015_window1_threshold1")
	{
	    #print STDERR "Got to fast NPP statement\n";
	    return 3512917705;
	}
	if($matrix eq "Matrices/1975_2015_window8")
	{
	    #print STDERR "Got to fast NPP statement\n";
	    return 27724843453;
	}
	if($matrix eq "Matrices/1975_2015_window8_threshold1")
	{
	    #print STDERR "Got to fast NPP statement\n";
	    return 27323428532;
	}
	if($matrix eq "Matrices/2000_2015_window1_ordered")
	{
	    #print STDERR "Got to fast NPP statement\n";
	    return 2369249353;
	}
	else
	{
	    open IN, $matrix or die "Cannot open $matrix for input: $!\n";
	    while (my $line = <IN>)
	    {
		chomp $line;
		my ($cui1, $cui2, $num) = split(/\t/,$line);
		
		$NPP += $num;
	    }
	    close IN;
	}

	print STDERR "Matrix: $matrix\nNPP: $NPP\n";
	return $NPP;
    }

    #  set up database
    my $db = $cuifinder->_getDB(); 
    
    my $nppRef = $db->selectcol_arrayref("select sum(N_11) from N_11"); 
    
    $NPP = shift @{$nppRef}; 

    if($NPP <= 0) { errorhandler->_error($pkg, $function, "", 5); } 
    
    return $NPP; 
}

#  Method to optimized data retrieval
#  input : $concept1 <- string containing a cui 1
#          $concept2 <- string containing a cui 2
#          $n11Ref   <- reference to the sorted array of lines containing co-occurrence information between all cuis (format CUI\tCUI\tNUM)
#          $n1pRef   <- reference to the sorted array of lines containing n1p values for all cuis (format CUI\tNUM)
#          $n1pRef   <- reference to the sorted array of lines containing np1 values for all cuis (format CUI\tNUM)
#  output: reference to @data = (n_11, n_1p, n_p1)

sub _getData{
    my $self = shift;
    
    my $concept1 = shift;
    my $concept2 = shift;
    
    my $function = "_getData"; 
    
    #  check self
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }
    
    #  check parameter exists
    if(!defined $concept1) {
        $errorhandler->_error($pkg, $function, "Error with input variable \$concept1.", 4);
    }
    if(!defined $concept2) {
        $errorhandler->_error($pkg, $function, "Error with input variable \$concept2.", 4);
    }
    
    #  check if concept exists
    if(! ($cuifinder->_exists($concept1)) ) {
        return -1; #$errorhandler->_error($pkg, $function, "Concept ($concept1) does not exist.", 6);
    }    
    #  check if concept exists
    if(! ($cuifinder->_exists($concept2)) ) {
        return -1; $errorhandler->_error($pkg, $function, "Concept ($concept2) does not exist.", 6);
    }   

    my $n11 = 0;
    my $n1p = 0;
    my $np1 = 0;

    #checks matrix file option
    if($matrix)
    {
	my $n11Ref = shift;
	my $n1pRef = shift;
	my $np1Ref = shift;

	my %n11s;
	my %n1ps;
	my %np1s;
	
	#checks if references to arrays are defined; if not, build the arrays
	if(!defined $n11Ref || !defined $n1pRef || !defined $np1Ref)
	{
	    print STDERR "Array references didn't get through :(\n";
	    open IN, $matrix or die "Cannot open $matrix for input: $!\n";
	    while (my $line = <IN>)
	    {
		chomp $line;
		my ($cui1, $cui2, $num) = split /\t/, $line;
		$n11s{$cui1}{$cui2} = $num;
	    }
	    close IN;

	    my $n1pFile = $matrix . "_N1P_sorted";
	    my $np1File = $matrix . "_NP1_sorted";

	    open IN, $n1pFile or die "Cannot open $n1pFile for input: $!\n\n";
	    while (my $line = <IN>)
	    {
		chomp $line;
		my ($cui1, $num) = split /\t/, $line;
		$n1ps{$cui1} = $num;
	    }
	    close IN;

	    open IN, $np1File or die "Cannot open $np1File for input: $!\n\n";
	    while (my $line = <IN>)
	    {
		chomp $line;
		my ($cui2, $num) = split /\t/, $line;
		$np1s{$cui2} = $num;
	    }
	    close IN;

	}
	#arrays references are defined, let's get the arrays
	else
	{
	    %n11s = %{$n11Ref};
	    %n1ps = %{$n1pRef};
	    %np1s = %{$np1Ref};
	}

	#get n11
        my $val = $n11s{$concept1}{$concept2};
	if(defined $val)
	{
	    $n11 += $val;
	}

	#get n1p
	$val = $n1ps{$concept1};
	if(defined $val)
	{
	    $n1p += $val;
	}

	#get np1
	$val =  $np1s{$concept2};
	if(defined $val)
	{
	    $np1 += $val;
	}

	#if order doesn't matter, we have to do more work
	if($noOrder)
	{
	    #get n11 again
	    $val = $n11s{$concept2}{$concept1};
	    if(defined $val)
	    {
		$n11 += $val;
	    }

	    #get n1p again
	    $val = $np1s{$concept1};
	    if(defined $val)
	    {
		$n1p += $val;
	    }

	    #you get the picture
	    $val =  $n1ps{$concept2};
	    if(defined $val)
	    {
		$np1 += $val;
	    }

	}
    }

    else
    {
	#  get frequency and marginal totals
	$n11 = $self->_getFrequency($concept1, $concept2); 
	$n1p = $self->_getN1p($concept1); 
	$np1 = $self->_getNp1($concept2); 
    }
    my @data = ($n11, $n1p, $np1); 
    #print STDERR "End of getData\n";
    
    return \@data;
    
}

#  Method to optimized data retrieval (using the descendants of each cui)
#  input : $concept1 <- string containing a cui 1
#          $concept2 <- string containing a cui 2
#          $n11Ref   <- reference to the hash of hashes containing co-occurrence information between all cuis (format $n11s{$cui1}{$cui2} = $num)
#          $n1pRef   <- reference to the hash containing n1p values for all cuis (format $n1ps{$cui} = $num)
#          $n1pRef   <- reference to the hash containing np1 values for all cuis (format $np1s{$cui} = $num)
#  output: reference to @data = (n_11, n_1p, n_p1)
sub _getDescendantData{
    
    my $self = shift;
    
    my $concept1 = shift;
    my $concept2 = shift;
    
    #  get descendants of each cui
    my @descendants1 =@{_findDescendants($concept1)};
    push @descendants1, $concept1;

    my @descendants2 = @{_findDescendants($concept2)};
    push @descendants2, $concept2;

    my $function = "_getDescendantData"; 
    
    #  check self
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }
    
    #  check parameter exists
    if(!defined $concept1) {
        $errorhandler->_error($pkg, $function, "Error with input variable \$concept1.", 4);
    }
    if(!defined $concept2) {
        $errorhandler->_error($pkg, $function, "Error with input variable \$concept2.", 4);
    }
    
    #  check if concept exists
    if(! ($cuifinder->_exists($concept1)) ) {
        return -1; #$errorhandler->_error($pkg, $function, "Concept ($concept1) does not exist.", 6);
    }    
    #  check if concept exists
    if(! ($cuifinder->_exists($concept2)) ) {
        return -1; $errorhandler->_error($pkg, $function, "Concept ($concept2) does not exist.", 6);
    }    
    #  set up database
    my $db = $cuifinder->_getDB(); 

    #get n11
    my $queryString;
    my $n11 = 0;
    my $n1p = 0;
    my $np1 = 0;

    #checks --matrix option
    if($matrix)
    {
	my $n11Ref = shift;
	my $n1pRef = shift;
	my $np1Ref = shift;

	my %n11s;
	my %n1ps;
	my %np1s;
	
	#checks if references to arrays are defined; if not, build the arrays
	if(!defined $n11Ref || !defined $n1pRef || !defined $np1Ref)
	{
	    print STDERR "Array references didn't get through :(\n";
	    open IN, $matrix or die "Cannot open $matrix for input: $!\n";
	    while (my $line = <IN>)
	    {
		chomp $line;
		my ($cui1, $cui2, $num) = split /\t/, $line;
		$n11s{$cui1}{$cui2} = $num;
	    }
	    close IN;

	    #sets up files from which we will obtain occurrence data
	    my $n1pFile = $matrix . "_N1P_sorted";
	    my $np1File = $matrix . "_NP1_sorted";

	    open IN, $n1pFile or die "Cannot open $n1pFile for input: $!\n\n";
	    while (my $line = <IN>)
	    {
		chomp $line;
		my ($cui1, $num) = split /\t/, $line;
		$n1ps{$cui1} = $num;
	    }
	    close IN;

	    open IN, $np1File or die "Cannot open $np1File for input: $!\n\n";
	    while (my $line = <IN>)
	    {
		chomp $line;
		my ($cui2, $num) = split /\t/, $line;
		$np1s{$cui2} = $num;
	    }
	    close IN;

	}
	#arrays references are defined, let's get the arrays
	else
	{
	    %n11s = %{$n11Ref};
	    %n1ps = %{$n1pRef};
	    %np1s = %{$np1Ref};
	}
	
	#add up n1ps of each descendant of concept1, as well as that of concept1 
	foreach my $desc1 (@descendants1)
	{
	    my $num = $n1ps{$desc1};
	    if(defined $num)
	    {
		$n1p += $num;
	    }
	}
	#same for concept2
	foreach my $desc2 (@descendants2)
	{
	    my $num = $np1s{$desc2};
	    if (defined $num)
	    {
		$np1 += $num;
	    }
	}
	#add up n11s for every permutation of descendants (order matters, concept 1 must be first)
	foreach my $desc1 (@descendants1)
	{
	    foreach my $desc2 (@descendants2)
	    {
		my $num = $n11s{$desc1}{$desc2};
		if(defined $num)
		{
		    $n11 += $num;
		}
	    }
	}

	if($noOrder)
	{
	    #add up occurrences of concept 1 as second cui to n1p
	    foreach my $desc1 (@descendants1)
	    {
		my $num = $np1s{$desc1};
		if(defined $num)
		{
		    $n1p += $num;
		}
	    }
	    #add up occurrences of concept 2 as first cui to np1
	    foreach my $desc2 (@descendants2)
	    {
		my $num = $n1ps{$desc2};
		if (defined $num)
		{
		    $np1 += $num;
		}
	    }
	    #add up n11s for every permutation of  descendants (order doesn't matter)
	    foreach my $desc1 (@descendants1)
	    {
		foreach my $desc2 (@descendants2)
		{
		    my $num = $n11s{$desc2}{$desc1};
		    if(defined $num)
		    {
			$n11 += $num;
		    }
		}
	    }
	}
    
    }

    else
    {

	if ($noOrder)
	{
	    $queryString = "select SUM(n_11) from N_11 where ((cui_2 = '$concept1' ";
	    
	    foreach my $desc (@descendants1)
	    {
		$queryString .= "or cui_2 = '$desc' ";
	    }

	    $queryString .= ") and (cui_1 = '$concept2' ";

	    foreach my $desc (@descendants2)
	    {
		$queryString .= "or cui_1 = '$desc' ";
	    }

	    $queryString .= ")) or ((cui_1 = '$concept1' ";
	    
	    foreach my $desc (@descendants1)
	    {
		$queryString .= "or cui_1 = '$desc' ";
	    }

	    $queryString .= ") and (cui_2 = '$concept2' ";

	    foreach my $desc (@descendants2)
	    {
		$queryString .= "or cui_2 = '$desc' ";
	    }

	    $queryString .= "));";

	    $n11 = shift @{$db->selectcol_arrayref($queryString)};
	}

	else
	{
	    $queryString = "select SUM(n_11) from N_11 where (cui_1 = '$concept1' ";
	    
	    foreach my $desc (@descendants1)
	    {
		$queryString .= "or cui_1 = '$desc' ";
	    }

	    $queryString .= ") and (cui_2 = '$concept2' ";

	    foreach my $desc (@descendants2)
	    {
		$queryString .= "or cui_2 = '$desc' ";
	    }

	    $queryString .= ");";

	    $n11 = shift @{$db->selectcol_arrayref($queryString)};
	}

	#get n1p

	if($noOrder)
	{
	    $queryString = "select SUM(n_11) from N_11 where (cui_2 = '$concept1' ";
	    
	    foreach my $desc (@descendants1)
	    {
		$queryString .= "or cui_2 = '$desc' ";
	    }

	    $queryString .= ") or (cui_1 = '$concept1' ";

	    foreach my $desc (@descendants1)
	    {
		$queryString .= "or cui_1 = '$desc' ";
	    }

	    $queryString .= ");";

	    $n1p = shift @{$db->selectcol_arrayref($queryString)};
	}

	else
	{
	    $queryString = "select SUM(n_11) from N_11 where (cui_1 = '$concept1' ";
	    
	    foreach my $desc (@descendants1)
	    {
		$queryString .= "or cui_1 = '$desc' ";
	    }

	    $queryString .= ");";

	    $n1p = shift @{$db->selectcol_arrayref($queryString)};
	}

	#get np1

	if($noOrder)
	{
	    $queryString = "select SUM(n_11) from N_11 where (cui_2 = '$concept2' ";
	    
	    foreach my $desc (@descendants2)
	    {
		$queryString .= "or cui_2 = '$desc' ";
	    }

	    $queryString .= ") or (cui_1 = '$concept2' ";

	    foreach my $desc (@descendants2)
	    {
		$queryString .= "or cui_1 = '$desc' ";
	    }

	    $queryString .= ");";

	    $np1 = shift @{$db->selectcol_arrayref($queryString)};
	}

	else
	{
	    $queryString = "select SUM(n_11) from N_11 where (cui_2 = '$concept2' ";
	    
	    foreach my $desc (@descendants2)
	    {
		$queryString .= "or cui_2 = '$desc' ";
	    }

	    $queryString .= ");";

	    $np1 = shift @{$db->selectcol_arrayref($queryString)};
	}
    }

    my @data = ($n11, $n1p, $np1);
    return \@data;
}

#  Method to optimized data retrieval (using implicit association)
#  input : $concept1 <- string containing a cui 1
#          $concept2 <- string containing a cui 2
#          $npp      <- Number of unique cuis (used for matrix option; can be used otherwise, but it's impractical to do so)
#          $n1pRef   <- reference to a hash containing implicit co-occurrence data for cui 1 (only for matrix option)
#          $np1Ref   <- reference to a hash containing implicit co-occurrence data for cui 2 (only for matrix option)
#  output: reference to @data = (n_11, n_1p, n_p1, n_pp)
sub _getImplicitData
{
    my $self = shift;

    my $concept1 = shift;
    my $concept2 = shift;

    #  set up database
    my $db = $cuifinder->_getDB(); 
    
    my %cui1Data;
    my %cui2Data;
    my @t1;
    my @t2;
    my @data;

    my $npp = shift;
    
    #check --matrix option
    if($matrix)
    {
	#get references
	my $n1pRef = shift;
	my $np1Ref = shift;

	my $n11 = 0;
	my $n1p = 0;
	my $np1 = 0;
	
	#get the hashes we need
	my %n1ps = %{$n1pRef};
	my %np1s = %{$np1Ref};

	#gets strings of cuis that occur with concept1 and concept2, separated by commas
	if (defined $n1ps{$concept1})
	{
	    foreach my $cui (split /,/,$n1ps{$concept1})
	    {
		if (! exists $cui1Data{$cui})
		{
		    $cui1Data{$cui} = 1;
		}
	    }
	}
	if (defined $np1s{$concept2})
	{
	    foreach my $cui (split /,/,$np1s{$concept2})
	    {
		if (! exists $cui2Data{$cui})
		{
		    $cui2Data{$cui} = 1;
		}
	    }
	}

	#gets other cuis that occur with concept1 and concept2, if order doesn't matter
	if ($noOrder)
	{
	    if (defined $np1s{$concept1})
	    {
		foreach my $cui (split /,/,$np1s{$concept1})
		{
		    if (! exists $cui1Data{$cui})
		    {
			$cui1Data{$cui} = 1;
		    }
		}
	    }
	    if (defined $n1ps{$concept2})
	    {
		foreach my $cui (split /,/,$n1ps{$concept2})
		{
		    if (! exists $cui2Data{$cui})
		    {
			$cui2Data{$cui} = 1;
		    }
		}
	    }
	}

	if(scalar keys %cui1Data > 0 && scalar keys %cui2Data > 0)
	{
	    #Find number of CUIs that co-occur with both CUI 1 and CUI 2
	    foreach my $c (keys %cui1Data)
	    {
		if (exists $cui2Data{$c})
		{
		    $n11++;
		}
	    }
	    
	    $n1p = scalar keys %cui1Data;
	    $np1 = scalar keys %cui2Data;
	}
	else
	{
	    if(! scalar keys %cui1Data > 0)
	    {
		$n1p = 0;
	    }
	    if(! scalar keys %cui2Data > 0)
	    {
		$np1 = 0;
	    }
	    $n11 = 0;
	}

	print STDERR "$concept1 $concept2 $n11 $n1p $np1 $npp\n";
	@data = ($n11, $n1p, $np1, $npp);
    }

    else
    {
	if($noOrder)
	{
	    #set up cui 1 hash
	    @t1 = @{$db->selectcol_arrayref("SELECT cui_2 FROM N_11 WHERE cui_1 = '$concept1' AND n_11 > 0;")};
	    @t2 = @{$db->selectcol_arrayref("SELECT cui_1 FROM N_11 WHERE cui_2 = '$concept1' AND n_11 > 0;")};
	    
	    foreach my $t (@t1)
	    {
		$cui1Data{$t} = shift @{$db->selectcol_arrayref("SELECT n_11 FROM N_11 WHERE (cui_1 = '$concept1' AND cui_2 = '$t') OR (cui_1 = '$t' AND cui_2 = '$concept1');")};
	    }

	    foreach my $t (@t2)
	    {
		if(!exists $cui1Data{$t})
		{
		    $cui1Data{$t} = shift @{$db->selectcol_arrayref("SELECT n_11 FROM N_11 WHERE (cui_1 = '$concept1' AND cui_2 = '$t') OR (cui_1 = '$t' AND cui_2 = '$concept1');")};
		}
	    }

	    #set up cui 2 hash
	    @t1 = @{$db->selectcol_arrayref("SELECT cui_2 FROM N_11 WHERE cui_1 = '$concept2' AND n_11 > 0;")};
	    @t2 = @{$db->selectcol_arrayref("SELECT cui_1 FROM N_11 WHERE cui_2 = '$concept2' AND n_11 > 0;")};

	    foreach my $t (@t1)
	    {
		$cui2Data{$t} = shift @{$db->selectcol_arrayref("SELECT n_11 FROM N_11 WHERE (cui_1 = '$concept2' AND cui_2 = '$t') OR (cui_1 = '$t' AND cui_2 = '$concept2');")};
	    }

	    foreach my $t (@t2)
	    {
		if(!exists $cui2Data{$t})
		{
		    $cui2Data{$t} = shift @{$db->selectcol_arrayref("SELECT n_11 FROM N_11 WHERE (cui_1 = '$concept2' AND cui_2 = '$t') OR (cui_1 = '$t' AND cui_2 = '$concept2');")};
		}
	    }
	}

	else
	{
	    #set up cui 1 hash
	    @t1 = @{$db->selectcol_arrayref("SELECT cui_2 FROM N_11 WHERE cui_1 = '$concept1' AND n_11 > 0;")};
	    
	    foreach my $t (@t1)
	    {
		$cui1Data{$t} = shift @{$db->selectcol_arrayref("SELECT n_11 FROM N_11 WHERE cui_1 = '$concept1' AND cui_2 = '$t';")};
	    }

	    #set up cui 2 hash
	    @t2 = @{$db->selectcol_arrayref("SELECT cui_1 FROM N_11 WHERE cui_2 = '$concept2' AND n_11 > 0;")};
	    
	    foreach my $t (@t2)
	    {
		$cui2Data{$t} = shift @{$db->selectcol_arrayref("SELECT n_11 FROM N_11 WHERE cui_2 = '$concept2' AND cui_1 = '$t' > 0;")};
	    }
	}

	#find number of vocabulary words (npp)

	my %cui1Counts;
	my %cui2Counts;

	my @cui1s = @{$db->selectcol_arrayref("SELECT cui_1 FROM N_11;")};
	foreach my $cui1 (@cui1s)
	{
	    if(!exists $cui1Counts{$cui1})
	    {
		$cui1Counts{$cui1} = shift @{$db->selectcol_arrayref("SELECT SUM(n_11) FROM N_11 WHERE cui_1 = '$cui1';")};
	    }
	}

	my @cui2s  = @{$db->selectcol_arrayref("SELECT cui_2 FROM N_11;")};
	foreach my $cui2 (@cui2s)
	{
	    if(!exists $cui2Counts{$cui2})
	    {
		$cui2Counts{$cui2} = shift @{$db->selectcol_arrayref("SELECT SUM(n_11) FROM N_11 WHERE cui_2 = '$cui2';")};
	    }
	}

	my $numUnigrams = scalar keys %cui1Counts;
	foreach my $cui2 (keys %cui2Counts)
	{
	    if(!exists $cui1Counts{$cui2})
	    {
		$numUnigrams++;
	    }
	}

	#set up intersecting data hash
	my %intersectingData;
	foreach my $c1 (keys %cui1Data)
	{
	    if(exists $cui2Data{$c1})
	    {
		my $min = $cui1Data{$c1};
		if($cui2Data{$c1} < $min)
		{
		    $min = $cui2Data{$c1};
		}

		$intersectingData{$c1} = $min;
	    }
	}

	#set up return data depending on method(n11, n1p, np1, npp)

	if($countMethod)
	{
	    my $n11 = scalar keys %intersectingData;
	    my $n1p = scalar keys %cui1Data;
	    my $np1 = scalar keys %cui2Data;
	    if(!defined $npp)
	    {
		$npp = $numUnigrams;
	    }

	    if($n11 > $n1p || $n11 > $np1 || $n11 > $npp)
	    {
		print STDERR "cui 1: $concept1\ncui 2: $concept2\n\nn11: $n11\nn1p: $n1p\nnp1: $np1\nnpp: $npp\n\n";
	    }
	    
	    @data = ($n11, $n1p, $np1, $npp);
	}

	else
	{
	    my $n11 = 0;
	    foreach my $cui (keys %intersectingData)
	    {
		$n11 += $intersectingData{$cui};
	    }

	    my $n1p = 0;
	    foreach my $cui (keys %cui1Data)
	    {
		$n1p += $cui1Data{$cui};
	    }

	    my $np1 = 0;
	    foreach my $cui (keys %cui2Data)
	    {
		$np1 += $cui2Data{$cui};
	    }
	    
	    if(!defined $npp)
	    {
		$npp = shift @{$db->selectcol_arrayref("SELECT n_pp FROM N_PP")};
	    }
	    
	    if($n11 > $n1p || $n11 > $np1 || $n11 > $npp)
	    {
		print STDERR "cui 1: $concept1\ncui 2: $concept2\n\nn11: $n11\nn1p: $n1p\nnp1: $np1\nnpp: $npp\n\n";
	    }

	    @data = ($n11, $n1p, $np1, $npp);
	}
    }

    return \@data;

}


#  Method to optimized data retrieval (using implicit association and descendants)
#  input : $concept1 <- string containing a cui 1
#          $concept2 <- string containing a cui 2
#          $npp      <- Number of unique cuis
#          $n1pRef   <- reference to a hash containing implicit co-occurrence data for cui 1 (only for matrix option)
#          $np1Ref   <- reference to a hash containing implicit co-occurrence data for cui 2 (only for matrix option)
#  output: reference to @data = (n_11, n_1p, n_p1, n_pp)
sub _getImplicitDataWithDescendants
{
    my $self = shift;

    my $concept1 = shift;
    my $concept2 = shift;

    my $npp = shift;

    #  get descendants of each cui, stored in strings $desc1 and $desc2
    my @descendants1 =@{_findDescendants($concept1)};
    push @descendants1, $concept1;

    my @descendants2 = @{_findDescendants($concept2)};
    push @descendants2, $concept2;

    my $function = "_getImplicitDataWithDescendants"; 
    
    #  check self
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }
    #  set up database
    my $db = $cuifinder->_getDB(); 
    
    my %cui1Data;
    my %cui2Data;
    my @t1;
    my @t2;
    my @data;
    
    #check --matrix option
    if($matrix)
    {
	#get references
	my $n1pRef = shift;
	my $np1Ref = shift;

	my $n11 = 0;
	my $n1p = 0;
	my $np1 = 0;
	
	#get the hashes we need
	my %n1ps = %{$n1pRef};
	my %np1s = %{$np1Ref};

	#for each descendant of concept 1, get cuis that occur with it and add them to the cui 1 data hash
	foreach my $desc1 (@descendants1)
	{
	    if (defined $n1ps{$desc1})
	    {
		foreach my $cui (split /,/,$n1ps{$desc1})
		{
		    if(! exists $cui1Data{$cui})
		    {
			$cui1Data{$cui} = 1;
		    }
		}
	    }
	}
	#same for cui 2
	foreach my $desc2 (@descendants2)
	{
	    if (defined $np1s{$desc2})
	    {
		foreach my $cui (split /,/,$np1s{$desc2})
		{
		    if(! exists $cui2Data{$cui})
		    {
			$cui2Data{$cui} = 1;
		    }
		}
	    }
	}

	#gets other cuis that occur with concept1 and concept2 and their descendants, if order doesn't matter
	if ($noOrder)
	{
	    foreach my $desc1 (@descendants1)
	    {
		if (defined $np1s{$desc1})
		{
		    foreach my $cui (split /,/,$np1s{$desc1})
		    {
			if(! exists $cui1Data{$cui})
			{
			    $cui1Data{$cui} = 1;
			}
		    }
		}
	    }
	    foreach my $desc2 (@descendants2)
	    {
		if(defined $n1ps{$desc2})
		{
		    foreach my $cui (split /,/,$n1ps{$desc2})
		    {
			if(! exists $cui2Data{$cui})
			{
			    $cui2Data{$cui} = 1;
			}
		    }
		}
	    }
	}

	#checks to see if concept1 and concept2 actually occur with anything; if so, calculate data
	if (scalar keys %cui1Data > 0 && scalar keys %cui2Data > 0)
	{
	    #Find number of CUIs that co-occur with both CUI 1 and CUI 2
	    foreach my $c (keys %cui1Data)
	    {
		if (exists $cui2Data{$c})
		{
		    $n11++;
		}
	    }
	    
	    $n1p = scalar keys %cui1Data;
	    $np1 = scalar keys %cui2Data;
	}
	#one of the two concepts doesn't occur with anything.  Adjust accordingly
	else
	{
	    if(scalar keys %cui1Data == 0)
	    {
		$n1p = 0;
	    }
	    if(scalar keys %cui2Data == 0)
	    {
		$np1 = 0;
	    }
	    $n11 = 0;
	}
	
	@data = ($n11, $n1p, $np1, $npp);
    }

    else
    {
	if($noOrder)
	{
	    #set up cui 1 hash
	    my $qs = "SELECT N_11.cui_2 FROM N_11 WHERE (N_11.cui_1 = '$concept1' ";
	    foreach my $desc (@descendants1)
	    {
		$qs .= "OR N_11.cui_1 = '$desc' ";
	    }
	    $qs .= ") AND N_11.n_11 > 0;";

	    @t1 = @{$db->selectcol_arrayref($qs)};

	    $qs =  "SELECT N_11.cui_1 FROM N_11 WHERE (N_11.cui_2 = '$concept1' ";
	    foreach my $desc (@descendants1)
	    {
		$qs .= "OR N_11.cui_2 = '$desc' ";
	    }
	    $qs .= ") AND N_11.n_11 > 0;";

	    @t2 = @{$db->selectcol_arrayref($qs)};

	    foreach my $t (@t1)
	    {
		$qs = "SELECT SUM(N_11.n_11) FROM N_11 WHERE (N_11.cui_1 = '$concept1' ";
		foreach my $desc (@descendants1)
		{
		    $qs .= "OR N_11.cui_1 = '$desc' ";
		}
		$qs .= ") AND N_11.cui_2 = '$t';";

		my $num =  shift @{$db->selectcol_arrayref($qs)};
		
		$qs = "SELECT SUM(N_11.n_11) FROM N_11 WHERE (N_11.cui_2 = '$concept1' ";
		foreach my $desc (@descendants1)
		{
		    $qs .= "OR N_11.cui_2 = '$desc' ";
		}
		$qs .= ") AND N_11.cui_1 = '$t';";

		$num += shift @{$db->selectcol_arrayref($qs)};

		$cui1Data{$t} = $num;
	    }

	    foreach my $t (@t2)
	    {
		if(!exists $cui1Data{$t})
		{
		    $qs = "SELECT SUM(N_11.n_11) FROM N_11 WHERE (N_11.cui_1 = '$concept1' ";
		    foreach my $desc (@descendants1)
		    {
			$qs .= "OR N_11.cui_1 = '$desc' ";
		    }
		    $qs .= ") AND N_11.cui_2 = '$t';";

		    my $num =  shift @{$db->selectcol_arrayref($qs)};
		    
		    $qs = "SELECT SUM(N_11.n_11) FROM N_11 WHERE (N_11.cui_2 = '$concept1' ";
		    foreach my $desc (@descendants1)
		    {
			$qs .= "OR N_11.cui_2 = '$desc' ";
		    }
		    $qs .= ") AND N_11.cui_1 = '$t';";

		    $num += shift @{$db->selectcol_arrayref($qs)};

		    $cui1Data{$t} = $num;
		}
	    }

	    #set up cui 2 hash
	    $qs = "SELECT N_11.cui_2 FROM N_11 WHERE (N_11.cui_1 = '$concept2' ";
	    foreach my $desc (@descendants2)
	    {
		$qs .= "OR N_11.cui_1 = '$desc' ";
	    }
	    $qs .= ") AND N_11.n_11 > 0;";

	    @t1 = @{$db->selectcol_arrayref($qs)};

	    $qs =  "SELECT N_11.cui_1 FROM N_11 WHERE (N_11.cui_2 = '$concept2' ";
	    foreach my $desc (@descendants2)
	    {
		$qs .= "OR N_11.cui_2 = '$desc' ";
	    }
	    $qs .= ") AND N_11.n_11 > 0;";

	    @t2 = @{$db->selectcol_arrayref($qs)};

	    foreach my $t (@t1)
	    {
		$qs = "SELECT SUM(N_11.n_11) FROM N_11 WHERE (N_11.cui_1 = '$concept2' ";
		foreach my $desc (@descendants2)
		{
		    $qs .= "OR N_11.cui_1 = '$desc' ";
		}
		$qs .= ") AND N_11.cui_2 = '$t';";

		my $num =  shift @{$db->selectcol_arrayref($qs)};
		
		$qs = "SELECT SUM(N_11.n_11) FROM N_11 WHERE (N_11.cui_2 = '$concept2' ";
		foreach my $desc (@descendants2)
		{
		    $qs .= "OR N_11.cui_2 = '$desc' ";
		}
		$qs .= ") AND N_11.cui_1 = '$t';";

		$num += shift @{$db->selectcol_arrayref($qs)};

		$cui2Data{$t} = $num;
	    }

	    foreach my $t (@t2)
	    {
		if(!exists $cui2Data{$t})
		{
		    $qs = "SELECT SUM(N_11.n_11) FROM N_11 WHERE (N_11.cui_1 = '$concept2' ";
		    foreach my $desc (@descendants2)
		    {
			$qs .= "OR N_11.cui_1 = '$desc' ";
		    }
		    $qs .= ") AND N_11.cui_2 = '$t';";

		    my $num =  shift @{$db->selectcol_arrayref($qs)};
		    
		    $qs = "SELECT SUM(N_11.n_11) FROM N_11 WHERE (N_11.cui_2 = '$concept2' ";
		    foreach my $desc (@descendants2)
		    {
			$qs .= "OR N_11.cui_2 = '$desc' ";
		    }
		    $qs .= ") AND N_11.cui_1 = '$t';";

		    $num += shift @{$db->selectcol_arrayref($qs)};

		    $cui2Data{$t} = $num;
		}
	    }
	}

	else
	{
	    #set up cui 1 hash
	    my $qs = "SELECT N_11.cui_2 FROM N_11 WHERE (N_11.cui_1 = '$concept1' ";
	    foreach my $desc (@descendants1)
	    {
		$qs .= "OR N_11.cui_1 = '$desc' ";
	    }
	    $qs .= ") AND N_11.n_11 > 0;";

	    @t1 = @{$db->selectcol_arrayref($qs)};
	    
	    foreach my $t (@t1)
	    {
		$qs = "SELECT SUM(N_11.n_11) FROM N_11 WHERE (N_11.cui_1 = '$concept1' ";
		foreach my $desc (@descendants1)
		{
		    $qs .= "OR N_11.cui_1 = '$desc' ";
		}
		$qs .= ") AND N_11.n_11 > 0;";

		$cui1Data{$t} =  shift @{$db->selectcol_arrayref($qs)};
	    }

	    #set up cui 2 hash
	    $qs = "SELECT N_11.cui_2 FROM N_11 WHERE (N_11.cui_1 = '$concept2' ";
	    foreach my $desc (@descendants2)
	    {
		$qs .= "OR N_11.cui_1 = '$desc' ";
	    }
	    $qs .= ") AND N_11.n_11 > 0;";

	    @t1 = @{$db->selectcol_arrayref($qs)};
	    
	    foreach my $t (@t1)
	    {
		$qs = "SELECT SUM(N_11.n_11) FROM N_11 WHERE (N_11.cui_1 = '$concept2' ";
		foreach my $desc (@descendants2)
		{
		    $qs .= "OR N_11.cui_1 = '$desc' ";
		}
		$qs .= ") AND N_11.n_11 > 0;";

		$cui2Data{$t} =  shift @{$db->selectcol_arrayref($qs)};
	    }

	}

	#find number of vocabulary words (npp)

	my %cui1Counts;
	my %cui2Counts;

	my @cui1s = @{$db->selectcol_arrayref("SELECT N_11.cui_1 FROM N_11;")};
	foreach my $cui1 (@cui1s)
	{
	    if(!exists $cui1Counts{$cui1})
	    {
		$cui1Counts{$cui1} = shift @{$db->selectcol_arrayref("SELECT SUM(N_11.n_11) FROM N_11 WHERE N_11.cui_1 = '$cui1';")};
	    }
	}

	my @cui2s  = @{$db->selectcol_arrayref("SELECT N_11.cui_2 FROM N_11;")};
	foreach my $cui2 (@cui2s)
	{
	    if(!exists $cui2Counts{$cui2})
	    {
		$cui2Counts{$cui2} = shift @{$db->selectcol_arrayref("SELECT SUM(N_11.n_11) FROM N_11 WHERE N_11.cui_2 = '$cui2';")};
	    }
	}

	my $numUnigrams = scalar keys %cui1Counts;
	foreach my $cui2 (keys %cui2Counts)
	{
	    if(!exists $cui1Counts{$cui2})
	    {
		$numUnigrams++;
	    }
	}

	#set up intersecting data hash
	my %intersectingData;
	foreach my $c1 (keys %cui1Data)
	{
	    if(exists $cui2Data{$c1})
	    {
		my $min = $cui1Data{$c1};
		if($cui2Data{$c1} < $min)
		{
		    $min = $cui2Data{$c1};
		}

		$intersectingData{$c1} = $min;
	    }
	}

	#set up return data depending on method(n11, n1p, np1, npp)

	if($countMethod)
	{
	    my $n11 = scalar keys %intersectingData;
	    my $n1p = scalar keys %cui1Data;
	    my $np1 = scalar keys %cui2Data;
	    $npp = $numUnigrams;

	    @data = ($n11, $n1p, $np1, $npp);
	}

	else
	{
	    my $n11 = 0;
	    foreach my $cui (keys %intersectingData)
	    {
		$n11 += $intersectingData{$cui};
	    }

	    my $n1p = 0;
	    foreach my $cui (keys %cui1Data)
	    {
		$n1p += $cui1Data{$cui};
	    }

	    my $np1 = 0;
	    foreach my $cui (keys %cui2Data)
	    {
		$np1 += $cui2Data{$cui};
	    }

	    $npp = shift @{$db->selectcol_arrayref("SELECT N_PP.n_pp FROM N_PP")};

	    @data = ($n11, $n1p, $np1, $npp);
	}
    }
    return \@data;

}

#  Method to retrieve descendants of a cui
#  input : $cui <- string containing a cui 
#  output: reference to @descendants, the descendants of the given cui
sub _findDescendants
{
    my $cui = shift;

    my $hashref = $umls->findDescendants($cui);
    my @descendants = (sort keys %{$hashref});
    return \@descendants;
}

# calculate a contingency table values, and a statistic (association score) 
# between two concepts
# input:  $concept1  <- the cui of the first concept 
#         $concept2  <- the cui of the second concept
#         $statistic <- the string specifying the stat to calc
# output: the statistic (Association score) between the two concepts
sub _calculateStatistic { 
    my $self = shift;
    my $concept1 = shift; 
    my $concept2 = shift; 
    my $statistic = shift; 
    
    my $function = "_calculateStatistic"; 

    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }
    
    # get frequency and marginal totals optimized
    my $valid = -1;

    if($getdescendants && !$implicit)
    {
	$valid = $self->_getDescendantData($concept1, $concept2);
    }
    elsif($implicit && !$getdescendants)
    {
	$valid = $self->_getImplicitData($concept1, $concept2);
    }
    elsif($implicit && $getdescendants)
    {
	$valid = $self->_getImplicitDataWithDescendants($concept1, $concept2);
    }
    else
    {
	$valid = $self->_getData($concept1, $concept2);
    }

    if($valid == -1){
        return -1;
    }
    my @data = @{$valid};
    
    my $n11 = $data[0];
    my $n1p = $data[1];
    my $np1 = $data[2];

    if(!defined $n11 || !defined $n1p || !defined $np1){
	return -1;
    }
    
    my $npp;
    if($implicit)
    {
	$npp = $data[3];
    }
    else
    {
	$npp = $self->_getNpp(); 
    }
    
    #calculate the statistic and return it
    my $stat = $self->_calculateStatisticFromContingencyTable(
	$n11, $n1p, $np1, $npp, $statistic);
    if ($debug)  {
	print "$n11,$n1p,$np1,$npp = $stat\n";
    }

    return $stat;
}

# calculate a statistic (association score) using just contingency table values
# input:  $n11 <- n11 for the cui pair
#         $n1p <- n1p for the cui pair
#         $np1 <- np1 for the cui pair
#         $npp <- npp for the cui pair
#         $statistic <- the string specifying the stat to calc
# output: the statistic (Association score) between the two concepts
sub _calculateStatisticFromContingencyTable {
    my $self = shift;
    my $n11 = shift;
    my $n1p = shift;
    my $np1 = shift;
    my $npp = shift;
    my $statistic = shift;

    # set frequency and marginal totals
    my %values = (n11=>$n11, 
		  n1p=>$n1p, 
		  np1=>$np1, 
		  npp=>$npp); 
    
    if($n11 < 0 || $n1p < 0 || $np1 < 0) { 	
	#print STDERR "N11: $n11 N1P: $n1p NP1: $np1 NPP: $npp Score: -1\n";
	return -1.000; 
    }
    
    if($n11 == 0) { 
	#print STDERR "N11: $n11 N1P: $n1p NP1: $np1 NPP: $npp Score: 0\n";
	return 0.000;
    }
    
    #  set default statistic
    if(! defined $statistic) { $statistic = "tscore";  }
    
    #  set statistic module
    my $includename = ""; my $usename = "";  my $ngram = 2; 
    if($statistic eq "ll")  { 
	$usename = 'Text::NSP::Measures::'.$ngram.'D::MI::'.$statistic;
	$includename = File::Spec->catfile('Text','NSP','Measures',$ngram.'D','MI',$statistic.'.pm');
    }
    elsif($statistic eq "pmi" || $statistic eq "tmi" || $statistic eq "ps") { 
	$usename = 'Text::NSP::Measures::'.$ngram.'D::MI::'.$statistic;
	$includename = File::Spec->catfile('Text','NSP','Measures',$ngram.'D','MI',$statistic.'.pm');
    }
    elsif($statistic eq "x2"||$statistic eq "phi") {
	$usename = 'Text::NSP::Measures::'.$ngram.'D::CHI::'.$statistic;
	$includename = File::Spec->catfile('Text','NSP','Measures',$ngram.'D','CHI',$statistic.'.pm');
    }
    elsif($statistic eq "leftFisher"||$statistic eq "rightFisher"||$statistic eq "twotailed") { 
	if($statistic eq "leftFisher")	       { $statistic = "left";  }
	elsif($statistic eq "rightFisher")  { $statistic = "right"; }
	$usename = 'Text::NSP::Measures::'.$ngram.'D::Fisher::'.$statistic;
	$includename = File::Spec->catfile('Text','NSP','Measures',$ngram.'D','Fisher',$statistic.'.pm');
    }
    elsif($statistic eq "dice" || $statistic eq "jaccard") {
	$usename = 'Text::NSP::Measures::'.$ngram.'D::Dice::'.$statistic;
	$includename = File::Spec->catfile('Text','NSP','Measures',$ngram.'D','Dice',$statistic.'.pm');
    }
    elsif($statistic eq "odds") { 
	$usename = 'Text::NSP::Measures::'.$ngram.'D::'.$statistic;
	$includename = File::Spec->catfile('Text','NSP','Measures',$ngram.'D',$statistic.'.pm');
    }
    elsif($statistic eq "tscore") { 
	$usename = 'Text::NSP::Measures::'.$ngram.'D::CHI::'.$statistic;
	$includename = File::Spec->catfile('Text','NSP','Measures',$ngram.'D','CHI',$statistic.'.pm');
    }
    
    #  import module
    require $includename;
    import $usename;
    
    #  get statistics
    my $statisticValue = calculateStatistic(%values); 
    
    # check for errors/warnings from statistics.pm     
    my $errorMessage=""; 
    my $errorCode = getErrorCode(); 
    if (defined $errorCode) { 
	if($errorCode =~ /^1/) { 
	    printf(STDERR "Error from statistic library!\n  Error code: %d\n", $errorCode);
	    $errorMessage = getErrorMessage();
	    print STDERR "  Error message: $errorMessage\n" if( $errorMessage ne "");
	    exit; # exit on error
	}
	if ($errorCode =~ /^2/)  { 
	    printf(STDERR "Warning from statistic library!\n  Warning code: %d\n", $errorCode);
	    $errorMessage = getErrorMessage();
	    print STDERR "  Warning message: $errorMessage\n" if( $errorMessage ne "");
	    print STDERR "Skipping ngram\n";
	    #print STDERR "Skipping ngram $concept1<>$concept2\n";
	    next; # if warning, dont save the statistic value just computed
	}
    }

    #return statistic to given precision.  if no precision given, default is 4
    my $floatFormat = join '', '%', '.', $precision, 'f';
    
    my $statScore = sprintf $floatFormat, $statisticValue;

    #print STDERR "N11: $n11 N1P: $n1p NP1: $np1 NPP: $npp Score: $statScore\n";
    return $statScore; 
}


1;

__END__

=head1 NAME

UMLS::Association::StatFinder - provides the statistical association information 
of the concept pairs in the UMLS 

=head1 DESCRIPTION
    For more information please see the UMLS::Association.pm documentation.

=head1 SYNOPSIS

use UMLS::Association::StatFinder;
use UMLS::Association::ErrorHandler;

%params = ();

$statfinder = UMLS::Association::StatFinder->new(\%params);
die "Unable to create UMLS::Association::StatFinder object.\n" if(!$statfinder);

my $cui1 = C0018563;   
my $cui2 = C0446516; 

#  get the frequecy
my $freq = $statfinder->_getFrequency($cui1, $cui2); 

#  get marginal totals
my $np1 = $statfinder->_getNp1($cui2);
my $n1p = $statfinder->_getN1p($cui1); 

# calculate measure assocation
my $measure = "ll"; 
my $score = $statfinder->_calculateStatistic($cui1, $cui2, $measure); 

=head1 INSTALL

    To install the module, run the following magic commands:

    perl Makefile.PL
    make
    make test
    make install

    This will install the module in the standard location. You will, most
    probably, require root privileges to install in standard system
    directories. To install in a non-standard directory, specify a prefix
    during the 'perl Makefile.PL' stage as:

    perl Makefile.PL PREFIX=/home/bridget

    It is possible to modify other parameters during installation. The
    details of these can be found in the ExtUtils::MakeMaker
    documentation. However, it is highly recommended not messing around
    with other parameters, unless you know what you're doing.

    =head1 SEE ALSO

    L<http://tech.groups.yahoo.com/group/umls-similarity/>

    =head1 AUTHOR

    Bridget T McInnes <bmcinnes@vcu.edu>
    Andriy Y. Mulyar  <andriy.mulyar@gmail.com>
    Alexander D. McQuilkin <alexmcq99@yahoo.com>

    =head1 COPYRIGHT

    Copyright (c) 2015
    Bridget T. McInnes, Virginia Commonwealth University
    btmcinnes at vcu.edu

    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the Free Software
    Foundation; either version 2 of the License, or (at your option) any later
    version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with
    this program; if not, write to

    The Free Software Foundation, Inc.,
    59 Temple Place - Suite 330,
    Boston, MA  02111-1307, USA.
