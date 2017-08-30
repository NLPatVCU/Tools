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
#local(*DEBUG_FILE);

#  global variables
my $debug     = 0;
my $umls = undef;
my $precision = 4;
my $conceptExpansion = 0;
my $lta = 0;
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
#  output: none
sub _initialize {

    my $self = shift;
    my $params = shift;
    my %params = %{$params};

    #set global variables using option hash
    $umls = $params{'umls'};
    $conceptExpansion = $params{'conceptExpansion'};
    $lta = $params{'lta'};
    $noOrder = $params{'noorder'};
    $matrix = $params{'matrix'};
    my $function = "_initialize";
    &_debug($function);

    #  check self
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    $params = {} if(!defined $params);
    if(defined $params{'precision'}) {
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
sub _getN11 {
    #grab parameters
    my $self = shift;
    my $concept1 = shift;
    my $concept2 = shift;

    my $function = "_getN11";

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
    #set up database
    my $db = $cuifinder->_getDB(); 
    
    #grab the frequency from the DB
    my $freqRef;
    if($noOrder) {
	$freqRef = $db->selectcol_arrayref("select sum(n_11) from N_11 where (cui_1='$concept1' and cui_2='$concept2') or (cui_1='$concept2' and cui_2='$concept1')"); 
    }
    else {
	$freqRef = $db->selectcol_arrayref("select n_11 from N_11 where cui_1='$concept1' and cui_2='$concept2'"); 
    }
    my $freq = shift @{$freqRef}; 
    
    #return the frequency if defined, else return 0
    if(defined $freq) { 
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

    #retrieve NP1
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
    
    #return np1 if defined, else 0
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
    
    #retreive n1p
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
sub getNpp {
    my $self = shift;
    my $function = "_getNpp"; 

    #  check self
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #calculate NPP from a matrix
    if($matrix) {
	my $npp = 0;
	if (!$lta) {
	    #TODO, weird so rather than using the pre-loaded values of the matrix I am reading it from file? Does that mean that I don't pre-load the matrix?  If so do I need to change npp as vocab size below?
	    #calculate npp as the number of co-occurrences
	    open IN, $matrix or die "Cannot open $matrix for input: $!\n";
	    while (my $line = <IN>) {
		chomp $line;
		my ($cui1, $cui2, $num) = split(/\t/,$line);
		$npp += $num;
	    }
	    close IN;
	} else {
	    #calculate npp as the vocabulary size
	    my %uniqueCuis = ();
	    foreach my $keyPair (keys $matrix) {
		(my $cui1, my $cui2) = split ("\t",$keyPair); #TODO, are the keys of $matrix tab seperated CUIs? 
		$uniqueCuis{$cui1} = 1;
		$uniqueCuis{$cui2} = 1;
	    }
	    $npp = scalar keys %uniqueCuis;
	}
	return $npp;
    }
    #calculate NPP from a database
    else {
	#set up database
	my $db = $cuifinder->_getDB(); 
	
	#get npp
	my $npp = 0;
	if (!$lta) {
	    #calculate npp as the number of co-occurrences
	    $npp = shift $db->selectcol_arrayref("select sum(N_11) from N_11"); 
	}
	else {
	    #calculate npp as the vocab size
	    $npp = shift $db->selectcol_arrayref("select count(N_11) from N_11"); 
	}
 
	#return npp
	if($npp <= 0) { errorhandler->_error($pkg, $function, "", 5); } 
	return $npp; 
    }
}



#TODO, these params are weird. Am I assuming stats are precomputed or not?
#  Method to compute n11, n1p, and np1 using precomputed hashes containing n11,np1,n1p values for all cui pairs and cuis.
#  input : $cui1     <- string containing a cui 1
#          $cui2     <- string containing a cui 2
#          $n11Ref   <- ref to a hash of n11 for all cui pairs (n11{"$cui1,$cui2"} = val)
#          $n1pRef   <- ref to a hash of n1p for all cuis (n1p{"$cui"} = val)
#          $n1pRef   <- reference to the sorted array of lines containing np1 values for all cuis (np1{"$cui"} = val)
#  output: \@data - a ref to an array containing n11, n1p, and np1
sub _getStats {
    #grab parameters
    my $self = shift;
    my $cui1 = shift;
    my $cui2 = shift;
    my $n11Ref = shift;
    my $n1pRef = shift;
    my $np1Ref = shift;
    my $function = "_getStats"; 

    #Do error checking
    #check self
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }
    #check parameter exists
    if(!defined $cui1) {
        $errorhandler->_error($pkg, $function, "Error with input variable \$cui1.", 4);
    }
    if(!defined $cui2) {
        $errorhandler->_error($pkg, $function, "Error with input variable \$cui2.", 4);
    }
    #check if cui1 exists
    if(! ($cuifinder->_exists($cui1)) ) {
        $errorhandler->_error($pkg, $function, "Concept ($cui1) does not exist.", 6);
    }    
    #check if cui2 exists
    if(! ($cuifinder->_exists($cui2)) ) {
        $errorhandler->_error($pkg, $function, "Concept ($cui2) does not exist.", 6);
    }   

    #calculate stats from a matrix or DB
    my $n11;
    my $n1p;
    my $np1;
    if (!$matrix) {
	#grab the data from a DB
	$n11 = $self->_getN11($cui1, $cui2); 
	$n1p = $self->_getN1p($cui1); 
	$np1 = $self->_getNp1($cui2); 
    } 
    else {	
        #Get values from the hashes. If the cui pair, or cui 
	# does not exist, set the value to 0
	#get n11
        $n11 = ${$n11Ref}{"$cui1,$cui2"};
	if(!defined $n11) {
	    $n11 = 0;
	}

	#get n1p
        $n1p = ${$n1pRef}{$cui1};
	if(!defined $n1p) {
	    $n1p = 0;
	}

	#get np1
        $np1 = ${$np1Ref}{$cui2};
	if(!defined $np1) {
	    $np1 = 0;
	}

	#add values to account for order not mattering
	if($noOrder) {
	    #sum n11 for $cui1 and $cui2
	    my $val;
	    $val = ${$n11Ref}{"$cui2,$cui1"};
	    if(defined $val) {
		$n11 += $val;
	    }

	    #sum, n1p = n1p + np1 
	    $val = ${$np1Ref}{$cui1};
	    if(defined $val) {
		$n1p += $val;
	    }

	    #sum, np1 = n1p + np1
	    $val = ${$n1pRef}{$cui2};
	    if(defined $val) {
		$np1 += $val;
	    }
	}
    }

    #set up and return the data
    my @data = ($n11, $n1p, $np1); 
    return \@data;
}

# calculates N11 for all cui pairs, and NP1 and N1P for all cuis
# in a single pass of the matrix file
# input: none
# output: \@stats <- an array containing three hash refs (\%n11, \%np1, \%n1p)
#                    \%n11 = a hash ref containing n11 values for all cui
#                            pairs (hash{"$cui1,$cui2"}=n11)
#                    \%n1p = a hash ref containing n1p values for all cuis
#                            hash{$cui}=n1p
#                    \%n1p = a hash ref containing np1 values for all cuis
#                            hash{$cui}=np1
#                    $npp = the value of npp
sub getStatsForAll {
    my $self = shiftl
    (defined $matrix) or die("Error, getStatsForAll requires a matrix file be defined\n");

    #calculate n11, n1p, and np1
    my %n11 = ();
    my %n1p = ();
    my %np1 = ();
    open IN, $matrix or die "Cannot open $matrix for input: $!\n";
    while (my $line = <IN>) {
	chomp $line;
	my ($cui1, $cui2, $num) = split /\t/, $line;
	$n11{"$cui1,$cui2"} = $num;
	$n1p{"$cui1"} += $num;
	$np1{"$cui2"} += $num;
    }
    close IN;

    my $npp = $self->getNpp();

    return (\%n11, \%n1p, \%np1, $npp);
}

#  Gets contingency table values with concept expansion
#  input : $cui1 <- string containing a cui 1
#          $cui2 <- string containing a cui 2
#          $n11Ref   <- ref to a hash of cui pairs containing n11  ($n11s{$cui1,$cui2} = $num)
#          $n1pRef   <- ref to a hash of n1p values for all cuis ($n1ps{$cui} = $num)
#          $n1pRef   <- ref to a hash of np1 values for all cuis ($np1s{$cui} = $num)
#  output: reference to @data = (n_11, n_1p, n_p1)
sub _getStats_withConceptExpansion{
    #grab parameters
    my $self = shift;
    my $cui1 = shift;
    my $cui2 = shift;
    my $n11Ref = shift;
    my $n1pRef = shift;
    my $np1Ref = shift;
    my $function = "_getStats_withConceptExpansion"; 
    
    #error checking
    #  check self
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }
    #check parameter exists
    if(!defined $cui1) {
        $errorhandler->_error($pkg, $function, "Error with input variable \$cui1.", 4);
    }
    if(!defined $cui2) {
        $errorhandler->_error($pkg, $function, "Error with input variable \$cui2.", 4);
    }
    #check if concept 1 exists
    if(! ($cuifinder->_exists($cui1)) ) {
        $errorhandler->_error($pkg, $function, "Concept ($cui1) does not exist.", 6);
    }    
    #check if concept 2 exists
    if(! ($cuifinder->_exists($cui2)) ) {
        $errorhandler->_error($pkg, $function, "Concept ($cui2) does not exist.", 6);
    }    

    #get descendants of each cui, and itself
    my @concepts1 =@{&_findDescendants($cui1)};
    push @concepts1, $cui1;

    my @concepts2 = @{&_findDescendants($cui2)};
    push @concepts2, $cui2;

    #get n11
    my $queryString;
    my $n11 = 0;
    my $n1p = 0;
    my $np1 = 0;

###
    #Get stats using a matrix or database
    if($matrix) {  #gets stats when using a matrix file
	#calculate n1p as the sum of n1p's for all concepts1
	foreach my $cui (@concepts1) {
	    my $num = ${$n1pRef}{$cui};
	    if(defined $num) {
		$n1p += $num;
	    }
	}
	#calculate np1 as the sum of np1's for all concepts2
	foreach my $cui (@concepts2) {
	    my $num = ${$np1Ref}{$cui};
	    if (defined $num) {
		$np1 += $num;
	    }
	}
	#calculate n11 as the sum n11s for all combinations of 
        # concepts1, concepts2 (order matters, cui1 must be first)
	foreach my $c1 (@concepts1) {
	    foreach my $c2 (@concepts2) {
		my $num = ${$n11Ref}{"$c1,$c2"};
		if(defined $num) {
		    $n11 += $num;
		}
	    }
	}

###
	#update values if order does not matter
	if($noOrder) {
	    #add all np1's to n1p
	    foreach my $cui (@concepts1) {
		my $num = ${$np1Ref}{$cui};
		if(defined $num) {
		    $n1p += $num;
		}
	    }
	    #add all n1p's to np1s
	    foreach my $cui (@concepts2) {
		my $num = ${$n1pRef}{$cui};
		if (defined $num) {
		    $np1 += $num;
		}
	    }
	    #add all n11's, now with the order reversed
	    foreach my $c1 (@concepts1) {
		foreach my $c2 (@concepts2) {
		    my $num = ${$n11Ref}{"$c2,$c1"};
		    if(defined $num) {
			$n11 += $num;
		    }
		}
	    }
	}
    
    }
###
    else  { #Get stats from a DB
	#set up database
	my $db = $cuifinder->_getDB(); 

	#build up query strings and query the DB
	if ($noOrder) {
	    #build up a query string for n11, where order doesn't matter
	    $queryString = "select SUM(n_11) from N_11 where ((cui_2 = '$cui1' ";

	    #set all cui_2 to be cui1's
	    foreach my $cui (@concepts1) {
		$queryString .= "or cui_2 = '$cui' ";
	    }

	    #set all cui_1 to be cui2's
	    $queryString .= ") and (cui_1 = '$cui2' ";
	    foreach my $cui (@concepts2) {
		$queryString .= "or cui_1 = '$cui' ";
	    }

	    #set all cui_1 to be cui1's
	    $queryString .= ")) or ((cui_1 = '$cui1' ";
	    foreach my $cui (@concepts1) {
		$queryString .= "or cui_1 = '$cui' ";
	    }

	    #set all cui_2 to be cui2's
	    $queryString .= ") and (cui_2 = '$cui2' ";
	    foreach my $cui (@concepts2) {
		$queryString .= "or cui_2 = '$cui' ";
	    }
	    $queryString .= "));";
	    #the end result is query of any combination of (any cui1,any cui2) OR (any cui2,any cui1)

	    #query the DB and get n11
	    $n11 = shift @{$db->selectcol_arrayref($queryString)};
	}
	else
	{
	    #build a query string for n11, where order does matter
	    $queryString = "select SUM(n_11) from N_11 where (cui_1 = '$cui1' ";
	    
            #set all cui1's
	    foreach my $cui (@concepts1) {
		$queryString .= "or cui_1 = '$cui' ";
	    }

	    #set all cui2's
	    $queryString .= ") and (cui_2 = '$cui2' ";
	    foreach my $cui (@concepts2) {
		$queryString .= "or cui_2 = '$cui' ";
	    }
	    $queryString .= ");";

	    #query the DB and get n11
	    $n11 = shift @{$db->selectcol_arrayref($queryString)};
	}

	#get n1p
	if($noOrder) {
	    #build a query where order doesn't matter
	    $queryString = "select SUM(n_11) from N_11 where (cui_2 = '$cui1' ";
	    foreach my $cui (@concepts1) {
		$queryString .= "or cui_2 = '$cui' ";
	    }

	    $queryString .= ") or (cui_1 = '$cui1' ";
	    foreach my $cui (@concepts1) {
		$queryString .= "or cui_1 = '$cui' ";
	    }
	    $queryString .= ");";

	    #query the db to retreive n1p
	    $n1p = shift @{$db->selectcol_arrayref($queryString)};
	}
	else {
	    #build a query where order does matter
	    $queryString = "select SUM(n_11) from N_11 where (cui_1 = '$cui1' ";
	    foreach my $cui (@concepts1) {
		$queryString .= "or cui_1 = '$cui' ";
	    }
	    $queryString .= ");";

	    #query the db to retrive n1p
	    $n1p = shift @{$db->selectcol_arrayref($queryString)};
	}

	#get np1
	if($noOrder) {
	    #build a query to get np1 where order does not matter
	    $queryString = "select SUM(n_11) from N_11 where (cui_2 = '$cui2' ";
	    foreach my $cui (@concepts2) {
		$queryString .= "or cui_2 = '$cui' ";
	    }

	    $queryString .= ") or (cui_1 = '$cui2' ";
	    foreach my $cui (@concepts2) {
		$queryString .= "or cui_1 = '$cui' ";
	    }
	    $queryString .= ");";

	    #query the DB to get np1
	    $np1 = shift @{$db->selectcol_arrayref($queryString)};
	}
	else {
	    #build a query to get np1 where order matters
	    $queryString = "select SUM(n_11) from N_11 where (cui_2 = '$cui2' ";
	    foreach my $cui (@concepts2) {
		$queryString .= "or cui_2 = '$cui' ";
	    }
	    $queryString .= ");";

	    #query the DB and get np1
	    $np1 = shift @{$db->selectcol_arrayref($queryString)};
	}
    }

    #return the values
    my @data = ($n11, $n1p, $np1);
    return \@data;
}

#TODO, I don't like these optional parameters that are only for the matrix option
#  Gets contingency table values for LTA
#  input : $concept1 <- string containing a cui 1
#          $concept2 <- string containing a cui 2
#          $npp      <- the value of npp, the vocabulary size
#          $n1pRef   <- reference to a hash containing implicit co-occurrence data for cui 1 (only for matrix option)
#          $np1Ref   <- reference to a hash containing implicit co-occurrence data for cui 2 (only for matrix option)
#  output: reference to @data = (n_11, n_1p, n_p1)
sub _getStats_LTA {
    #grab input
    my $self = shift;
    my $concept1 = shift;
    my $concept2 = shift;
    my $npp = shift;
    my $n1pRef = shift;
    my $np1Ref = shift;

    #get values with or without a matrix input
    if($matrix) {
	#find the number of linking terms by creating a list of 
	# cui-cooccurrences for each input cui (concept1, concept2)
	# then use count of the overlap between those terms, and 
	# their counts to find n11,n1p,np1 values.

	#TODO, so n1pRef is a list of CUIs that co-occurr with cui1?
	#create cui1 and cui2 co-occurring cui lists
	my %cui1Data = ();
	if (defined ${$n1pRef}{$concept1}) {
	    foreach my $cui (split /,/,${$n1pRef}{$concept1}) {
		$cui1Data{$cui} = 1;
	    }
	}
	my %cui2Data = ();
	if (defined ${$np1Ref}{$concept2}) {
	    foreach my $cui (split /,/,${$np1Ref}{$concept2}) {
		$cui2Data{$cui} = 1;
	    }
	}

	#add extra CUIs to lists if order doesn't matter
	if ($noOrder) {
	    if (defined ${$np1Ref}{$concept1}) {
		foreach my $cui (split /,/,${$np1Ref}{$concept1}){
		    $cui1Data{$cui} = 1;
		}
	    }
	    if (defined ${$n1pRef}{$concept2}) {
		foreach my $cui (split /,/,${$n1pRef}{$concept2}) {
		    $cui2Data{$cui} = 1;
		}
	    }
	}

	#calculate n11 as the overlap of co-occurring cuis
	my $n11 = 0;
	if(scalar keys %cui1Data > 0 && scalar keys %cui2Data > 0) {
	    #Find number of CUIs that co-occur with both CUI 1 and CUI 2
	    foreach my $cui (keys %cui1Data) {
		if (exists $cui2Data{$cui}) {
		    $n11++;
		}
	    }
	}

	#set and return values
	my $n1p = scalar keys %cui1Data;
	my $np1 = scalar keys %cui2Data;
	
	my @data = ($n11, $n1p, $np1);
	return \@data;
    }
    else {
	#A matrix was not passed in, so find n11, n1p, np1 (ordered or not ordered)
	# using DB queries.

	#set up database
	my $db = $cuifinder->_getDB(); 

	#find cuis that co-occur with cui1 and cui2 
	my %cui1Data = ();
	my %cui2Data = ();
	if($noOrder)
	{
	    #get a list of co-occurring CUIs for concepts 1&2, where order doesn't matter

	    #Create a list of CUIs that co-occurr with $concept1
	    #get a list of cuiA's and cuiB's in a (cuiA,cuiB) co-occurrence pair $concept1 is in the pair
	    my @cuisB = @{$db->selectcol_arrayref("SELECT cui_2 FROM N_11 WHERE cui_1 = '$concept1' AND n_11 > 0;")};
	    my @cuisA = @{$db->selectcol_arrayref("SELECT cui_1 FROM N_11 WHERE cui_2 = '$concept1' AND n_11 > 0;")};
	    
	    #create cui1Data from cuisA, and cuisB = the co-occurring CUIs with $concept1
	    foreach my $key(@cuisA) {
		$cui1Data{$key} = 1;
	    }
	    foreach my $key(@cuisB) {
		$cui1Data{$key} = 1;
	    }

	    #Create a list of CUIs that co-occurr with $concept1
	    #get a list of cuiA's and cuiB's in a (cuiA,cuiB) co-occurrence pair $concept1 is in the pair
	    @cuisB = @{$db->selectcol_arrayref("SELECT cui_2 FROM N_11 WHERE cui_1 = '$concept2' AND n_11 > 0;")};
	    @cuisA = @{$db->selectcol_arrayref("SELECT cui_1 FROM N_11 WHERE cui_2 = '$concept2' AND n_11 > 0;")};

	    #create cui1Data from cuisA, and cuisB = the co-occurring CUIs with $concept1
	    foreach my $key(@cuisA) {
		$cui2Data{$key} = 1;
	    }
	    foreach my $key(@cuisB) {
		$cui2Data{$key} = 1;
	    }
	}
	else {	    
	    #get a list of co-occurring CUIs for concepts 1&2, where order matters
	    #get cui1Data
	    my @cuis = @{$db->selectcol_arrayref("SELECT cui_2 FROM N_11 WHERE cui_1 = '$concept1' AND n_11 > 0;")}; #same as above
	    foreach my $cui (@cuis) {
		$cui1Data{$cui} = 1;
	    }

	    #get cui2Data
	    my @cuis2 = @{$db->selectcol_arrayref("SELECT cui_1 FROM N_11 WHERE cui_2 = '$concept2' AND n_11 > 0;")};  #same as above
	    foreach my $cui (@cuis) {
		$cui2Data{$cui} = 1;
	    }
	}
	my $n1p = scalar keys %cui1Data;
	my $np1 = scalar keys %cui2Data;

	#find n11 = the number of cuis in common
	my $n11 = 0;
	foreach my $cui (keys %cui1Data) {
	    if(exists $cui2Data{$cui}){
		$n11++;
	    }
	}
	
	#return the values
	my @data = ($n11, $n1p, $np1);
	return \@data;
    }
}

#  Gets contingency table values for LTA with concept expansion
#  input : $concept1 <- string containing a cui 1
#          $concept2 <- string containing a cui 2
#          $n1pRef   <- reference to a hash containing implicit co-occurrence data for cui1 (only for matrix option)
#          $np1Ref   <- reference to a hash containing implicit co-occurrence data for cui2 (only for matrix option)
#  output: reference to @data = (n_11, n_1p, n_p1)
sub _getStats_LTAWithConceptExpansion {
    #grab parameters
    my $self = shift;
    my $cui1 = shift;
    my $cui2 = shift;
    my $n1pRef = shift;
    my $np1Ref = shift;
    my $function = "_getStats_LTAWithConceptExpansion"; 
    
    #check self
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #get the expanded concepts cuis 1 and 2
    my @concepts1 =@{_findDescendants($cui1)};
    push @concepts1, $cui1;
    my @concepts2 = @{_findDescendants($cui2)};
    push @concepts2, $cui2;

    #get values with or without a matrix input
    if($matrix) {
	#Get values from a matrix
	
	#get lists of explicitly co-occurring CUIs for each concept
	#add trailing cui co-occurrences to cui1Data
	my %cui1Data;
	foreach my $desc1 (@concepts1) {
	    if (defined ${$n1pRef}{$desc1}) {
		foreach my $cui (split /,/,${$n1pRef}{$desc1}) {
		    $cui1Data{$cui} = 1;
		}
	    }
	}
	#add leading cui co-occurrences to cui2Data
	my %cui2Data;
	foreach my $desc2 (@concepts2) {
	    if (defined ${$np1Ref}{$desc2}) {
		foreach my $cui (split /,/,${$np1Ref}{$desc2}) {
		    $cui2Data{$cui} = 1;
		}
	    }
	}

	#add more CUIs if order doesn't matter
	if ($noOrder) {
	    #add leading co-occurring cuis to cui1Data
	    foreach my $desc1 (@concepts1) {
		if (defined ${$np1Ref}{$desc1}) {
		    foreach my $cui (split /,/,${$np1Ref}{$desc1}) {
			$cui1Data{$cui} = 1;
		    }
		}
	    }
	    #add trailling co-occurring cuis to cui2Data
	    foreach my $desc2 (@concepts2) {
		if(defined ${$n1pRef}{$desc2}) {
		    foreach my $cui (split /,/,${$n1pRef}{$desc2}) {
			$cui2Data{$cui} = 1;
		    }
		}
	    }
	}

	#checks to see if concept1 and concept2 actually occur with anything; if so, calculate data
	my $n11 = 0;
	my $n1p = 0;
	my $np1 = 0;
	if (scalar keys %cui1Data > 0 && scalar keys %cui2Data > 0) {
	    #Find number of CUIs that co-occur with both CUI 1 and CUI 2
	    foreach my $c (keys %cui1Data) {
		if (exists $cui2Data{$c}) {
		    $n11++;
		}
	    }
	    $n1p = scalar keys %cui1Data;
	    $np1 = scalar keys %cui2Data;
	}
	#one of the two concepts doesn't occur with anything.  Adjust accordingly
	else
	{
	    if(scalar keys %cui1Data == 0) {
		$n1p = 0;
	    }
	    if(scalar keys %cui2Data == 0) {
		$np1 = 0;
	    }
	    $n11 = 0;
	}
	
	#return the values
	return ($n11, $n1p, $np1)
    }
    else {
	#get values from a DB (rather than matrix)
	#set up database
	my $db = $cuifinder->_getDB(); 

	#get co-occurring CUIs with or without order mattering
	my %cui1Data = ();
	my %cui2Data = ();
	if($noOrder) {
	    #Build query strings and query DB for co-occurring CUIs with CUI 1
	    #get n11 tables, where concept 1 is the leading cui
	    my $query = "SELECT N_11.cui_2 FROM N_11 WHERE (N_11.cui_1 = '$cui1' ";
	    foreach my $desc (@concepts1) {
		$query .= "OR N_11.cui_1 = '$desc' ";
	    }
	    $query .= ") AND N_11.n_11 > 0;";
	    my @trailingCUIs = @{$db->selectcol_arrayref($query)};

	    #get n11 tables, where concept 1 is the trailing cui
	    $query =  "SELECT N_11.cui_1 FROM N_11 WHERE (N_11.cui_2 = '$cui1' ";
	    foreach my $desc (@concepts1) {
		$query .= "OR N_11.cui_2 = '$desc' ";
	    }
	    $query .= ") AND N_11.n_11 > 0;";
	    my @leadingCUIs = @{$db->selectcol_arrayref($query)};
	    
	    #combine the trailing and leading CUIs to a list of co-occurring CUIs with cui1
	    foreach my $cui (@trailingCUIs) {
		$cui1Data{$cui} = 1;
	    }
	    foreach my $cui (@leadingCUIs) {
		$cui1Data{$cui} = 1;
	    }

	    #Build query strings and query DB for co-occurring CUIs with CUI 1
	    #get n11 tables, where concept 1 is the leading cui
	    $query = "SELECT N_11.cui_2 FROM N_11 WHERE (N_11.cui_1 = '$cui2' ";
	    foreach my $desc (@concepts1) {
		$query .= "OR N_11.cui_1 = '$desc' ";
	    }
	    $query .= ") AND N_11.n_11 > 0;";
	    @trailingCUIs = @{$db->selectcol_arrayref($query)};

	    #get n11 tables, where concept 1 is the trailing cui
	    $query =  "SELECT N_11.cui_1 FROM N_11 WHERE (N_11.cui_2 = '$cui2' ";
	    foreach my $desc (@concepts1) {
		$query .= "OR N_11.cui_2 = '$desc' ";
	    }
	    $query .= ") AND N_11.n_11 > 0;";
	    @leadingCUIs = @{$db->selectcol_arrayref($query)};
	    
	    #combine the trailing and leading CUIs to a list of co-occurring CUIs with cui1
	    foreach my $cui (@trailingCUIs) {
		$cui2Data{$cui} = 1;
	    }
	    foreach my $cui (@leadingCUIs) {
		$cui2Data{$cui} = 1;
	    }
	}
	else { #order does matter
	    #Build query strings and query DB for co-occurring CUIs with CUI 1
	    #get n11 tables, where concept 1 is the leading cuir 
	    my $query = "SELECT N_11.cui_2 FROM N_11 WHERE (N_11.cui_1 = '$cui1' ";
	    foreach my $desc (@concepts1) {
		$query .= "OR N_11.cui_1 = '$desc' ";
	    }
	    $query .= ") AND N_11.n_11 > 0;";
	    my @trailingCUIs = @{$db->selectcol_arrayref($query)};

	    #make trailing CUIs a hash (co-occurring CUIs with cui1)
	    foreach my $cui (@trailingCUIs) {
		$cui1Data{$cui} = 1;
	    }
	    
            #get n11 tables, where concept 2 is the trailing cui
	    $query =  "SELECT N_11.cui_1 FROM N_11 WHERE (N_11.cui_2 = '$cui2' ";
	    foreach my $desc (@concepts1) {
		$query .= "OR N_11.cui_2 = '$desc' ";
	    }
	    $query .= ") AND N_11.n_11 > 0;";
	    my @leadingCUIs = @{$db->selectcol_arrayref($query)};

	    #make leading CUIs a hash (co-occurring CUIs with cui2)
	    foreach my $cui (@leadingCUIs) {
		$cui2Data{$cui} = 1;
	    }

	}

	#find n11 = the intersecting CUIS (cuis that co-occur with cui1 and cui2)
	my $n11 = 0;
	foreach my $cui (keys %cui1Data) {
	    if(exists $cui2Data{$cui}) {
		$n11++;
	    }
	}
	#set up and return data
	my $n1p = scalar keys %cui1Data;
	my $np1 = scalar keys %cui2Data;

	#return the values
	return ($n11, $n1p, $np1)
    }
}

#  Method to retrieve descendants of a cui
#  input : $cui <- string containing a cui 
#  output: reference to @descendants, the descendants of the given cui
sub _findDescendants
{
    my $cui = shift;
    my $hashref = $umls->findDescendants($cui);
    my @descendants = (keys %{$hashref});
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
    
    #TODO, these functions are problematic because I think they assume n1p, np1, are pre-computed for all, which isnt necassarily the case
    # get frequency and marginal totals
    my $valid = -1;
    if($conceptExpansion && !$lta) {
	$valid = $self->_getStats_withConceptExpansion($concept1, $concept2);
    }
    elsif($lta && !$conceptExpansion) {
	$valid = $self->_getStats_LTA($concept1, $concept2);
    }
    elsif($lta && $conceptExpansion) {
	$valid = $self->_getStats_LTAWithConceptExpansion($concept1, $concept2);
    }
    else {
	$valid = $self->_getStats($concept1, $concept2);
    }

    #error checking for n11, n1p, np1
    if($valid == -1){
        return -1;
    }
    my @data = @{$valid};
    
    my $n11 = $data[0];
    my $n1p = $data[1];
    my $np1 = $data[2];
    my $npp = $self->_getNPP();
    if(!defined $n11 || !defined $n1p || !defined $np1){
	return -1;
    }
    
    #calculate the statistic and return it
    my $stat = $self->_calculateStatisticFromContingencyTable(
	$n11, $n1p, $np1, $npp, $statistic);
    return $stat;
}

# calculate a statistic (association score) using pre-computed contingency table values
# input:  $n11 <- n11 for the cui pair
#         $npp <- npp for the dataset
#         $n1p <- n1p for the cui pair
#         $np1 <- np1 for the cui pair
#         $statistic <- the string specifying the stat to calc
# output: the statistic (Association score) between the two concepts
sub calculateStatisticFromContingencyTable {
    my $self = shift;
    my $n11 = shift;
    my $npp = shift;
    my $n1p = shift;
    my $np1 = shift;
    my $statistic = shift;

    # set frequency and marginal totals
    my %values = (n11=>$n11, 
		  n1p=>$n1p, 
		  np1=>$np1, 
		  npp=>$npp); 
    
    #return cannot compute, or 0
    if($n11 < 0 || $n1p < 0 || $np1 < 0) { 	
	return -1.000; 
    }
    if($n11 == 0) { 
	return 0.000;
    }
    
    #set default statistic
    if(! defined $statistic) { $statistic = "tscore";  }
    
    #set statistic module
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
	    next; # if warning, dont save the statistic value just computed
	}
    }

    #return statistic to given precision.  if no precision given, default is 4
    my $floatFormat = join '', '%', '.', $precision, 'f';
    my $statScore = sprintf $floatFormat, $statisticValue;

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
my $freq = $statfinder->_getN11($cui1, $cui2); 

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
