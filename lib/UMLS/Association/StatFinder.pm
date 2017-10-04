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

my $pkg = "UMLS::Association::StatFinder";

#  debug variables
#local(*DEBUG_FILE);

#NOTE: every global variable is followed by a _G with the 
# exception of debugm error handler, and constants which are all caps
#  global variables
my $DEFAULT_STATISTIC = "tscore"; #association measure to use
my $debug     = 0; #in debug mode or not
my $umls_G = undef; #UMLS interface instance
my $precision_G = 4; #precision of the output
my $cuifinder_G      = ""; 

#global options variables
my $conceptExpansion_G = 0; #1 or 0 if using concept expansion or not
my $lta_G = 0; #1 or 0 is using lta or not
my $noOrder_G = 0; #1 or 0 if noOrder is enabled or not
my $matrix_G = 0; #matrix file name is using a matrix file rather than DB

######################################################################
#                 Initialization Functions
######################################################################
#  method to create a new UMLS::Association::StatFinder object
#  input : $params <- reference to hash of database parameters
#          $handler <- reference to cuifinder object 
#  output: $self
sub new {
    #grab params and create self
    my $self = {};
    my $className = shift;
    my $params = shift;
    my $handler = shift; 

    #bless the object.
    bless($self, $className);

    #initialize error handler
    $errorhandler = UMLS::Association::ErrorHandler->new();
    if(! defined $errorhandler) {
        print STDERR "The error handler did not get passed properly.\n";
        exit;
    }

    # initialize the object.
    $cuifinder_G = $handler; 
    $debug = 0; 
    $self->_initialize($params);
    return $self;
}

#  method to initialize the UMLS::Association::StatFinder object.
#  input : $parameters <- reference to a hash of database parameters
#  output: none
sub _initialize {
    #grab parameters
    my $self = shift;
    my $params = shift;
    my %params = %{$params};

    #set global variables using option hash
    $umls_G = $params{'umls'};
    $conceptExpansion_G = $params{'conceptexpansion'};
    $lta_G = $params{'lta'};
    $noOrder_G = $params{'noorder'};
    $matrix_G = $params{'matrix'};

     #set precision
    if(defined $params{'precision'}) {
	$precision_G = $params{'precision'};
    }

    #error checking
    my $function = "_initialize";
    &_debug($function);
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #require UMLS::Interface instance with concept expansion
    if ($conceptexpansion_G && !(defined $umls_G)) {
	$umls_
    }
}

sub _debug {
    my $function = shift;
    if($debug) { print STDERR "In UMLS::Association::StatFinder::$function\n"; }
}

######################################################################
#      Public Functions to get Association Scores
######################################################################

# calculate a contingency table values and an associationScore
# for the cui pair
# input:  $cui1  <- the cui of the first concept 
#         $cui2  <- the cui of the second concept
#         $statistic <- the string specifying the stat to calc
# output: the statistic (Association score) between the two concepts
sub calculateAssociation { 
    #grab parameters
    my $self = shift;
    my $cui1 = shift; 
    my $cui2 = shift; 
    my $statistic = shift; 

    #error checking
    my $function = "calculateAssociation"; 
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #create the sets of cuis, default is just a single cui
    my @cuis1 = ();
    push @cuis1, $cui1;
    my @cuis2 = ();
    push @cuis2, $cui2;

    #add descendants to the sets if needed
    if ($conceptExpansion_G) {
	push @cuis1, @{&_findDescendants($cui1)};
	push @cuis2, @{&_findDescendants($cui2)};	
    }

    #calculate n11, n1p, np1, npp using a matrix or DB
    # and according to the method of various other options
    my $valid = -1; #if invalid, valid stays -1, else if becomes a hash ref
    if ($lta_G) {
	$valid = $self->_getStats_LTA(\@cuis1, \@cuis2);
    }
    else {
	if ($matrix_G) {
	    $valid = $self->_getStats_matrix(\@cuis1, \@cuis2);
	}
	else {
	    $valid = $self->_getStats_DB(\@cuis1, \@cuis2);
	}
    }
    
    #error checking for n11, n1p, np1, npp
    if($valid == -1){
	return -1;
    }
    my @data = @{$valid};

    my $n11 = $data[0];
    my $n1p = $data[1];
    my $np1 = $data[2];
    my $npp = $data[3];
    if(!defined $n11 || !defined $n1p || !defined $np1){
	return -1;
    }
    
    #calculate the statistic and return it
    my $stat = $self->calculateAssociationFromValues(
	$n11, $n1p, $np1, $npp, $statistic);

    return $stat;
}

# calculates an association score from the provided values
# NOTE: Please be careful when writing code that uses this
# method. Results may become inconsistent if you don't check
# that CUIs occur in the hierarchy before calling
# e.g. C0009951 does not occur in the SNOMEDCT Hierarchy but
# it likely occurs in the association database so if not check
# is made an association score will be calculate for it, but it has not
# been done in reported results from this application
# input:  $n11 <- n11 for the cui pair
#         $npp <- npp for the dataset
#         $n1p <- n1p for the cui pair
#         $np1 <- np1 for the cui pair
#         $statistic <- the string specifying the stat to calc
# output: the statistic (association score) between the two concepts
sub calculateAssociationFromValues {
    #grab parameters
    my $self = shift;
    my $n11 = shift;
    my $n1p = shift;
    my $np1 = shift;
    my $npp = shift;
    my $statistic = shift;

    #set frequency and marginal totals
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
    if(!defined $statistic) { 
	$statistic = $DEFAULT_STATISTIC;  
    }

    #set statistic module (Text::NSP)
    my $includename = ""; my $usename = "";  my $ngram = 2; #TODO, what is this ngram parameter
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
    
    # import module
    require $includename;
    import $usename;
    
    # get statistics (From NSP package)
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
    my $floatFormat = join '', '%', '.', $precision_G, 'f';
    my $statScore = sprintf $floatFormat, $statisticValue;

    return $statScore; 
}

######################################################################
# functions to get statistical information about the cuis using a DB
######################################################################

# gets N11, N1P, NP1, NPP for the two sets of CUIs using a database
#  input : $cuis1Ref <- ref to an array of the first cuis in a set of cui pairs
#          $cuis2Ref <- ref to an array of the second cuis in a set of cui pairs
#  output: $\@data  <- array ref containing four values: $n11, $n1p, $np1, and 
#                      $npp for the sets of cui pairs
sub _getStats_DB {
    #grab parameters
    my $self = shift;
    my $cuis1Ref = shift;
    my $cuis2Ref = shift;
    
    #error checking
    my $function = "_getStats_DB"; 
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #grab the data from a DB
    my $n11 = $self->_getN11_DB($cuis1Ref, $cuis2Ref); 
    my $n1p = $self->_getN1p_DB($cuis1Ref); 
    my $np1 = $self->_getNp1_DB($cuis2Ref); 
    my $npp = $self->_getNpp_DB();

    #return the data
    my @data = ($n11, $n1p, $np1, $npp);
    return  \@data;
}

#  Gets N11 of the cui pair using a database
#  input:  $cuis1Ref <- ref to an array of the first cuis in a set of cui pairs
#          $cuis2Ref <- ref to an array of the second cuis in a set of cui pairs
#  output: $n11  <- n11 of cui sets 
sub _getN11_DB {
    #grab parameters
    my $self = shift;
    my $cuis1Ref = shift;
    my $cuis2Ref = shift;

    #error checking
    my $function = "_getN11";
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #set up database
    my $db = $cuifinder_G->_getDB(); 
    
    #build a query string for n11
    my $firstCui = shift @{$cuis1Ref};
    my $queryString = "select SUM(n_11) from N_11 where ((cui_1 = '$firstCui' ";
    foreach my $cui (@{$cuis1Ref}) {
	$queryString .= "or cui_1 = '$cui' ";
    }
    unshift @{$cuis1Ref}, $firstCui;

    #set all cui2's
    $firstCui = shift @{$cuis2Ref};
    $queryString .= ") and (cui_2 = '$firstCui' ";
    foreach my $cui (@{$cuis2Ref}) {
	$queryString .= "or cui_2 = '$cui' ";
    }
    unshift @{$cuis2Ref}, $firstCui;

    #finalize the query string
    if ($noOrder_G) {
	#swap the positions of the cuis
	$firstCui = shift @{$cuis2Ref};
	$queryString .= ")) or ((cui_1 = '$firstCui' ";
	foreach my $cui (@{$cuis2Ref}) {
	    $queryString .= "or cui_1 = '$cui' ";
	}
	unshift @{$cuis2Ref}, $firstCui;

	$firstCui = shift @{$cuis1Ref};
	$queryString .= ") and (cui_2 = '$firstCui' ";
	foreach my $cui (@{$cuis1Ref}) {
	    $queryString .= "or cui_2 = '$cui' ";
	}
	unshift @{$cuis1Ref}, $firstCui;
    }
    $queryString .= "));";
    
    #query the DB and return n11
    my $n11 = shift @{$db->selectcol_arrayref($queryString)};
    if (!defined $n11) {
	$n11 = 0;
    }
    return $n11;
}

#  Method to return the np1 of a concept using a database
#  input : $cuis2Ref <- ref to an array of the second cuis in a set of cui pairs
#  output: $np1 <- number of times the cuis2Ref set occurs in second bigram position
sub _getNp1_DB {
    my $self = shift;
    my $cuis2Ref = shift; 
    
    #error checking
    my $function = "_getNp1_DB";
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #set up database
    my $db = $cuifinder_G->_getDB(); 

    #build a query string for all where cui2's are in the second position
    my $firstCui = shift @{$cuis2Ref};
    my $queryString = "select SUM(n_11) from N_11 where (cui_2 = '$firstCui' ";
    foreach my $cui (@{$cuis2Ref}) {
	$queryString .= "or cui_2 = '$cui' ";
    }
    unshift @{$cuis2Ref}, $firstCui;

    #finalize the query string
    if ($noOrder_G) {
	#add where cui2 is in the first position
	$firstCui = shift @{$cuis2Ref};
	$queryString .= ") or (cui_1 = '$firstCui' ";
	foreach my $cui (@{$cuis2Ref}) {
	    $queryString .= "or cui_1 = '$cui' ";
	}
	unshift @{$cuis2Ref}, $firstCui;
    }
    $queryString .= ");";
  
    #query the db to retrive np1
    my $np1 = shift @{$db->selectcol_arrayref($queryString)};
    if (!defined $np1) {
	$np1 = -1;
    }
    return $np1;
}

#  Method to return the n1p of a concept from a database
#  input : $cuis1Ref <- ref to an array of the first cuis in a set of cui pairs
#  output: $n1p <- number of times cuis in cuis1 set occurs in first bigram position
sub _getN1p_DB {
    my $self = shift;
    my $cuis1Ref = shift; 

    #error checking
    my $function = "_getN1p";
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #set up database
    my $db = $cuifinder_G->_getDB(); 
    
    #build the query string for all where cui1's are in the first position
    my $firstCui = shift @{$cuis1Ref};
    my $queryString = "select SUM(n_11) from N_11 where (cui_1 = '$firstCui' ";
    foreach my $cui (@{$cuis1Ref}) {
	$queryString .= "or cui_1 = '$cui' ";
    }
    unshift @{$cuis1Ref}, $firstCui;

    #finalize the query string
    if ($noOrder_G) {
	#add where cui1 is in the second position
	$firstCui = shift @{$cuis1Ref};
	$queryString .= ") or (cui_2 = '$firstCui' ";
	foreach my $cui (@{$cuis1Ref}) {
	    $queryString .= "or cui_2 = '$cui' ";
	}
	unshift @{$cuis1Ref}, $firstCui;
    }
    $queryString .= ");";

    #query the db to retrive n1p
    my $n1p = shift @{$db->selectcol_arrayref($queryString)};
    if (!defined $n1p) {
        $n1p = -1;
    }
    return $n1p;
}

#  Method to calculate npp from a DB
#  input : none
#  output: $npp
sub _getNpp_DB {
    my $self = shift;
    
    #error checking
    my $function = "getNpp_DB";
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #set up database
    my $db = $cuifinder_G->_getDB(); 
    
    #get npp, the number of co-occurrences
    my $npp = shift $db->selectcol_arrayref("select sum(N_11) from N_11"); 

    #update $npp for noOrder, since Cuis can be trailing or leading its 2x ordered npp
    if ($noOrder_G) {
	$npp *= 2;
    }

    #return npp
    if($npp <= 0) { errorhandler->_error($pkg, $function, "", 5); } 
    return $npp; 
}

########################################################################
# functions to get statistical information about the cuis using a matrix 
########################################################################

# gets N11, N1P, NP1, NPP for the two sets of CUIs using a matrix
#  input : $cuis1Ref <- ref to an array of the first cuis in a set of cui pairs
#          $cuis2Ref <- ref to an array of the second cuis in a set of cui pairs
#  output: $\@data  <- array ref containing four values: $n11, $n1p, $np1, and 
#                      $npp for the sets of cui pairs
sub _getStats_matrix {
    #grab parameters
    my $self = shift;
    my $cuis1Ref = shift;
    my $cuis2Ref = shift;
    
    #error checking
    my $function = "_getStats_DB"; 
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #get all observed counts for all cuis in the term pairs
    my $countsRef = $self->_getObservedCounts_matrix($cuis1Ref, $cuis2Ref);

    #calculate stats for the term pair
    my $n11 = $self->_getN11_matrix($cuis1Ref, $cuis2Ref, ${$countsRef}[0]); 
    my $n1p = $self->_getN1p_matrix($cuis1Ref, ${$countsRef}[1], ${$countsRef}[2]); 
    my $np1 = $self->_getNp1_matrix($cuis2Ref, ${$countsRef}[1], ${$countsRef}[2]); 
    my $npp = ${$countsRef}[3];

    #update $npp for noOrder, since Cuis can be trailing or leading its 2x ordered npp
    if ($noOrder_G) {
	$npp *= 2;
    }

    #return the data
    my @data = ($n11, $n1p, $np1, $npp);
    return \@data;
}

#computes the observed counts for all combinations of the cuis passed in
#doing this in a single function makes it so all values can be computed with a 
#single pass of the input file, making execution time much faster
#  input : $cuis1Ref <- ref to an array of the first cuis in a set of cui pairs
#          $cuis2Ref <- ref to an array of the second cuis in a set of cui pairs
#  output: $\@counts  <- array ref containing four sets of values: 
#                      \%n11, \%n1p, \%np1, and $npp for the cui pairs
#                      hashes are indexed: $n11{"$cui1,$cui2"}, $n1p{$cui},
#                                          $np1{$cui}
sub _getObservedCounts_matrix {
    my $self = shift;
    my $cuis1Ref = shift;
    my $cuis2Ref = shift;

    #convert cui arrays to hashes, makes looping thru
    # the file faster
    my %cuis1 = ();
    foreach my $cui(@{$cuis1Ref}) {
	$cuis1{$cui} = 1;
    }
    my %cuis2 = ();
    foreach my $cui(@{$cuis2Ref}) {
	$cuis2{$cui} = 1;
    }

    #precalculate values for all cuis and cui pairs
    my %n11 = ();
    my %n1p = ();
    my %np1 = ();
    my $npp = 0;
    open IN, $matrix_G or die "Cannot open $matrix_G for input: $!\n";
    while (my $line = <IN>) {
	#get cuis and value from the line
	chomp $line;
	my ($cui1, $cui2, $num) = split /\t/, $line;
	
	#update n11
	if (exists $cuis1{$cui1} && exists $cuis2{$cui2}) {
	    $n11{"$cui1,$cui2"} = $num;
	}

	#udpate cui1 stats
	if (exists $cuis1{$cui1}) {
	    $n1p{$cui1} += $num;
	}

	#update cui2 stats
	if (exists $cuis2{$cui2}) {
	    $np1{$cui2} += $num;
	}

	#update npp
	$npp += $num;
    }
    close IN;

    #return counts
    my @counts = (\%n11, \%n1p, \%np1, $npp);
    return \@counts;
}

#  Gets N11 of the cui pair using a matrix
#  input : $cuis1Ref <- ref to an array of the first cuis in a set of cui pairs
#          $cuis2Ref <- ref to an array of the second cuis in a set of cui pairs
#          $n11AllRef <- ref to an array containing n11 values for all possible
#                        cui pairs of the cuis1 and cuis2, of the form
#                        n11All{"$cui1,$cui2"}=value. See _getObservedCounts_matrix
#  output: $n11      <- frequency of co-occurrences of the cuis in the cui sets
sub _getN11_matrix {
    #grab parameters
    my $self = shift;
    my $cuis1Ref = shift;
    my $cuis2Ref = shift;
    my $n11AllRef = shift;

    #error checking
    my $function = "_getN11_matrix"; 
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #calculate n11 as the sum n11s for all combinations of 
    # cuis1, cuis2 (order matters, cui1 must be first)
    my $n11 = 0;
    foreach my $cui1 (@{$cuis1Ref}) {
	foreach my $cui2 (@{$cuis2Ref}) {
	    my $num = ${$n11AllRef}{"$cui1,$cui2"};
	    if(defined $num) {
		$n11 += $num;
	    }
	}
    }

    #update values if ignoring word order
    if($noOrder_G) {
	#add all n11's, now with the order reversed
	foreach my $cui1 (@{$cuis1Ref}) {
	    foreach my $cui2 (@{$cuis2Ref}) {
		my $num = ${$n11AllRef}{"$cui2,$cui1"};
		if(defined $num) {
		    $n11 += $num;
		}
	    }
	}
    }

    return $n11;
}

#  gets N1P for a concept using a matrix
#  input : $cuis1Ref <- reference to an array containing the first cuis in a set of cui pairs
#          $countsRef <- ref to an array containing n11, n1p, np1, and npp counts
#                        for the cui combinations. See _getObservedCounts_matrix()
#          $n1pAllRef <- ref to an array containing n1p values for all cuis of cuis1 and cuis2, 
#                        of the form n1pAll{$cui} = value. See _getObservedCounts_matrix
#          $np1AllRef <- ref to an array containing n1p values for all cuis of cuis1 and cuis2, 
#                        of the form np1All{$cui} = value. See _getObservedCounts_matrix
#  output: $n1p      <- the number of times the set of concepts occurs in first position
sub _getN1p_matrix {
    #grab parameters
    my $self = shift;
    my $cuis1Ref = shift;
    my $n1pAllRef = shift;
    my $np1AllRef = shift;

    #error checking
    my $function = "_getN1P_matrix"; 
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }
    
    #calculate n1p as the sum of n1p's for all cuis1
    my $n1p = 0;
    foreach my $cui (@{$cuis1Ref}) {
	my $num = ${$n1pAllRef}{$cui};
	if(defined $num) {
	    $n1p += $num;
	}
    }

    #update values if ignoring word order
    if ($noOrder_G) {
	#add all np1's to n1p
	foreach my $cui (@{$cuis1Ref}) {
	    my $num = ${$np1AllRef}{$cui};
	    if(defined $num) {
		$n1p += $num;
	    }
	}
    }

    #set n1p to -1 if there are no values for it since this indicates
    # there is not enough information to calculate the score
    if ($n1p == 0) {
	$n1p = -1;
    }

    #return the value
    return $n1p;
}

#  gets NP1 for a concept using a matrix
#  input : $cuis2Ref <- reference to an array containing the first cuis in a set of cui pairs
#          $countsRef <- ref to an array containing n11, n1p, np1, and npp counts
#                        for the cui combinations. See _getObservedCounts_matrix()
#          $n1pAllRef <- ref to an array containing n1p values for all cuis of cuis1 and cuis2, 
#                        of the form n1pAll{$cui} = value. See _getObservedCounts_matrix
#          $np1AllRef <- ref to an array containing n1p values for all cuis of cuis1 and cuis2, 
#                        of the form np1All{$cui} = value. See _getObservedCounts_matrix
#  output: $np1      <- the number of times the set of concepts occurs in second position
sub _getNp1_matrix {
    #grab parameters
    my $self = shift;
    my $cuis2Ref = shift;
    my $n1pAllRef = shift;
    my $np1AllRef = shift;

    #calculate np1 as the sum of np1's for all cuis2
    my $np1 = 0;
    foreach my $cui (@{$cuis2Ref}) {
	my $num = ${$np1AllRef}{$cui};
	if (defined $num) {
	    $np1 += $num;
	}
    }

    #update values if ignoring word order
    if ($noOrder_G) {
	#add all n1p's to np1s
	foreach my $cui (@{$cuis2Ref}) {
	    my $num = ${$n1pAllRef}{$cui};
	    if (defined $num) {
		$np1 += $num;
	    }
	}
    }

    #set n1p to -1 if there are no values for it since this indicates
    # there is not enough information to calculate the score
    if ($np1 == 0) {
	$np1 = -1;
    }

    #return the value
    return $np1;
}

########################################################################
# functions to get statistical information about the cuis LTA
########################################################################
#  Gets contingency table values for LTA using a matrix
#  input : $cuis1Ref <- ref to an array of the first cuis in a set of cui pairs
#          $cuis2Ref <- ref to an array of the second cuis in a set of cui pairs
#  output: $\@data  <- array ref containing four values: $n11, $n1p, $np1, and 
#                      $npp for the sets of cui pairs
sub _getStats_LTA {
    #grab parameters
    my $self = shift;
    my $cuis1Ref = shift;
    my $cuis2Ref = shift;
   
    #error checking
    my $function = "_getStats_LTA"; 
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #Get co-occurrences with each set of CUIs
    my $cooccurrences1Ref;
    my $cooccurrences2Ref;
    my $npp = 0;
    if ($matrix_G) {
	#get observed counts
	my $observedRef = self->_getObserved_matrix_LTA($cuis1Ref, $cuis2Ref);

	#get co-occurrence data
	($cooccurrences1Ref, $cooccurrences2Ref) = $self
	    ->_getCUICooccurrences_matrix($cuis1Ref, $cuis2Ref, 
					  ${$observedRef}[0], ${$observedRef}[1]);

	#get npp
	$npp = ${$observedRef}[2];
    }
    else {
	#get co-occurrence data
	($cooccurrences1Ref, $cooccurrences2Ref) = $self
	    ->_getCUICooccurrences_DB($cuis1Ref, $cuis2Ref);

	#get npp
	my $db = $cuifinder_G->_getDB(); 
	$npp = shift $db->selectcol_arrayref("select count(cui_1) from N_11"); 
    }

    #calculate n1p and np1 as the number of co-occurring terms
    my $n1p = scalar keys %{$cooccurrences1Ref};
    my $np1 = scalar keys %{$cooccurrences2Ref};

    #calculate n11
    my $n11 = 0;
    #Find number of CUIs that co-occur with both CUI 1 and CUI 2
    foreach my $cui (keys %{$cooccurrences1Ref}) {
	if (exists ${$cooccurrences2Ref}{$cui}) {
	    $n11++;
	}
    }

    #return the values
    my @data = ($n11, $n1p, $np1, $npp);
    return  \@data;
}

#computes the observed co-occurrences for all combinations of the cuis passed in
#doing this in a single function makes it so all values can be computed with a 
#single pass of the input file, making execution time much faster
#  input : $cuis1Ref <- ref to an array of the first cuis in a set of cui pairs
#          $cuis2Ref <- ref to an array of the second cuis in a set of cui pairs
#  output: $\@counts  <- array ref containing three sets of values: 
#                      \%n1p, \%np1, and $npp for the cui pairs.
#                      n1p and np1 are hashes where the key is a cui, and 
#                      the value is a comma seperated list of cuis it co-occurs 
#                      with. Npp is the number of unique cuis in the vocabular
#                      which is the vocabulary size
sub _getObserved_matrix_LTA {
    #grab parameters
    my $self = shift;
    my $cuis1Ref = shift;
    my $cuis2Ref = shift;

    #convert cui arrays to hashes, makes looping thru
    # the file faster
    my %cuis1 = ();
    foreach my $cui(@{$cuis1Ref}) {
	$cuis1{$cui} = 1;
    }
    my %cuis2 = ();
    foreach my $cui(@{$cuis2Ref}) {
	$cuis2{$cui} = 1;
    }

    #get stats
    my %uniqueCuis = ();
    my %n1pAll = ();
    my %np1All = ();
    open IN, $matrix_G or die "Cannot open matrix_G for input: $matrix_G\n";
    while (my $line = <IN>) {
	#get cuis and value fro mthe line
	chomp $line;
	my ($cui1, $cui2, $num) = split /\t/, $line;

	#update n1p and np1
	if (exists $cuis1{$cui1}) {
	    $n1pAll{$cui1} .= ",$cui2";
	}
	if (exists $cuis2{$cui2}) {
	    $np1All{$cui2} .= ",$cui1";
	}
	
	#update unique cui lists to calculate npp
	$uniqueCuis{$cui1} = 1;
	$uniqueCuis{$cui2} = 1;
    }
    close IN;

    #npp is the number of unique cuis (vocab size)
    my $npp = scalar keys %uniqueCuis;
   
    #rest up and the observed
    my @observed = (\%n1pAll, \%np1All, $npp);
    return \@observed;
}

# Gets hashes of CUIs that co-occurr with the sets of cuis1 and cuis 2 using
# a matrix. This is the first step in computing linking term associations
#  input : $cuis1Ref <- ref to an array of the first cuis in a set of cui pairs
#          $cuis2Ref <- ref to an array of the second cuis in a set of cui pairs
# output: \%cooccurrences1 <- hash ref, keys are co-occurring cuis with cui 1, 
#                             values are 1
#         \%cooccurrences1 <- hash ref, keys are co-occurring cuis with cui 2, 
#                             values are 1
sub _getCUICooccurrences_matrix {
    #grab parameters
    my $self = shift;
    my $cuis1Ref = shift;
    my $cuis2Ref = shift;
    my $n1pAllRef = shift;
    my $np1AllRef = shift;
    
    #error checking
    my $function = "_getCUICooccurrences"; 
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #get lists of explicitly co-occurring CUIs for each concept
    #add trailing cui co-occurrences to cui1Data
    my %cooccurrences1;
    foreach my $cui1 (@{$cuis1Ref}){
	if (defined ${$n1pAllRef}{$cui1}) {
	    foreach my $cui (split /,/,${$n1pAllRef}{$cui1}) {
		$cooccurrences1{$cui} = 1;
	    }
	}
    }
    #add leading cui co-occurrences to cui2Data
    my %cooccurrences2;
    foreach my $cui2 (@{$cuis2Ref}) {
	if (defined ${$np1AllRef}{$cui2}) {
	    foreach my $cui (split /,/,${$np1AllRef}{$cui2}) {
		$cooccurrences2{$cui} = 1;
	    }
	}
    }

    #add more CUIs if order doesn't matter
    if ($noOrder_G) {
	#add leading co-occurring cuis to cui1Data
	foreach my $cui1 (@{$cuis1Ref}) {
	    if (defined ${$np1AllRef}{$cui1}) {
		foreach my $cui (split /,/,${$np1AllRef}{$cui1}) {
		    $cooccurrences1{$cui} = 1;
		}
	    }
	}
	#add trailling co-occurring cuis to cui2Data
	foreach my $cui2 (@{$cuis2Ref}) {
	    if(defined ${$n1pAllRef}{$cui2}) {
		foreach my $cui (split /,/,${$n1pAllRef}{$cui2}) {
		    $cooccurrences2{$cui} = 1;
		}
	    }
	}
    }

    return (\%cooccurrences1, \%cooccurrences2);
}


# Gets hashes of CUIs that co-occurr with the sets of cuis1 and cuis 2 using
# a database. This is the first step in computing linking term associations
#  input : $cuis1Ref <- ref to an array of the first cuis in a set of cui pairs
#          $cuis2Ref <- ref to an array of the second cuis in a set of cui pairs
# output: \%cooccurrences1 <- hash ref, keys are co-occurring cuis with cui 1, 
#                             values are 1
#         \%cooccurrences1 <- hash ref, keys are co-occurring cuis with cui 2, 
#                             values are 1
sub _getCUICooccurrences_DB {
    #grab parameters
    my $self = shift;
    my $cuis1Ref = shift;
    my $cuis2Ref = shift;
    
    #error checking
    my $function = "_getStats_LTA_DB"; 
    if(!defined $self || !ref $self) {
        $errorhandler->_error($pkg, $function, "", 2);
    }

    #set up database
    my $db = $cuifinder_G->_getDB(); 

    #get hashes of co-occurring CUIs
    my %cooccurrences1 = ();
    my %cooccurrences2 = ();

    #query DB to get cuis, where concept 1 is the leading cui
    my $firstCui = shift @{$cuis1Ref};
    my $query = "SELECT N_11.cui_2 FROM N_11 WHERE (N_11.cui_1 = '$firstCui' ";
    foreach my $cui (@{$cuis1Ref}) {
	$query .= "OR N_11.cui_1 = '$cui' ";
    }
    $query .= ") AND N_11.n_11 > 0;";
    my @cuis = @{$db->selectcol_arrayref($query)};
    unshift @{$cuis1Ref}, $firstCui;

    #turn CUIs into a hash of cui1's cooccurrences
    foreach my $cui (@cuis) {
	$cooccurrences1{$cui} = 1;
    }

    #query DB to get cuis, where concept 2 is the trailing cui
    $firstCui = shift @{$cuis2Ref};
    $query =  "SELECT N_11.cui_1 FROM N_11 WHERE (N_11.cui_2 = '$firstCui' ";
    foreach my $cui (@{$cuis2Ref}) {
	$query .= "OR N_11.cui_2 = '$cui' ";
    }
    $query .= ") AND N_11.n_11 > 0;";
    @cuis = @{$db->selectcol_arrayref($query)};
    unshift @{$cuis2Ref}, $firstCui;

    #turn CUIs into a hash of cui2's co-occurrences
    foreach my $cui (@cuis) {
	$cooccurrences2{$cui} = 1;
    }

    #add additional cuis if order doesn't matter
    if($noOrder_G) {
	#get cuis, where concept 1 is the trailing cui
	$firstCui = shift @{$cuis1Ref};
	my $query = "SELECT N_11.cui_1 FROM N_11 WHERE (N_11.cui_2 = '$firstCui' ";
	foreach my $cui (@{$cuis1Ref}) {
	    $query .= "OR N_11.cui_2 = '$cui' ";
	}
	$query .= ") AND N_11.n_11 > 0;";
	@cuis = @{$db->selectcol_arrayref($query)};
	unshift @{$cuis1Ref}, $firstCui;

	#add cuis to the hash of cui1's co-occurrences
	foreach my $cui (@cuis) {
	    $cooccurrences1{$cui} = 1;
	}

	#get cuis, where concept 2 is the leading cui
	$firstCui = shift @{$cuis2Ref};
	$query =  "SELECT N_11.cui_2 FROM N_11 WHERE (N_11.cui_1 = '$firstCui' ";
	foreach my $cui (@{$cuis2Ref}) {
	    $query .= "OR N_11.cui_1 = '$cui' ";
	}
	$query .= ") AND N_11.n_11 > 0;";
	@cuis = @{$db->selectcol_arrayref($query)};
	unshift @{$cuis2Ref}, $firstCui;

	#add cuis to the hash of cui2's co-occurrences
	foreach my $cui (@cuis) {
	    $cooccurrences2{$cui} = 1;
	}
    }

    #return the cui co-occurrences
    return (\%cooccurrences1, \%cooccurrences2);
}

######################################################################
#                  Utility Functions
######################################################################

#  Method to retrieve descendants of a cui
#  input : $cui <- string containing a cui 
#  output: reference to @descendants, the descendants of the given cui
sub _findDescendants {
    my $cui = shift;

    print "CUI in desc = $cui\n";
    print "umls_G = $umls_G\n";

    my $hashref = $umls_G->findDescendants($cui);
    my @descendants = (keys %{$hashref});
    return \@descendants;
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

# calculate measure assocation
my $measure = "ll"; 
my $score = $statfinder->calculateStatistic($cui1, $cui2, $measure); 

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

    <http://tech.groups.yahoo.com/group/umls-similarity/>

    =head1 AUTHOR

    Bridget T McInnes <bmcinnes@vcu.edu>
    Andriy Y. Mulyar  <andriy.mulyar@gmail.com>
    Alexander D. McQuilkin <alexmcq99@yahoo.com>
    Alex McQuilken <alexmcq99@yahoo.com>
    Sam Henry <henryst@vcu.edu>

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
