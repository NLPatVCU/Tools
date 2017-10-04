use strict;
use warnings;
use UMLS::Association::StatFinder;
use UMLS::Association::CuiFinder;
use UMLS::Association;

#user parameters
my $cuisFileName = '/home/sam/UMLS-Assoc_new/UMLS-Association/Demos/DataSets/MiniMayoSRS.snomedct.cuis';
my $n11MatrixFileName = '/home/sam/semmeddb';
my $outputFileName = 'scores';
my $noOrder =1;
my $conceptExpansion = 0;
my $measure = 'll';


#####################################
#####         BEGIN CODE        #####
#####################################

#initialize the statFinder
my %params = ();
my $cuifinder = UMLS::Association::CuiFinder->new(\%params);
my $statFinder = UMLS::Association::StatFinder->new(\%params, $cuifinder);


#read in all the first and second cuis
open IN, $cuisFileName 
    or die ("Error: unable to open cui list file: $cuisFileName");
my %cuis1 = ();
my %cuis2 = ();
my %cuiPair = ();
my %n11 = ();
foreach my $line (<IN>) {
    chomp $line;
    (my $cui1, my $cui2) = split('<>',$line);
    $cuis1{$cui1} = 0;
    $cuis2{$cui2} = 0;
    $cuiPair{"$cui1<>$cui2"} = 0;
    $n11{"$cui1<>$cui2"}=0;
}
close IN;
 
#apply concept expansion to all cuis if needed
my %originalCuiPair = %cuiPair;
if($conceptExpansion) {
    &applyConceptExpansion();
}

#read in all vaules need with a single pass of the file
open IN, $n11MatrixFileName 
    or die ("Error: unable to open input matrix file: $n11MatrixFileName\n");
my %n1p = ();
my %np1 = ();
my $npp = 0;
foreach my $line(<IN>) {
    #grab line values
    chomp $line;
    (my $cui1, my $cui2, my $val) = split(/\t/,$line);

    #update counts (n11, n1p, np1, npp)
    &updateCounts($cui1, $cui2, $val);
    if ($noOrder) {
	print "NO ORDER\n";
	&updateCounts($cui2,$cui1,$val);
    }
}
close IN;

#update original Pair counts with expanded counts
if ($conceptExpansion) {
    &updateCountsWithExpansions();
}


#calucate the association scores for all term pairs
foreach my $pairKey(keys %cuiPair) {
    (my $cui1, my $cui2) = split('<>',$pairKey);
    $cuiPair{$pairKey} = $statFinder->calculateAssociationFromValues(
	$n11{"$cui1<>$cui2"}, $n1p{$cui1}, $np1{$cui2}, $npp, $measure);
}

#output the results
open OUT, ">$outputFileName" 
    or die ("Error: Unable to open output file: $outputFileName");
foreach my $pairKey(keys %cuiPair) {
    print OUT "$cuiPair{$pairKey}<>$pairKey\n";
} 
close OUT;






#####################
######################


# updates the n11, n1p, np1, and npp counts
# Input:
# Output:
sub updateCounts {
    my $cui1 = shift;
    my $cui2 = shift;
    my $val = shift;

    #update needed values
    if (defined $cuiPair{"$cui1<>$cui2"}) {
	$n11{"$cui1<>$cui2"}+=$val;
    }
    if (defined $cuis1{$cui1}) {
	$n1p{$cui1} += $val;
    }
    if (defined $cuis2{$cui2}) {
	$np1{$cui2} += $val;
    }
    $npp += $val;
}




# updates $cuis1, $cuis2, and $cuiPairs with expanded concepts
# Input: 
# Ouptut: 
sub conceptExpansion {
   
    #expanded each pair
    foreach my $pairKey (keys %cuiPair) {
	#get cuis from each pair key
	(my $cui1, my $cui2) = split('<>',$pairKey);

	#find all descendants
	my $expandedCuis1Ref = $statFinder->_findDescendants($cui1);
	my $expandedCuis2Ref = $statFinder->_findDescendants($cui2);

	#add all cuis1 to cuis1 list
	foreach my $key (@{$expandedCuis1Ref}) {
	    $cuis1{$key} = 0;
	}

	#add all cuis2 to cuis2 list
	foreach my $key (@{$expandedCuis2Ref}) {
	    $cuis2{$key} = 0;
	}

	#add cui pairs that need to be found
	push @{$expandedCuis1Ref}, $cui1;
	push @{$expandedCuis2Ref}, $cui2;
	foreach my $key1 (@{$expandedCuis1Ref}) {
	    foreach my $key2 (@{$expandedCuis2Ref}) {
		$cuiPair{"$key1<>$key2"} = 0;
	    }
	}
    }
}


#updates the original cui pair counts with their expansions
sub updateCountsWithExpansion {
    
    #expanded each pair
    foreach my $pairKey (keys %originalCuiPair) {
	#get cuis from each pair key
	(my $cui1, my $cui2) = split('<>',$pairKey);

	#find all descendants
	my $expandedCuis1Ref = $statFinder->_findDescendants($cui1);
	my $expandedCuis2Ref = $statFinder->_findDescendants($cui2);

	#add all cuis1 to cuis1 list
	foreach my $key (@{$expandedCuis1Ref}) {
	    $n1p{$cui1} += $n1p{$key};
	}

	#add all cuis2 to cuis2 list
	foreach my $key (@{$expandedCuis2Ref}) {
	    $np1{$cui2} += $np1{$key};
	}

	#add cui pairs that need to be found
	foreach my $key1 (@{$expandedCuis1Ref}) {
	    foreach my $key2 (@{$expandedCuis2Ref}) {
		 $cuiPair{"$cui1<>$cui2"} += $n11{"$key1<>$key2"};
	    }
	}
    }
}



