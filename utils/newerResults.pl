#!usr/bin/perl
use strict;
use warnings;
use UMLS::Interface;
use UMLS::Association;
use UMLS::Association::StatFinder;
use UMLS::Association::CuiFinder;

#User configurable parameters
my @databases = qw(1975_2015_window8);
my @datasets = qw(coders physicians sim rel);
my @aMeasures = qw(ll);
my @hierarchies = qw(SNOMEDCT_US);
my @types = qw(reg conceptExpansion lta ltaWithConceptExpansion);
my @order = (0,1);

my %option_hash;
my %assoc_option_hash;
$option_hash{'t'} = 1;
$assoc_option_hash{'precision'} = 100;
$assoc_option_hash{'database'} = "1975onward";

my %n11s;
my %n1ps;
my %np1s;

if (! -d "NewerResults") {`mkdir NewerResults`;}
if (! -d "NewerResults/Data") {`mkdir NewerResults/Data`;}
if (! -d "NewerResults/Scores") {`mkdir NewerResults/Scores`;}

#iterate through different hierarchies (SNOMEDCT, MSH, UMLS_ALL)
foreach my $hierarchy (@hierarchies) {
    my $d1 = "NewerResults/Data/$hierarchy";
    my $s1 = $d1;
    $s1 =~ s/Data/Scores/;

    if (! -d $d1){`mkdir $d1`;}
    if (! -d $s1) {`mkdir $s1`;}
    

    #load the configuration and instantiate the UMLS
    open OUT, ">configuration" or die "Cannot open configuration for output: $!\n\n";
    print OUT "SAB :: include $hierarchy\nREL :: include PAR, CHD";
    close OUT;

    $option_hash{'config'} = 'configuration';
    my $umls = UMLS::Interface->new(\%option_hash);
    $assoc_option_hash{"umls"} = $umls;

    #iterate over each co-occurrence database (1809+, 1975+)
    foreach my $database (@databases) {
	my $temp = $database; 
	$temp =~ s/Matrices\///;
	my $d2 = $d1;
	$d2 .= "/$temp";
	my $s2 = $d2;
	$s2 =~ s/Data/Scores/;
	
	print STDERR "In $database\n";
	$assoc_option_hash{"matrix"} = $database;
	#iterate over the association score type (reg, conceptExpansion, lta, ltaWithConceptExpansion)
	foreach my $type (@types) {
	    print STDERR "Setting type options\n";
	    
	    my $d3 = $d2;
	    $d3 .= "/$type";
	    my $s3 = $d3;
	    $s3 =~ s/Data/Scores/;
	   
	    if (! -d $d3) {`mkdir $d3`;} 
	    if (! -d $s3) {`mkdir $s3`;}  

	    print STDERR "In $type\n";

	    #set parameters
	    if($type eq "lta" || $type eq "ltaWithConceptExpansion") {
		$assoc_option_hash{'lta'} = 1;
	    }
	    else {
		$assoc_option_hash{'lta'} = 0;
	    }

	    if($type eq "conceptExpansion" || $type eq "ltaWithConceptExpansion") {
		$assoc_option_hash{'conceptExpansion'} = 1;
	    }
	    else {
		$assoc_option_hash{'conceptExpansion'} = 0;
	    }

	    #iterate over whether order matters or not
	    foreach my $notOrdered (@order) {
		my $d4 = $d3;
		if($notOrdered) {
		    $d4 .= "/notOrdered";
		}
		else {
		    $d4 .= "/ordered";
		}

		my $s4 = $d4;
		$s4 =~ s/Data/Scores/;

		if (! -d $d4){`mkdir $d4`;}
		if (! -d $s4){`mkdir $s4`;}  

		print STDERR "      Not Ordered: $notOrdered\n";
	    
		#set parameters
		$assoc_option_hash{'noorder'} = $notOrdered;
		my $params = \%assoc_option_hash;
		
		print STDERR "Creating statfinder object\n";
		#create association objects with the now fully defined parameters
		my $cuifinder = UMLS::Association::CuiFinder->new($params);
		my $statfinder = UMLS::Association::StatFinder->new($params, $cuifinder);

		#iterate over each dataset to get all the results with this configuration
		foreach my $dataset (@datasets) {
		    my $d5 = $d4;
		    $d5 .= "/$dataset";
		    my $s5 = $d5;
		    $s5 =~ s/Data/Scores/;
		  
		    if (-d $s5) 
		    {
			`rm -rf $s5`;
		    }
		    `mkdir $s5`;

		    print STDERR "         In $dataset\n";
		    print STDERR "            Data file: $d5\n";

		    #set the CUI input and gold standard files
		    my $cuifile = "Data/";
		    my $goldfile = "Data/";
		    if($dataset eq "coders" || $dataset eq "physicians") {
			$cuifile .= "MiniMayoSRS.snomedct.cuis";
			$goldfile .= "MiniMayoSRS.snomedct." . $dataset;
		    }
		    else {
			$cuifile .= "UMNSRS_reduced_" . $dataset . ".cuis";
			$goldfile .= "UMNSRS_reduced_" . $dataset . ".gold";
		    }

		    print STDERR "Building n11 array\n";
		    my $timeBuildN11 = time;
		    my %cui1s = ();
		    my %cui2s = (); 
		    open (IN, $cuifile) or die "Can't open $cuifile for input: $!\n";
		    if($type eq "conceptExpansion" || $type eq 'ltaWithConceptExpansion') {
			while(my $line = <IN>) {
			    #format line and add to hash
			    chomp $line;
			    my ($cui1, $cui2) = split(/<>/, $line);
			    $cui1s{$cui1} = 1;
			    $cui2s{$cui2} = 1;

			    #print STDERR "$cui1\n";
			    my $hashref1 = $umls->findDescendants($cui1);
			    my @descendants1 = (sort keys %{$hashref1}); 
			    #print STDERR "$cui2\n";
			    my $hashref2 = $umls->findDescendants($cui2);
			    my @descendants2 = (sort keys %{$hashref2});
			    
			    foreach my $desc (@descendants1) {
				$cui1s{$desc} = 1;
			    }
			    foreach my $desc (@descendants2) {
				$cui2s{$desc} = 1;
			    }

			    if($notOrdered) {
				 $cui1s{$cui2} = 1;
				 $cui2s{$cui1} = 1;
				 foreach my $desc (@descendants1)
				 {
				     $cui2s{$desc} = 1;
				 }
				 foreach my $desc (@descendants2)
				 {
				     $cui1s{$desc} = 1;
				 }
			    }
			}
		    }
		    
		    #TODO this may not be the most efficient method
		    #get a list of CUIs for which N11 values matter
		    else {
			while(my $line = <IN>)
			{
			    #format line and add to hash
			    chomp $line;
			    my ($cui1, $cui2) = split(/<>/,$line);
			    $cui1s{$cui1} = 1;
			    $cui2s{$cui2} = 1;
			    
			    if($notOrdered)
			    {
				$cui1s{$cui2} = 1;
				$cui2s{$cui1} = 1;
			    }
			}
		    }
		    close IN;

		    #read in N11 for the CUI pairs that matter
		    my %cuis = ();
		    my %cui1Data = ();
		    my %cui2Data = ();
		    open IN, $database or die "Cannot open $database for input: $!\n";
		    my ($cui1,$cui2,$val);
		    while(my $line = <IN>) { 
			chomp $line;
			($cui1,$cui2,$val) = split(/\t/,$line);
			if ($type eq 'desc' || $type eq 'reg') 
			{
			    if(exists $cui1s{$cui1} && exists $cui2s{$cui2})
			    {
				$n11s{$cui1}{$cui2} = $val;
			    }
			}
			else
			{
			    if(!exists $cuis{$cui1})
			    {
				$cuis{$cui1} = 1;
			    }
			    if(!exists $cuis{$cui2})
			    {
				$cuis{$cui2} = 1;
			    }
			}
		    }  
		    close IN;
		    my $numCuis = scalar keys %cuis;
		    print STDERR "$numCuis\n";
		    

		    if($type eq 'implicit' || $type eq 'implicitDesc') {
			my $cui1File = $database . "_ImplicitCUI1Data_sorted";
			my $cui2File = $database . "_ImplicitCUI2Data_sorted";

			if (! -e $cui1File || ! -e $cui2File) {`perl getImplicitData.pl $database`;}

			#read in the cui1 file
			open IN, $cui1File or die "Cannot open $cui1File for input: $!\n\n";
			while (my $line = <IN>)
			{
			    chomp $line;
			    my ($cui1, $data) = split /\t/, $line;
			    if(exists $cui1s{$cui1})
			    {
				$cui1Data{$cui1} = $data;
			    }
			}
			close IN;

			#read in the cui2 file
			open IN, $cui2File or die "Cannot open $cui2File for input: $!\n\n";
			while (my $line = <IN>) {
			    chomp $line;
			    my ($cui2, $data) = split /\t/, $line;
			    if(exists $cui2s{$cui2})
			    {
				$cui2Data{$cui2} = $data;
			    }
			}
			close IN;
		    }

		    $timeBuildN11 = time - $timeBuildN11;
		    print STDERR "Time to build n11 array: $timeBuildN11\n";
		    
		    #get association scores for each CUI pair in the cuiFile
		    open (IN, $cuifile) or die "Can't open $cuifile for input: $!\n";
		    while(my $line = <IN>)
		    {
			my $timeOuter = time;
			#get the cui pairs from the line
			chomp $line;
			my ($concept1,$concept2) = split(/<>/,$line);
			print STDERR "On $concept1, $concept2\n";
			#get extra data (if needed)
			my $valid = -1;
			if($assoc_option_hash{'getdescendants'} && !$assoc_option_hash{'implicit'})
			{
			    print STDERR "Going into descendant data\n";
			    $valid = $statfinder->_getDescendantData($concept1, $concept2, \%n11s, \%n1ps, \%np1s);
			    print STDERR "Done with descendant data\n";
			}
			elsif($assoc_option_hash{'implicit'} && !$assoc_option_hash{'getdescendants'})
			{
			    print STDERR "Going into implicit data\n";
			    $valid = $statfinder->_getImplicitData($concept1, $concept2, $numCuis, \%cui1Data, \%cui2Data);
			}
			elsif($assoc_option_hash{'implicit'} && $assoc_option_hash{'getdescendants'})
			{
			    print STDERR "Going into implicit with descendants data\n";
			    my $time = time;
			    $valid = $statfinder->_getImplicitDataWithDescendants($concept1, $concept2, $numCuis, \%cui1Data, \%cui2Data);
			    $time = time - $time;
			    #print STDERR "Time to get values: $time\n";
			}
			else
			{
			    print STDERR "Going into get data\n";
			    my $time = time;
			    $valid = $statfinder->_getData($concept1, $concept2, \%n11s, \%n1ps, \%np1s);
			    $time = time - $time; 
			    #print STDERR "Time to get values: $time\n";
			}

			#calculate n11, n1p, np1, npp
			my $n11 = -1;
			my $n1p = -1;
			my $np1 = -1;
			my $npp = -1;
			my @data;
			if ($valid != -1)
			{
			    @data = @{$valid};
			    $n11 = $data[0];
 			    $n1p = $data[1];
			    $np1 = $data[2];
			    if($assoc_option_hash{'implicit'})
			    {
				$npp = $data[3];
			    }
			    else
			    {
				$npp = $statfinder->_getNpp();
			    }
			}

			
			#calculate the stat for each measure (ll, chi, etc...)
			foreach my $aM (@aMeasures)
			{
			    my $s6 = $d5;
			    $s6 =~ s/Data/Scores/;
			    $s6 .= "/$aM";
			    
			    my $score = -1;
			    if ($valid != -1) 
			    {
				$score = $n11;
				$score = $statfinder->_calculateStatisticFromContingencyTable($n11, $n1p, $np1, $npp, $aM);
			    }
			    
			    open SCORES, ">>$s6" or die "Error: Can't open $s6 for output: $!\n";
			    print SCORES "$score<>$line\n";
			    close SCORES;
			}
		        $timeOuter = time - $timeOuter;
			print STDERR "Time for $concept1,$concept2 : $timeOuter\n";
		    }
		    close IN;
		    
		    #generate spearmans correlations for each of the association measures
		    open (DATA, ">$d5") or die "Can't open $d5 for input: $!\n";
		    foreach my $aM (@aMeasures)
		    {
			my $s6 = $d5;
			$s6 =~ s/Data/Scores/;
			$s6 .= "/$aM";

			my @outputs = `perl spearmans.pl $goldfile $s6`;
			
			$outputs[0] =~ /(\d+)/g;
			my $n = $1;
			$outputs[1] =~ /(\d+\.\d+)/g;
			my $score = $1;			
		
			print DATA _cutString($aM) . "\t" . _cutString($score) . "\t$n\n";
		    }
		    close DATA;
		}
	    }
	}
	
    }
    
}

sub _cutString
{
    my $str = shift;
    my $max = 7; #maximum string length such that the formatting won't mess up
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
sub _genScores
{
    my ($inFile, $outFile, $database, $aM,$desc,$implicit,$cooccurrence,$notOrdered) = @_;

    open (INFILE, $inFile) or die "Can't open $inFile for input: $!\n";
    open OUTFILE, ">$outFile" or die "Error: Can't open $outFile for output: $!\n";
    while(my $line = <INFILE>)
    {
	chomp $line;
	my ($cui1,$cui2) = split(/<>/,$line);
	
	my $commStr = "umls-association.pl --assocdatabase $database --config configuration --measure $aM --precision 100 ";

	if($desc)
	{
	    $commStr .= "--getdescendants ";
	}
	if($implicit)
	{
	    $commStr .= "--implicit ";
	}
	if($cooccurrence)
	{
	    $commStr .= "--cooccurrence ";
	}
	if($notOrdered)
	{
	    $commStr .= "--noorder ";
	}
	
	$commStr .= "$cui1 $cui2";
	my $output = (`$commStr`)[0];

	my $score = (split(/<>/,$output))[0];
	print OUTFILE "$score<>$line\n";
    }
    close INFILE;
    close OUTFILE;

}	 


