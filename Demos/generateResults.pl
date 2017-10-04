#sample program to generate results using the sample matrix
use strict;
use warnings;
use lib '../lib/';
use UMLS::Association;

#things to loop over
my @cuiFiles = qw(DataSets/UMNSRS_reduced_sim.cuis);
my @goldFiles =  qw(DataSets/UMNSRS_reduced_sim.gold);

my @assocMeasures = qw(ll);
my @assocTypes = qw(reg conceptexpansion);
my @orderOptions = (0,1);
#my $dataMatrix = 'sampleMatrix';
my $dataMatrix = '/home/sam/semmeddb';
my $assocDB = 'CUI_Bigram';

#output parameters
my $outputFile = 'results_matrixSemMedDB_simReg.txt';
my $tempResultsOutFile = 'tempResultsOut_matrixSemMedDB_simReg.txt';

######################################################################
#  Begin Code to loop over files and parameters and generate scores
######################################################################

#open the output file
open RESULTS_OUT, ">$outputFile" or die("Error: cannot open outputFile = $outputFile\n");

#create options hash
my %assocOptions;
#$option_hash{'t'} = 1;
$assocOptions{'precision'} = 100;
$assocOptions{'database'} = $assocDB;
if ($dataMatrix ne '') {
    $assocOptions{'matrix'} = $dataMatrix;
}
$assocOptions{'hostname'} = '192.168.24.89';
$assocOptions{'socket'} = '/var/run/mysqld.sock';
$assocOptions{'username'} = 'henryst';
$assocOptions{'password'} = 'OhFaht3eique';

#generate scores over association types
foreach my $assocType (@assocTypes) {
    print "assocType = $assocType\n";

    #generate scores over whether or not order matters
    foreach my $orderMatter (@orderOptions) {
	print "      orderMatter = $orderMatter\n";
	
	#set parameters for this run
	delete $assocOptions{'lta'};
	if ($assocType =~ /lta/) {
	    $assocOptions{'lta'} = 1;
	}
	delete $assocOptions{'conceptexpansion'};
	if ($assocType =~ /conceptexpansion/) {
	    $assocOptions{'conceptexpansion'} = 1;
	}
	delete $assocOptions{'noorder'};
	if ($orderMatter) {
	    $assocOptions{'noorder'} = 1;
	}
	
        #create UMLS::Association with parameters for this run
	print "         initalizing UMLS::Association\n";
	my $assoc = UMLS::Association->new(\%assocOptions);

	#loop over each input file (dataset)
	for (my $fileIndex = 0; $fileIndex < scalar @cuiFiles; $fileIndex++) {
	    print "         cuiFile = $cuiFiles[$fileIndex]\n";

	    #loop over each measure
	    foreach my $measure (@assocMeasures) {
		print "            measure = $measure\n";

		#calculate stats for each cui pair, and write to file
		open IN, $cuiFiles[$fileIndex] or die ("Error: cannot open cuiFiles[$fileIndex] = $cuiFiles[$fileIndex]\n");
		open SCORES_OUT, ">$tempResultsOutFile" or die ("Error: cannot open tempResultsOutFile = $tempResultsOutFile\n");
		while (my $line = <IN>) {
		    chomp $line;
		    my ($cui1, $cui2) = split(/<>/, $line);
		    my $score = $assoc->calculateAssociation(
			$cui1, $cui2, $measure);
		    print SCORES_OUT "$score<>$line\n";
		}
		close IN;
		close SCORES_OUT;

		#ensure gold can be opened
		open GOLD, $goldFiles[$fileIndex] or die ("Error: cannot open goldFiles[$fileIndex] = $goldFiles[$fileIndex]\n");
		close GOLD;

		#get the spearmans correlation
		my @outputs = `perl spearmans.pl $goldFiles[$fileIndex] $tempResultsOutFile`;
	        $outputs[0] =~ /(\d+)/g;
		my $n = $1;
		$outputs[1] =~ /(\d+\.\d+)/g;
		my $score = $1;		
		
		#output the results
		print RESULTS_OUT "$goldFiles[$fileIndex]\t$measure\t$score\t$n\n";
	    }
	}
    }
}
close RESULTS_OUT;

print "DONE!\n";
