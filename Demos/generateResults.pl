#function to generate results using a series of datasets 
# and parameters. This can be used to replicate or generate
# published results.
# Please hard code the parameters before the 'Begin Code'
# section.

use strict;
use warnings;

#the output file for all results
my $outputFile = 'dateRangeComparison.csv';

#a list of matrices to run
my $matrixDirectory = '/home/sam/matrices/'
my @matrices = ();
push @matrices, '1809_2015_window8';
push @matrices, '1975_2015_window8';
push @matrices, '2000_2015_window8';


#a list of measures to test
my @measures = ('ll','x2','dice','odds','leftFisher');

#A hash of parameters combinations to test
#regO, regU, regOExp, regUExp, ltaO, ltaU, ltaOExp, ltaUExp
my @paramTags = ('regO','lta','regOExp','ltaOExp');
my @paramCommands = ('', '--lta', '--conceptexpansion', '--lta --conceptexpansion');
my @orderCommands = ('','--noorder');

#dataset params
my @dataTags = ('mmCoders','mmPhys','sim','rel');
my @goldFileNames = ('DataSets/MiniMayoSRS.snomedct.coders', 'DataSets/MiniMayoSRS.snomedct.physicians', 'DataSets/UMNSRS_reduced_sim.gold','DataSets/UMNSRS_reduced_rel.gold');
my @cuisFileNames = ('DataSets/MiniMayoSRS.snomedct.cuis', 'DataSets/MiniMayoSRS.snomedct.cuis', 'DataSets/UMNSRS_reduced_sim.cuis','DataSets/UMNSRS_reduced_rel.cuis');


###################################################
#                  Begin Code
###################################################

#open the output file
open OUT, ">$outputFile" or die ("Error: unable to open outputFile: $outputFile\n");

#loop over each co-occurrence matrix
foreach my $matrix(@matrices) {
    #loop over ordered or no order
    foreach my $orderCommand(@orderCommands) {
	#set up the order string, part of the table title
	my $orderString = 'ordered';
	if ($orderCommand eq "--noorder") {
	    $orderString = "No Order";
	}

	#print titles to output file
	print OUT "$matrix $orderString\n";
	foreach my $label (@paramTags) {
	    print OUT "$label,,,,";
	}
	print OUT "\n";

	#loop over each dataset (mmCod, mmPhys, sim, reg)
	for (my $datasetNum = 0; $datasetNum < scalar @dataTags; $datasetNum++) {
	    #grab values for this dataset
	    my $dataTag = $dataTags[$datasetNum];
	    my $goldFile = $goldFileNames[$datasetNum];
	    my $cuisFileName = $cuisFileNames[$datasetNum];

	    #print titles to output file
	    #print titles to output file
	    foreach my $label (@paramTags) {
		print OUT "$dataTag,,,,";
	    }
	    print OUT "\n";\

	    #generate the line for this measure
	    foreach my $measure(@measures) {	
		#generate scores using each set of parameters
		for (my $paramNum = 0; $paramNum < scalar @paramTags; $paramNum++) {
		    #grab values for these params
		    my $paramsTag = $paramTags[$paramNum];
		    my $paramsCommand = $paramCommands[$paramNum];

		    #generate and output the correlation
		    (my $correlation, my $n) = &_getCorrelation($matrix, $matrixDirectory, $orderCommand, $dataTag, $goldFile, $cuisFileName, $measure, $paramsTag, $paramsCommand);
		    print OUT "$measure,$correlation,$n,,";
		}
		print OUT "\n";
	    }
	    print OUT "\n";
	}
	print OUT "\n";
    }
    print OUT "\n";
}


#gets spearnmans correlation using the passed in params
sub _getCorrelation {
    my $matrixName = shift;
    my $matrixPath = shift;
    my $orderCommand = shift;
    my $dataTag = shift;
    my $goldFile = shift;
    my $cuisFileName = shift;
    my $measure = shift;
    my $paramsTag = shift;
    my $paramsCommand = shift;

    #generate the results file
    my $resultsFileName = "results/$paramsTag\_$matrixName\_$dataTag\_$measure";
    my $command = "umls-association-runDataSet.pl $cuisFileName $resultsFileName --matrix $matrixPath$matrixName --measure $measure $paramsCommand $orderCommand";
    print "$command\n";
    `$command`;
    
    #find correlation with results and gold
     print "perl spearmans.pl $goldFile $resultsFileName\n";
    my @outputs = `perl spearmans.pl $goldFile $resultsFileName`;
    $outputs[0] =~ /(\d+)/g;
    my $n = $1;
    $outputs[1] =~ /(\d+\.\d+)/g;
    my $correlation = $1;
    
    print "          $correlation, $n\n";
    return ($correlation, $n);
}
