
my $matrixName = '1975_2015_window8';
my $matrixPath = "../../../$matrixName";

&runForMeasure('ll'); print "ll\n";
&runForMeasure('x2'); print "x2\n";
&runForMeasure('dice'); print "dice\n";
&runForMeasure('odds'); print "odds\n";
&runForMeasure('leftFisher'); print "leftFisher\n";

sub runForMeasure {
    my $measure = shift;

    `umls-association_runDataSet.pl DataSets/MiniMayoSRS.snomedct.cuis results/minimayo_$matrixName_reg_$measure --matrix $matrixPath --measure $measure`;

    `umls-association_runDataSet.pl DataSets/UMNSRS_reduced_rel.cuis results/rel_$matrixName_reg_$measure --matrix $matrixPath --measure $measure`;

    `umls-association_runDataSet.pl DataSets/UMNSRS_reduced_sim.cuis results/sim_1975_8_reg_ll --matrix /home/sam/1975_2015_window8 --measure $measure`;
}
