#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use UMLS::Interface;
use UMLS::Association;
use Getopt::Long;
#############################################
#  Get Options and params
#############################################
eval(GetOptions( "version", "help", "debug", "username=s", "password=s", "hostname=s", "umlsdatabase=s", "assocdatabase=s", "socket=s", "measure=s", "conceptexpansion", "noorder", "lta", "matrix=s", "config=s","precision=s")) or die ("Please check the above mentioned option(s).\n");

#get required input
my $cuisFileName = shift;
my $outputFileName = shift;

#set default measure
my $measure = $opt_measure;
if (!$measure) {
    $measure = "tscore";
}

#############################################
#  Check help, version, minimal usage notes
#############################################
#  if help is defined, print out help
if( defined $opt_help ) {
    $opt_help = 1;
    &showHelp();
    exit;
}

#  if version is requested, show version
if( defined $opt_version ) {
    $opt_version = 1;
    &showVersion();
    exit;
}

# a single input file and output file must be passed in 
if(!(defined $cuisFileName)) {
    print STDERR "No CUI Pair Input File provided\n";
    &minimalUsageNotes();
    exit;
}
if(!(defined $outputFileName)) {
    print STDERR "No Output File provided\n";
    &minimalUsageNotes();
    exit;
}


#############################################
#  Set Up UMLS::Interface
#############################################
#  set UMLS-Interface options
my %umls_option_hash = ();

$umls_option_hash{"t"} = 1; 

if(defined $opt_debug) {
    $umls_option_hash{"debug"} = $opt_debug;
}
if(defined $opt_verbose) {
    $umls_option_hash{"verbose"} = $opt_verbose;
}
if(defined $opt_username) {
    $umls_option_hash{"username"} = $opt_username;
}
if(defined $opt_driver) {
    $umls_option_hash{"driver"}   = $opt_driver;
}
if(defined $opt_umlsdatabase) {
    $umls_option_hash{"database"} = $opt_umlsdatabase;
}
if(defined $opt_password) {
    $umls_option_hash{"password"} = $opt_password;
}
if(defined $opt_hostname) {
    $umls_option_hash{"hostname"} = $opt_hostname;
}
if(defined $opt_socket) {
    $umls_option_hash{"socket"}   = $opt_socket;
}
if(defined $opt_config){
    $umls_option_hash{"config"} = $opt_config;
}

#  instantiate instance of UMLS-Interface
my $umls = UMLS::Interface->new(\%umls_option_hash); 
die "Unable to create UMLS::Interface object.\n" if(!$umls);

#############################################
#  Set Up UMLS::Association
#############################################
#  set UMLS-Association option hash
my %assoc_option_hash = ();
$assoc_option_hash{'umls'} = $umls;

if(defined $opt_debug) {
    $assoc_option_hash{"debug"} = $opt_debug;
}
if(defined $opt_verbose) {
    $assoc_option_hash{"verbose"} = $opt_verbose;
}
if(defined $opt_username) {
    $assoc_option_hash{"username"} = $opt_username;
}
if(defined $opt_driver) {
    $assoc_option_hash{"driver"}   = $opt_driver;
}
if(defined $opt_assocdatabase) {
    $assoc_option_hash{"database"} = $opt_assocdatabase;
}
if(defined $opt_password) {
    $assoc_option_hash{"password"} = $opt_password;
}
if(defined $opt_hostname) {
    $assoc_option_hash{"hostname"} = $opt_hostname;
}
if(defined $opt_socket) {
    $assoc_option_hash{"socket"}   = $opt_socket;
}
if(defined $opt_conceptexpansion) {
    $assoc_option_hash{"conceptexpansion"}   = $opt_conceptexpansion;
}
if(defined $opt_precision){
    $assoc_option_hash{"precision"} = $opt_precision;
}
if(defined $opt_lta){
    $assoc_option_hash{"lta"} = $opt_lta;
}
if(defined $opt_noorder){
    $assoc_option_hash{"noorder"} = $opt_noorder;
}
if(defined $opt_matrix){
    $assoc_option_hash{"matrix"} = $opt_matrix;
}

#  instantiate instance of UMLS-Assocation
my $association = UMLS::Association->new(\%assoc_option_hash); 
die "Unable to create UMLS::Association object.\n" if(!$association);

#############################################
#  Calculate Association
#############################################

#read in all the first and second cuis
open IN, $cuisFileName 
    or die ("Error: unable to open cui list file: $cuisFileName");
my @cuiPairs = ();
foreach my $line (<IN>) {
    chomp $line;
    (my $cui1, my $cui2) = split('<>',$line);
    push @cuiPairs, "$cui1,$cui2";
}
close IN;

#calculate association scores for each term pair
my $scoresRef = $association->calculateAssociation_termPairList(\@cuiPairs, $params->{'measure'});

#output the results
open OUT, ">$outputFileName" 
    or die ("Error: Unable to open output file: $outputFileName");
for (my $i = 0; $i < scalar @cuiPairs; $i++) {
    (my $cui1, my $cui2) = split(',',$cuiPairs[$i]);
    print OUT "${$scoresRef}[$i]<>$cui1<>$cui2\n";
} 
close OUT;



#########################
###########################
sub minimalUsageNotes {
#TODO
}

sub showHelp {
#TODO
}

sub showVersion {
#TODO
}
