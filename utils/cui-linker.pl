#!/usr/bin/perl

=head1 cui-linker.pl

Finds shortest path between any two Content Unique Identifiers (CUIs). 

=head1 SYNOPSIS

    $ cui-linker.pl C000001 C999999
    C000001 C012345 C654321 C999999

=head1 USAGE
    cui-linker.pl [DATABASE OPTIONS] [OTHER OPTIONS] CUI_1 CUI_2

=head1 INPUT

=cut


###############################################################################
#                               THE CODE STARTS HERE
###############################################################################

use warnings;
use strict;

use UMLS::Association;
use Getopt::Long;
use DBI;
use Data::Dumper;

use feature qw(say);

###############################################################################
#                               CONSTANT STRINGS
###############################################################################

my $version = "0.001";
my $header =
"CUI-Linker $version - (c) 2015 Keith Herbert and Bridget McInnes, PhD\n"
."Released under the GNU GPL";

my $usage = $header."\n"
."FLAGS\n"
."--debug       Print EVERYTHING to STDERR.\n"
."--verbose     Print log of files processed to STDOUT (DEFAULT)\n"
."--quiet       Print NOTHING to STDERR or STDOUT.\n"
."--help        Print this help screen.\n"
."DATABASE OPTIONS\n"
."--database    Name of MySQL database to store CUI data. (Default=CUI_DB)\n"
."--hostname    Hostname for MySQL database. (Default=localhost)\n"
."--username    Username for MySQL access. (Optional, but will be prompted)\n"
."--password    Password for MySQL access. (Optional, but will be prompted)\n"
;

###############################################################################
#                           Parse command line options 
###############################################################################

my $DEBUG = 0;      # Prints EVERYTHING. Use with small testing files.
my $VERBOSE = 1;    #           
my $HELP = '';      # Prints usage and exits if true.


my $database = "CUI_Bigrams";        # Values needed to connect to MySQL dbase.
my $hostname = "localhost";
my $port     = "3306";
my $socket   = "/home/mysql/mysql.sock";
my $username;
my $password;

my $cui_1;
my $cui_2;

GetOptions( 'debug'         => \$DEBUG, 
            'help'          => \$HELP,
            'verbose!'      => \$VERBOSE,
            'quiet'         => sub { $VERBOSE = 0 }, 
         
            'database=s'    => \$database,
            'hostname=s'    => \$hostname,
            'port=s'        => \$port,
            'username=s'    => \$username,
            'password=s'    => \$password,
            
            'cui_1=s'       => \$cui_1,
            'cui_2=s'       => \$cui_2
            );

die $usage unless $#ARGV;    
die $usage if $HELP;               

say $header if $VERBOSE;

## Prompt for username/pass if they weren't provided.
if (not defined $username){
    print "Enter username for MySQL server on $hostname: ";
    $username = <STDIN>;
    chomp $username;
}
if (not defined $password){     
    print "Enter password for $username: ";
    $password = <STDIN>;
    chomp $password;
}
  

###############################################################################
#                                   Main 
###############################################################################
{
    # Initialize our CUI Associator
    my $assoc = UMLS::Association->new({
            "database" => $database, 
            "hostname" => $hostname,
            "port"     => $port,
            "socket"   => $socket,
            "username" => $username,  
            "password" => $password
        });


    say "cui_1: $cui_1" if $DEBUG;
    say "cui_2: $cui_2" if $DEBUG;
    

    # Check if the specified CUIs exist in the database
    $assoc->exists($cui_1) or die "$cui_1 not found in $database";
    $assoc->exists($cui_2) or die "$cui_2 not found in $database";

    # Breadth-first search from CUI_1 to CUI_2
#     my $shortest_path_is_found;
#     my @nodes_visited;
#     my @shortest_path;
#     my $current = $cui_1;
#     
#     until ($shortest_path_is_found) {
#         foreach my $child_node (sort @{ $assoc->getChildren($current) }) {
#             next if $child_node ~~ @nodes_visited;
#             push @nodes_visited, $child_node;
#             
#         }
#     }

    my @cui_queue;
    my @dist_queue;
    
    my @shortest_paths;
    my $shortest_dist = 2**53;
    
    push @cui_queue, $cui_1;
    push @dist_queue, 0;
    
    while ($#cui_queue >= 0) {
        
        my $curr_path = shift @cui_queue;
        my $curr_dist = shift @dist_queue;
        
        say $curr_dist." ".$curr_path if $DEBUG;
        
        my @path = split /\s+/, $curr_path;
        my $current = $path[$#path];
        
        
        
        if($current eq $cui_2){
            push @shortest_paths, $curr_path;
            $shortest_dist = $curr_dist;
        }
        last if($curr_dist > $shortest_dist);
        
        my $dist = $curr_dist + 1;
        
        foreach my $child_cui (sort @{ $assoc->getChildren($current) }) {
            if ($curr_path =~ /$child_cui/) { next; }
            
            push @cui_queue, "$curr_path $child_cui";
            push @dist_queue, $dist;
        }
       
    }
    
    print Dumper(@shortest_paths);
    
}