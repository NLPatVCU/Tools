#!/usr/bin/perl

=head1 NAME

ImportHadoopBigrams.pl - Reads bigrams from Hadoop output files and loads into specified MySQL database.

=head1 USAGE

ImportHadoopBigrams.pl [DATABASE OPTIONS] [OTHER OPTIONS] [FILES | DIRECTORIES]

=head1 INPUT

=head2 Required Arguments:

=head3 [FILE]

Specify the file Hadoop results are saved in 
    --files text.out_01.txt.gz 


=head2 Optional Arguments:


=head3 --database STRING        

Database to contain the CUI bigram scores. DEFAULT: CUI_Bigrams

If the database is not found in MySQL, CUICollector will create it for you. 

=head3 --username STRING

Username is required to access the CUI bigram database on MySql. You will be prompted
for it if it is not supplied as an argument. 

=head3 --password STRING

Password is required to access the CUI bigram database on MySql. You will be prompted
for it if it is not supplied as an argument. 

=head3 --hostname STRING

Hostname where mysql is located. DEFAULT: localhost

=head3 --socket STRING

Socket where the mysql.sock or mysqld.sock is located. 
DEFAULT: mysql.sock 

=head3 --port STRING

The port your mysql is using. DEFAULT: 3306

=head3 --file_step INTEGER

How many MetaMap files to read between writes to the database. 
DEFAULT: 5

MMO files can be rather large so setting a low file_step reduces the memory footprint of the script. However, setting a higher file_step reduces the number of write operations to the database.

=head3 --debug 

Sets the debug flag for testing. NOTE: extremely verbose.

=head3 --verbose 

Print the current status of the program to STDOUT. This indicates the files being processed and when the program is writing to the database. This is the default output setting.

=head3 --quiet 

Don't print anything to STDOUT.

=head3 --help

Displays the quick summary of program options.

=head1 OUTPUT

By default, CUICollector prints he current status of the program as it works through the Metamapped Medline Output files (disable with `--quiet`). It creates a database (or connects to an existing one) and adds bigram scores of the CUIs it encounters in the MMO files. 

The resulting database will have four tables:

=over

=item N_11

    cui_1   cui_2   n_11
    
This shows the count (n_11) for every time a particular CUI (cui_1) is immediately followed by another particular CUI (cui_2) in an utterance. 

=back

=head1 AUTHOR

 Amy Olex, Virginia Commonwealth University
 Keith Herbert, Virginia Commonwealth University
 Bridget McInnes, Virginia Commonwealth University

=head1 COPYRIGHT

Copyright (c) 2017
Amy Olex, Virginia Commonwealth University
alolex at vcu dot edu

Bridget McInnes, Virginia Commonwealth University
btmcinnes at vcu dot edu


This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to:

 The Free Software Foundation, Inc.,
 59 Temple Place - Suite 330,
 Boston, MA  02111-1307, USA.

=cut

###############################################################################
#                               THE CODE STARTS HERE
###############################################################################

use strict;
use warnings;

use Data::Dumper;       # Helps with debugging
use DBI;                # Database stuff
use Getopt::Long;       # Parse command line options

use feature qw(say);

$|=1;   # toggle buffering to STDOUT. Essential when CPU is fully worked.

###############################################################################
# CONSTANT STRINGS
###############################################################################

my $version = "0.20";
my $header = 
"ImportHadoopBigrams $version - (C) 2017 Amy Olex and Bridget McInnes, PhD\n"
."Released under the GNU GPL.";

my $usage = $header."\n"
."FLAGS\n"
."--debug       Print EVERYTHING to STDERR.\n"
."--verbose     Print log of files processed to STDOUT (DEFAULT)\n"
."--quiet       Print NOTHING to STDERR or STDOUT.\n"
."--help        Print this help screen.\n"
."DATABASE OPTIONS\n"
."--database    Name of MySQL database to store CUI data. (Default=CUI_Bigrams)\n"
."--hostname    Hostname for MySQL database. (Default=localhost)\n"
."--socket      Socket where the mysql.sock or mysqld.sock is located.\n"
."--username    Username for MySQL access. (Optional, but will be prompted)\n"
."--password    Password for MySQL access. (Optional, but will be prompted)\n"
."Hadoop File Options\n"
."--file       Name and path to Hadoop output file.\n"
."\nUSAGE EXAMPLES\n"
."Open the Hadoop output file and write to default database:\n"
."\tperl CUICollector.pl --file hadoop-output.txt\n"
;

###############################################################################
#                           Parse command line options 
###############################################################################
my $DEBUG = 0;      # Prints EVERYTHING. Use with small testing files.
my $VERBOSE = 1;    # Only print for reading from files, writing to database          
my $HELP = '';      # Prints usage and exits if true.


my $database = "Hadoop_CUI_Bigrams";        # Values needed to connect to MySQL dbase.
my $hostname = "localhost";
my $socket = "/var/run/mysqld/mysqld.sock";
my $port     = "3306";
my $username;
my $password;

my $file;  # Name of hadoop output file

GetOptions( 'debug'         => \$DEBUG, 
            'help'          => \$HELP,
            'verbose!'      => \$VERBOSE,
            'quiet'         => sub { $VERBOSE = 0 },
            'database=s'    => \$database,
            'hostname=s'    => \$hostname,
	    'socket=s'      => \$socket,
            'port=s'        => \$port,
            'username=s'    => \$username,
            'password=s'    => \$password,
            'file=s'   => \$file, 
);

die $usage unless $#ARGV;    
die $usage if $HELP;               
die "*** No input files ***\n$usage" unless $file;

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
 

## Test if file is readable to avoid a nasty surprise later.
die "Cannot read $file" unless -r $file;

###############################################################################
#                                   Main 
###############################################################################
{
	
say "Connecting to $database on $hostname" if $VERBOSE;
my $dbh = open_mysql_database($database, $hostname, $port, $username, $password, $socket);

update_database($dbh, $file);

# Close connection to database                       
$dbh->disconnect;

say "Finished." if $VERBOSE;
}

###############################################################################
#                           Database Subroutines
###############################################################################
sub open_mysql_database {
    my ($dbase, $host, $port, $user, $pass, $socket) = (@_);
    
       # See if database exists in the specified DBMS                        
    my @dbases = DBI->data_sources("mysql",
      {"host" => $host, "user" => $user, 
       password => $pass, "socket"=> $socket});

    my $dbase_exists = grep /DBI:mysql:$dbase/, @dbases;

    # Connect to the database if it exists. Otherwise create it from scratch.
    my $dbh;                
    if ($dbase_exists) {

        $dbh =  DBI->connect("DBI:mysql:database=$dbase;host=$host;"
			     ."mysql_socket=$socket",
			     $user, $pass,
			     {'RaiseError' => 1, 'AutoCommit' => 0});
    } 
    else {
	#connect to the DB and turn AutoCommit Off (value of 1)
        $dbh = DBI->connect("DBI:mysql:host=$host", $user, $pass,
                            {'RaiseError' => 1,'AutoCommit' => 0});
        create_database($dbh, $dbase);
    }

    return $dbh;        # Return the handler to keep connection alive.
}

###############################################################################
sub create_database {
    my $dbh = shift;
    my $dbase = shift;
    
    $dbh->do("CREATE DATABASE $dbase");
    $dbh->do("USE $dbase");
    
     $dbh->do("CREATE TABLE N_11 (   
                    cui_1   CHAR(10)    NOT NULL,
                    cui_2   CHAR(10)    NOT NULL, 
                    n_11    BIGINT      NOT NULL, 
                    PRIMARY KEY (cui_1, cui_2) )"   );
                    
}

##############################################################################
sub update_database {
    my($dbh, $file) = (@_);
    
    $dbh->do("LOAD DATA LOCAL INFILE \"$file\" INTO TABLE N_11 COLUMNS TERMINATED BY ' '");

    $dbh->commit or die $dbh->errstr;
}

