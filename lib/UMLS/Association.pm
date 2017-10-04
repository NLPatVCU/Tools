# UMLS::Association 
#
# Perl module for scoring the semantic association of terms in the Unified
# Medical Language System (UMLS).
#
# This module borrows heavily from the UMLS::Interface package so you will 
# see similarities
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

=head1 NAME

UMLS::Association -  A suite of Perl modules that implement a number of semantic
association measures in order to calculate the semantic association between two
concepts in the UMLS. 

=head1 SYNOPSIS


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

  perl Makefile.PL PREFIX=/home/sid

It is possible to modify other parameters during installation. The
details of these can be found in the ExtUtils::MakeMaker
documentation. However, it is highly recommended not messing around
with other parameters, unless you know what you're doing.

=head1 DESCRIPTION

This package provides a Perl interface to 

=head1 DATABASE SETUP

The interface assumes that the CUI network extracted from the MetaMapped 
Medline Baseline is present in a mysql database. The name of the database 
can be passed as configuration options at initialization. However, if the 
names of the databases are not provided at initialization, then default 
value is used -- the database is called 'CUI_BIGRAMS'.

The CUI_BIGRAMS database must contain four? tables: 
	1. N11
	2. N1P
	3. NP1
	4. NPP

All other tables in the databases will be ignored, and any of these
tables missing would raise an error.

A script explaining how to create the CUI network and the mysql database 
are in the INSTALL file.

If the files that are being parsed are large, "ERROR 1206: The total number
of locks exceeds the lock table size" may occur. This can be corrected by increasing 
the lock table size of mysql. This is done by increasing the innodb_buffer_pool_size
variable in your my.cnf file. If the variable does not exist in the my.cnf file simply
add a line such as:
"innodb_buffer_pool_size=1G"
which sets the size to 1 GB. Once updated mysql must be restarted for the changes to 
take effect.

=head1 INITIALIZING THE MODULE

To create an instance of the interface object, using default values
for all configuration options:

  use UMLS::Association;
  my $associaton = UMLS::Association->new();

Database connection options can be passed through the my.cnf file. For 
example: 
           [client]
	    user            = <username>
	    password    = <password>
	    port	      = 3306
	    socket        = /tmp/mysql.sock
	    database     = mmb

Or through the by passing the connection information when first 
instantiating an instance. For example:

    $associaton = UMLS::Association->new({"driver" => "mysql", 
				  "database" => "$database", 
				  "username" => "$username",  
				  "password" => "$password", 
				  "hostname" => "$hostname", 
				  "socket"   => "$socket"}); 

  'driver'       -> Default value 'mysql'. This option specifies the Perl 
                    DBD driver that should be used to access the
                    database. This implies that the some other DBMS
                    system (such as PostgresSQL) could also be used,
                    as long as there exist Perl DBD drivers to
                    access the database.
  'database'     -> Default value 'CUI_BIGRAM'. This option specifies the name
                    of the database.
  'hostname'     -> Default value 'localhost'. The name or the IP address
                    of the machine on which the database server is
                    running.
  'socket'       -> Default value '/tmp/mysql.sock'. The socket on which 
                    the database server is using.
  'port'         -> The port number on which the database server accepts
                    connections.
  'username'     -> Username to use to connect to the database server. If
                    not provided, the module attempts to connect as an
                    anonymous user.
  'password'     -> Password for access to the database server. If not
                    provided, the module attempts to access the server
                    without a password.

More information is provided in the INSTALL file. 

=head1 PARAMETERS

You can also pass other parameters which controls the functionality 
of the Association.pm module. 

    $assoc = UMLS::Association->new({"measure"     => "lch"});

   'measure'    -> This modifies the association measure 

=head1 FUNCTION DESCRIPTIONS

=cut

package UMLS::Association;

use Fcntl;
use strict;
use warnings;
use DBI;
use bytes;

use UMLS::Association::StatFinder;
use UMLS::Association::ErrorHandler; 

my $errorhandler     = ""; 
my $statfinder = ""; 

my $pkg = "UMLS::Association";

use vars qw($VERSION);

$VERSION = '0.13';

my $debug = 0;

# UMLS-specific stuff ends ----------

# -------------------- Class methods start here --------------------

#  method to create a new UMLS::Association object
#  input : $params <- reference to hash containing the parameters 
#  output:
sub new {
    my $self      = {};
    my $className = shift;
    my $params    = shift;

    # bless the object.
    bless($self, $className);

    # initialize error handler
    $errorhandler = UMLS::Association::ErrorHandler->new();
    if(! defined $errorhandler) {
	print STDERR "The error handler did not get passed properly.\n";
	exit;
    }
    
    #  check options
    $self->_checkOptions($params);

    # Initialize the object.
    $self->_initialize($params);

    return $self;
}

#  initialize the variables and set the parameters
#  input : $params <- reference to hash containing the parameters 
#  output:
sub _initialize {
    
    my $self = shift;
    my $params = shift;

    my $function = "_initialize";

    #  check self
    if(!defined $self || !ref $self) {
    $errorhandler->_error($pkg, $function, "", 2);
    }
    
    #  set the statfinder
    $statfinder = UMLS::Association::StatFinder->new($params);
    if(! defined $statfinder) { 
	my $str = "The UMLS::Association::StatFinder object was not created.";
	$errorhandler->_error($pkg, $function, $str, 8);
    }
}

#  method checks the parameters based to the UMLS::Association package
#  input : $params <- reference to hash containing the parameters 
#  output:
sub _checkOptions {

    my $self = shift;
    my $params = shift;

    my $function = "_checkOptions";

    #  check self
    if(!defined $self || !ref $self) {
	$errorhandler->_error($pkg, $function, "", 2);
    }

    #  database options
    my $database     = $params->{'database'};
    my $hostname     = $params->{'hostname'};
    my $socket       = $params->{'socket'};
    my $port         = $params->{'port'};
    my $username     = $params->{'username'};
    my $password     = $params->{'password'};
    my $conceptExpansion = $params->{'conceptExpansion'};
   
    #  cuifinder options
    my $measure = $params->{'config'}; 
    
    #  general options
    my $debugoption  = $params->{'debug'};
    my $verbose      = $params->{'verbose'};

    if( (defined $username) && (!defined $password) ) {
	my $str = "The --password option must be defined when using --username.";
	$errorhandler->_error($pkg, $function, $str, 10);
    }

    if( (!defined $username) && (defined $password) ) {
	my $str = "The --username option must be defined when using --password.";
	$errorhandler->_error($pkg, $function, $str, 10);
    }
}

 
=head3 getFrequency

description:
 
 function returns the frequency of a given concept pair

input:   

 $concept1 <- cui
 $concept2 <- cui

output:

$frequency <- number

example:

 use UMLS::Association;
 my $associator = UMLS::Association->new(); 
 my $freq = $mmb->getFrequency($concept1, $concept2)

=cut
sub getFrequency { 
    my $self = shift;
    my $c1 = shift; 
    my $c2 = shift; 
    
    return $statfinder->_getFrequency($c1, $c2); 
}

=head3 calculateAssociation

description:
 
 function returns the given association score of a given concept pair

input:   

 $concept1 <- cui
 $concept2 <- cui 
 $measure <- statistical measure
output:

$score <- float

example:

 use UMLS::Association;
 my $associator = UMLS::Association->new(); 
 my $stat = $associator->calculateAssociation($concept1, $concept2, $measure)

=cut
sub calculateAssociation { 
    my $self = shift;
    my $c1 = shift; 
    my $c2 = shift; 
    my $meas = shift; 

    return $statfinder->calculateAssociation($c1, $c2, $meas); 
}

=head3 calculateAssociationFromValues

description:
 
 function returns the given association score of a given set of contingency 
 table values.

input:   

 $n11 <- n11 for the cui pair
 $n1p <- n1p for the cui pair
 $np1 <- np1 for the cui pair
 $npp <- npp for the cui pair
 $measure <- statistical measure
output:

$score <- float

example:

 use UMLS::Association;
 my $associator = UMLS::Association->new(); 
 my $stat = $associator->calculateStatisticFromContingencyTable(
    $n11, $n1p, $np1, $npp, $measure)

=cut
sub calculateAssociationFromValues {
    my $self = shift;
    my $n11 = shift;
    my $n1p = shift;
    my $np1 = shift;
    my $npp = shift;
    my $meas = shift;

    return $statfinder->calculateAssociationFromValues(
	$n11, $n1p, $np1, $npp, $meas); 
}

1;

__END__

=head1 REFERENCING

If you write a paper that has used UMLS-Association in some way, we'd 
certainly be grateful if you sent us a copy. Currently we have no paper
referrencing the package hopefully we will soon. 

=head1 SEE ALSO

http://search.cpan.org/dist/UMLS-Association

=head1 AUTHOR

Bridget T McInnes <btmcinnes@vcu.edu>
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

=cut
