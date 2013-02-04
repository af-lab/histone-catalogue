#!/usr/bin/perl
## Copyright (C) 2011 Carnë Draug <carandraug+dev@gmail.com>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>.

## all this script does is run bp_genbank_ref_extractor (now part of bioperl) and
## save the sequences in the sequences directory

use 5.010;                      # Use Perl 5.10
use strict;                     # Enforce some good programming rules
use warnings;                   # Replacement for the -w flag, but lexically scoped
use File::Path;                 # Create or remove directory trees
use FindBin;                    # Locate directory of original perl script
use Getopt::Long;               # Parse program arguments
use Email::Valid;               # Validate e-mail address

use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVar;                      # Load variables

## check the e-mail is valid to avoid problems with NCBI
my $email;
GetOptions(
            'email=s'  => sub {
                                $email = Email::Valid->address($_[1])
                                  or die "Invalid e-mail adress $_[1]: $Email::Valid::Details";
                              },
          ) or die "Error processing options";
die "No email specified. Use the --email option." unless $email;

## remove old files to avoid problems with previous results
File::Path::remove_tree($MyVar::sequences_dir, {verbose => 1});

## create search string
## note that "Right side truncation with wild card does work for gene symbol" <-- from NCBI helpdesk in September 2011
my $search = '"homo sapiens"[organism] ';
$search   .= '(';
$search   .= "$_*[gene name] OR " for (@MyVar::histones);                     # get all variants
$search   .= "HIST$_*[gene name] OR " for (1 .. $MyVar::cluster_number + 1);  # all clusters and try +1 to see if there's a new one
$search   .= 'CENPA[gene name]';                                              # CENP-A name is special
$search   .= ')';

## run sequence extractor
my @extractor_args = (
                      '--assembly',     'primary assembly',
                      '--genes',        'uid',
                      '--pseudo',
                      '--non-coding',
                      '--upstream',     '0',
                      '--downstream',   '0',
                      '--transcripts',  'accession',
                      '--proteins',     'accession',
                      '--limit',        '300',
                      '--format',       'genbank',
                      '--save',         $MyVar::sequences_dir,
                      '--save-data',    'csv',
                      '--email',        $email,
                      );
unshift (@extractor_args, $MyVar::seq_extractor);
push    (@extractor_args, $search);
system  (@extractor_args) == 0 or die "Running @extractor_args failed: $?";
