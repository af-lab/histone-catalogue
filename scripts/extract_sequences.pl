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
## save the sequences in the sequences directory. Usage is:
##
## extract_sequences --email you@there.eu path_to_save_sequences

use 5.010;                      # Use Perl 5.10
use strict;                     # Enforce some good programming rules
use warnings;                   # Replacement for the -w flag, but lexically scoped
use File::Path;                 # Create or remove directory trees
use Email::Valid;               # Check validity of Internet email addresses

use FindBin;                    # Locate directory of original perl script
use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVar;                      # Load variables
use MyLib;                      # Load functions

## Make sure the email is valid. It is important that the email is correct
## since this script allows one to abuse (even if by accident), the NCBI
## servers who may block access. With an email address they will contact
## the user first.
my %opts = MyLib::parse_argv ("email");
$opts{email} = Email::Valid->address($opts{email})
  or die "Invalid e-mail adress $opts{email}: $Email::Valid::Details";

## Path to save the donwloaded sequences
my $seq_dir = $ARGV[0];

## remove old files to avoid problems with previous results
File::Path::remove_tree($seq_dir, {verbose => 1});

## create search string
## note that "Right side truncation with wild card does work for gene symbol" <-- from NCBI helpdesk in September 2011
my $search = '"homo sapiens"[organism] ';
$search   .= '(';
$search   .= "$_*[gene name] OR " for (@MyVar::histones, "H1");               # get all variants
$search   .= "HIST$_*[gene name] OR " for (1 .. $MyVar::cluster_number + 1);  # all clusters and try +1 to see if there's a new one
$search   .= 'CENPA[gene name]';                                              # CENP-A name is special
$search   .= ')';

## run sequence extractor
my @args = (
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
  '--save',         $seq_dir,
  '--save-data',    'csv',
  '--email',        $opts{'email'},
);
my @call = ($MyVar::seq_extractor, @args, $search);
system (@call) == 0 or die "Running @call failed: $?";
