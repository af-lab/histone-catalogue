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

use 5.010;                      # Use Perl 5.10
use strict;                     # Enforce some good programming rules
use warnings;                   # Replacement for the -w flag, but lexically scoped
use File::Spec;                 # Perform operation on file names
use File::Path;                 # Create or remove directory trees

## all this does is to run the bp_genbank_ref_extractor (now part of bioperl) and
## save the sequences in the sequences directory

my $results_dir   = 'results';
my $sequences_dir = 'sequences';
my $seq_extractor = 'bp_genbank_ref_extractor';

## set path to save sequences
my @dirs           = File::Spec->splitdir( $0 );
my $sequences_path = File::Spec->catdir(@dirs[0 .. ($#dirs - 2)], $results_dir, $sequences_dir);

## remove old files to avoid problems with previous results
File::Path::remove_tree($sequences_path, {verbose => 1});

## create search string
## note that "Right side truncation with wild card does work for gene symbol" <-- from NCBI helpdesk in September 2011
my $search =  '"homo sapiens"[organism] ';
$search .= "(";
$search .= "$_*[gene name] OR " for ("H2A", "H2B", "H3", "H4", "HIST2", "HIST2", "HIST3", "HIST4", "HIST5");
$search .= "CENPA[gene name])";

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
                      '--save',         $sequences_path,
                      '--save-data',    'csv',
                      );

unshift (@extractor_args, $seq_extractor);
push    (@extractor_args, $search);
system  (@extractor_args) == 0 or die "Running @extractor_args failed: $?";
