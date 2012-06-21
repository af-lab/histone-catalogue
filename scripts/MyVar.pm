package MyVar;
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

## This 'package' has the stuff that is common to more than one script

use 5.010;                      # Use Perl 5.10
use strict;                     # Enforce some good programming rules
use warnings;                   # Replacement for the -w flag, but lexically scoped
use File::Spec;                 # Perform operation on file names
use FindBin;                    # Locate directory of original perl script

## histones that we care about (in case one day we start caring about H1)
our @histones       = qw(H2A H2B H3 H4);
our $histone_regexp = join ("|", @histones);
## how to call bp_genbank_ref_extractor
our $seq_extractor  = 'bp_genbank_ref_extractor';
## directory where the results are saved
our $results_dir    = 'results';
## directory where bp_genbank_ref_extractor saves the sequences
our $sequences_dir  = 'sequences';
## current number of known clusters
our $cluster_number = 4;
## max number of significant figures (digits) for sizes (cluster length)
our $size_precision = 2;
## file (inside $results_dir) to store the results for the LaTeX compiler
our $results_clust  = 'variables-cluster_stats.tex';

################################################################################
## calculate complete relative paths
################################################################################

my @dirs            = File::Spec->splitdir($0);

$results_dir        = File::Spec->catdir(@dirs[0 .. ($#dirs - 2)], $results_dir);
$sequences_dir      = File::Spec->catdir($results_dir, $sequences_dir);
$results_clust      = File::Spec->catdir($results_dir, $results_clust);

## path for data saved by sequence extractor
our $data_path      = File::Spec->catdir($sequences_dir, 'data.csv');

1; # a package must return true
