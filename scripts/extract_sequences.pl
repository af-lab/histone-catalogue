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

## This script runs bp_genbank_ref_extractor (now part of bioperl), saves
## the sequences in the sequences directory, as well as CSV files for
## some sets.  It will generate the following files:
##
##    * canonical.csv
##    * variant.csv
##    * h1.csv
##
## Usage is:
##
## extract_sequences path_to_save_sequences

use 5.010;
use strict;
use warnings;
use File::Spec;
use Storable;

use MyLib; # FIXME: we should stop using this.

## Path to save the donwloaded sequences
my $seq_dir = $ARGV[0];

my %genes = MyLib::load_csv (File::Spec->catdir ($seq_dir, "data.csv"));
my @canon = MyLib::select_canonical (%genes);
my @variants = MyLib::select_variant (%genes);
my @h1 = MyLib::select_H1 (%genes);

for ((["canonical", \@canon], ["variant", \@variants], ["h1", \@h1])) {
  Storable::store ($_->[1], File::Spec->catdir ($seq_dir, $_->[0] . ".store"));
}
