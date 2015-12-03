#!/usr/bin/perl
## Copyright (C) 2014 Carnë Draug <carandraug+dev@gmail.com>
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
use List::Util;                 # Includes min and max

use MyLib;

## This script will deal with histone variants. They are not really the subject
## of this paper but a table listing all of the variant histones is useful and
## easy. It will generate the following files:
##
##    * variables-variants.tex
##
## Usage is:
##
## variants.pl --sequences path/for/sequences --results path/for/results

## Check input options
my %path = MyLib::parse_argv("sequences", "results");

my @variants = MyLib::load_variant ($path{sequences});

my $var_path = File::Spec->catdir($path{results}, "variables-variants.tex");
open (my $var_file, ">", $var_path)
  or die "Could not open $var_path for writing: $!";

say {$var_file} HistoneCatalogue::latex_newcommand (
  "TotalVariantGenes",
  scalar @variants,
  "Total number of histone variants genes"
);

close ($var_file)
  or die "Couldn't close $var_path after writing: $!";
