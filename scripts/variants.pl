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
use Number::Format;             # pretty format of long numbers with SI prefix and precision

use FindBin;                    # Locate directory of original perl script
use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVar;                      # Load variables
use MyLib;                      # Load functions

## This script will deal with histone variants. They are not really the subject
## of this paper but a table listing all of the variant histones is useful and
## easy. It will generate the following files:
##
##    * table-variant_catalogue.tex (a very long LaTeX table with all variant
##      histones, their UIDs, and transcript and protein accession numbers)
##
## Usage is:
##
## variants.pl --sequences path/for/sequences --results path/for/results

## Check input options
my %path = MyLib::parse_argv("sequences", "results");

## Get the variant genes and order them by gene symbol
my @variants = sort {$$a{'symbol'} cmp $$b{'symbol'}}
  MyLib::load_variants ($path{sequences});

my $table_path = File::Spec->catdir($path{results}, "table-variant_catalogue.tex");
MyLib::make_catalogue ($table_path, @variants);
