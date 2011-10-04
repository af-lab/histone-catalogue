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
use Lingua::EN::Numbers;        # Turn a number into its 'english' form
use FindBin;                    # Locate directory of original perl script

use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVars;                     # Load variables
use MyLib;                      # Load functions

## this returns an array of references to an hash whose keys are the column names
my @data = MyLib::load_csv;
my %canon;
foreach my $gene (@data) {

  my $symbol = $$gene{'gene symbol'};

  ## skip genes unless they are on cluster numbered 1-$limit
  next unless $symbol =~ m/^(HIST\d*)/;
  my $cluster = $1;

  ## warn if a gene is found who nomeclature suggests belongs to a cluster whose
  ## number is higher than the currently known
  if ( $cluster =~ m/(\d*)$/ && $1 > $MyVars::cluster_number) {
    warn ("Update/Check the code, found possible NEW histone cluster $1 with gene '$symbol'");
  }

  ## this will skip genes that do not have genomic information
  if ( !$$gene{'chromosome accession'}) {
    warn ("Gene '$symbol' has no genomic information. Skipping it!");
    next;
  }

  $canon{$cluster}{$symbol}{'start'} = $$gene{'chromosome start coordinates'};
  $canon{$cluster}{$symbol}{'end'}   = $$gene{'chromosome stop coordinates'};

}

my $number = Lingua::EN::Numbers::num2en($cluster =~ m/(\d*)$/);
$nn =~ s/^(\w)/\U$1/g;
say $nn;

