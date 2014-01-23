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

use FindBin;                    # Locate directory of original perl script
use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVar;                      # Load variables
use MyLib;                      # Load functions

## This script will calculate stats for the proteins. It will
## create the following files:
##    * variables-protein_stats.tex (LaTeX variables with the arg by lys ratio)
##
## Usage is:
##
## protein_stats.pl --sequences path/for/sequences --results path/for/results

my %path = MyLib::parse_argv ("sequences", "results");
my $stats_path = File::Spec->catdir($path{results}, "variables-protein_stats.tex");
open (my $stats, ">", $stats_path) or die "Could not open $stats_path for writing: $!";

my @data = MyLib::load_canonical ($path{sequences});

## Measure the arginine by lysine ratio in both the core and linker histones
my $core_ratio   = arg_lys_ratio(@data);
my $linker_ratio = arg_lys_ratio(MyLib::load_H1 ($path{sequences}));
say {$stats} MyLib::latex_newcommand ("CoreArgLysRatio", $core_ratio);
say {$stats} MyLib::latex_newcommand ("LinkerArgLysRatio", $linker_ratio);

close($stats) or die "Couldn't close $stats_path after writing: $!";

sub arg_lys_ratio {
  my $arg = 0; # count of arginine residues
  my $lys = 0; # count of lysine residues

  foreach my $gene (@_) {
    my $access = (keys $gene->{'proteins'})[0];
    next unless $access; # skip entries with no protein acession such as pseudogenes
    my $seq = MyLib::load_seq("protein", $access, $path{sequences})->seq;
    ## we know that some genes will encode proteins with the same sequence. We
    ## are counting those again on purpose. This gives us the Arg/Lys ratio for
    ## the proteins being expressed. Of course, we do not know the expression
    ## levels of each gene so we have to assume they are equal
    $arg++ while $seq =~ m/R/ig;
    $lys++ while $seq =~ m/K/ig;
  }
  if ($arg == 0 || $lys == 0) {
    die "Found 0 arginine or lysine while calculating arg/lys ratio";
  }

  ## the actual fractions we get will be a bit unwildy (827/977) and can't be
  ## reduced or we would use Math::BigRat->new("$arg/$lys")->bnorm;
  return sprintf ("%.2f", $arg/$lys);
}
