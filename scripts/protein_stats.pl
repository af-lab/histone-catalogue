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
use Getopt::Long;               # Parse program arguments

use FindBin;                    # Locate directory of original perl script
use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVar;                      # Load variables
use MyLib;                      # Load functions

## Check input options
my %path = ("sequences" => "",
            "results"   => "");
GetOptions(
            "sequences=s" => \$path{sequences},
            "results=s"   => \$path{results},
          ) or die "Error processing options. Paths must be strings";
for (keys %path) {
  die "No path for $_ specified. Use the --$_ option." unless $path{$_};
}


my $core_ratio   = arg_lys_ratio(MyLib::load_canonical ($path{sequences}));
my $linker_ratio = arg_lys_ratio(MyLib::load_H1 ($path{sequences}));

my $filepath = File::Spec->catdir($path{results}, "variables-protein_stats.tex");
open (my $fh, ">", $filepath) or die "Could not open $filepath for writing: $!";
say {$fh} MyLib::latex_newcommand ("CoreArgLysRatio", $core_ratio);
say {$fh} MyLib::latex_newcommand ("LinkerArgLysRatio", $linker_ratio);
close($fh) or die "Couldn't close $filepath after writing: $!";


sub arg_lys_ratio {
  my $arg = 0; # count of arginine residues
  my $lys = 0; # count of lysine residues
  foreach my $gene (@_) {
    my $access = $$gene{'protein accession'};
    next unless $access; # skip entries with no protein acession such as pseudogenes
    my $seq = MyLib::load_protein($path{sequences}, $access)->seq;
    ## we know that some genes will encode proteins with the same sequence. We
    ## are counting those again on purpose. This gives us the Arg/Lys ratio for
    ## the proteins being expressed. Of course, we do not know the expression
    ## levels of each gene so we have to assume they are equal
    $arg++ while $seq =~ m/R/ig;
    $lys++ while $seq =~ m/K/ig;
  }
  ## the actual fractions we get will be a bit unwildy (827/977) and can can't be
  ## reduced or we would use Math::BigRat->new("$arg/$lys")->bnorm;
  return sprintf ("%.2f", $arg/$lys);
}
