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
use File::Path;                 # Create or remove directory trees
use Text::CSV;                  # Comma-separated values manipulator
use FindBin;                    # Locate directory of original perl script

use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVars;                     # Load variables

## To cover the widest range of parsing options, you will always want to set binary.
my $csv = Text::CSV->new ({ binary => 1, eol => $/ });
open (my $data, "<", $MyVars::data_path) or die "Could not open $MyVars::data_path for reading: $!";
while (my $row = $csv->getline ($data)) {
   my @fields = @$row;
  say @$row;
}
#my @rows;
#my $csv = Text::CSV->new ( { binary => 1 } )  # should set binary attribute.
#               or die "Cannot use CSV: ".Text::CSV->error_diag ();

#open my $fh, "<:encoding(utf8)", "test.csv" or die "test.csv: $!";
#while ( my $row = $csv->getline( $fh ) ) {
#   $row->[2] =~ m/pattern/ or next; # 3rd field should match
#   push @rows, $row;
#}
#$csv->eof or $csv->error_diag();
#close $fh;

#$csv->eol ("\r\n");

#open $fh, ">:encoding(utf8)", "new.csv" or die "new.csv: $!";
#$csv->print ($fh, $_) for @rows;
#close $fh or die "new.csv: $!";

