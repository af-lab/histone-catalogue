package MyLib;
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

use 5.010;                                  # use Perl 5.10
use strict;                                 # enforce some good programming rules
use warnings;                               # replacement for the -w flag, but lexically scoped
use Carp;                                   # alternative warn and die for modules
use Text::CSV 1.21;                         # Comma-separated values manipulator (require 1.21 for getline_hr_all
use FindBin;                                # Locate directory of original perl script

use lib $FindBin::Bin;                      # Add script directory to @INC to find 'package'
use MyVar;                                  # Load variables

sub load_csv {
  ## To cover the widest range of parsing options, you will always want to set binary
  my $csv = Text::CSV->new ({
                              binary => 1,
                              eol    => $/,
                              }) or die "Cannot use Text::CSV: ". Text::CSV->error_diag ();
  open (my $file, "<", $MyVar::data_path) or die "Could not open $MyVar::data_path for reading: $!";

  $csv->column_names ($csv->getline ($file));   # read first line and sets it as the column name

  ## note that get_line_hr_all was only implemented on 1.21. If using 1.18, would
  ## need a while loop and use get_line_hr
  my $data_ref = $csv->getline_hr_all ($file);  # reads all lines of file into an array of hashes (returns ref to array)
  close $file;                                  # close file
  return @$data_ref;                            # dereference the array and return it
}

1; # a package must return true
