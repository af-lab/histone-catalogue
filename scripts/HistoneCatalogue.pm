package HistoneCatalogue;
use utf8;

## Copyright (C) 2011-2015 Carnë Draug <carandraug+dev@gmail.com>
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

use 5.010;
use strict;
use warnings;
use Carp;


=func get_sequences_date

Returns the date when sequences where retrieved.  This is parsed from
the log file of bp_genbank_ref_extractor.  Note that this date may differ
from the date when the sequences are analysed.

Params:
  fpath - path to the extractor.log by bp_genbank_ref_extractor

Returns:
  String with date.

=cut

sub get_sequences_date
{
  my $fpath = "../results/sequences/extractor.log";
  open (my $data_log, "<", $fpath)
    or croak "Could not open $fpath for reading: $!";
  my $data_header = <$data_log>; # read the first line only
  close $data_log;

  $data_header =~ m/(?<=\[)([\d\-: ]+)(?=\])/;
  return $1;
}


1;
