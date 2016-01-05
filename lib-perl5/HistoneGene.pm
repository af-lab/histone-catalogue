package HistoneGene;
use utf8;

## Copyright (C) 2015 Carnë Draug <carandraug+dev@gmail.com>
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

use strict;
use warnings;
use Carp;

use Moose;
use Moose::Util::TypeConstraints qw(enum);

use namespace::autoclean;

extends 'Gene';

enum 'HistoneType', [qw(H1 H2A H2B H3 H4 H5)];

has ['histone_type']
  => (is => 'ro', isa => 'HistoneType', required => 1);


=method is_core_histone
True if is a core histone (H2A, H2B, H3, or H4)
=cut
sub is_core_histone
{
  return not shift ()->is_linker_histone ();
}

=method is_linker_histone
True if is a linker histone (H1 or H5)
=cut
sub is_linker_histone
{
  my $type = shift()->histone_type ();
  return $type eq 'H1' || $type eq 'H5';
}


=func symbol2type
Convert gene symbol to histone type.

Args:
  symbol - string with gene symbol.

Returns:
  string or undef (if symbol does not look like histone).
=cut
sub symbol2type
{
  my $symbol = shift;

  my $type;
  if ($symbol =~ m/^HIST(\d+)(H1|H2A|H2B|H3|H4)/i) # canonical
    { $type = uc $2; }
  elsif ($symbol =~ m/^((H1|H2A|H2B|H3|H4)F|CENPA$)/i) # variant
    {
      ## $2 will be the histone if followed by F. If it's empty,
      ## then $1 will be CENPA which is a H3 variant.
      if ($2)
        { $type = uc $2; }
      elsif ($1 eq "CENPA")
        { $type = "H3"; }
    }
  return $type;
}

__PACKAGE__->meta->make_immutable;

1;
