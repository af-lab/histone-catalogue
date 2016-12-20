package CanonicalHistoneGene;
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
use Moose::Util::TypeConstraints qw(subtype as where);
use MooseX::StrictConstructor;

use namespace::autoclean;

extends 'HistoneGene';

subtype 'ClusterNumber',
  as 'Int',
  where { $_ >= 1 and $_ <= 4};

## 'init_arg => undef' will prevent from setting the value in the constructor
## (we want to be able to guess the values from the gene symbol only)

## 'is => "ro", writer => "_foo"' will create a private writer so that we
## can still make use of the isa type checks (otherwise we would have to
## set the value with 'self->{prop} = val'.

has ['cluster'] =>
(
  is => 'ro',
  isa => 'ClusterNumber',
  writer => '_cluster',
  init_arg => undef
);

has ['histone_type'] =>
(
  is => 'ro',
  isa => 'HistoneType',
  writer => '_histone_type',
  init_arg => undef
);


sub BUILD
{
  my $self = shift;

  $self->symbol() =~ m/^HIST(\d+)(H1|H2A|H2B|H3|H4|H5)/i;
  if (! $1)
    { croak 'unable to find cluster number from symbol ' . $self->symbol (); }
  if (! $2)
    { croak 'unable to find histone type from symbol ' . $self->symbol (); }

  $self->_cluster ($1);
  $self->_histone_type (uc ($2));
}


__PACKAGE__->meta->make_immutable;

1;
