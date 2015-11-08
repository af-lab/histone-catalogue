package Gene;
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
use Moose::Util::TypeConstraints qw(enum subtype as where);

use namespace::autoclean;

## Entrez Gene types (http://www.ncbi.nlm.nih.gov/books/NBK3841/#EntrezGene.Properties)
enum 'GeneType', [qw(tRNA rRNA snRNA scRNA snoRNA miscRNA ncRNA
                     coding pseudo other unknown)];

subtype 'NonEmptyStr',
  as 'Str',
  where { length ($_) > 0 };

subtype 'PositiveInt',
  as 'Int',
  where { $_ > 0 };


## NCBI gene ID
has ['uid']
  => (is => 'ro', isa => 'PositiveInt', required => 1);

has ['symbol']
  => (is => 'ro', isa => 'NonEmptyStr', required => 1);

has ['type']
  => (is => 'ro', isa => 'GeneType', required => 1);

has ['ensembl_id', 'species']
  => (is => 'ro', isa => 'NonEmptyStr', required => 0);

has ['description']
  => (is => 'ro', isa => 'Str', required => 0);

## If we chr_start or chr_end exist, so must this (handled by BUILD)
has ['chr_acc']
  => (is => 'ro', isa => 'NonEmptyStr', required => 0);

## If one of these exist, so must the other (handled by BUILD).
has ['chr_start', 'chr_end']
  => (is => 'ro', isa => 'PositiveInt', required => 0);

## Keys and values transcription and proteins accession numbers
has ['products']
  => (is => 'ro', isa => 'HashRef', lazy => 1, default => sub { {} });


sub BUILD
{
  my $self = shift;
  if ($self->chr_start xor $self->chr_end)
    { croak "attempt to create a Gene with only chromosome start or end coordinates"; }
  elsif (($self->chr_start or $self->chr_start) and not $self->chr_acc)
    { croak "attempt to create a Gene with chromosome coordinates but no accession"; }

  if ($self->is_coding())
    {
      if (scalar (keys %{$self->products}) == 0)
        { croak "attempt to create a Gene of coding type without products"; }
    }
  elsif (scalar (keys %{$self->products}) > 0)
    { croak "attempt to create a Gene of non coding type with products"; }
}

sub is_coding
{
  my $self = shift;
  return $self->type eq 'coding';
}

sub transcripts
{
  my $self = shift;
  return keys %{$self->products()};
}

sub proteins
{
  my $self = shift;
  return values %{$self->products()};
}

__PACKAGE__->meta->make_immutable;

1;
