#!/usr/bin/perl
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

use Test::More;
use Test::Exception;

use HistoneGene;
my $ctor = sub { HistoneGene->new (@_) };

{
  my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'coding',
                  histone_type => 'H3');
  ok ($g->symbol eq 'CENPA');
  ok ($g->histone_type eq 'H3');
}

dies_ok {my $g = &$ctor ()}
  'dies on constructor without arguments';

dies_ok {my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'coding',
                         histone_type => 'H2AX')}
  'dies when setting inavlid histone type';

{
  my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'coding',
                  histone_type => 'H3');
  ok ($g->is_core_histone ());
  ok (not $g->is_linker_histone ());
}
{
  my $g = &$ctor (uid => 42, symbol => 'HIST1H1A', type => 'coding',
                  histone_type => 'H1');
  ok (not $g->is_core_histone ());
  ok ($g->is_linker_histone ());
}
{
  ## H5 is a linker histone in avian erythrocytes
  my $g = &$ctor (uid => 42, symbol => 'LOCfoo', type => 'coding',
                  histone_type => 'H5');
  ok (not $g->is_core_histone ());
  ok ($g->is_linker_histone ());
}

done_testing;
