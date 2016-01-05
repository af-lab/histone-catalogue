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
                  histone_type => 'H3', products => {'NM_1' => 'NP_1'});
  is ($g->uid, 42, "retrieval of UID");
  is ($g->type, 'coding', "retrieval of gene type");
  is ($g->symbol, 'CENPA', "retrieval of gene symbol");
  is ($g->histone_type, 'H3', "retrieval of histone type");
}

dies_ok {my $g = &$ctor ()}
  'dies on constructor without arguments';

dies_ok {my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'coding',
                         histone_type => 'H2AX'), products => {'NM_1' => 'NP_1'}}
  'dies when setting inavlid histone type';

{
  my $g = &$ctor (uid => 42, symbol => 'HIST1H2AA4', type => 'coding',
                  histone_type => 'H2A', products => {'NM_1' => 'NP_1'});
  ok ($g->is_core_histone (), 'H2A canonical histone is core histone');
  ok (! $g->is_linker_histone (), 'H2A canonical histone is not linker histone');
}
{
  my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'coding',
                  histone_type => 'H3', products => {'NM_1' => 'NP_1'});
  ok ($g->is_core_histone (), 'H3 histone variant is core histone');
  ok (! $g->is_linker_histone (), 'H3 histone variant is not linker histone');
}
{
  my $g = &$ctor (uid => 42, symbol => 'HIST1H1A', type => 'coding',
                  histone_type => 'H1', products => {'NM_1' => 'NP_1'});
  ok (! $g->is_core_histone (), 'H1 histone is not core histone');
  ok ($g->is_linker_histone (), 'H1 histone is linker histone');
}
{
  ## H5 is a linker histone in avian erythrocytes
  my $g = &$ctor (uid => 42, symbol => 'LOCfoo', type => 'coding',
                  histone_type => 'H5', products => {'NM_1' => 'NP_1'});
  ok (! $g->is_core_histone (), 'H5 histone is not core histone');
  ok ($g->is_linker_histone (), 'H5 histone is core histone');
}

is(HistoneGene::symbol2type("HIST1H1A"), "H1",
   "find H1 (linker) histone type from gene symbol");
is(HistoneGene::symbol2type("HIST3H2AA4"), "H2A",
   "find H2A (canonical) histone type from gene symbol from another cluster");
is(HistoneGene::symbol2type("CENPA"), "H3",
   "recognize CEPNA as H3 histone type");
ok(! defined HistoneGene::symbol2type("FOOH2AFX"),
   "return undef for non histone gene symbol");
ok(! defined HistoneGene::symbol2type("CENPAA"),
   "return undef for something similar but different from CENPA");

done_testing;
