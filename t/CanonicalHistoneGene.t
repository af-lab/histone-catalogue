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

use CanonicalHistoneGene;
my $ctor = sub { CanonicalHistoneGene->new (@_) };

{
  my $g = &$ctor (uid => 42, symbol => 'HIST2H2BA4', type => 'coding',
                  products => {'NM_1' => 'NP_1'});
  is ($g->uid, 42, "retrieval of UID");
  is ($g->type, 'coding', "retrieval of gene type");
  is ($g->symbol, 'HIST2H2BA4', "retrieval of gene symbol");
  is ($g->histone_type, 'H2B', "retrieval of histone type");
  is ($g->cluster, 2, "retrieval of histone cluster");
}

dies_ok {my $g = &$ctor ()} 'constructor dies without arguments';

throws_ok {my $g = &$ctor (uid => 42, symbol => 'HIST2H2BA4', type => 'coding',
                           cluster => 2, products => {'NM_1' => 'NP_1'})}
  qr/unknown attribute\(s\) passed to the constructor: cluster/,
  'setting cluster manually should throw an error';

throws_ok {my $g = &$ctor (uid => 42, symbol => 'HIST2H2BA4', type => 'coding',
                           histone_type => 'H2B', products => {'NM_1' => 'NP_1'})}
  qr/unknown attribute\(s\) passed to the constructor: histone_type/,
  'setting histone type manually should throw an error';

throws_ok {my $g = &$ctor (uid => 42, symbol => 'HIST5H2BA4', type => 'coding',
                           products => {'NM_1' => 'NP_1'})}
  qr/Validation failed for 'ClusterNumber' with value /,
  'attempt to use an histone with suggested cluster above 4';

throws_ok {my $g = &$ctor (uid => 42, symbol => 'HISTH2BA4', type => 'coding',
                           products => {'NM_1' => 'NP_1'})}
  qr/unable to find cluster number from symbol HISTH2BA4/,
  'input check of cluster number from symbol';

throws_ok {my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'coding',
                           products => {'NM_1' => 'NP_1'})}
  qr/unable to find cluster number from symbol CENPA/,
  'input check of cluster number from symbol';

throws_ok {my $g = &$ctor (uid => 42, symbol => 'HIST2A4', type => 'coding',
                           products => {'NM_1' => 'NP_1'})}
  qr/unable to find cluster number from symbol HIST2A4/,
  'input check of histone type from symbol';

done_testing;
