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

use Gene;
my $ctor = sub { Gene->new (@_) };

{
  my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'coding',
                  chr_acc => 'NC_007', chr_start => 500, chr_end => 1000,
                  ensembl_id => 'foo78');
  ok ($g->uid       == 42);
  ok ($g->symbol    eq 'CENPA');
  ok ($g->type      eq 'coding');
  ok ($g->chr_acc   eq 'NC_007');
  ok ($g->chr_start == 500);
  ok ($g->chr_end   == 1000);
  ok ($g->ensembl_id eq 'foo78');
}

dies_ok {my $g = &$ctor ()};
dies_ok {my $g = &$ctor (uid => 42, symbol => 'foo')};
dies_ok {my $g = &$ctor (uid => 42, symbol => '', type => 'coding')};
dies_ok {my $g = &$ctor (uid => -42, symbol => 'CENPA', type => 'coding')};
dies_ok {my $g = &$ctor (uid => -42, symbol => 'CENPA', type => 'garbage')};
dies_ok {my $g = &$ctor (uid => -42, symbol => 'CENPA', type => 'protein coding')};

throws_ok {my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'coding',
                           chr_acc => 'NC_007', chr_start => 67)}
  qr/only chromosome start or end coordinates/,
  'chr_start and chr_end must not exist without the other';

throws_ok {my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'coding',
                           chr_acc => 'NC_007', chr_end => 67)}
  qr/only chromosome start or end coordinates/,
  'chr_start and chr_end must not exist without the other';

throws_ok {my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'coding',
                           chr_start => 100, chr_end => 67)}
  qr/chromosome coordinates but no accession/,
  'chr_start and chr_end must not exist without chr_acc';

dies_ok {my $g = &$ctor (uid => -42, symbol => 'CENPA', type => 'coding',
                         species => '')};
dies_ok {my $g = &$ctor (uid => -42, symbol => 'CENPA', type => 'coding',
                         species => undef)};

{
  my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'coding',
                  chr_acc => 'NC_007');
  ok (not defined ($g->chr_start));
}

ok (not &$ctor (uid => 5, symbol => 'A', type => 'pseudo')->is_coding ());
ok (not &$ctor (uid => 5, symbol => 'A', type => 'unknown')->is_coding ());
ok (not &$ctor (uid => 5, symbol => 'A', type => 'ncRNA')->is_coding ());
ok (&$ctor (uid => 5, symbol => 'A', type => 'coding')->is_coding ());

done_testing;
