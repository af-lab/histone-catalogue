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
                  ensembl_id => 'foo78', products => {'NM_1' => 'NP_1'});
  is ($g->uid, 42, "retrieval of UID");
  is ($g->symbol, 'CENPA', "retrieval of gene symbol");
  is ($g->type, 'coding', "retrieval of gene type");
  is ($g->chr_acc, 'NC_007', "retrieval of chr acc");
  is ($g->chr_start, 500, "retrieval of chr start");
  is ($g->chr_end, 1000, "retrieval of chr end");
  is ($g->ensembl_id, 'foo78', "retrieval of ensembl id");
}

dies_ok {my $g = &$ctor ()}
  'constructor dies without arguments';
dies_ok {my $g = &$ctor (uid => 42, symbol => 'foo')}
  'dies without type';
dies_ok {my $g = &$ctor (uid => 42, symbol => '', type => 'coding')}
  'dies with empty symbol';
dies_ok {my $g = &$ctor (uid => -42, symbol => 'CENPA', type => 'coding')}
  'dies with non positive uid';
dies_ok {my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'garbage')}
  'dies with invalid, made up, type';
dies_ok {my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'protein coding')}
  "dies with 'protein coding' as type (very similar to coding)";

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

dies_ok {my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'coding',
                         species => '')}
  'dies with invalid species (empty string)';
dies_ok {my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'coding',
                         species => undef)}
  'dies with invalid species (undef)';

{
  my $g = &$ctor (uid => 42, symbol => 'CENPA', type => 'coding',
                  chr_acc => 'NC_007', products => {'NM_1' => 'NP_1'});
  ok (! defined ($g->chr_start),
      'chr start defaulting to undef');
}

ok (! &$ctor (uid => 5, symbol => 'A', type => 'pseudo')->is_coding (),
    'type pseudo not is_coding()');
ok (! &$ctor (uid => 5, symbol => 'A', type => 'unknown')->is_coding (),
     'type unknown not is_coding()');
ok (! &$ctor (uid => 5, symbol => 'A', type => 'ncRNA')->is_coding (),
    'type ncRNA not is_coding()');
ok (&$ctor (uid => 5, symbol => 'A', type => 'coding',
            products => {'NM_1' => 'NP_1'})->is_coding (),
    'type coding is_coding()');

{
  my $g = &$ctor (uid => 5, symbol => 'A', type => 'coding',
                  products => {'NM_1' => 'NP_1'});
  is_deeply ($g->products, {'NM_1' => 'NP_1'}, 'retrieve hash ref of products');
  is_deeply ([$g->transcripts], ['NM_1'], 'retrieve array of transcripts');
  is_deeply ([$g->proteins], ['NP_1'], 'retrieve array of proteins');
  is_deeply ($g->coding_products, {'NM_1' => 'NP_1'},
             'retrieve hash ref of coding products');
}

{
  my $g = &$ctor (uid => 5, symbol => 'A', type => 'coding',
                  products => {'NM_1' => 'NP_1', 'NM_4' => 'NP_6'});
  is_deeply ($g->products, {'NM_1' => 'NP_1', 'NM_4' => 'NP_6'},
    'retrieve hash ref of products');
  is_deeply ([sort $g->transcripts], ['NM_1', 'NM_4'],
    'retrieve array of transcripts');
  is_deeply ([sort $g->proteins], ['NP_1', 'NP_6'],
    'retrieve array of proteins');
  is_deeply ($g->coding_products, {'NM_1' => 'NP_1', 'NM_4' => 'NP_6'},
             'retrieve hash ref of coding products');
}

{
  my $g = &$ctor (uid => 5, symbol => 'A', type => 'coding',
                  products => {'NM_1' => 'NP_1', 'NM_4' => ''});
  is_deeply ($g->products, {'NM_1' => 'NP_1', 'NM_4' => ''},
    'retrieve hash ref of products with non-coding transcripts');
  is_deeply ($g->coding_products, {'NM_1' => 'NP_1'},
             'retrieve hash ref of coding products with non-coding transcripts');
}

throws_ok
  { my $g = &$ctor (uid => 5, symbol => 'A', type => 'coding',
                    products => {'' => ''})}
  qr/Validation failed for 'GeneProducts'/,
  "Fails with empty products pair";

throws_ok
  { my $g = &$ctor (uid => 5, symbol => 'A', type => 'coding',
                    products => {'NM_1' => 'NP_1', '' => ''})}
  qr/Validation failed for 'GeneProducts'/,
  "Fails with a empty products pair mixed with 'good' pairs";

throws_ok
  { my $g = &$ctor (uid => 5, symbol => 'A', type => 'coding',
                    products => {'NM_1' => '', 'NM_4' => ''})}
  qr/Gene of coding type without products/,
  "Fails with coding gene with all non-coding transcripts";

throws_ok
  { my $g = &$ctor (uid => 5, symbol => 'A', type => 'coding',
                    products => {'NM_1' => 'NP_1', '' => 'NP_4'})}
  qr/Validation failed for 'GeneProducts'/,
  "Fails with an empty transcript";

done_testing;
