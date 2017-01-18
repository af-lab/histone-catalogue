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

use FindBin;
use File::Temp;
use File::Spec;
use File::Compare;
use File::Copy;

use Test::More;
use Test::Exception;

use WebLogo;


throws_ok { WebLogo->new(bin_path => '/bin') }
  qr/Attribute \(bin_path\) does not pass the type constraint/,
  'checks validity of bin_path';

## We can't really test this well, but we can check that at least
## we only get numbers and dots.
like (WebLogo->new()->version, qr/^\d[\d\.]+\z/, 'get a version number');

my $data_dir    = File::Spec->catfile($FindBin::Bin, "test-data");
my $input_test  = File::Spec->catfile($data_dir, "seqlogo-test_input.eps");
my $expected    = File::Spec->catfile($data_dir, "seqlogo-test_output.eps");

{
  my $fh = File::Temp->new();
  WebLogo::seqlogo_2_fancy_align($input_test, "$fh");
  ok(File::Compare::compare("$fh", $expected) == 0,
    "convert weblogo eps into a simpler alignment eps");
}

{
  my $fh = File::Temp->new();
  File::Copy::copy($input_test, "$fh");
  WebLogo::seqlogo_2_fancy_align("$fh", "$fh");
  ok(File::Compare::compare("$fh", $expected) == 0,
    "Test seqlogo_2_fancy_align with inplace modify");
}

done_testing;
