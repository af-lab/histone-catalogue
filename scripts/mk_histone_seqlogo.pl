#!/usr/bin/perl
use utf8;

## Copyright (C) 2010-2015 Carnë Draug <carandraug+dev@gmail.com>
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

## SYNOPSIS
##
##   mk_histone_seqlogo.pl align_file seqlogo
##
## DESCRIPTION
##
## Takes an alignment file and calls weblogo on it.  It sets our settings
## for histone proteins and then modifies it for our prefereed visualization
## where only the changes are in black.

use strict;
use warnings;

use WebLogo;

if (@ARGV != 2)
  {
    print "Usage error -- no input arguments.\n";
    print "Correct usage is:\n";
    print "\n";
    print "  \$ mk_histone_seqlogo.pl align_file seqlogo\n";
    exit (1);
  }

my $align_path = $ARGV[0];
my $logo_path = $ARGV[1];

my $weblogo = WebLogo->new();
$weblogo->call(
  $align_path,
  $logo_path,
  {
    "--units",          "probability",
    "--show-yaxis",     "no",
    "--stacks-per-line", 50,
    "--datatype",       "fasta",
    "--errorbars",      "no",
  },
);
WebLogo::seqlogo_2_fancy_align($logo_path, $logo_path);
