#!/usr/bin/perl
use utf8;

## Copyright (C) 2011-2016 Carnë Draug <carandraug+dev@gmail.com>
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
##   describe_isoforms.pl path/for/dbfile path/for/alignment
##
## DESCRIPTION
##
## It will read the store file of an HistoneSequencesDB object, and
## print to stdout the LaTeX table describing the isoforms.

use strict;
use warnings;

use Bio::AlignIO;

use HistoneCatalogue;
use HistoneSequencesDB;

if (@ARGV != 2)
  {
    print "Usage error -- no input arguments.\n";
    print "Correct usage is:\n";
    print "\n";
    print "  \$ describe_isoforms.pl path/for/dbfile path/for/alignment\n";
    exit (1);
  }

my $db = HistoneSequencesDB::read_db($ARGV[0]);
my $align = Bio::AlignIO->new(-file => $ARGV[1])->next_aln();

my %pacc2gsym; # Protein ACCession 2 Gene SYMbol
for my $gene (@{ $db->genes })
  {
    my @proteins;
    for my $acc (values %{ $gene->coding_products })
      { push (@proteins, $acc); }
    next unless @proteins;

    my $symbol = $gene->symbol;
    if (@proteins == 1)
      { $pacc2gsym{$proteins[0]} = $symbol; }
    elsif (@proteins > 1)
      {
        ## We don't really use the gene symbol.  We need something useful
        ## for identification, so in case of genes with multiple proteins,
        ## that can't be the gene symbol.  So we append a number to it.

        ## Also, we sort the proteins because we want reproducible results.
        ## Their numbers can't change each time we rune the code.  And we
        ## use the accession number (not ascii, because they may be NM_001
        ## with NM_00004)

        my @proteins = map { $_->[0] }
                       sort { $a->[1] <=> $b->[1] }
                       map { [$_, /_(\d+)$/] } @proteins;

        for (0..$#proteins)
          { $pacc2gsym{$proteins[$_]} = "$symbol." . ($_+1); }
      }
  }
HistoneCatalogue::say_table_isoforms_description($align, %pacc2gsym);
