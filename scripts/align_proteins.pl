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
##   align_proteins.pl path/dbfile histone_type path/output/transcript/alignment
##
## DESCRIPTION
##
## Simply calls TCoffee (actually uses the bioperl to do it) to align
## all histone sequences of a specific type.
##
## Since we are using TCoffee http://www.tcoffee.org/ we should reference
## it with:
##
## Notredame, C, Higgins, DG, Heringa, J (2000). "T-Coffee: A novel method
## for fast and accurate multiple sequence alignment" Journal of molecular
## biology, 302(1):205-218
##
## @article{tcoffee2000,
##   title={{T-Coffee}: a novel method for fast and accurate multiple sequence alignment},
##   author={Notredame, C. and Higgins, D.G. and Heringa, J. and others},
##   journal={Journal of molecular biology},
##   volume={302},
##   number={1},
##   pages={205--218},
##   year={2000},
## }

use strict;
use warnings;

use Bio::Tools::Run::Alignment::TCoffee;

use HistoneSequencesDB;

if (@ARGV != 3)
  {
    print "Usage error -- no input arguments.\n";
    print "Correct usage is:\n";
    print "\n";
    print "  \$ align_proteins.pl path/dbfile histone_type path/output/transcript/alignment\n";
    exit (1);
  }

my $db = HistoneSequencesDB::read_db($ARGV[0]);
my $histone_type = $ARGV[1];
my $out_fpath = $ARGV[2];

my @genes = grep {$_->histone_type eq $histone_type} $db->canonical_core();

my @proteins;
for my $g (@genes)
  {
    my @p = map { $db->get_protein($_) } values %{ $g->coding_products };
    push (@proteins, @p);
  }

## This works but goes against the documentation.  We can't fix the
## problem upstream because we don't know what's really wrong.  Should
## we fix the documentation or should we fix the code?  If the later
## ever happens, we will need to fix ours.
## See https://redmine.open-bio.org/issues/3406
my $tcoffee = Bio::Tools::Run::Alignment::TCoffee->new(
  'aformat' => 'fasta',
  'output'  => 'fasta', # do not be fooled by documentation
  'quiet'   => 1,       # do not be fooled by documentation
  'gapopen' => -10,     # tuned as discussed on issue #32
);

$tcoffee->outfile($out_fpath);
$tcoffee->align(\@proteins);
