#!/usr/bin/perl
## Copyright (C) 2010, 2011, 2013 Carnë Draug <carandraug+dev@gmail.com>
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

use 5.010;                      # Use Perl 5.10
use strict;                     # Enforce some good programming rules
use warnings;                   # Replacement for the -w flag, but lexically scoped

use FindBin;                    # Locate directory of original perl script
use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVar;                      # Load variables
use MyLib;                      # Load functions

## this script performs a list of tests on the sequences, on what we expect
## from an histone gene, and give warnings about weird things
my %path = MyLib::input_check ("sequences", "results");

my @data = MyLib::load_canonical ($path{sequences});

foreach my $gene (@data) {
  my $symbol = $gene->{'symbol'};

  ## check if gene has multiple products
  my $nP = keys ($gene->{'transcripts'});
  if (! $gene->{'pseudo'} && $nP != 1) {
    say "Gene $symbol has $nP transcripts.";
  }

  ## check if we have possibly discovered a new cluster
  if ($gene->{'cluster'} > $MyVar::cluster_number) {
    say "Gene $symbol belongs to unknown cluster $gene->{'cluster'}.";
  }

  foreach my $acc (keys $gene->{'transcripts'}) {
    my $seq = MyLib::load_seq ("transcripts", $acc, $path{sequences});

    my $exon_count = 0;
    my $polyA_tail = 0; # did we found one?
    my $stem_loop  = 0;
    my $cds        = 0;

    my @feats = $seq->get_SeqFeatures;
    foreach my $feat (@feats) {
      my $tag = $feat->primary_tag;
      $exon_count++      if $tag eq "exon";
      $polyA_tail++      if $tag eq "polyA_signal";
      $stem_loop = $feat if $tag eq "stem_loop";
      $cds       = $feat if $tag eq "CDS";
    }

    if ($exon_count != 1) {
      say "Gene $symbol has $exon_count exons on transcript $acc.";
    }
    if ($polyA_tail) {
      say "Gene $symbol has a polyA signal on transcript $acc.";
    }
    if (! $stem_loop) {
      say "Gene $symbol has no annotated stem-loop on transcript $acc.";
      ## it's not annotated, but can we find it somewhere?
      my $str = $seq->seq;
      if ($str =~ m/($MyVar::stlp_seq)/gi) {
        my $start = pos ($str) - length ($1) +1; # start of *last* match
        say "Gene $symbol has possible stem loop starting at position $start";
      }
    } else {
      ## the stem-loop is never too far away from the stop codon
      my $dist = $stem_loop->start - $cds->end;
      if ($dist > $MyVar::stlp_dist) {
        say "Gene $symbol has stem-loop $dist bp away from end of CDS on transcripts $acc.";
      }
      ## has a specific length
      if ($stem_loop->length != $MyVar::stlp_length) {
        say "Gene $symbol has stem-loop ".$stem_loop->length ." bp long on transcripts $acc.";
      }
      ## and a specific sequence
      if ($stem_loop->seq->seq !~ m/^$MyVar::stlp_seq$/i) {
        say "Gene $symbol has unmatched stem-loop sequence ".$stem_loop->seq->seq." on transcript $acc.";
      }
    }
  }
}
