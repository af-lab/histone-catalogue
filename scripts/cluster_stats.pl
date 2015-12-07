#!/usr/bin/perl
## Copyright (C) 2011-2015 Carnë Draug <carandraug+dev@gmail.com>
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
use List::Util;                 # Includes min and max

use HistoneCatalogue;
use MyLib;

## This script will print stats for each of the histone clusters. It will
## include:
##    (LaTeX variables with the number of histones
##    in each cluster, how many are pseudo and coding, the length of each
##    each cluster, and the location in the genome of each cluster)
##
## Usage is:
##
## cluster_stats.pl --sequences path/for/sequences

## Check input options
my %path = MyLib::parse_argv("sequences");

## Read the data and create a data structure for the analysis
my %canon;   # organized by cluster, with counts of histones and other info
my @data = MyLib::load_canonical ($path{sequences});
foreach my $gene (@data) {
  my $symbol  = $$gene{'symbol'};
  my $cluster = "HIST" . $$gene{'cluster'};
  my $histone = $$gene{'histone'};

  ## to find the start and end of each cluster, we just list all the start and
  ## end coordinates. At the end, we get the min and max of them.
  push (@{$canon{$cluster}{"coordinates"}}, $$gene{'start'}, $$gene{'end'});

  ## Get the locus.
  ## It is not possible to calculate it from the genomic coordinates (but
  ## should be possible to make bp_genbank_ref_extractor do it). However, we
  ## can find the value in the features of the transcripts files. Because of
  ## that, this only works on coding genes.
  my @transcripts = keys %{$$gene{"transcripts"}};
  if (@transcripts) {
    my  $seq      = MyLib::load_seq("transcript", $transcripts[0], $path{sequences});
    my ($feature) = $seq->get_SeqFeatures("source");
    my ($locus)   = $feature->get_tag_values("map");
    if ($locus =~ m/[\d]+[qp][\d\.]+/) {
      push (@{$canon{$cluster}{'locus'}}, $locus);
    } else {
      warn ("Could not find locus for $symbol in $transcripts[0]");
    }
  }
}

## Get the counts and stats for each of the clusters
foreach my $cluster_k (keys %canon) {
  my $cluster = $canon{$cluster_k};
  ## Calculate the length (in bp) of each cluster
  my $coord_start  = List::Util::min (@{$$cluster{'coordinates'}});
  my $coord_end    = List::Util::max (@{$$cluster{'coordinates'}});
  my $coord_length = abs ($coord_start - $coord_end);
  say HistoneCatalogue::latex_newcommand(
    $cluster_k."Span",
    $coord_length,
    "Span, in bp, of the histone cluster $cluster_k"
  );

  ## Get a nice LaTeX string showing the range of locus for each cluster,
  ## e.g, 6p21.3--6p22.2. The problem is that some locus have less precision
  ## than others. For example, 1q21 does not mean 1q21.0, it only means
  ## somewhere in 1q21. While it is tempting to just ignore it, it is not
  ## correct, we should use the one with lowest precision. Luckily, sorting
  ## seems to already do that for us. We are left with the problems of
  ## having some genes in the short and long arm of a chromosome but a cluster
  ## should not be spanning the two arms.
  my @locus = @{$$cluster{'locus'}};
  if (@locus == 0) {
    warn ("No locus information for cluster $cluster_k.");
  } else {
    my $locus_start = List::Util::minstr (@locus);
    my $locus_end   = List::Util::maxstr (@locus);
    my $locus = $locus_start eq $locus_end ?
                $locus_start : "$locus_start--$locus_end";
    say HistoneCatalogue::latex_newcommand(
      $cluster_k."Locus",
      $locus,
      "Locus of the histone cluster $cluster_k"
    );
  }
}
