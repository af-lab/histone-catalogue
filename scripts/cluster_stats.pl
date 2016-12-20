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

## SYNOPSIS
##
##   cluster_stats.pl path/for/dbfile
##
## DESCRIPTION
##
## It will read the store file of an HistoneSequencesDB object, and
## print to stdout length of in base pairs, and the chromosome locus,
## of each cluster.

use 5.010;
use strict;
use warnings;
use List::Util; # for min and max

use HistoneSequencesDB;
use HistoneCatalogue;

=func get_locus
Return locus of a Gene object (reads it from the metadata on its transcript).

Implementation detail:

It is not possible to calculate it from the genomic coordinates (but
should be possible to make bp_genbank_ref_extractor do it). However, we
can find the value in the features of the transcripts files. Because of
that, this only works on coding genes.

Args:
  db (HistoneSequencesDB)
  gene (HistoneGene)

Returns:
  string with locus (arm letter and region and band)
=cut
sub get_locus
{
  my $db = shift;
  my $gene = shift;

  ## XXX this is a safe assumption in histone genes but may not be so for
  ##      other cases.  We can use the locus of a single transcript, and
  ##      assume it applies to the whole gene.
  for my $acc ($gene->transcripts)
    {
      my $seq = $db->get_transcript($acc);
      my ($feature) = $seq->get_SeqFeatures("source");
      my ($locus)   = $feature->get_tag_values("map");
      if ($locus)
        { return $locus; }
      ## else, try the next transcript
    }
  die ("Could not find locus for " . $gene->symbol);
}

sub main
{
  if (@_ != 1)
    {
      print "Usage error -- no input arguments.\n";
      print "Correct usage is:\n";
      print "\n";
      print "  \$ cluster_stats.pl path/for/dbfile\n";
      exit (1);
    }

  my $db = HistoneSequencesDB::read_db ($_[0]);
  my @canonical_core = $db->canonical_core;

  my %clusters = map {$_->cluster => 1} @canonical_core;
  my @clusters = sort keys %clusters;
  for my $cluster_k (@clusters)
    {
      my @this_cluster = grep {$_->cluster == $cluster_k} @canonical_core;


      ##
      ## Compute length in bp of each cluster
      ##

      ## to find the start and end of each cluster, we just list all the
      ## start and end coordinates. Then we only need the min and max.
      my @coordinates;
      push (@coordinates, $_->chr_start, $_->chr_end) for (@this_cluster);
      my $coord_start  = List::Util::min (@coordinates);
      my $coord_end    = List::Util::max (@coordinates);
      my $coord_length = abs ($coord_start - $coord_end);
      say HistoneCatalogue::latex_newcommand(
        "HIST${cluster_k}Span",
        $coord_length,
        "Span, in bp, of the histone cluster HIST$cluster_k"
      );


      ##
      ## Compute cluster genomic locus
      ##

      ## FIXME get_locus only works with coding genes...
      my @coding = grep {$_->type eq "coding"} @this_cluster;
      my @this_locus = map {get_locus($db, $_)} @coding;

      ## Get a nice LaTeX string showing the range of locus for each cluster,
      ## e.g, 6p21.3--6p22.2. The problem is that some locus have less precision
      ## than others. For example, 1q21 does not mean 1q21.0, it only means
      ## somewhere in 1q21. While it is tempting to just ignore it, it is not
      ## correct, we should use the one with lowest precision. Luckily, sorting
      ## seems to already do that for us. We are left with the problems of
      ## having some genes in the short and long arm of a chromosome but a cluster
      ## should not be spanning the two arms.
      if (@this_locus == 0) {
        warn ("No locus information for cluster HIST$cluster_k.");
      } else {
        my $locus_start = List::Util::minstr (@this_locus);
        my $locus_end   = List::Util::maxstr (@this_locus);
        my $locus = $locus_start eq $locus_end ?
                    $locus_start : "$locus_start--$locus_end";
        say HistoneCatalogue::latex_newcommand(
          "HIST${cluster_k}Locus",
          $locus,
          "Locus of the histone cluster HIST$cluster_k"
        );
      }
    }
  return 0;
}

main (@ARGV);
