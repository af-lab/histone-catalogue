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
##   histone_sanity_checks.pl path/for/dbfile
##
## DESCRIPTION
##
## It will read the store file of an HistoneSequencesDB object, and
## perform a series of tests on the histone genes, on what we expect
## from an histone gene, and give warnings about weird things.  It
## prints to stdout, a list of LaTeX \\item for use in an itemize
## envioronment.

use 5.010;
use strict;
use warnings;

use HistoneCatalogue;
use HistoneSequencesDB;
use MyLib; # FIXME: we should stop using this.

if (@ARGV != 1)
  {
    print "Usage error -- no input arguments.\n";
    print "Correct usage is:\n";
    print "\n";
    print "  \$ histone_sanity_checks.pl path/for/dbfile\n";
    exit (1);
  }

my $db = HistoneSequencesDB::read_db ($ARGV[0]);
my @genes = HistoneSequencesDB::sort_histones($db->canonical_core());

my @weirds;
foreach my $gene (@genes) {
  my $symbol = $gene->symbol;
  my $products = $gene->products;

  ## check if gene has multiple products
  my $nP = keys %{$products};
  if ($gene->is_coding && $nP != 1)
    { push @weirds, "$symbol has $nP transcripts."; };

  ## check if we have possibly discovered a new cluster
  if ($gene->cluster > $HistoneCatalogue::cluster_number)
    { push @weirds, "$symbol belongs to unknown cluster ".$gene->cluster; }

  foreach my $acc (keys %$products)
    {
      my $seq = $db->get_transcript($acc);

      my $exon_count = 0;
      my $polyA_tail = 0; # did we found one?
      my $stem_loop  = 0;
      my $cds        = 0;

      my @feats = $seq->get_SeqFeatures;
      foreach my $feat (@feats)
        {
          my $tag = $feat->primary_tag;
          $exon_count++      if $tag eq "exon";
          $polyA_tail++      if $tag eq "polyA_signal";
          $stem_loop = $feat if $tag eq "stem_loop";
          $cds       = $feat if $tag eq "CDS";
        }

      ## Canonical histone genes should have only 1 exon
      if ($exon_count != 1)
        { push @weirds, "$symbol has $exon_count exons on transcript $acc."; }
      ## Canonical histone genes should not have polyA tails
      if ($polyA_tail)
        { push @weirds, "$symbol has a polyA signal on transcript $acc."; }

      ## Canonical histone genes should have a stem loop
      if (! $stem_loop)
        {
          push @weirds, "Gene $symbol has no annotated stem-loop on transcript $acc.";

          ## it's not annotated, but can we find it somewhere?
          my $str = $seq->seq;
          if ($str =~ m/($HistoneCatalogue::stlp_seq)/gi)
            {
              my $start = pos ($str) - length ($1) +1; # start of *last* match
              push @weirds, "$symbol has possible stem loop starting at position $start of $acc";
            }
          else
            {
              ## if we can't find it the stem-loop on the transcript, could it
              ## be that it's actually on the genome, but whoever made the curation
              ## thinks it's good to remove it?
              my $gseq = $db->get_genomic($gene->uid);
              ## find the CDS for this specific gene
              foreach my $feat ($gseq->get_SeqFeatures)
                {
                  next unless $feat->primary_tag eq "CDS";
                  next unless scalar (grep {lc ($_) eq lc ($symbol)} $feat->get_tag_values("gene"));

                  ## There should only be one (also, ideally we would refer
                  ## the transcript but we don't have that information from
                  ## the CDS feature (and we don't want to use the mRNA
                  ## feature because that includes the UTR which complicates
                  ## to find the actual CDS).
                  my $protein_acc = ($feat->get_tag_values("protein_id"))[0];

                  my $start = $feat->end;
                  my $end   = $start + $HistoneCatalogue::stlp_dist + $HistoneCatalogue::stlp_length;
                  if ($end > $seq->length)
                    { $end = $gseq->length ; }
                  my $subseq = $gseq->subseq($start, $end);
                  if ($subseq =~ m/($HistoneCatalogue::stlp_seq)/gi)
                    {
                      my $sl = pos ($subseq) - length ($1) +1; # start of *last* match
                      push @weirds, "$symbol has possible stem loop in genomic sequence starting at $sl bp from the end of protein $protein_acc CDS";
                    }
                }
            }
        }
        ## Confirm that the annotated stem-loop is correct
      else
        {
          ## the stem-loop is never too far away from the stop codon
          my $dist = $stem_loop->start - $cds->end;
          if ($dist > $HistoneCatalogue::stlp_dist)
            { push @weirds, "$symbol has stem-loop $dist bp away from end of CDS on transcripts $acc."; }
          ## has a specific length
          if ($stem_loop->length != $HistoneCatalogue::stlp_length)
            { push @weirds, "$symbol has stem-loop ".$stem_loop->length ." bp long on transcripts $acc."; }
          ## and a specific sequence
          if ($stem_loop->seq->seq !~ m/^$HistoneCatalogue::stlp_seq$/i)
            { push @weirds, "$symbol has unmatched stem-loop sequence ".$stem_loop->seq->seq." on transcript $acc."; }
        }
    }
}

say "\\item ". HistoneCatalogue::mk_latex_string ($_) foreach (@weirds);
