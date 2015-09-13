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
use HistoneCatalogue;
use MyLib;

## This script performs a list of tests on the sequences, on what we expect
## from an histone gene, and give warnings about weird things. It will write
## about weird things it finds to:
##    * histone_insanities.tex
##
## Usage is:
##
## histone_sanity_checks.pl --sequences path/for/sequences --results path/for/log_file

my %path = MyLib::parse_argv ("sequences", "results");


my @weirds;
my @data = MyLib::load_canonical ($path{sequences});
foreach my $gene (@data) {
  my $symbol = $gene->{'symbol'};

  ## check if gene has multiple products
  my $nP = keys ($gene->{'transcripts'});
  if (! $gene->{'pseudo'} && $nP != 1) {
    push @weirds, "Gene $symbol has $nP transcripts.";
  }

  ## check if we have possibly discovered a new cluster
  if ($gene->{'cluster'} > $HistoneCatalogue::cluster_number) {
    push @weirds, "Gene $symbol belongs to unknown cluster $gene->{'cluster'}.";
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

    ## Canonical histone genes should have only 1 exon
    if ($exon_count != 1) {
      push @weirds, "Gene $symbol has $exon_count exons on transcript $acc.";
    }
    ## Canonical histone genes should not have polyA tails
    if ($polyA_tail) {
      push @weirds, "Gene $symbol has a polyA signal on transcript $acc.";
    }

    ## Canonical histone genes should have a stem loop
    if (! $stem_loop) {
      push @weirds, "Gene $symbol has no annotated stem-loop on transcript $acc.";

      ## it's not annotated, but can we find it somewhere?
      my $str = $seq->seq;
      if ($str =~ m/($HistoneCatalogue::stlp_seq)/gi) {
        my $start = pos ($str) - length ($1) +1; # start of *last* match
        push @weirds, "Gene $symbol has possible stem loop starting at position $start of $acc";
      } else {
        ## if we can't find it the stem-loop on the transcript, could it
        ## be that it's actually on the genome, but whoever made the curation
        ## thinks it's good to remove it?
        my $gseq = MyLib::load_seq ("gene", $gene->{'uid'}, $path{sequences});
        ## find the CDS for this specific gene
        foreach my $feat ($gseq->get_SeqFeatures) {
          next unless $feat->primary_tag eq "CDS";
          next unless scalar (grep {lc ($_) eq lc ($symbol)} $feat->get_tag_values("gene"));
          my $start = $feat->end;
          my $end   = $start + $HistoneCatalogue::stlp_dist + $HistoneCatalogue::stlp_length;
          if ($end > $seq->length) {
            $end = $gseq->length ;
          }
          my $subseq = $gseq->subseq($start, $end);
          if ($subseq =~ m/($HistoneCatalogue::stlp_seq)/gi) {
            my $sl = pos ($subseq) - length ($1) +1; # start of *last* match
            push @weirds, "Gene $symbol has possible stem loop in genomic starting at $sl bp from the end of its CDS";
          }
        }
      }
    ## Confirm that the annotated stem-loop is correct
    } else {
      ## the stem-loop is never too far away from the stop codon
      my $dist = $stem_loop->start - $cds->end;
      if ($dist > $HistoneCatalogue::stlp_dist) {
        push @weirds, "Gene $symbol has stem-loop $dist bp away from end of CDS on transcripts $acc.";
      }
      ## has a specific length
      if ($stem_loop->length != $HistoneCatalogue::stlp_length) {
        push @weirds, "Gene $symbol has stem-loop ".$stem_loop->length ." bp long on transcripts $acc.";
      }
      ## and a specific sequence
      if ($stem_loop->seq->seq !~ m/^$HistoneCatalogue::stlp_seq$/i) {
        push @weirds, "Gene $symbol has unmatched stem-loop sequence ".$stem_loop->seq->seq." on transcript $acc.";
      }
    }
  }
}

my $log_path = File::Spec->catdir($path{results}, "histone_insanities.tex");
open (my $log, ">", $log_path)
  or die "Could not open $log_path for writing: $!";

say {$log} "\\begin{itemize}";
say {$log} "  \\item ". HistoneCatalogue::mk_latex_string ($_) foreach (@weirds);
say {$log} "\\end{itemize}";

close ($log) or die "Couldn't close $log_path after writing: $!";
