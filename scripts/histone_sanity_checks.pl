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
use POSIX ();                   # We want to use strftime

use FindBin;                    # Locate directory of original perl script
use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVar;                      # Load variables
use MyLib;                      # Load functions

## This script performs a list of tests on the sequences, on what we expect
## from an histone gene, and give warnings about weird things. It will write
## about weird things it finds to:
##    * sanity_checks.log
##
## Usage is:
##
## histone_sanity_checks.pl --sequences path/for/sequences --results path/for/log_file

my %path = MyLib::parse_argv ("sequences", "results");

my $log_path = File::Spec->catdir($path{results}, "sanity_checks.log");
open (my $log, ">", $log_path) or die "Could not open $log_path for writing: $!";

## Get the first line of the extractor log which is the date when the
## sequences were actually retrieved
my $data_log_path = File::Spec->catdir($path{sequences}, "extractor.log");
open (my $data_log, "<", $data_log_path) or die "Could not open $data_log_path for reading: $!";
chomp (my $data_header = <$data_log>); # read the first line only
close $data_log;

my $time = POSIX::strftime ("%F %T", localtime $^T);
say {$log} <<"END_HEADER";
Running histone_sanity_checks on $time using data from:
$data_header
--------------------------------------------------------------------------------
END_HEADER

my @data = MyLib::load_canonical ($path{sequences});
foreach my $gene (@data) {
  my $symbol = $gene->{'symbol'};

  ## check if gene has multiple products
  my $nP = keys ($gene->{'transcripts'});
  if (! $gene->{'pseudo'} && $nP != 1) {
    say {$log} "Gene $symbol has $nP transcripts.";
  }

  ## check if we have possibly discovered a new cluster
  if ($gene->{'cluster'} > $MyVar::cluster_number) {
    say {$log} "Gene $symbol belongs to unknown cluster $gene->{'cluster'}.";
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
    say {$log} "Gene $symbol has $exon_count exons on transcript $acc."
      if ($exon_count != 1);
    ## Canonical histone genes should not have polyA tails
    say {$log} "Gene $symbol has a polyA signal on transcript $acc."
      if ($polyA_tail);

    ## Canonical histone genes should have a stem loop
    if (! $stem_loop) {
      say {$log} "Gene $symbol has no annotated stem-loop on transcript $acc.";
      ## it's not annotated, but can we find it somewhere?
      my $str = $seq->seq;
      if ($str =~ m/($MyVar::stlp_seq)/gi) {
        my $start = pos ($str) - length ($1) +1; # start of *last* match
        say {$log} "Gene $symbol has possible stem loop starting at position $start";
      }
    ## Confirm that the annotated stem-loop is correct
    } else {
      ## the stem-loop is never too far away from the stop codon
      my $dist = $stem_loop->start - $cds->end;
      if ($dist > $MyVar::stlp_dist) {
        say {$log} "Gene $symbol has stem-loop $dist bp away from end of CDS on transcripts $acc.";
      }
      ## has a specific length
      if ($stem_loop->length != $MyVar::stlp_length) {
        say {$log} "Gene $symbol has stem-loop ".$stem_loop->length ." bp long on transcripts $acc.";
      }
      ## and a specific sequence
      if ($stem_loop->seq->seq !~ m/^$MyVar::stlp_seq$/i) {
        say {$log} "Gene $symbol has unmatched stem-loop sequence ".$stem_loop->seq->seq." on transcript $acc.";
      }
    }
  }
}
close ($log) or die "Couldn't close $log_path after writing: $!";
