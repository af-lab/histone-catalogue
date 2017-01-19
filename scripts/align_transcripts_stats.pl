#!/usr/bin/perl
use utf8;

## This script was highly based on bp_pairwise_kaks, by Jason Stajich,
## which is under the Perl 5 license (Artistic License 1.0 or GNU GPL)
##
## Copyright (C) 2003-2015 Jason Stajich <jason@bioperl.org>
## Copyright (C) 2015-2016 Carnë Draug <carandraug+dev@gmail.com>
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
##   align_transcripts_stats.pl path/for/transcript/align_1 path/for/transcript/align_2 path/for/transcript/align_X
##
## DESCRIPTION
##
## It read several alignment of transcripts and prints to stdout stats about
## them.  At the moment, that's only the mean dN, dS, and dN/dS.

use 5.010;
use strict;
use warnings;

use File::Spec;
use List::Util;

## Codeml implements the method of Goldman and Yang (1994)
## Yn00 implements the method of Yang and Nielsen (2000)
use Bio::Tools::Run::Phylo::PAML::Codeml;
use Bio::AlignIO;

use Statistics::Basic;

use HistoneCatalogue;

=func get_dNdS

Computes the pairwise dN, dS, and omega (dN/dS) values from an alignment.

Since the purpose is to estimate synonymous and nonsynonymous substitutions,
the alignment should have been performed in multiples of 3.  This can be
done by aligning the protein sequence and the convert the aligned protein
sequences back to their CDS sequences.

Args:
  align (Bio::Align::AlignI)

Returns:
  {'dN' => float, 'dS' => float, 'omega' => float}
=cut
sub get_dNdS_stats
{
  my $align = shift;

  my $dsdn_factory = Bio::Tools::Run::Phylo::PAML::Codeml->new(
    -params => {'runmode' => -2, 'seqtype' => 1},
  );
  $dsdn_factory->alignment($align);
  my ($rc, $parser) = $dsdn_factory->run();
  if ($rc <= 0)
    { die "Failed on phylogenetic analysis: ${dsdn_factory->error_string}"; }

  my $result = $parser->next_result();
  my $MLmatrix = $result->get_MLmatrix();

  ## For inspection of individual values, we can get the sequences
  ## names by using $otus[$i]->display_id
#  my @otus = $result->get_seqs();

  my $stats = {'dN' => [], 'dS' => [], 'omega' => []};

  ## Iterate over the lower triangular elements below the main diagonal.
  my $n_seqs = $align->num_sequences;
  for my $i (0 .. ($n_seqs -2))
    {
      for my $j (($i +1) .. ($n_seqs -1))
        {
          my $dn = $MLmatrix->[$i]->[$j]->{'dN'};
          my $ds = $MLmatrix->[$i]->[$j]->{'dS'};

          ## Some sequences are exactly the same and will have a dN and
          ## dS of zero, which then takes a dn/ds of 99.  This is all
          ## non-sense so we just skip those cases.  See issue #23
          next if $dn == 0 && $ds == 0;

          push @{$stats->{'dN'}}, $dn;
          push @{$stats->{'dS'}}, $ds;
          push @{$stats->{'omega'}}, $MLmatrix->[$i]->[$j]->{'omega'};
#          say join("\t", $otus[$i]->display_id,
#                         $otus[$j]->display_id,
#                         $MLmatrix->[$i]->[$j]->{'dN'},
#                         $MLmatrix->[$i]->[$j]->{'dS'},
#                         $MLmatrix->[$i]->[$j]->{'omega'});
        }
    }
  return $stats;
}

=func min_decimal_places

Returns the smallest number of decimal places in an array of numbers.

Args:
  s (ArrayRef)

Returns:
  integer
=cut
sub min_decimal_places
{
    my $s = shift;
    my @nd = map {$_ =~ /^\d+\.(\d+)$/; length $1} @{$s};
    return List::Util::min (@nd);
}


sub main
{
  if (@_ < 1)
    {
      print "Usage error -- no input arguments.\n";
      print "Correct usage is:\n";
      print "\n";
      print "  \$ align_transcripts_stats.pl path/for/transcript/align_1 path/for/transcript/align_2 path/for/transcript/align_X \n";
      exit (1);
    }

  for my $aln_file (@_)
    {
      if ((File::Spec->splitpath($aln_file))[2]
           !~ /aligned_($HistoneCatalogue::histone_regexp)_cds.fasta/)
        { die "unable to identify histone type from alignment filename '$aln_file'"; }
      my $histone = $1;

      my $aln = Bio::AlignIO->new(-file => $aln_file)->next_aln();

      ##
      ## dN/dS stats per histone type
      ##

      my $dNdS_stats = get_dNdS_stats ($aln);

      ## Use the same number of decimal places as the numbers we are
      ## averaging.  The rules of thumb for significant figures don't
      ## seem to apply when computing the average (see issue #24).

      local $Statistics::Basic::IPRES = min_decimal_places ($dNdS_stats->{'dN'});
      say HistoneCatalogue::latex_newcommand (
        "Mean".$histone."dN",
        Statistics::Basic::mean ($dNdS_stats->{'dN'}),
        "Mean of estimates of non-synonymous substitutions per "
        . "non-synonymous site (dN or Ka) between all core $histone pairs"
      );
      local $Statistics::Basic::IPRES = min_decimal_places ($dNdS_stats->{'dS'});
      say HistoneCatalogue::latex_newcommand (
        "Mean".$histone."dS",
        Statistics::Basic::mean ($dNdS_stats->{'dS'}),
        "Mean of estimates of synonymous substitutions per "
        . "synonymous site (dS or Ks) between all core $histone pairs"
      );

      local $Statistics::Basic::IPRES = min_decimal_places ($dNdS_stats->{'omega'});
      say HistoneCatalogue::latex_newcommand (
        "Mean".$histone."dNdS",
        Statistics::Basic::mean ($dNdS_stats->{'omega'}),
        "Mean of dN/dS (or Ka/Ks) ratios between all core $histone pairs"
      );
    }
  return 0;
}

main (@ARGV);
