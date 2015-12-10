#!/usr/bin/perl
use utf8;

## This script was highly based on bp_pairwise_kaks, by Jason Stajich,
## which is under the Perl 5 license (Artistic License 1.0 or GNU GPL)
##
## Copyright (C) 2003-2015 Jason Stajich <jason@bioperl.org>
## Copyright (C) 2015 Carnë Draug <carandraug+dev@gmail.com>
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
##   pairwise_dnds.pl path/for/dbfile
##
## DESCRIPTION
##
## It will read the store file of an HistoneSequencesDB object, and
## print to stdout Tex newcommands that define the mean dN, dS, and
## dN/dS, one for each core histone type.

use 5.010;
use strict;
use warnings;

use File::Spec;

## Codeml implements the method of Goldman and Yang (1994)
## Yn00 implements the method of Yang and Nielsen (2000)
use Bio::Tools::Run::Phylo::PAML::Codeml;

use Bio::Tools::Run::Alignment::TCoffee;    # alternatively, Clustalw
use Bio::Align::Utilities;
use Statistics::Basic;

use HistoneCatalogue;
use HistoneSequencesDB;


=func get_dNdS

From a list of Genes, it computes the pairwise dN, dS, and omega (dN/dS)
values, and returns their mean.

Since the purpose is to estimate synonymous and nonsynonymous substitutions,
we need to make sure the alignment is done as to perform multiples of 3.
So instead of aligning the CDS, we align the protein sequences, and then
convert the aligned proteins sequences back to their CDS sequences.

Args:
  $db (HistoneSequencesDB)
  @genes ([Gene])

Returns:
  {'dN' => float, 'dS' => float, 'omega' => float}
=cut
sub get_dNdS_stats
{
  my $db = shift;
  my @genes = @_;

  ## Get all the protein as Bio::Seq objects for the alignment.
  ## Also get all cds with protein accession as key.
  my %seqs;
  my @prots;
  for my $g (@genes)
    {
      my $products = $g->coding_products();
      for my $t_acc (keys %{$products})
        {
          my $p_acc = $products->{$t_acc};
          $seqs{$p_acc} = $db->get_transcript_cds($t_acc);
          my $protein = $db->get_protein($p_acc);
          push @prots, $protein;
        }
    }

  ## Align the protein sequences, and then convert the aligned protein
  ## sequence back to their original cDNA.
  my $aln_factory = Bio::Tools::Run::Alignment::TCoffee->new(-verbose => -1);
  my $aa_aln = $aln_factory->align(\@prots);
  my $dna_aln = Bio::Align::Utilities::aa_to_dna_aln($aa_aln, \%seqs);

  my $dsdn_factory = Bio::Tools::Run::Phylo::PAML::Codeml->new(
    -params => {'runmode' => -2, 'seqtype' => 1},
  );
  $dsdn_factory->alignment($dna_aln);
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
  my $n_seqs = scalar keys %seqs;
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


sub main
{
  if (@_ != 1)
    {
      print "Usage error -- no input arguments.\n";
      print "Correct usage is:\n";
      print "\n";
      print "  \$ pairwise_dnds.pl path/for/dbfile\n";
      exit (1);
    }

  my $db = HistoneSequencesDB::read_db ($_[0]);

  my @genes = $db->canonical_core();
  foreach my $histone (@HistoneCatalogue::histones)
    {
      my @this_histones = grep {$_->histone_type eq $histone
                                && $_->type eq "coding"} @genes;
      my $dNdS_stats = get_dNdS_stats ($db, @this_histones);

      ## We set IPRES to 4 for the dN and dS values, and 5 for omega,
      ## because that's the precision we get back from codeml

      local $Statistics::Basic::IPRES = 4;
      say HistoneCatalogue::latex_newcommand (
        "Mean".$histone."dN",
        Statistics::Basic::mean ($dNdS_stats->{'dN'}),
        "Mean of estimates of non-synonymous substitutions per "
        . "non-synonymous site (dN or Ka) between all core $histone pairs"
      );
      say HistoneCatalogue::latex_newcommand (
        "Mean".$histone."dS",
        Statistics::Basic::mean ($dNdS_stats->{'dS'}),
        "Mean of estimates of synonymous substitutions per "
        . "synonymous site (dS or Ks) between all core $histone pairs"
      );

      local $Statistics::Basic::IPRES = 5;
      say HistoneCatalogue::latex_newcommand (
        "Mean".$histone."dNdS",
        Statistics::Basic::mean ($dNdS_stats->{'omega'}),
        "Mean of dN/dS (or Ka/Ks) ratios between all core $histone pairs"
      );
    }
  return 0;
}

main (@ARGV);
