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
##   condon_usage.pl path/for/dbfile
##
## DESCRIPTION
##
## It will read the store file of an HistoneSequencesDB object, and
## print to stdout a very long Latex table with the relative frequency
## of each codon.

use 5.010;
use strict;
use warnings;

use List::Util;

use Bio::CodonUsage::Table;
use Bio::Tools::CodonTable;
use Bio::Tools::SeqStats;
use Bio::SeqUtils;

use HistoneCatalogue;
use HistoneSequencesDB;


=func compute_cds_codon_counts

Gets the counts for each codon in the cds of each gene.  Note that
HistoneSequencesDB removes the start and stop codon from the transcript
sequences.

Args:
  db (HistoneSequencesDB)
  Array of [Gene]

Returns:
  HashRef - codons are keys, values are conts

=cut
sub compute_cds_codon_counts
{
  my $db = shift;

  my $codonstats = {};
  for my $gene (@_)
    {
      for my $transcript ($gene->transcripts)
        {
          my $seq = $db->get_transcript_cds($transcript);
          my $t_codonstats = Bio::Tools::SeqStats->count_codons($seq);
          while (my ($k, $v) = each (%$t_codonstats))
            { $codonstats->{$k} += $v; }
        }
    }
  return $codonstats;
}

=func say_codon_table
Print the actual laex table with the relative frequency of each codon.

Args:
  usage - hash with table headers as keys, and Bio::CodonUsage::Table
    as values

=cut
sub say_codon_table
{
  my %usage = @_;
  my @headers = sort keys %usage;

  my $codon_table = Bio::Tools::CodonTable->new();

  ## Do not include Asx, Glx, Sel, and Xaa.
  ## Also do not include Ter because we are not including the stop codon.
  my @aa = qw(Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys
              Met Phe Pro Ser Thr Trp Tyr Val);

  say "\\begin{ctabular}{l l " . join(" ", ("l") x @headers) . "}";
  say "  \\toprule";
  say "  \\null & \\null & " . join (" & ", @headers). " \\\\";
  say "  \\midrule";

  foreach my $aa (@aa)
    {
      my @codons = map { uc $_ } $codon_table->revtranslate($aa);

      foreach my $idx (0 .. $#codons)
        {
          my $codon = $codons[$idx];
          my @freqs; # frequency of this codon for each histone.

          ## Some codons have count 0 and codon_rel_frequency()
          ## returns undef.  If all of the codons for this amino acid
          ## are undefined, then this amino acid is not present at all
          ## in which case it's NA instead of zero.  If only some are
          ## undefined, then that codon frequency is zero (TODO: fix
          ## that upstream).
          foreach my $histone (@headers)
            {
              my $table = $usage{$histone};
              my $val = $table->codon_rel_frequency ($codon);
              if (List::Util::none {defined $table->codon_rel_frequency ($_)} @codons)
                {
                  $val = "NA";
                }
              elsif (! defined $val)
                {
                  $val = 0;
                }
              push @freqs, $val;
            }

          my $tex_codon = "\\texttt{$codon}";
          ## Only print the 3 letter amino acid on its first codon
          if ($idx == 0)
            { say "  " . join (" & ", $aa, $tex_codon, @freqs) . " \\\\"; }
          else
            { say "  \\null & " . join (" & ", $tex_codon, @freqs) . " \\\\"; }
        }
    }

  say "  \\bottomrule";
  say "\\end{ctabular}";
}


sub main
{
  if (@_ != 1)
    {
      print "Usage error -- no input arguments.\n";
      print "Correct usage is:\n";
      print "\n";
      print "  \$ codon_usage.pl path/for/dbfile\n";
      exit (1);
    }

  my $db = HistoneSequencesDB::read_db ($_[0]);
  my @genes = $db->canonical_core();
  my @histone_types = @HistoneCatalogue::histones;

  my %usage;
  foreach my $histone (@histone_types)
    {
      my @this_histones = grep {$_->histone_type eq $histone
                                && $_->type eq "coding"} @genes;

      my $codon_counts = compute_cds_codon_counts ($db, @this_histones);
      $usage{$histone} = Bio::CodonUsage::Table->new(-data => $codon_counts);
    }

  say_codon_table(%usage);
  return 0;
}

main (@ARGV);
