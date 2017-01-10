#!/usr/bin/perl
use utf8;

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
##   align_transcripts.pl path/dbfile path/input/protein/alignment path/output/transcript/alignment
##
## DESCRIPTION
##
## To align histone transcript sequences, we actually want it to
## perform the alignment in multiples of 3.  So we use the protein
## alignment, and simply map it back to its CDS sequences.

## Takes an alignment file and calls weblogo on it.  It sets our settings
## for histone proteins and then modifies it for our prefereed visualization
## where only the changes are in black.

use strict;
use warnings;

use Bio::AlignIO;
use Bio::Align::Utilities;

use HistoneSequencesDB;

if (@ARGV != 3)
  {
    print "Usage error -- no input arguments.\n";
    print "Correct usage is:\n";
    print "\n";
    print "  \$ align_transcripts.pl path/dbfile path/input/protein/alignment path/output/transcript/alignment\n";
    exit (1);
  }

my $db = HistoneSequencesDB::read_db($ARGV[0]);
my $aa_aln = Bio::AlignIO->new(-file => $ARGV[1])->next_aln();
my $out_fpath = $ARGV[2];

my @acc = map { $_->display_id } $aa_aln->each_seq;

## Required to convert the protein alignment to mRNA, maps the protein
## accessions to the corresponding transcript Bio::Seq objects.
my %seqs;
for my $gene (@{ $db->genes })
  {
    my $products = $gene->coding_products();
    for my $t_acc (keys %{$products})
      {
        my $p_acc = $products->{$t_acc};
        $seqs{$p_acc} = $db->get_transcript_cds($t_acc);
      }
  }
my $mrna_aln = Bio::Align::Utilities::aa_to_dna_aln($aa_aln, \%seqs);

## In older versions of bioperl, aa_to_dna_aln() only changed the
## actual sequences from protein to dna but kept the protein
## display_id on the Bio::LocatableSeq objects.  When we later write
## the transcript alignment, we want to have the transcript accession
## numbers so we needed that fixed.
## See https://github.com/bioperl/bioperl-live/issues/137
##
## We worked around that by changing the display_id ourselves but the
## fix on bioperl was backwards incompatible.  So we check the output
## of aa_to_dna_aln and apply our workaround only if needed.
## See https://github.com/af-lab/histone-catalogue/issues/36
if ($mrna_aln->get_seq_by_id ((keys %seqs)[0]))
  {
    for my $aln_seq ($mrna_aln->each_seq)
      {
        $mrna_aln->remove_seq($aln_seq);

        my $mrna = $seqs{$aln_seq->display_id};
        $aln_seq->display_id($mrna->display_id);
        $aln_seq->desc($mrna->desc);

        $mrna_aln->add_seq($aln_seq);
      }
  }


my $out = Bio::AlignIO->new(
  -file   => ">$out_fpath",
  -format => 'fasta',
  -displayname_flat => 1,
);
$out->write_aln($mrna_aln);
