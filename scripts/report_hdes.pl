#!/usr/bin/perl
use utf8;

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

## This is the script we are using now to find the HDE sequences.  The idea
## is to print this in a nicely formatted table, and provide it to the NCBI
## curators so it becomes annotatted and we can use it on our paper.
##
## usage: report_hdes.pl --sequences path/for/sequences

use 5.010;
use strict;
use warnings;

use Bio::DB::EUtilities;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Seq;

use HistoneCatalogue;
use MyLib; # FIXME: we should stop using this.

my %path = MyLib::parse_argv ("sequences");

## to perform the pair sequence alignment and find the location of the HDE
my $clustalw = Bio::Tools::Run::Alignment::Clustalw->new(
  -matrix => 'BLOSUM',
  -quiet  => 1,
);

## to confirm the sequences match and we did not mess up the maths
my $fetcher = Bio::DB::EUtilities->new(
  -eutil    => 'efetch',
  -db       => 'nucleotide',
  -retmode  => 'text',
  -rettype  => 'gb',
  -email    => 'david.pinto@nuigalway.ie',
);

foreach my $gene (MyLib::load_canonical ($path{sequences})) {

  foreach my $acc (keys %{$gene->{transcripts}}) {
    my $seq = MyLib::load_seq ("transcripts", $acc, $path{sequences});

    my $sl  = ($seq->get_SeqFeatures ("stem_loop"))[0];
    my $cds = ($seq->get_SeqFeatures ("CDS"))[0];
    next if (! $sl || ! $cds);

    ## currently this is a problem on the model. Histone genes are not
    ## supposed to have two transcripts, so why is this one annotated?
    next if $seq->display_id eq "NM_001161334";

    ## This end is relative to the end of the CDS
    my $stem_loop_end = $sl->end - $cds->end -1;

    my $gseq = MyLib::load_seq ("gene", $gene->{uid}, $path{sequences});
    my $genomic_cds = (grep {($_->get_tag_values("gene"))[0] eq $gene->{symbol}}
                            $gseq->get_SeqFeatures ("CDS"))[0];

    ## index for the first nucleotide after the stem loop
    my $genomic_post_sl = $genomic_cds->end + $stem_loop_end +2;

    my $aln = $clustalw->align([
      $HistoneCatalogue::u7_srna,
      Bio::Seq->new(
        -seq => $gseq->subseq($genomic_post_sl, $genomic_post_sl + 40),
        -id  => $seq->display_id,
      )
    ]);

#  say ">" . $gene->{symbol};
#  print $_->seq . "\n" foreach $aln->each_seq();

  my $aln_seq_snrna = $aln->get_seq_by_pos(1)->seq;
  my @match_idx;
  push (@match_idx, pos ($aln_seq_snrna) -1) while ($aln_seq_snrna =~ m/[^\.]/g);


  ## The strand information has been lost because the start and end of a
  ## gene are relative to the start and end of the contig (and independent
  ## of transcription direction).  The genomic sequence files we got are
  ## the reverse complement when appropriate so the start and end properties
  ## of the CDS feature are relative to start of the sequence in file.  The
  ## only way to get the strand is to check if anywhere in the accession
  ## it is mentioned the genomic sequence is the complement sequence.
  my $strand = scalar (grep {$_->value =~ m/complement/i}
                                $gseq->get_Annotations ('secondary_accession'))
                      ? 1 : 0; # complement is strand true

  my $genomic_hde_start = $genomic_post_sl + $match_idx[ 0];
  my $genomic_hde_end   = $genomic_post_sl + $match_idx[-1];
  my $hde_seq = $gseq->subseq($genomic_hde_start, $genomic_hde_end);
  if ($strand) {
    $genomic_hde_start  = $gene->{end} - $genomic_hde_start +1;
    $genomic_hde_end    = $gene->{end} - $genomic_hde_end +1;
  } else {
    $genomic_hde_start  = $gene->{start} + $genomic_hde_start -1;
    $genomic_hde_end    = $gene->{start} + $genomic_hde_end -1;
  }

  ## if we make several requests too quickly, some of them will fail.
  ## If so, we sleep for a while before trying again.
  my $tries = 0;
  my $online_seq;
  do {
    sleep 10 if $tries > 0;
    $fetcher->set_parameters (
      -id         => $gene->{chr_acc},
      -seq_start  => $genomic_hde_start,
      -seq_stop   => $genomic_hde_end,
      -strand     => $strand ? 2 : 1,
    );
    open(my $fh, "<", \$fetcher->get_Response->content)
      or die "Could not open response content string for reading: $!";
    $online_seq = Bio::SeqIO->new(
      -fh      => $fh,
      -format  => "genbank",
    )->next_seq()->seq;
    close ($fh);
  } until ($online_seq || $tries++ > 3);

  warn "wrong sequence for " . $gene->{uid} if ($hde_seq ne $online_seq);

  say join (" ", ($gene->{uid}, $gene->{symbol}, $gene->{chr_acc},
                  $strand ? 2 : 1, $genomic_hde_start, $genomic_hde_end,
                  $online_seq));

  }
}
