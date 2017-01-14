#!/usr/bin/perl
## Copyright (C) 2014 Carnë Draug <carandraug+dev@gmail.com>
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

use 5.010;
use strict;
use warnings;
use File::Spec;
use File::Temp;
use List::Util;

use Bio::Tools::Run::Alignment::Clustalw;
use Statistics::Basic;
use Bio::Seq;

use WebLogo;
use MyLib; # FIXME: we should stop using this.

## This script will look at the UTR (currently, only the stem loop)
##
## It will create the following files:
##    * results/aligned_stem_loops.fasta
##    * figs/seqlogo_stem_loops.eps
##    * results/variables-utr.tex
##
## Usage is:
##
## utr_analysis.pl --sequences path/for/sequences --results path/for/results --figures path/for/figures
##
## If we are using this script, we will be using TCoffee http://www.tcoffee.org/
## and weblogo http://weblogo.threeplusone.com/ and they should be referenced
## (see the align_sequences script)

my %path = MyLib::parse_argv ("sequences", "figures", "results");

my @hdes;
my @stem_loops;
my @inits;  # start position of the stem loop, in bps after the stop codon
my @ends;   # end position of the stem loop, in bps after the stop codon
my @sl2hde; # distance between end of stem loop and start of HDE

## TODO the HDEs are currently not annotated, so we try to identify them
##      by alignment.  We are submitting them the annotations and waiting.

## to perform the pair sequence alignment and find the location of the HDE
my $clustalw = Bio::Tools::Run::Alignment::Clustalw->new(
  -matrix => 'BLOSUM',
  -quiet  => 1,
);

foreach my $gene (MyLib::load_canonical ($path{sequences})) {
  foreach my $acc (keys %{$gene->{transcripts}}) {
    my $seq = MyLib::load_seq ("transcripts", $acc, $path{sequences});

    my $sl  = ($seq->get_SeqFeatures ("stem_loop"))[0];
    my $cds = ($seq->get_SeqFeatures ("CDS"))[0];
    next if (! $sl || ! $cds);

    push (@stem_loops, $sl->seq);
    push (@inits, $sl->start - $cds->end -1);
    push (@ends,  $sl->end - $cds->end -1);

    my $gseq = MyLib::load_seq ("gene", $gene->{uid}, $path{sequences});
    my $genomic_cds = (grep {($_->get_tag_values("gene"))[0] eq $gene->{symbol}}
                            $gseq->get_SeqFeatures ("CDS"))[0];

    ## FIXME this is temporary because at the moment, the HDEs are not
    ##       annotated features.  Once they are, we should simply fetch
    ##       its details.

    ## index for the first nucleotide after the stem loop
    my $genomic_post_sl = $genomic_cds->end + $ends[-1] +2;

    my $aln = $clustalw->align([
      $HistoneCatalogue::u7_srna,
      Bio::Seq->new(
        -seq => $gseq->subseq($genomic_post_sl, $genomic_post_sl + 80),
        -id  => $seq->display_id,
      )
    ]);

    my $aln_seq_snrna = $aln->get_seq_by_pos(1)->seq;
    my @match_idx;
    push (@match_idx, pos ($aln_seq_snrna) -1) while ($aln_seq_snrna =~ m/[^\.]/g);
    push (@sl2hde, $match_idx[0]);

    push (@hdes,Bio::Seq->new(
      -seq => $gseq->subseq ($genomic_post_sl + $match_idx[ 0],
                             $genomic_post_sl + $match_idx[-1]),
      -id => $seq->display_id,
      -accession_number => $seq->accession_number,
    ));
  }
}

my $var_path = File::Spec->catdir($path{results}, "variables-utr.tex");
open (my $var_file, ">", $var_path)
  or die "Could not open $var_path for writing: $!";

say {$var_file} HistoneCatalogue::latex_newcommand (
  "StemLoopStartMode",
  Statistics::Basic::mode (@inits),
  "Mode of distances, in bp, between the end of the CDS and the start of the stem loop."
);
say {$var_file} HistoneCatalogue::latex_newcommand (
  "StemLoopEndMode",
  Statistics::Basic::mode (@ends),
  "Mode of distances, in bp, between the end of the CDS and the end of the stem loop."
);
say {$var_file} HistoneCatalogue::latex_newcommand (
  "StemLoopStartMin",
  List::Util::min (@inits),
  "Shortest distances, in bp, between the end of a CDS and the start of a stem loop."
);
say {$var_file} HistoneCatalogue::latex_newcommand (
  "StemLoopStartMax",
  List::Util::max (@inits),
  "Longest distances, in bp, between the end of a CDS and the start of a stem loop."
);
say {$var_file} HistoneCatalogue::latex_newcommand (
  "HDEsDistanceFromStemLoopMode",
  Statistics::Basic::mode (@sl2hde),
  "Mode of distances, in bp, between the HDE and the stem loop."
);

close ($var_file)
  or die "Couldn't close $var_path after writing: $!";

sub align_and_seqlogo {
  my $align_path = File::Spec->catdir($path{results}, shift);
  my $slogo_path = File::Spec->catdir($path{figures}, shift);

  my $align = Bio::Tools::Run::Alignment::Clustalw->new(
    'aformat' => 'fasta',
    'output'  => 'fasta', # do not be fooled by documentation
    'quiet'   => 1,       # do not be fooled by documentation
    'outfile' => $align_path,
  )->align(\@_);

  my $weblogo = WebLogo->new();
  $weblogo->call(
    $align_path,
    $slogo_path,
    {
      "--datatype",       "fasta",
      "--sequence-type",  "dna",
    },
  );
}

align_and_seqlogo ("aligned_stem_loops.fasta", "seqlogo_stem_loops.eps", @stem_loops);
align_and_seqlogo ("aligned_HDEs.fasta", "seqlogo_HDEs.eps", @hdes);
