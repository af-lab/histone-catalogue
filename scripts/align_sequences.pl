#!/usr/bin/perl
## Copyright (C) 2011,2013 Carnë Draug <carandraug+dev@gmail.com>
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
use File::Spec;                 # Perform operation on file names
use File::Temp;                 # Create temporary files
use Bio::Tools::Run::Alignment::TCoffee;  # Multiple sequence alignment with TCoffee

use Bio::Seq;

use HistoneCatalogue;
use MyLib;

## This script will get the downloaded sequences from all canonical histones, align
## them, and create a LaTeX table with
## the differences between each histone.
##
## It will create the following files:
##    * results/aligned_H2A_proteins.fasta (one for each histone)
##    * results/aligned_H2A_cds.fasta (one for each histone)
##    * results/table-H2A-proteins-align.tex (one per histone, with the differences
##      between the different histone proteins in tabular form)
##    * results/variables-align_results.tex (LaTeX variables for the numbers
##      of unique histone proteins)
##
## Usage is:
##
## align_sequences --sequences path/for/sequences --results path/for/results --figures path/for/figures
##
## where:
##    * path/for/sequences - directory used by bp_genbank_ref_extractor to save sequences
##    * path/for/results   - directory where this script will save the LaTeX tables and aligned sequences
##    * path/for/figures   - directory where this script will save the seqlogo eps files
##
## If we are using this script, we will be using TCoffee http://www.tcoffee.org/
## We should reference it with:
##
## Notredame, C, Higgins, DG, Heringa, J (2000). "T-Coffee: A novel method for fast and
## accurate multiple sequence alignment" Journal of molecular biology, 302(1):205-218
##
##@article{tcoffee2000,
##  title={{T-Coffee}: a novel method for fast and accurate multiple sequence alignment},
##  author={Notredame, C. and Higgins, D.G. and Heringa, J. and others},
##  journal={Journal of molecular biology},
##  volume={302},
##  number={1},
##  pages={205--218},
##  year={2000},
##}
##
##
## Flow of this script:
##    1 - read in all protein sequences of all canonical histone, creating a Bio::Seq
##        object for each of them
##    2 - multiple sequence alignment for each
##    3 - compare each protein to the most common sequence, listing each difference
##    4 - make pretty LaTeX table to display it

my %path = MyLib::parse_argv("sequences", "figures", "results");

## Two hashes whose keys are the histone types (H2A, H2B, etc), and values
## are Bio::Seq objects of:
my %proteins;
my %cds;

## maps protein accession number to a gene symbol (because the
## sequence we get from the alignment only has the accession
## number and we want to use the gene symbol for the tables
my %pacc2gsym; # Protein ACCession 2 Gene SYMbol

foreach my $gene (MyLib::load_canonical ($path{sequences})) {
  my $symbol = $$gene{'symbol'};

  my @access = keys %{$$gene{proteins}};
  my $access = $access[0];
  next unless $access; # skip entries with no protein acession such as pseudogenes
  if (@access > 1) {
    ## Thank the Flying Spaghetti Monster we are working with canonical histones
    ## where each gene should have only one transcript and one protein. Any gene
    ## with reference to multiple proteins must be fixed on the databases
    warn ("Gene $symbol has more than one protein. Will use the first one ($access) only!");
  }

  $pacc2gsym{$access} = $symbol;
  push (@{$proteins{$$gene{histone}}},
        MyLib::load_seq("protein", $access, $path{sequences}));

  my $mrna = MyLib::load_seq ("transcript", $$gene{proteins}{$access}, $path{sequences});
  my $cds = ($mrna->get_SeqFeatures ("CDS"))[0];
  if (! $cds) {
    warn ("Unable to find CDS for $symbol");
    next;
  }
  push (@{$cds{$$gene{histone}}}, $cds->seq ());
}

my $var_path = File::Spec->catdir($path{results}, "variables-align_results.tex");
open (my $var_file, ">", $var_path) or die "Could not open $var_path for writing: $!";

foreach my $histone (keys %proteins) {
  my $cds_align = align(
    File::Spec->catdir($path{results}, "aligned_${histone}_cds.fasta"),
    @{$cds{$histone}}
  );

  my $protein_align = align(
    File::Spec->catdir($path{results}, "aligned_${histone}_proteins.fasta"),
    @{$proteins{$histone}}
  );
  tex_compare_histone_proteins ($var_file, $histone, $protein_align);
}
close ($var_file) or die "Couldn't close $var_path after writing: $!";


## Align all sequences with TCoffee, saving the alignment as a fasta file.
sub align
{
  my $align_path  = shift;

  ## this works but goes against the documentation. We can't fix the problem upstream because
  ## we don't know what's really wrong. Should we fix the documentation or should we fix the
  ## code? If the later ever happens, we will need to fix ours.
  ## See https://redmine.open-bio.org/issues/3406
  my $tcoffee = Bio::Tools::Run::Alignment::TCoffee->new(
    'aformat' => 'fasta',
    'output'  => 'fasta', # do not be fooled by documentation
    'quiet'   => 1,       # do not be fooled by documentation
  );

  $tcoffee->outfile($align_path);
  return $tcoffee->align(\@_);
}


sub tex_compare_histone_proteins {
  my $var_file  = shift;
  my $histone   = shift;
  my $align     = shift;

  say {$var_file} HistoneCatalogue::latex_newcommand (
    $histone."PID",
    $align->overall_percentage_identity,
    "Overall percentage identity between all histone $histone proteins"
  );

  ## Why we do not get the consensus sequence from the align object, why is it wrong to use
  ## the consensus sequence, and what did Marzluff used on the paper then?
  ##
  ## A consensus sequence is the most frequent residue at _each_ position, not the most
  ## frequent sequence. So, how should we act with respect to insertions and deletions?
  ## The aligned sequences will have "-" (nothing) at such locations, which means that
  ## even if only one of the proteins has a residue at a certain location, the consensus
  ## sequence will keep it. For example:
  ##
  ## SIHK----K
  ## SKHKAKGLK <-- the only that is different
  ## SIHK----K
  ## SIHK----K
  ## SIHK----K
  ##
  ## SIHKAKGLK <-- consensus sequence (different from all of them)
  ##
  ## In this case, the consensus sequence does not actually exist. Even considering the
  ## empty positions (-) if it was a residue and count it on the frequency. To work around
  ## this cases, programs to define a consensus sequence often have tuning parameters such
  ## as threshold. Marzluff's paper, says that the consensus sequence was calculated with
  ## the PRETTYBOX program, part of the GCG package http://www.csd.hku.hk/bruhk/gcgdoc/prettybox.html
  ## which indeed does have such paremeters. He does not mention what parameteres were used
  ## but should be safe to assume he used the default values.
  ##
  ## Anyway, the consensus can still lead to a new sequence, one that is different from all
  ## the sequences used in the alignment and I'm surprised that he did not. What we actually
  ## want to use in the tables describing the variants is the most common sequence, not the
  ## consensus.

  my $most_common = HistoneCatalogue::most_common_seq_in_alignment($align);

  ## Get a list of the genes whose sequence is equal to the most common,
  ## and the text describing the difference against it for the others.
  my @eq2common;
  my %description;
  foreach my $seq ($align->each_seq) {
    ## FIXME we use display_id to get the accession number. The accession_number
    ## method is not working, and this is likely a bug on the TCoffee method
    ## which is not creating the Bio::Seq object properly for the alignment.
    my $symbol = $pacc2gsym{$seq->display_id};
    if ($seq->seq ne $most_common->seq) {
      $description{$symbol} = HistoneCatalogue::describe_protein_variant($most_common, $seq);
    } else {
      push @eq2common, $symbol;
    }
  }

  (my $most_common_seq = $most_common->seq) =~ tr/-//d; # remove the gaps from the alignment

  my $filepath = File::Spec->catdir($path{results}, "table-${histone}-proteins-align.tex");
  open (my $table, ">", $filepath)
    or die "Couldn't open $filepath for writing: $!";

  say {$table} "\\begin{tabular}{F p{\\dimexpr(\\textwidth-\\eqboxwidth{firstentry}-4\\tabcolsep)}}";
  say {$table} "  \\toprule";
  say {$table} "  \\multicolumn{2}{p{\\dimexpr\\textwidth-2\\tabcolsep\\relax}}{Most common $histone isoform (" .
                length ($most_common_seq) . " amino acids; " .
                HistoneCatalogue::mk_latex_list_name_isoforms (@eq2common) . ")}\\\\";
  say {$table} "  \\multicolumn{2}{p{\\dimexpr\\textwidth-2\\tabcolsep\\relax}}{\\texttt{\\seqsplit{$most_common_seq}}} \\\\";
  say {$table} "  \\midrule";

  foreach my $symbol (sort keys %description) {
    ## Having each equal proteins that are different from the most common
    ## in a single row could be handy (easy to see the groups) but it would
    ## look horrible. Just image: the first column taking 70% of the table
    ## width because one of the different sequences has 5 gene names on it.
    say {$table} "  $symbol & " . HistoneCatalogue::mk_latex_string ($description{$symbol}) ." \\\\";
  }

  say {$table} "  \\bottomrule";
  say {$table} "\\end{tabular}";
  close($table) or die "Couldn't close $filepath after writing: $!";
}
