package HistoneCatalogue;
use utf8;

## Copyright (C) 2011-2015 Carnë Draug <carandraug+dev@gmail.com>
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
use Carp;

use Bio::Root::Version;
use Bio::Tools::EUtilities;
use Bio::Tools::Run::Alignment::TCoffee;

use WebLogo;
use HistoneSequencesDB;

=var cluster_number
An integer value with the current number of known clusters.
=cut
our $cluster_number = 4;

## histones that we care about (in case one day we start caring about H1)
our @histones       = qw(H2A H2B H3 H4);
our $histone_regexp = join ("|", @histones);

## LaTeX code to mark genes as pseudo on tables
our $pseudo_mark    = '($\psi$)';

## max distance since end of CDS and start of stem-loop
our $stlp_dist      = 70;

=var stlp_seq
A regexp to identify stem-loop.  Currently, according to PMID:17531405.
=cut
our $stlp_seq = 'GG[CT][CT]CTT[CT]T[CTA]AG[GA]GCC';

=var stlp_length
Known length of the stem-loop
=cut
our $stlp_length    = 16;

=var tex_macro_name
The name of the LaTeX command that we should use for our automatically
generated values.  The idea is to control the appearance of these values
automatically from the LaTeX sources.
=cut
our $tex_macro_name = "ScriptValue";

=var u7_srna
Bio::Seq object for the human U7 sRNA sequence that binds to the HDE.

Reverse complement subsequence of human U7 small nuclear 1 (RNU7-1), snRNA
between position 8 and 17 (http://www.ncbi.nlm.nih.gov/nuccore/NR_023317)
as described on Figure 1 of Marzluff, William F., Eric J. Wagner, and Robert
J. Duronio. "Metabolism and regulation of canonical histone mRNAs: life
without a poly (A) tail." Nature Reviews Genetics 9.11 (2008): 843--854.
=cut
our $u7_srna = Bio::Seq->new(-seq =>"AAAGAGCTGT", -id => 'U7');


=func print_config_variables

Print the "configuration" variables such as BioPerl and TCoffee versions
used, date of the sequences, number of clusters, etc...

Params:
  seq_log_path - path for the bp_genbank_ref_extractor log file.

=cut

sub write_config_variables
{
  my $seq_log_path = shift;

  say latex_newcommand(
    "SequencesDate",
    get_sequences_date ($seq_log_path),
    "Date when the sequences were obtained and RefSeq was queried"
  );
  say latex_newcommand(
    "BioPerlVersion",
    $Bio::Root::Version::VERSION,
    "Version of BioPerl used"
  );
  say latex_newcommand(
    "BioEUtilitiesVersion",
    $Bio::Tools::EUtilities::VERSION,
    "Version of Bio-EUtilities used"
  );
  say latex_newcommand(
    "TCoffeVersion",
    Bio::Tools::Run::Alignment::TCoffee->new()->version(),
    "Version of TCoffee used"
  );
  say latex_newcommand(
    "WebLogoVersion",
    WebLogo->new()->version,
    "Version of WebLogo used"
  );
  say latex_newcommand(
    "NumberOfClusters",
    $cluster_number,
    "Number of histone clusters assumed"
  );
}


=func get_sequences_date

Returns the date when sequences where retrieved.  This is parsed from
the log file of bp_genbank_ref_extractor.  Note that this date may differ
from the date when the sequences are analysed.

Params:
  fpath - path to the extractor.log by bp_genbank_ref_extractor

Returns:
  String with date.
=cut

sub get_sequences_date
{
  my $fpath = shift;
  open (my $data_log, "<", $fpath)
    or croak "Could not open $fpath for reading: $!";
  my $data_header = <$data_log>; # read the first line only
  close $data_log;

  $data_header =~ m/(?<=\[)(\d\d\d\d\-\d\d\-\d\d)[\d: ]+(?=\])/;
  return $1;
}


=func say_histone_catalogue

Prints a very long LaTeX table listing all input histone genes, sorted by
histone type and symbol, their UIDs, and transcript and protein accession
numbers.

Args:
  genes (array of HistoneGene)

Returns:
  void
=cut
sub say_histone_catalogue
{
  my @genes = HistoneSequencesDB::sort_histones(@_);

  say "\\begin{ctabular}{l l l l l}";
  say "  \\toprule";
  say "  Type & Gene name & Gene UID & Transcript accession & Protein accession \\\\";
  say "  \\midrule";
  foreach my $gene (@genes)
    {
      my @cols = ($gene->histone_type(), $gene->symbol(), $gene->uid());
      if (not $gene->is_coding)
        { say "  " . join (" & ", @cols, "pseudogene", "pseudogene") . "\\\\"; }
      else
        {
          my $products = $gene->products();
          my %tex_products = map
            {
              if ($products->{$_})
                { mk_latex_string($_) => mk_latex_string($products->{$_}); }
              else
                ## Non-coding trancripts are an issue:
                ##  https://github.com/af-lab/histone-catalogue/issues/21
                { mk_latex_string($_) => mk_latex_string("non-coding"); }
            } keys %$products;

          my @transcripts = sort keys %tex_products;
          foreach my $transcript (@transcripts)
            {
              say "  " . join (" & ", @cols, $transcript,
                                      $tex_products{$transcript}) . "\\\\";
            }
        }
    }

  say "  \\bottomrule";
  say "\\end{ctabular}";
  return;
}


=func say_histone_counts

Prints a series of LaTeX newcommands with total number of histone
genes, total per cluster and histone type, variants and canonical.

Args:
  db (HistoneSequencesDB)

Returns:
  void
=cut
sub say_histone_counts
{
  my $db = shift;

  ##
  ## Canonical core counts
  ##
  my @canonical_core = $db->canonical_core;

  say HistoneCatalogue::latex_newcommand(
    "TotalCoreGenes",
    scalar (@canonical_core),
    "Total number of canonical core histone genes in the genome"
  );
  say HistoneCatalogue::latex_newcommand(
    "TotalCoreCodingGenes",
    scalar grep ($_->type eq "coding", @canonical_core),
    "Total number of canonical core histone protein coding genes in the genome"
  );
  say HistoneCatalogue::latex_newcommand(
    "TotalCorePseudoGenes",
    scalar grep ($_->type eq "pseudo", @canonical_core),
    "Total number of canonical core histone protein pseudogenes in the genome"
  );

  ##
  ## Count by histone type
  ##
  for my $histone (@HistoneCatalogue::histones)
    {
      my @this_histone = grep {$_->histone_type eq $histone} @canonical_core;
      say HistoneCatalogue::latex_newcommand(
        $histone."TotalGenes",
        scalar (@this_histone),
        "Total number of histone $histone genes"
      );

      say HistoneCatalogue::latex_newcommand(
        $histone."CodingGenes",
        scalar grep ($_->type eq "coding", @this_histone),
        "Number of histone $histone coding genes"
      );
      say HistoneCatalogue::latex_newcommand(
        $histone."PseudoGenes",
        scalar grep ($_->type eq "pseudo", @this_histone),
        "Number of histone $histone pseudogenes"
      );
    }

  ##
  ## Count by cluster number
  ##
  my %clusters = map {$_->cluster => 1} @canonical_core;
  my @clusters = sort keys %clusters;
  for my $cluster_k (@clusters)
    {
      my @this_cluster = grep {$_->cluster == $cluster_k} @canonical_core;
      say HistoneCatalogue::latex_newcommand(
        "CoreCodingGenesInHIST$cluster_k",
        scalar grep ($_->type eq "coding", @this_cluster),
        "Number of core histone coding genes in the histone cluster $cluster_k"
      );
      say HistoneCatalogue::latex_newcommand(
        "CorePseudoGenesInHIST$cluster_k",
        scalar grep ($_->type eq "pseudo", @this_cluster),
        "Number of core histone pseudogenes in the histone cluster $cluster_k"
      );
      say HistoneCatalogue::latex_newcommand(
        "TotalCoreGenesInHIST$cluster_k",
        scalar @this_cluster,
        "Total number of core histone genes in the histone cluster $cluster_k"
      );
    }

  ##
  ## Count by cluster and histone type
  ##
  for my $histone (@HistoneCatalogue::histones)
    {
      my @this_histone = grep {$_->histone_type eq $histone} @canonical_core;
      for my $cluster_k (@clusters)
        {
          my @this_cluster = grep {$_->cluster == $cluster_k} @this_histone;
          say HistoneCatalogue::latex_newcommand(
            $histone."CodingInHIST$cluster_k",
            scalar grep ($_->type eq "coding", @this_cluster),
            "Number of $histone coding genes in the histone cluster $cluster_k"
          );
          say HistoneCatalogue::latex_newcommand(
            $histone."PseudoInHIST$cluster_k",
            scalar grep ($_->type eq "pseudo", @this_cluster),
            "Number of $histone pseudogenes in the histone cluster $cluster_k"
          );
          say HistoneCatalogue::latex_newcommand(
            $histone."TotalInHIST$cluster_k",
            scalar @this_cluster,
            "Total Number of $histone genes in the histone cluster $cluster_k"
          );
        }
    }

  ##
  ## Histone variant counts
  ##
  say HistoneCatalogue::latex_newcommand(
    "TotalCoreVariantGenes",
    scalar (() = $db->variants_core),
    "Total number of core histone variants genes"
  );

}


=func mk_latex_string

Convert a string into a LaTeX usable string.  Basically, this means
escaping characters for LaTeX.  We have really basic needs for now.
Hopefully it will remain that way.

Params:
  string - A string to be converted for LaTeX.

Returns:
  The input string with escaped characters for use in LaTeX sources.
=cut

sub mk_latex_string
{
  my $string = shift;
  my %replace = (
    '&'   => '\&',
    '%'   => '\%',
    '$'   => '\$',
    '#'   => '\#',
    '_'   => '\_',
    '{'   => '\{',
    '}'   => '\}',
    '~'   => '\\textasciitilde{}',
    '^'  => '\\textasciicircum{}',
    '\\'  => '\\textbackslash{}',
  );
  $string =~ s/(&|%|\$|#|_|\{|\}|~|\^|\\)/$replace{$1}/g;
  return $string;
}


=func latex_newcommand

Return LaTeX code for new LaTeX command with documentation as a LaTeX
comment.

Creates the LaTeX code that generates a new LaTeX command with that
value as argument to $HistoneCatalogue::tex_macro_name. Because of
TeX limitations on valid identifiers, commands must be alphabetic only.
Numeric characters are replaced by their capitalized English name.  Any
other character will cause it to fail.

Params:
  name - name of the command.  Numbers get replaced, e.g., the
    name "HIST1" will be converted to "HISTOne"
  value - value corresponding to the new command.
  docs - description of the command.  It will be added as a comment
    on the previous line for documentation of the generated sources.
    Defaults to "Not documented"

Returns:
  string - with documentation as comments and \newcommand macro

Exception:
  dies if name has any non-alphanumeric character.
=cut

sub latex_newcommand
{
  my $name  = shift;
  my $value = shift;
  my $docs  = shift || "Not documented";

  if ($name =~ m/[^a-z0-9]/i)
    { croak "Unable to form LaTeX command '$name': must be alphanumeric"; }
  $value = mk_latex_string ($value);

  ## Replace numeric characters with words.  This is not meant complete,
  ## there is Lingua::EN::Nums2Words, Lingua::EN::Numbers or Number::Spell
  ## for that.  Simply replace 0-9 is enough for us.
  my %trans = (
    0 => "Zero",
    1 => "One",
    2 => "Two",
    3 => "Three",
    4 => "Four",
    5 => "Five",
    6 => "Six",
    7 => "Seven",
    8 => "Eight",
    9 => "Nine",
  );
  $name =~ s/([0-9])/$trans{$1}/g;

  ## we need to set the height of the color box manually otherwise \colorbox
  ## will change the height of the line
  my $command = "\\newcommand{\\$name}{\\${tex_macro_name}{$value}}";

  ## Append a comment character for each line of the docs
  $docs = join ("\n", map { "%% $_" } split ("\n", $docs));

  return $docs . "\n" . $command;
}


=func mk_latex_list_name_isoforms

Creates a string for use in LaTeX that describes a list of histone
isoforms gene symbols.  This is meant for use in the table of
isoforms to refer to all isoforms that encode the most common
isoform.  For example, it will create a string such as:

  HIST1H4 --B, --C, --E, --F, --I, --L; HIST2H4 --A, --B; HIST4H4

Note how the last one is not followed by any -- isoform letter, because
there isn't any.

Params:
  symbols - list of gene symbols to be used.  Must all belong to the
    same histone.

Returns:
  A string describing all the isoforms.

Exception:
  If regexp fails for some reason (probably wrong input).
=cut

sub mk_latex_list_name_isoforms
{
  my @symbols = @_;

  my %clusters;
  foreach my $symbol (@symbols)
    {
      ## Do not forget cases such as HIST4H4 where there is no isoform
      ## letter.  But also don't skip those, we want the cluster name
      ## anyway to generate the final string, we will filter out undefs
      ## later.
      if ($symbol !~ m/^((.+)($HistoneCatalogue::histone_regexp))(.*)$/i)
        { croak "Unable to list most common sequence"; }
      push (@{$clusters{$1}}, $4);
    }

  my @cluster_strs;
  foreach my $cluster (sort keys %clusters)
    {
      my @symbols = grep {$_} @{$clusters{$cluster}};
      my $str;
      if (@symbols == 0)
        { $str = $cluster; }
      elsif (@symbols == 1)
        { $str = $cluster . $symbols[0]; }
      else
        {
          @symbols = map {"--$_"} sort @symbols;
          $str = $cluster . " " . join (", ", @symbols);
        }
      push (@cluster_strs, $str);
    }
  return join ("; ", @cluster_strs);
}


=func describe_protein_variant
Create HGVS recommended description of difference between two protein sequences.

Create string with sequence difference according to nomenclature set by HGVS
at http://www.hgvs.org/mutnomen/recs-prot.html

TODO this should be added to Bioperl instead

Args:
  standard (Bio::Seq)
  variant (Bio::Seq)

Returns:
  String with difference.
=cut
sub describe_protein_variant
{
  my $standard = (shift)->seq;
  my $variant = (shift)->seq;

  if (length ($standard) != length ($variant))
    { croak "length of common and variant sequences is different"; }

  ## It is possible for the two sequences to have gaps ("-") on the same
  ## locations (they may be a pair from a multiple sequence alignment).
  ## So we start by removing those "gaps".
  while ($standard =~ m/-/g)
    {
      if (substr ($variant, $-[0], 1) eq "-")
        {
          substr ($standard, $-[0], 1) = "";
          substr ($variant, $-[0], 1) = "";
        }
    }
  my $seq_len = length ($standard);

  ## We will need this later to adjust the position numbers.
  my @gaps;
  push (@gaps, $-[0]) while ($standard =~ m/-/g);

  ## Find differences between the two sequences and store the start
  ## and end position of each **interval** of differences.
  my $mask = $standard ^ $variant;
  my @diffs;
  push (@diffs, [@-, @+]) while ($mask =~ /[^\0]+/g);

  my @desc; # the list of differences description
  foreach my $idx (@diffs)
    {
      my $start = $idx->[0];
      my $end   = $idx->[1];

      my $pre  = substr ($standard,  $start, $end - $start);
      my $post = substr ($variant, $start, $end - $start);

      ## adjust the position number due to the missing residues in $standard
      my $spos = $start +1 - (grep {$_ < $start} @gaps);  # start position
      my $epos = $end      - (grep {$_ < $end  } @gaps);  # end position

      my $str;
      if ($start == 0 && substr ($standard, 0, 1) eq "-")
        {
          ## extension of N terminus: Met1ext-5 or Met1Valext-12

          ## If variant retains the methionine/start codon of the original
          ## sequence, the extension was caused by a new start codon upstream,
          if (substr ($variant, $end, 1) eq "M")
            { $str = "M1ext-$end"; }
          ## otherwise it was forced because the original start codon is gone.
          else
            {
              my $ext_len = $end - 1;
              my $new_aa = substr ($variant, $ext_len, 1);
              $str = "M1${new_aa}ext-$ext_len";
            }
        }
      elsif ($end == $seq_len && substr ($standard, -1, 1) eq "-")
        {
          ## extension of C terminus: *110Glnext*17
          my $ext_len = $end - $spos;
          my $new_aa = substr ($variant, $start, 1);
          $str = "*${spos}${new_aa}ext*${ext_len}";
        }
      elsif ($pre !~ m/-/ && $post !~ m/-/)
        {
          ## substitution: Gly10Ser or Gly10_Met13LysCysHisVal
          $str = substr ($pre, 0, 1) . $spos;
          if ($spos != $epos)
            { $str .= "_" . substr ($pre, -1) . $epos; }
          $str .= $post;
        }
      elsif ($pre !~ m/[^-]/)
        {
          ## all insertion: Lys2_met3insGln
          ## FIXME we should differentiate with duplication
          ## FIXME we should check that before and after it's not a --
          $str = substr ($standard, $start -1, 1) . ($spos -1)
                 . "_"
                 . substr ($standard, $end, 1) . $spos
                 . "ins" . $post;
        }
      elsif ($post !~ m/[^-]/)
        {
          ## all deletion: Cys28del or Cys28_Met32del
          $str = substr ($pre, 0, 1) . $spos;
          if ($spos != $epos)
            { $str .= "_" . substr ($pre, -1) . $epos; }
          $str .= "del";
        }
      else
        {
          ## FIXME a mix of deletion and insertion and substitution... we should
          ##       have something better for this but what?

          ## remove the - from the sequence for display
          $pre  =~ s/-//g;
          $post =~ s/-//g;

          $str = substr ($pre, 0, 1) . $spos
                 . "_" . substr ($pre, -1) . $epos
                 . $post;
        }
      push (@desc, $str);
    }
  return join (" ", @desc);
}


=func most_common_seq_in_alignment

Returns the most common sequence (actual Bio::Seq->seq) from an
alignment.  If there is more than one sequence with equal
frequencies, returns the first one after sorting.

Args:
  align (Bio::Align::AlignI)

Returns:
  Bio::Seq object with most common sequence on $align
=cut
sub most_common_seq_in_alignment
{
  my $align = shift;

  my $max = 0;
  my %seqs;
  foreach ($align->each_seq)
    {
      if (++$seqs{$_->seq} > $max)
        { $max++; }
    }
  ## the sort is important for reproducibility.  If there is more than
  ## one sequence with same frequency, we want to always return the same.
  my @seqs = sort grep {$seqs{$_} == $max} keys %seqs;

  if (! @seqs) # probably only happens if $align was empty
    { croak "Unable to find most common sequence (maybe empty alignment)"; }

  return Bio::Seq->new(-seq => $seqs[0]);
}

=func say_table_isoforms_description

Prints a LaTeX table with two columns describing the differences from
each sequence against the most common sequence.  The most common sequence
is displayed at the top in a single row.  Note that we use the most common
sequence.  We do not use the consensus sequence.

* Why do not use the consensus sequence?

Using the consensus sequences is wrong because it can be a sequence that
doesn't actually exist.

A consensus sequence is the most frequent residue at _each_ position, not
the most frequent sequence.  So, how should we act with respect to insertions
and deletions?  The aligned sequences will have "-" (nothing) at such
locations, which means that even if only one of the proteins has a residue
at a certain location, the consensus sequence will keep it.  For example:

    SIHK----K
    SKHKAKGLK <-- the only that is different
    SIHK----K
    SIHK----K
    SIHK----K

    SIHKAKGLK <-- consensus sequence (different from all of them)

In this case, the consensus sequence does not actually exist.  Even
considering the empty positions (-) if it was a residue and count it
on the frequency.

* What did Marzluff used on the paper then?

To work around this cases, programs to define a consensus sequence often
have tuning parameters such as threshold.  Marzluff's paper, says that the
consensus sequence was calculated with the PRETTYBOX program, part of the
GCG package http://www.csd.hku.hk/bruhk/gcgdoc/prettybox.html which indeed
does have such paremeters.  He does not mention what parameteres were used
but should be safe to assume he used the default values.

Anyway, the consensus can still lead to a new sequence, one that is different
from all the sequences used in the alignment and I'm surprised that didn't
happen on Marzluff's paper.  And if we think about it, what we really want
to compare against is the most common sequence, so why should we be tuning
algorithms to compute a consensus that is equal to the most common sequence?

Args:
  align (Bio::Align)
  pacc2gsym (Hash) Protein ACCession 2 Gene SYMbol - a hash that maps
    the protein accession (display_id of the sequences in $align) to
    gene symbols to be used in the table.  May not actually be the gene
    symbol.  For example, genes with multiple transcripts must handle
    it somehow.

Returns:
  void
=cut
sub say_table_isoforms_description
{
  my $align = shift;
  my %pacc2gsym = @_;

  my $most_common = HistoneCatalogue::most_common_seq_in_alignment($align);

  my @eq2common;
  my %desc;
  foreach my $seq ($align->each_seq)
    {
      ## We use display_id to get the accession number because
      ## accession_number is not defined on the align object.
      my $symbol = $pacc2gsym{$seq->display_id};
      if ($seq->seq eq $most_common->seq)
        { push (@eq2common, $symbol); }
      else
        { $desc{$symbol} = HistoneCatalogue::describe_protein_variant($most_common, $seq); }
    }

  (my $most_common_seq = $most_common->seq) =~ tr/-//d; # remove the gaps

  ## tabularx has a column option "X" to set a column width to whatever
  ## space is available after the others.  We use it for the second column,
  ## the one describing the differences.

  ## The "\dimexpr\linewidth-2\tabcolsep\relax" below refers to the
  ## available space (width) for text in a table with 2 columns.  The
  ## actual space for the table is \linewidth
  ## See http://tex.stackexchange.com/a/99330/24374

  say "\\begin{tabularx}{\\linewidth}{l >{\\raggedright\\arraybackslash}X}";
  say "  \\toprule";
  say "  \\multicolumn{2}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}{Most common isoform (" .
                length ($most_common_seq) . " amino acids; " .
                HistoneCatalogue::mk_latex_list_name_isoforms (@eq2common) . ")}\\\\";
  say "  \\multicolumn{2}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}{\\texttt{\\seqsplit{$most_common_seq}}} \\\\";
  say "  \\midrule";

  foreach my $symbol (sort keys %desc)
    {
      ## Having each equal proteins that are different from the most common
      ## in a single row could be handy (easy to see the groups) but it would
      ## look horrible. Just image: the first column taking 70% of the table
      ## width because one of the different sequences has 5 gene names on it.
      say "  $symbol & " . HistoneCatalogue::mk_latex_string ($desc{$symbol}) ." \\\\";
    }

  say "  \\bottomrule";
  say "\\end{tabularx}";
}

1;
