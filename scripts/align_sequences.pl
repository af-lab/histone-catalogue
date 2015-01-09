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

use FindBin;                    # Locate directory of original perl script
use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVar;                      # Load variables
use MyLib;                      # Load functions

## This script will get the downloaded sequences from all canonical histones, align
## them, use the alignment to create a sequence logo and create a LaTeX table with
## the differences between each histone.
##
## It will create the following files:
##    * results/aligned_H2A.fasta (one for each histone)
##    * figs/seqlogo_H2A.eps (one for each histone)
##    * results/table-H2A-align.tex (one per histone, with the differences
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
## and weblogo http://weblogo.threeplusone.com/ We should reference them with:
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
## Crooks, GE, Hon, G, Chandonia, JM, Brenner SE (2004). "WebLogo: a sequence logo
## generator", Genome Research, 14:1188-1190
##
##@article{weblogo2004,
##  title={{WebLogo}: a sequence logo generator},
##  author={Crooks, G.E. and Hon, G. and Chandonia, J.M. and Brenner, S.E.},
##  journal={Genome research},
##  volume={14},
##  number={6},
##  pages={1188--1190},
##  year={2004},
##}
##
## Flow of this script:
##    1 - read in all protein sequences of all canonical histone, creating a Bio::Seq
##        object for each of them
##    2 - multiple sequence alignment for each
##    3 - use weblogo to create a sequence logo for each
##    4 - modify the sequence logo postscript to leave in dark only the positions
##        where the sequence differs
##    5 - compare each protein to the most common sequence, listing each difference
##    6 - make pretty LaTeX table to display it

my %path = MyLib::parse_argv("sequences", "figures", "results");

my @weblogo_params = (
  "--units",          "probability",
  "--format",         "eps",
  "--show-yaxis",     "no",
  "--stacks-per-line", 50,
  "--datatype",       "fasta",
  "--sequence-type",  "protein",
  "--fineprint",      "",   # empty fineprint
  "--errorbars",      "no",
  "--color-scheme",   "monochrome",
);

my %multi_seq; # an hash for arrays of Bio::Seq objects, one for each histone type

my %pacc2gsym; # maps protein accession number to a gene symbol (because the
               # sequence we get from the alignment only has the accession
               # number and we want to use the gene symbol for the tables

foreach my $gene (MyLib::load_canonical ($path{sequences})) {
  my $symbol = $$gene{'symbol'};

  ## Get protein accessions
  my @access = keys $$gene{'proteins'};
  my $access = $access[0];
  next unless $access; # skip entries with no protein acession such as pseudogenes
  if (@access > 1) {
    ## Thank the Flying Spaghetti Monster we are working with canonical histones
    ## where each gene should have only one transcript and one protein. Any gene
    ## with reference to multiple proteins must be fixed on the databases
    warn ("Gene $symbol has more than one protein. Will use the first one ($access) only!");
  }

  $pacc2gsym{$access} = $symbol;
  ## Add the protein Bio::Seq to the array of histone proteins of that type
  push (@{$multi_seq{$$gene{'histone'}}},
        MyLib::load_seq("protein", $access, $path{sequences}));
}

## this works but goes against the documentation. We can't fix the problem upstream because
## we don't know what's really wrong. Should we fix the documentation or should we fix the
## code? If the later ever happens, we will need to fix ours.
## See https://redmine.open-bio.org/issues/3406
my $tcoffee = Bio::Tools::Run::Alignment::TCoffee->new(
  'aformat' => 'fasta',
  'output'  => 'fasta', # do not be fooled by documentation
  'quiet'   => 1,       # do not be fooled by documentation
);

## Align all sequences with TCoffee, saving the alignment as a fasta file.
## Then use that file to create a logo (with WebLogo) as an eps file.

my $var_path = File::Spec->catdir($path{results}, "variables-align_results.tex");
open (my $var_file, ">", $var_path) or die "Could not open $var_path for writing: $!";

my %aligned; # an hash of hashs, of aligned sequences (just strings)
foreach my $histone (keys %multi_seq) {
  my $align_path = File::Spec->catdir($path{results}, "aligned_$histone.fasta");
  my $slogo_path = File::Spec->catdir($path{figures}, "seqlogo_$histone.eps");
  $tcoffee->outfile($align_path);
  my $align = $tcoffee->align(\@{$multi_seq{$histone}});

  say {$var_file} MyLib::latex_newcommand ($histone."PID" ,
     sprintf ("%.${MyVar::size_precision}f", $align->overall_percentage_identity));

  system (
    "weblogo", @weblogo_params,
    "--fin",   $align_path,
    "--fout",  $slogo_path,
  ) == 0 or die "Call to weblogo failed: $?";
  mod_weblogo_eps ($slogo_path);

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


  ## Getting the value of $max here, even though the values for each key is
  ## already the number of times for each sequence, saves us having to loop
  ## through the values later to find it.
  my $max = 0;
  my %seqs; # keys will be aligned sequences
  foreach ($align->each_seq) {
    ## increment the value in the hash everytime and also increment the
    ## value of $max if we go above it
    $max++ if (++$seqs{$_->seq} > $max);
  }
  my @common = grep ($seqs{$_} == $max, keys %seqs);
  if (@common > 1) {
    my $n = @common;
    warn ("Found $n `most common' sequences for $histone. Only the first will be used.");
  }
  my $most_common = $common[0];

  say {$var_file} MyLib::latex_newcommand ($histone."UniqueProteins" , scalar keys %seqs);

  ## Get a list of the genes whose sequence is equal to the most common,
  ## and the text describing the difference against it for the others.
  my @eq2common;
  my %description;
  foreach my $seq ($align->each_seq) {
    ## FIXME we use display_id to get the accession number. The accession_number
    ## method is not working, and this is likely a bug on the TCoffee method
    ## which is not creating the Bio::Seq object properly for the alignment.
    my $symbol = $pacc2gsym{$seq->display_id};
    if ($seq->seq ne $most_common) {
      $description{$symbol} = seq_diff_str (\$most_common, \$seq->seq);
    } else {
      push @eq2common, $symbol;
    }
  }

  $most_common =~ tr/-//d; # remove the gaps from the alignment

  ## Print LaTeX table
  my $filepath = File::Spec->catdir($path{results}, "table-$histone-align.tex");
  open (my $table, ">", $filepath) or die "Couldn't open $filepath for writing: $!";

  say {$table} "\\begin{tabular}{F p{\\dimexpr(\\textwidth-\\eqboxwidth{firstentry}-4\\tabcolsep)}}";
  say {$table} "  \\toprule";
  say {$table} "  \\multicolumn{2}{p{\\dimexpr\\textwidth-2\\tabcolsep\\relax}}{Most common $histone isoform (" .
                length ($most_common) . " amino acids" .
                names_for_most_common ($histone, @eq2common) . ")}\\\\";
  say {$table} "  \\multicolumn{2}{p{\\dimexpr\\textwidth-2\\tabcolsep\\relax}}{\\texttt{\\seqsplit{$most_common}}} \\\\";
  say {$table} "  \\midrule";

  foreach my $symbol (sort keys %description) {
    ## Having each equal proteins that are different from the most common
    ## in a single row could be handy (easy to see the groups) but it would
    ## look horrible. Just image: the first column taking 70% of the table
    ## width because one of the different sequences has 5 gene names on it.
    say {$table} "  $symbol & " . MyLib::latex_string ($description{$symbol}) ." \\\\";
  }

  say {$table} "  \\bottomrule";
  say {$table} "\\end{tabular}";
  close($table) or die "Couldn't close $filepath after writing: $!";
}
close ($var_file) or die "Couldn't close $var_path after writing: $!";

## function that will modify the weblogo eps file
##    usage: mod_weblogo_eps (path_to_eps_file)
sub mod_weblogo_eps {
  ## It seems to me that the postscript files generated by weblogo are: 1) options set as variable
  ## at the top; 2) followed by a piece of code common to all usage; and 3) the stacks at the end
  ## where all characters and height values are. Technically, we could create the file ourselves,
  ## no need for weblogo. We would just copy the start of one, and use it as template. The stacks
  ## should be trivial to create. Or we could even use weblogo to create the stacks and replace
  ## everything that comes before before. But that's just duplicating code which should be avoided.
  ## It means that changes on their side can break ours but the right solution is to use autoconf
  ## and change our code depending on the weblogo version.
  ##
  ## So we process the postscript file and change the color used for the letter depending on its
  ## height. This is done in the operator DrawChar (called from ShowSymbol). We change the color
  ## right before the call to show() to minimize problems if something changes in the future. There
  ## is only line with "tc show" in the whole file. It would be easy to make sure we are in the
  ## DrawChar operator but if the code changes so much that there's two "tc show" in the file, then
  ## it's time for us to update this.
  ##
  ## FIXME: in the case of amino acids that are only present in some of the proteins, their height
  ## will be the height of the stack so not in black. We probably should highlight those as well.

  open (my $read, "<", $_[0]) or die "Couldn't open $_[0] for reading: $!\n";
  my @weblogo = <$read>;
  close($read) or die "Couldn't close $_[0] after reading: $!";
  ## we are overwriting the file
  open (my $save, ">", $_[0]) or die "Couldn't open $_[0] for writing: $!\n";
  my $fixed = 0;
  foreach (@weblogo) {
    if ($_ =~ m/tc show/) {
      if ($fixed) {
        die "Trying to fix previously fixed sequence Logo. Probably new version of weblogo and we need to update our script.\n"
      }
      print {$save} "        ysize stack_height stack_margin sub ge {\n" .
                    "            0.7 0.7 0.7 setrgbcolor\n" .
                    "        } {\n" .
                    "            0.0 0.0 0.0 setrgbcolor\n" .
                    "        }ifelse\n";
      $fixed = 1;
    }
    print {$save} $_;
  }
  close($save) or die "Couldn't close $_[0] after writing: $!";
}

## Create string with sequence difference according to nomenclature set by HGVS
## at http://www.hgvs.org/mutnomen/recs-prot.html
##    usage: $string = seq_diff_str ($normal, $variant)
sub seq_diff_str {
  ## TODO this should be made into a bioperl module

  my ($common, $variant) = @_;
  my $str;

  ## finding differences between the sequences and store start and end position
  ## of each interval of differences
  my $mask = $$common ^ $$variant;
  my @pos;
  push (@pos, [@-, @+]) while ($mask =~ /[^\0]+/g);

  ## positions of all - in $$common to adjust position number
  my @pos_adj;
  push (@pos_adj, @-) while ($$common =~ m/-/g);

  foreach my $idx (@pos) {
    my $start = ${$idx}[0];
    my $end   = ${$idx}[1];

    my $pre  = substr ($$common,  $start, $end - $start);
    my $post = substr ($$variant, $start, $end - $start);

    ## adjust the position number due to the missing residues in $ori
    my $spos = $start +1 - (grep {$_ < $start} @pos_adj);  # start position
    my $epos = $end      - (grep {$_ < $end  } @pos_adj);  # end position

    if ($pre !~ m/-/ && $post !~ m/-/) {
      ## substitution: Gly10Ser or Gly10_Met13LysCysHisVal
      $str .= substr ($pre, 0, 1) . $spos;
      $str .= "_" . substr ($pre, -1) . $epos if ($spos != $epos);
      $str .= $post;
    } elsif ($pre !~ m/[^-]/) {
      ## all insertion: Lys2_met3insGln
      ## FIXME we should differentiate with duplication
      ## FIXME we should check that before and after it's not a --
      $str .= substr ($$common, $start -1, 1) . ($spos -1) .
              "_" .
              substr ($$common, $end, 1) . $spos .
              "ins" . $post;
    } elsif ($post !~ m/[^-]/) {
      ## all deletion: Cys28del or Cys28_Met32del
      $str .= substr ($pre, 0, 1) . $spos;
      $str .= "_" . substr ($pre, -1) . $epos if ($spos != $epos);
      $str .= "del";
    } else {
      ## FIXME a mix of deletion and insertion and substitution... we should
      ##       have something better for this but what?

      ## remove the - from the sequence for display
      $pre  =~ s/-//g;
      $post =~ s/-//g;

      $str .= substr ($pre, 0, 1) . $spos;
      $str .= "_" . substr ($pre, -1) . $epos;
      $str .= $post;
    }
    $str .= " ";
  }
  return $str;
}

## From a list of histone symbols, creates a string listing them for display
## of the type "HIST1H2A; --A, --B, --D; HIST2H2A: --B". Goes on top of the
## table of differences between isoforms.
##    usage: $str = most_comon_str ($histone_type, @list_of_gene_names)
sub names_for_most_common {
  my $histone = shift (@_);
  my @genes   = sort (@_);
  my ($str, $cluster) = ("", "");
  foreach my $symbol (@genes) {
    ## in some cases, there may be no isoform letter, for example, HIST4H4 so we
    ## use (.*) instead of (.+) at the end
    die "Unable to list most common sequence." unless $symbol =~ m/^(.+)$histone(.*)$/;
    if ($cluster ne $1) {
      $cluster = $1;
      $str    .= "; $cluster$histone";
      $str    .= " " if $2;
    }
    $str .= ", --$2" if $2;
  }
  return $str;
}
