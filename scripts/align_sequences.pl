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
use Getopt::Long;               # Parse program arguments
use Bio::SeqIO;                 # Handler for SeqIO formats
use Bio::Tools::Run::Alignment::TCoffee;  # Multiple sequence alignment with TCoffee

use FindBin;                    # Locate directory of original perl script
use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVar;                      # Load variables
use MyLib;                      # Load functions

## This script will get the downloaded sequences from all canonical histones, align
## them, use the alignment to create a sequence logo and create a LaTeX table with
## the differences between each histone. Usage is:
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
## object for each of them
##    2 - multiple sequence alignment for each
##    3 - use weblogo to create a sequence logo for each
##    4 - modify the sequence logo postscript to leave in dark only the positions
## where the sequence differs
##    5 - compare each protein to the most common sequence, listing each difference
##    6 - make pretty LaTeX table to display it

## Check input options
my %path = ("sequences" => "",
            "figures"   => "",
            "results"   => "");
GetOptions(
            "sequences=s" => \$path{sequences},
            "figures=s"   => \$path{figures},
            "results=s"   => \$path{results},
          ) or die "Error processing options. Paths must be strings";
for (keys %path) {
  die "No path for $_ specified. Use the --$_ option." unless $path{$_};
}


my @weblogo_params = ("--units",          "probability",
                      "--format",         "eps",
                      "--show-yaxis",     "no",
                      "--stacks-per-line", 50,
                      ## bug on weblogo 3.3. We can't set datatype but default is good, as long as we
                      ## have correct file extension. See http://code.google.com/p/weblogo/issues/detail?id=32
#                      "--datatype",       "fasta",
                      "--sequence-type",  "protein",
                      "--fineprint",      "",   # empty fineprint
                      "--errorbars",      "no",
                      "--color-scheme",   "monochrome",);

my %multi_seq; # an hash of arrays of Bio::Seq objects, one for each histone
my %pacc2gsym; # maps each protein accession number to a gene symbol

foreach my $gene (MyLib::load_canonical ($path{sequences})) {
  my $symbol = $$gene{'gene symbol'};
  my $access = $$gene{'protein accession'};
  next unless $access; # skip entries with no protein acession such as pseudogenes
  next unless $symbol =~ m/^HIST(\d+)($MyVar::histone_regexp)/;
  my $histone = $2;
  $pacc2gsym{$access} = $symbol;
  my $path = File::Spec->catdir($path{sequences}, "proteins", "$access.gb");
  ## We make no next_seq loop because we know there's only one sequence in those genbank file
  my $seq = Bio::SeqIO->new(-file => $path)->next_seq;
  ## in histones, the first amino-acid, methionine, since it's cleaved off
  $seq =~ s/^M//i;
  push (@{$multi_seq{$histone}}, $seq);
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

my %aligned; # an hash of hashs, of aligned sequences (justs strings)

foreach my $histone (keys %multi_seq) {
  my $align_path = File::Spec->catdir($path{results}, "aligned_$histone.fasta");
  my $slogo_path = File::Spec->catdir($path{figures}, "seqlogo_$histone.eps");
  $tcoffee->outfile($align_path);
  my $align = $tcoffee->align(\@{$multi_seq{$histone}});
  system ('weblogo', @weblogo_params,
          "--fin",   $align_path,
          "--fout",  $slogo_path,
          ) == 0 or die "Call to weblogo failed: $?";
  mod_weblogo_eps ($slogo_path);

  ## we use display_id to get the accession number. The accession_number method
  ## is not working, and this is likely a bug on the TCoffee method which is not
  ## crateing the Bio::Seq object properly for the alignment
  $aligned{$histone}{$pacc2gsym{$_->display_id}} = $_->seq foreach ($align->each_seq);
}

## Why we did not store the consensus sequence after the alignment, and why Marzluff
## notation was wrong on the paper.
##
## Because he was calling `consensus sequence' to `most common isoform sequence'.
## A consensus sequence is the most frequent residue at _each_ position, not the most
## frequent sequence. Compared to the most frequent sequence, some isoforms will have
## inserts (in our case, specially longer C-terminal). The consensus sequences includes
## all these inserts since after the alignment, there are no other residues at those
## locations which makes it the most frequent. If indeed we use the consensus sequence
## (instead of the most common sequence), we will end up listing "mutation" for all
## proteins that have inserts or deletions.
##
## SHHK----K
## SKHKAKGLK <-- the only that is different
## SHHK----K
## SHHK----K
## SHHK----K
## SHHKAKGLK <-- consensus sequence (different from all of them)

foreach my $histone (keys %aligned) {
  my %seqs;
  ## the protein sequence is the key for an array of histone gene names that encode it
  foreach my $symbol (keys %{$aligned{$histone}}) {
    push (@{$seqs{$aligned{$histone}{$symbol}}}, $symbol);
  }
  my $max = 0;
  my $common;
  foreach (keys %seqs) {
    if (@{$seqs{$_}} > $max) {
      $max    = @{$seqs{$_}};  # number of proteins that have the sequence currently in $common
      $common = $_;            # most frequent protein sequence thus far
    }
  }

  ## Print LaTeX table
  ##
  ## On table design (from memoir class manual)
  ## In the simplest of cases a table begins with a top rule, has a single row of column
  ## headings, then a dividing rule, and after the columns of data it is finished off with
  ## a final rule. The top and bottom rules are normally set heavier (i.e., thicker or
  ## darker) than any intermediate rules.

  ## the actual sequence (without the - generated by the alignement)
  (my $common_seq = $common) =~ s/-//g;

  my $filepath = File::Spec->catdir($path{results}, "table-$histone-align.tex");
  open (my $table, ">", $filepath) or die "Couldn't open $filepath for writing: $!";

  say {$table} "\\begin{tabular}{F p{\\dimexpr(\\textwidth-\\eqboxwidth{firstentry}-4\\tabcolsep)}}";
  say {$table} "  \\toprule";
  say {$table} "  \\multicolumn{2}{p{\\dimexpr\\textwidth-2\\tabcolsep\\relax}}{Most common $histone isoform (" .
                length ($common_seq) .
                " amino acids" .
                most_common_str ($histone, @{$seqs{$common}}) .
                ")}\\\\";
  say {$table} "  \\multicolumn{2}{p{\\dimexpr\\textwidth-2\\tabcolsep\\relax}}{\\texttt{\\seqsplit{$common_seq}}} \\\\";
  say {$table} "  \\midrule";

  ## write the variance description for each sequence
  my %muts;
  foreach my $sequence (keys %seqs) {
    next if $sequence eq $common;
    my $str = MyLib::latex_string (seq_diff_str (\$common, \$sequence));
    ## not the most memory efficient but makes it easier to sort the table alphabeticaly
    $muts{$_} = $str foreach (@{$seqs{$sequence}});
  }
  foreach (sort keys %muts) {
    say {$table} "  $_ & $muts{$_} \\\\";
  }

  ## close table
  say {$table} "  \\bottomrule";
  say {$table} "\\end{tabular}";
  close($table) or die "Couldn't close $filepath after writing: $!";
}


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
  ## will be the height of the stack so not in black. We might want to highlight those.

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
sub most_common_str {
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
