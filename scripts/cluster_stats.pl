#!/usr/bin/perl
## Copyright (C) 2011 Carnë Draug <carandraug+dev@gmail.com>
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
use Getopt::Long;               # Parse program arguments
use List::Util;                 # Includes min and max

use FindBin;                    # Locate directory of original perl script
use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVar;                      # Load variables
use MyLib;                      # Load functions

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


## Read the data and create a data structure for the analysis
my %canon;          # will store the gene info
my %id_tables;      # will store data for the ID tables
my @data = MyLib::load_canonical ($path{sequences});
foreach my $gene (@data) {
  my $symbol = $$gene{'gene symbol'};

  my $cluster = $1 if $symbol =~ m/^(HIST\d+)/;
  die "Unable to identify cluster of $symbol" unless defined $cluster;
  my $histone = $1 if $symbol =~ m/^${cluster}($MyVar::histone_regexp)/;
  die "Unable to identify histone type of $symbol" unless defined $histone;

  ## to find the start and end of each cluster, we just list all the start and
  ## end coordinates. At the end, we get the min and max of them.
  push (@{$canon{$cluster}{"coordinates"}},
        $$gene{'chromosome start coordinates'},
        $$gene{'chromosome stop coordinates'});

  ## count histones (by cluster, type and coding/pseudo)
  $canon{$cluster}{"total"}++;
  $canon{$cluster}{$histone}{"total"}++;
  if ($$gene{'pseudo'}) {
    $canon{$cluster}{"pseudo"}++;
    $canon{$cluster}{$histone}{"pseudo"}++;
  } else {
    $canon{$cluster}{"coding"}++;
    $canon{$cluster}{$histone}{"coding"}++;
  }

  ## get the locus. It is not possible to get it from genomic coordinates (it should
  ## be possible to make bp_genbank_ref_extractor do it though) but we can get it
  ## from the features of the transcript source. Because of that, it only works for
  ## coding genes
  if ($$gene{"transcript accession"}) {
    my  $seq      = MyLib::load_seq("transcript", $$gene{"transcript accession"}, $path{sequences});
    my ($feature) = $seq->get_SeqFeatures("source");
    my ($locus)   = $feature->get_tag_values("map");
    push (@{$canon{$cluster}{'locus'}}, $locus);
  }

  ## it's a pain to create a good data structure here and then open it up again
  ## later. Specially considering there's genes with multiple transcripts,
  ## pseudo-genes, and non-coding transcripts. We could create a nice class for
  ## it, like we did for bp_genbank_ref_extractor, but that's overkill. Instead,
  ## we prepare the latex text in advance
  if (! exists $id_tables{$histone}{$symbol}{"name"}) {
    $id_tables{$histone}{$symbol}{"name"}  = $symbol;
    $id_tables{$histone}{$symbol}{"name"} .= " $MyVar::pseudo_mark" if $$gene{'pseudo'};
    $id_tables{$histone}{$symbol}{"uid"}   = $$gene{'gene UID'};
  }
  my $accession .= MyLib::latex_string ($$gene{"transcript accession"} || "n/a");
     $accession .= " & ";
     $accession .= MyLib::latex_string ($$gene{"protein accession"} || "n/a");
  push (@{$id_tables{$histone}{$symbol}{"accessions"}}, $accession);
}

## Make a LaTeX table with all canonical histones
my $table_path = File::Spec->catdir($path{results}, "table-histone_catalogue.tex");
open (my $table, ">", $table_path) or die "Could not open $table_path for writing: $!";

say {$table} "\\begin{ctabular}{l l l l}";
say {$table} "  \\toprule";
say {$table} "  Gene name & Gene UID & Transcript accession & Protein accession \\\\";
say {$table} "  \\midrule";
foreach my $histone (keys %id_tables) {
  my $spaced = 1; # have we left a space yet?
  print {$table} "  \\addlinespace\n" unless $spaced;
  foreach my $symbol (sort keys %{$id_tables{$histone}}) {
    print {$table} "  $id_tables{$histone}{$symbol}{'name'} & $id_tables{$histone}{$symbol}{'uid'} & ";
    ## in case of multiple transcripts and proteins, the first two columns are
    ## empty for the other rows
    my @accessions = sort @{$id_tables{$histone}{$symbol}{'accessions'}};
    print {$table} (shift (@accessions)) . " \\\\\n";
    foreach (@accessions) {
      say {$table} "      & & $_ \\\\";
    }
  }
  $spaced = 0;
}
say {$table} "  \\bottomrule";
say {$table} "\\end{ctabular}";
close($table) or die "Couldn't close $table_path after writing: $!";

## Write down results
my $stats_path = File::Spec->catdir($path{results}, "variables-cluster_stats.tex");
open (my $stats, ">", $stats_path) or die "Could not open $stats_path for writing: $!";

foreach my $cluster (keys %canon) {
  my $coord_start  = List::Util::min (@{$canon{$cluster}{"coordinates"}});
  my $coord_end    = List::Util::max (@{$canon{$cluster}{"coordinates"}});
  my $coord_length = MyLib::pretty_length (abs ($coord_start - $coord_end));

  say {$stats} MyLib::latex_newcommand ($cluster."Span", $coord_length);

  ## some genes do not have the locus well defined (will have 1q21 instead of 1q21.2)
  ## so we filter the ones without enough precision or
  @{$canon{$cluster}{"locus"}} = grep (m/\./, @{$canon{$cluster}{"locus"}});
  my $locus_start = List::Util::minstr (@{$canon{$cluster}{"locus"}});
  my $locus_end   = List::Util::maxstr (@{$canon{$cluster}{"locus"}});
  my $locus;
  if ($locus_start eq $locus_end) {
    $locus = $locus_start;
  } else {
    $locus = "$locus_start--$locus_end";
  }
  say {$stats} MyLib::latex_newcommand ($cluster."Locus", $locus);

  ## because some clusters may not have coding or pseudo genes in which case
  ## these were never initialized
  $canon{$cluster}{'coding'} //= 0;
  $canon{$cluster}{'pseudo'} //= 0;


  say {$stats} MyLib::latex_newcommand ("CodingGenesIn$cluster", $canon{$cluster}{'coding'});
  say {$stats} MyLib::latex_newcommand ("PseudoGenesIn$cluster", $canon{$cluster}{'pseudo'});
  say {$stats} MyLib::latex_newcommand ("TotalGenesIn$cluster",  $canon{$cluster}{'total'});
}

close ($stats) or die "Couldn't close $stats_path after writing: $!";
