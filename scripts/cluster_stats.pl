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
  $canon{$cluster}{$histone}{"pseudo"}++ if $$gene{'pseudo'};
  $canon{$cluster}{$histone}{"coding"}++ if !$$gene{'pseudo'};

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


## Make a LaTeX table with all canonical histones (one per type)
foreach my $histone (keys %id_tables) {
  my $ids_path = File::Spec->catdir($path{results}, "table-$histone-ids.tex");
  open (my $table, ">", $ids_path) or die "Could not open $ids_path for writing: $!";

  say {$table} "\\begin{ctabular}{l l l l}";
  say {$table} "  \\toprule";
  say {$table} "  Gene name & Gene UID & Transcript accession & Protein accession \\\\";
  say {$table} "  \\midrule";
  foreach my $symbol (sort keys %{$id_tables{$histone}}) {
    print {$table} "  $id_tables{$histone}{$symbol}{'name'} & $id_tables{$histone}{$symbol}{'uid'} &";
    foreach (sort @{$id_tables{$histone}{$symbol}{'accessions'}}) {
      print {$table} " $_ \\\\\n"
    }
  }
  say {$table} "  \\bottomrule";
  say {$table} "\\end{ctabular}";
  close($table) or die "Couldn't close $ids_path after writing: $!";
}

## Calculate:
##    * start and end coordinates of each cluster
##    * length of each cluster
for (keys %canon) {
  $canon{$_}{"start"} = List::Util::min (@{$canon{$_}{"coordinates"}});
  $canon{$_}{"end"}   = List::Util::max (@{$canon{$_}{"coordinates"}});
  my $length  = abs ($canon{$_}{'start'} - $canon{$_}{'end'});
  $canon{$_}{'length'} = MyLib::pretty_length($length);
}

### write to results file
#open($file, ">", $MyVar::results_clust) or die "Could not open $MyVar::results_clust for reading: $!";

#close($file);
#say "$_ is $canon{$_}{'length'}" for (keys %canon);
#for (keys %canon) {

#  say $canon{$_}{'start'};
#  my $remainder = length($length) % 3;
#  my $etc;
#  given ($remainder) {
#    when (0) {  $etc = sprintf("%1.3g", $length); }

#    default {  $etc = sprintf("%2.2g", $length); }
#  }
##  my $etc = sprintf("%1.${MyVar::size_precision}g", $length);
#  say "$_ has pot $pot from $length";
#say length( "23.4535433" );

#printf("%3.1f", value / pow(10, pot) )
#$pot = 0;
#while ( 10**($pot +3) ) {
#  $pot += 3;
#}
#printf("%3.1f", value / pow(10, pot) )
#use Data::Dumper;
#print Dumper %canon;
#my $number = $cluster =~ m/(\d*)$/;
#$nn =~ s/^(\w)/\U$1/g;
#say $nn;
