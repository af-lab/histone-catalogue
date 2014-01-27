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

use 5.010;                      # Use Perl 5.10
use strict;                     # Enforce some good programming rules
use warnings;                   # Replacement for the -w flag, but lexically scoped
use File::Spec;                 # Perform operation on file names
use File::Temp;                 # Create temporary files
use Bio::Tools::Run::Alignment::TCoffee;  # Multiple sequence alignment with TCoffee
use Statistics::Basic;          # we want to calculate mode

use FindBin;                    # Locate directory of original perl script
use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVar;                      # Load variables
use MyLib;                      # Load functions

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

my @weblogo_params = (
  "--format",         "eps",
  ## bug on weblogo 3.3. We can't set datatype but default is good, as long as we
  ## have correct file extension. See http://code.google.com/p/weblogo/issues/detail?id=32
#  "--datatype",       "fasta",
  "--sequence-type",  "dna",
  "--fineprint",      "",   # empty fineprint
  "--color-scheme",   "monochrome",
);

my @stem_loops;
my @inits;  # start position of the stem loop, in bps after the stop codon
my @ends;   # end position of the stem loop, in bps after the stop codon

foreach my $gene (MyLib::load_canonical ($path{sequences})) {
  foreach my $acc (keys $gene->{'transcripts'}) {
    my $seq = MyLib::load_seq ("transcripts", $acc, $path{sequences});

    my $sl  = ($seq->get_SeqFeatures ("stem_loop"))[0];
    my $cds = ($seq->get_SeqFeatures ("CDS"))[0];
    next if (! $sl || ! $cds);

    push (@stem_loops, $sl->seq);
    push (@inits, $sl->start - $cds->end -1);
    push (@ends,  $sl->end - $cds->end -1);
  }
}

my $var_path = File::Spec->catdir($path{results}, "variables-utr.tex");
open (my $var_file, ">", $var_path)
  or die "Could not open $var_path for writing: $!";

say {$var_file} MyLib::latex_newcommand ("StemLoopStart", Statistics::Basic::mode (@inits));
say {$var_file} MyLib::latex_newcommand ("StemLoopEnd", Statistics::Basic::mode (@ends));

close ($var_file)
  or die "Couldn't close $var_path after writing: $!";

my $align_path = File::Spec->catdir($path{results}, "aligned_stem_loops.fasta");
my $slogo_path = File::Spec->catdir($path{figures}, "seqlogo_stem_loops.eps");

my $align = Bio::Tools::Run::Alignment::TCoffee->new(
  'aformat' => 'fasta',
  'output'  => 'fasta', # do not be fooled by documentation
  'quiet'   => 1,       # do not be fooled by documentation
  'outfile' => $align_path,
)->align(\@stem_loops);

system (
  "weblogo", @weblogo_params,
  "--fin",   $align_path,
  "--fout",  $slogo_path,
) == 0 or die "Call to weblogo failed: $?";

