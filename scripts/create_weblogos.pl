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
use FindBin;                    # Locate directory of original perl script
use File::Spec;                 # Perform operation on file names
use File::Temp;                 # Create temporary files
use Bio::SeqIO;                 # Handler for SeqIO formats
use Bio::Tools::Run::Alignment::TCoffee;  # Multiple sequence alignment with TCoffee

use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVar;                      # Load variables
use MyLib;                      # Load functions

## This script will get the downloaded sequences from all canonical histones, align
## them, use the alignment to create a sequence logo and create a LaTeX table with
## the differences between each histone.
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

## Flow of this script:
##    1 - read in all protein sequences of all canonical histone, creating a Bio::Seq
## object for each of them
##    2 - multiple sequence alignment for each
##    3 - use weblogo to create a sequence logo for each
##    4 - modify the sequence logo postscript to leave in dark only the positions
## where the sequence differs

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

foreach my $gene (MyLib::load_canonical) {
  my $symbol = $$gene{'gene symbol'};
  next unless $$gene{'protein accession'}; # skip entries with no protein acession such as pseudogenes
  next unless $symbol =~ m/^HIST(\d+)($MyVar::histone_regexp)/;
  my $histone = $2;
  my $path = File::Spec->catdir($MyVar::sequences_dir, "proteins", "$$gene{'protein accession'}.gb");
  ## We make no next_seq loop because we know there's only one sequence in those genbank file
  push (@{$multi_seq{$histone}}, Bio::SeqIO->new(-file => $path)->next_seq);
}

my $tcoffee = Bio::Tools::Run::Alignment::TCoffee->new(
                                                        'aformat' => 'fasta',
                                                        'output'  => 'fasta', # do not be fooled by documentation
                                                        'quiet'   => 1,       # do not be fooled by documentation
                                                        );
foreach (keys %multi_seq) {
  my $align_path = File::Spec->catdir($MyVar::results_dir, "aligned_$_.fasta");
  my $slogo_path = File::Spec->catdir($MyVar::figs_dir, "seqlogo_$_.eps");
  $tcoffee->outfile($align_path);
  $tcoffee->align(\@{$multi_seq{$_}});
  system ('weblogo', @weblogo_params,
          "--fin",   $align_path,
          "--fout",  $slogo_path,
          ) == 0 or die "Call to weblogo failed: $?";
  mod_weblogo_eps($slogo_path);
}

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

sub mod_weblogo_eps {
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
