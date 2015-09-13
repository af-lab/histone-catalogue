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

## This script runs bp_genbank_ref_extractor (now part of bioperl), saves
## the sequences in the sequences directory, as well as CSV files for
## some sets.  It will generate the following files:
##
##    * canonical.csv
##    * variant.csv
##    * h1.csv
##
## Usage is:
##
## extract_sequences --email you@there.eu path_to_save_sequences

use 5.010;                      # Use Perl 5.10
use strict;                     # Enforce some good programming rules
use warnings;                   # Replacement for the -w flag, but lexically scoped
use File::Path;                 # Create or remove directory trees
use File::Spec;                 # Perform operation on file names
use Email::Valid;               # Check validity of Internet email addresses
use Text::CSV 1.21;             # Comma-separated values manipulator
use Storable;                   # persistence for Perl data structures

use FindBin;                    # Locate directory of original perl script
use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use HistoneCatalogue;
use MyLib;

## Make sure the email is valid. It is important that the email is correct
## since this script allows one to abuse (even if by accident), the NCBI
## servers who may block access. With an email address they will contact
## the user first.
my %opts = MyLib::parse_argv ("email");
$opts{email} = Email::Valid->address($opts{email})
  or die "Invalid e-mail adress $opts{email}: $Email::Valid::Details";

## Path to save the donwloaded sequences
my $seq_dir = $ARGV[0];

## remove old files to avoid problems with previous results
File::Path::remove_tree($seq_dir, {verbose => 0});

## create search string
## note that "Right side truncation with wild card does work for gene symbol" <-- from NCBI helpdesk in September 2011
my $search = '"homo sapiens"[organism] ';
$search   .= '(';
$search   .= "$_*[gene name] OR " for (@HistoneCatalogue::histones, "H1");               # get all variants
$search   .= "HIST$_*[gene name] OR " for (1 .. $HistoneCatalogue::cluster_number + 1);  # all clusters and try +1 to see if there's a new one
$search   .= 'CENPA[gene name]';                                              # CENP-A name is special
$search   .= ')';

## run sequence extractor
my @args = (
  '--assembly',     'Reference GRCh',
  '--genes',        'uid',
  '--pseudo',
  '--non-coding',
  '--upstream',     '500',
  '--downstream',   '500',
  '--transcripts',  'accession',
  '--proteins',     'accession',
  '--limit',        '300',
  '--format',       'genbank',
  '--save',         $seq_dir,
  '--save-data',    'csv',
  '--email',        $opts{'email'},
);
my @call = ($HistoneCatalogue::seq_extractor, @args, $search);
system (@call) == 0 or die "Running @call failed: $?";

my %genes = MyLib::load_csv (File::Spec->catdir ($seq_dir, "data.csv"));
my @canon = MyLib::select_canonical (%genes);
my @variants = MyLib::select_variant (%genes);
my @h1 = MyLib::select_H1 (%genes);

## I wish there was a file format for gene information but there is not.
## Because histones are simple, we can get away with a CSV file.
sub gene2csv {
  my $fpath = shift;
  my $csv = Text::CSV->new ({
    binary => 1,
    eol    => $/,
  }) or die "Cannot use Text::CSV: ". Text::CSV->error_diag ();
  open (my $fh, ">:encoding(utf8)", $fpath)
    or die "Could not open $fpath for writing: $!";

  $csv->print ($fh, ["Histone type", "Symbol", "NCBI UID", "EnsEMBL ID",
                     "chr_acc", "chr_start", "chr_end",
                     "Transcript accession", "Protein accession"]);

  foreach my $gene (@_) {
    my @line = ($gene->{histone}, $gene->{symbol}, $gene->{uid},
                $gene->{chr_acc}, $gene->{start}, $gene->{end},
                "", "");
    if ($gene->{pseudo}) {
      $csv->print ($fh, \@line);
    } else {
      while (my ($mrna, $prot) = each %{$gene->{transcripts}}) {
        @line[-2,-1] = ($mrna, $prot);
        $csv->print ($fh, \@line);
      }
    }
  }
  close $fh or die "Could not close $fpath after writing: $!";
}

for ((["canonical", \@canon], ["variant", \@variants], ["h1", \@h1])) {
  gene2csv (File::Spec->catdir ($seq_dir, $_->[0] . ".csv"), @{$_->[1]});
  Storable::store ($_->[1], File::Spec->catdir ($seq_dir, $_->[0] . ".store"));
}

