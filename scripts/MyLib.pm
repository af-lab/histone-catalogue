package MyLib;
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

use 5.010;                                  # use Perl 5.10
use strict;                                 # enforce some good programming rules
use warnings;                               # replacement for the -w flag, but lexically scoped
use File::Spec;                             # Perform operation on file names
use Text::CSV 1.21;                         # Comma-separated values manipulator (require 1.21 for getline_hr_all)
use POSIX;                                  # Perl interface to IEEE Std 1003.1
use Bio::SeqIO;                             # Handler for SeqIO formats
use Getopt::Long;                           # Parse program arguments
use Email::Valid;                           # Validate e-mail address

use FindBin;                                # Locate directory of original perl script
use lib $FindBin::Bin;                      # Add script directory to @INC to find 'package'
use MyVar;                                  # Load variables

## performs input check. Takes an array with the name of the required options
sub input_check {
  my %opts = (
    "sequences" => "",
    "figures"   => "",
    "results"   => "",
    "email"     => "",
  );
  GetOptions (
    "sequences=s" => \$opts{sequences},
    "figures=s"   => \$opts{figures},
    "results=s"   => \$opts{results},
    "email=s"     => sub {
      $opts{email} = Email::Valid->address($_[1])
        or die "Invalid e-mail adress $_[1]: $Email::Valid::Details";
    },
  ) or die "Error processing options. Paths must be strings";

  my %args;
  foreach (@_) {
    die "No value for $_ specified. Use the --$_ option." unless $opts{$_};
    $args{$_} = $opts{$_};
  }
  return %args;
}

## load the gene information from all genes found
sub load_csv {
  my $data_path = File::Spec->catdir($_[0], 'data.csv');
  ## To cover the widest range of parsing options, you will always want to set binary
  my $csv = Text::CSV->new ({
                              binary => 1,
                              eol    => $/,
                              }) or die "Cannot use Text::CSV: ". Text::CSV->error_diag ();
  open (my $file, "<", $data_path) or die "Could not open $data_path for reading: $!";

  $csv->column_names ($csv->getline ($file)); # read first line and sets it as the column name
  my $data = $csv->getline_hr_all ($file);    # reads all lines of file into an array of hashes (returns ref to array)
  close $file;                                  # close file

  my %genes;
  foreach my $entry (@$data) {
    my $uid = $$entry{'gene UID'};
    ## skip genes without genomic information
    if (! $$entry{'chromosome accession'}) {
      warn ("Gene with UID '$uid' has no genomic information. Skipping it!");
      next;
    }
    ## having a field for UID is not duplicating data because later we will
    ## use this to make sets of each entry, and won't have access to the key
    $genes{$uid}{'uid'}     //= $uid;
    $genes{$uid}{'symbol'}  //= $$entry{'gene symbol'};
    $genes{$uid}{'desc'}    //= $$entry{'gene name'};
    $genes{$uid}{'species'} //= $$entry{'species'};
    $genes{$uid}{'pseudo'}  //= $$entry{'pseudo'};
    $genes{$uid}{'ensembl'} //= $$entry{'EnsEMBL ID'};
    $genes{$uid}{'chr_acc'} //= $$entry{'chromosome accession'};
    $genes{$uid}{'start'}   //= $$entry{'chromosome start coordinates'};
    $genes{$uid}{'end'}     //= $$entry{'chromosome stop coordinates'};
    unless ($genes{$uid}{'pseudo'}) {
      $genes{$uid}{'transcripts'}{$$entry{'transcript accession'}} = $$entry{'protein accession'};
      $genes{$uid}{'proteins'   }{$$entry{'protein accession'   }} = $$entry{'transcript accession'};
    }
    ## if they are pseudo genes, we create an empty hash. If we don't
    ## then this would not exist and we'd have to keep cheking if it's
    ## a pseudo gene
    $genes{$uid}{'transcripts'} //= {};
    $genes{$uid}{'proteins'   } //= {};
  }
  return %genes;
}

## rather than load information from all genes found and extracted, get only
## the canonical histones
sub load_canonical {
  my %genes = load_csv (@_);
  my @canon;
  foreach my $uid (keys %genes) {
    my $symbol = $genes{$uid}{'symbol'};

    ## skip genes that don't look canonical and get cluster number
    next unless $symbol =~ m/^HIST(\d+)($MyVar::histone_regexp)/;

    $genes{$uid}{'cluster'} = $1;
    $genes{$uid}{'histone'} = $2;

    ## warn if a gene is found whose nomeclature mentions an unknown cluster
    if ($genes{$uid}{'cluster'} > $MyVar::cluster_number) {
      warn ("Update/Check the code, found possible NEW histone cluster $1 with gene '$symbol'");
    }
    push (@canon, $genes{$uid}); # pushing references to that struct into an array
  }
  return @canon;
}

## load csv but return only the H1 genes
sub load_H1 {
  my %genes = load_csv (@_);
  my @h1;
  foreach my $uid (keys %genes) {
    next unless $genes{$uid}{'symbol'} =~ m/^HIST\dH1/;
    push (@h1, $genes{$uid});
  }
  return @h1;
}

## loads a sequence file returning a Bio::Seq object
sub load_seq {
  my ($type, $access, $path) = @_;
  given ($type) {
    when (/^gene/)       {$type = "genes";}
    when (/^transcrip/)  {$type = "transcripts";}
    when (/^protein/)    {$type = "proteins";}
  }
  $path = File::Spec->catdir($path, $type, "$access.gb");
  ## we make no next_seq loop because we know there's only one sequence in
  ## those genbank files
  my $seq = Bio::SeqIO->new(-file => $path)->next_seq;
  if ($type eq "protein") {
    ## we remove the first amino acid since in histones it's cleaved off
    $seq = $seq->trunc(2, $seq->length);
  }
  return $seq;
}

## turn a distance in base pairs (a positive integer) into a string appropriate
## for text (in LaTeX format), e.g. 1000000 -> 1\,Mbp; 1540 -> 1.54\,kbp (the
## precision is defined in MyVar.pm
sub pretty_length {
  my $length   = $_[0];
  ## (ceil() -1 ) because when /3, if it's 3, floor will give 1 when we want 0
  my $power    = ( ceil(length($length) / 3) -1) * 3;
  my $dec_case = 0;
  my $number   = sprintf("%1.${dec_case}f", $length / (10 ** $power) );
  if ( length($length) < $MyVar::size_precision) {
    ## nothing, no decimal cases at all in these cases
  } elsif (length($number) < $MyVar::size_precision) {
    my $dec_case = $MyVar::size_precision - length ($number);
    $number      = sprintf("%1.${dec_case}f", $length / (10 ** $power) );
  }
  my $prefix;
  given ($power) {
    when  (0) { $prefix = ''  }
    when  (3) { $prefix = 'k' }
    when  (6) { $prefix = 'M' }
    when  (9) { $prefix = 'G' }
    when (12) { $prefix = 'T' }
    when (15) { $prefix = 'P' }
    when (18) { $prefix = 'E' }
    when (21) { $prefix = 'Z' }
    when (24) { $prefix = 'Y' }
    default   { $power -= 24; $prefix = "e+${power}Y" }
  }
  return "$number\\,${prefix}bp";
}

## escape necessary characters for latex (we will have really basic needs,
## probably only the underscore)
sub latex_string {
  ## the list of characters to escape in a look-ahead operator so they are not
  ## captured. This way the substution actually only adds a \ before them
  (my $fixed = $_[0]) =~ s/(?=[_])/\\/g;
  return $fixed;
}

## gets a string that should be printed to define a new latex command
sub latex_newcommand {
  my $command = latex_string (num2en ($_[0]));
  my $value   = $_[1];
  return "\\newcommand{\\$command}{$value}";
}


## Replaces numbers in a string by their english word, and capitalizes the
## first character. This is not meant to be correct, there's perl modules for
## that (Lingua::EN::Nums2Words, Lingua::EN::Numbers or Number::Spell). We
## really just want this to create valid LaTeX commands
sub num2en {
  my $string = shift;
  my %trans = (
               1 => "One",
               2 => "Two",
               3 => "Three",
               4 => "Four",
               5 => "Five",
               6 => "Six",
               7 => "Seven",
               8 => "Eight",
               9 => "Nine",
               0 => "Zero",
               );
  my $keys = join ('', keys %trans);
  $string =~ s/([$keys])/$trans{$1}/g;
  return $string;
}

1; # a package must return true
