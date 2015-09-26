package MyLib;
## Copyright (C) 2011-2014 Carnë Draug <carandraug+dev@gmail.com>
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
use Storable;                               # persistence for Perl data structures

use HistoneCatalogue;

## This function will parse the arguments used to call the script. It
## takes an array with the name of the required options, take them
## out from the caller $ARGV, and return the leftovers. It
## will give an error if any of the requested options is missing.
##
## NOTE: this function will access and modify @ARGV on the caller.
sub parse_argv {
  ## List of all possible options accross all our scripts. Any function
  ## trying to look for an option not listed here will cause an error.
  my @possible_opt = qw(sequences figures results email reference);

  ## Initialize the hash options for all requested options
  my %req_opt = map { $_ => "" } @_;

  ## Build the parsing instructions for GetOptions
  my @parsing;
  foreach my $opt (@_) {
    if (! grep {$_ eq $opt} @possible_opt) {
      die "Unrecognized option to parse $_";
    }
    ## We don't bother doing checking of the values. That will be
    ## up to the caller.
    push (@parsing, "$opt=s" => \$req_opt{$opt});
  }
  GetOptions (@parsing) or die "Error processing options. Paths must be strings";

  ## Make sure we have values for all the requested options
  foreach (keys (%req_opt)) {
    if (! $req_opt{$_}) {
      die "No value for $_ specified. Use the --$_ option.";
    }
  }
  return %req_opt;
}

## Load the gene information from all genes found
## Returns an hash with UID values as keys. Each of the values is an hash ref
## with the following keys:
##    uid           -> UID of the gene
##    symbol        -> official gene symbol (short name)
##    desc          -> gene description (long name)
##    species       -> species name
##    pseudo        -> true if a pseudo gene, false otherwise
##    ensembl       -> EnsEMBL ID
##    chr_acc       -> chromosome accession number
##    start         -> chromosome start coordinates
##    end           -> chromosome end coordinates
##    transcripts   -> hash with accession number for transcripts as keys and
##                     its products (proteins) accession numbers as value.
##    proteins      -> hash with accession number for proteins as keys and
##                     its origins (transcripts) accession numbers as value.
##
## Why do we return an hash of hash refs rather than some class to identify
## a gene? Because we don't want to edit its values. We really only want to
## access the values after reading the file, so anything other than this is
## overkill.
sub load_csv {
  my $data_path = shift;
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
    $genes{$uid}{'symbol'}  //= uc ($$entry{'gene symbol'});
    $genes{$uid}{'desc'}    //= $$entry{'gene name'};
    $genes{$uid}{'species'} //= $$entry{'species'};
    $genes{$uid}{'pseudo'}  //= $$entry{'pseudo'};
    $genes{$uid}{'ensembl'} //= $$entry{'EnsEMBL ID'};
    $genes{$uid}{'chr_acc'} //= $$entry{'chromosome accession'};
    $genes{$uid}{'start'}   //= $$entry{'chromosome start coordinates'};
    $genes{$uid}{'end'}     //= $$entry{'chromosome stop coordinates'};

    my $nm_acc = $$entry{'transcript accession'};
    my $np_acc = $$entry{'protein accession'};
    if ($nm_acc && $np_acc) {
      $genes{$uid}{'transcripts'}{$$entry{'transcript accession'}} = $$entry{'protein accession'};
      $genes{$uid}{'proteins'   }{$$entry{'protein accession'   }} = $$entry{'transcript accession'};
    } else {
      warn ("Coding gene $uid named $genes{$uid}{'symbol'} has no protein and mRNA accession number")
        unless $genes{$uid}{'pseudo'};

      ## if they are pseudo genes, we create an empty hash. If we don't
      ## then this would not exist and we'd have to keep cheking if it's
      ## a pseudo gene
      $genes{$uid}{'transcripts'} //= {};
      $genes{$uid}{'proteins'   } //= {};
    }
  }
  return %genes;
}

sub select_canonical {
  my %genes = @_;

  my @canon;
  foreach my $uid (keys %genes) {
    my $symbol = $genes{$uid}{'symbol'};

    ## skip genes that don't look canonical and get cluster number
    next unless $symbol =~ m/^HIST(\d+)($HistoneCatalogue::histone_regexp)/i;

    $genes{$uid}{'cluster'} = $1;
    $genes{$uid}{'histone'} = $2;

    ## warn if a gene is found whose nomeclature mentions an unknown cluster
    if ($genes{$uid}{'cluster'} > $HistoneCatalogue::cluster_number) {
      warn ("Update/Check the code, found possible NEW histone cluster $1 with gene '$symbol'");
    }
    push (@canon, $genes{$uid}); # pushing references to that struct into an array
  }
  return @canon;
}

## rather than load information from all genes found and extracted, get only
## the canonical histones. Check load_csv() for the fieldnames. In addition
## the following two fieldnames are added:
##    cluster   -> cluster number (only the number. Does NOT include HIST)
##    histone   -> histone type (H2A, H2B, H3, H4)
sub load_canonical {
  my $fpath = File::Spec->catdir (shift, "canonical.store");
  my $canon_ref = Storable::retrieve ($fpath);
  return @$canon_ref;
}

sub select_H1 {
  my %genes = @_;
  my @h1;
  foreach my $uid (keys %genes) {
    next unless $genes{$uid}{'symbol'} =~ m/^HIST(\d+)H1/i;
    $genes{$uid}{'cluster'} = $1;
    $genes{$uid}{'histone'} = "H1";
    push (@h1, $genes{$uid});
  }
  return @h1;
}

## load csv but return only the H1 genes
sub load_H1 {
  my $fpath = File::Spec->catdir (shift, "h1.store");
  my $h1_ref = Storable::retrieve ($fpath);
  return @$h1_ref;
}

sub select_variant {
  my %genes = @_;
  my @variants;
  foreach my $uid (keys %genes) {
    my $symbol = $genes{$uid}{'symbol'};

    ## skip genes that don't look canonical and get cluster number
    next unless $symbol =~ m/^(($HistoneCatalogue::histone_regexp)F|CENPA)/i;

    ## $2 will be the histone if followed by F. If it's empty, then $1 will
    ## be CENPA which is a H3 variant
    if ($2) {
      $genes{$uid}{'histone'} = $2;
    } elsif ($1 eq "CENPA") {
      $genes{$uid}{'histone'} = "H3";
    } else {
      warn ("Regexp to guess histone variant type needs fixing");
      next;
    }

    push @variants, $genes{$uid};
  }
  return @variants;
}

sub load_variant {
  my $fpath = File::Spec->catdir (shift, "variant.store");
  my $variant_ref = Storable::retrieve ($fpath);
  return @$variant_ref;
}

## loads a sequence file returning a Bio::Seq object. First argument
## must be the type (gene, transcript, or protein), second argument
## its access number, and third argument the directory where to look
## for it
sub load_seq {
  my ($type, $access, $path) = @_;
  for ($type) {
    if    (/^gene/)       {$type = "genes";}
    elsif (/^transcript/) {$type = "transcripts";}
    elsif (/^protein/)    {$type = "proteins";}
  }
  $path = File::Spec->catdir($path, $type, "$access.gb");
  ## we make no next_seq loop because we know there's only one sequence in
  ## those genbank files
  my $seq = Bio::SeqIO->new(-file => $path)->next_seq;
  ## XXX  it is standard to remove the initial methionine when dealing
  ##      with histone proteins (it does not even count when giving position
  ##      in the protein). This is really really not recommended by the HGVS
  ##      (one would think is also common sense), but the number of the
  ##      amino acids is so ingrained in the field that we can't start
  ##      now to give them other numbers (everyone knows H3 K4Me3, we can't
  ##      start to correct them and call it H3 K5Me3 because we are not
  ##      important enough to break the convention).  Comment and Uncomment
  ##      the following block, to keep or remove the N-terminal methionine
  if ($type eq "proteins") {
    ## we remove the first amino acid since in histones it's cleaved off
    $seq = $seq->trunc(2, $seq->length);
  }
  return $seq;
}

## turn a distance in base pairs (a positive integer) into a string appropriate
## for text (in LaTeX format), e.g. 1000000 -> 1\,Mbp; 1540 -> 1.54\,kbp (the
## precision is defined in HistoneCatalogue.pm
sub pretty_length {
  my $length   = $_[0];
  ## (ceil() -1 ) because when /3, if it's 3, floor will give 1 when we want 0
  my $power    = ( ceil(length($length) / 3) -1) * 3;
  my $dec_case = 0;
  my $number   = sprintf("%1.${dec_case}f", $length / (10 ** $power) );
  if ( length($length) < $HistoneCatalogue::size_precision) {
    ## nothing, no decimal cases at all in these cases
  } elsif (length($number) < $HistoneCatalogue::size_precision) {
    my $dec_case = $HistoneCatalogue::size_precision - length ($number);
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

## fill the catalogue. First argument is the file path for the
## table, while the rest is an array of genes
## fill_catalogue ($path, @genes)
sub make_tex_catalogue {
  my $path = shift;
  open (my $table, ">", $path)
    or die "Could not open $path for writing: $!";

  say {$table} "\\begin{ctabular}{l l l l l}";
  say {$table} "  \\toprule";
  say {$table} "  Type & Gene name & Gene UID & Transcript accession & Protein accession \\\\";
  say {$table} "  \\midrule";

  foreach my $gene (@_) {

    print {$table} "  $$gene{'histone'} & $$gene{'symbol'} & $$gene{'uid'} & ";
    if ($$gene{'pseudo'}) {
      print {$table} "pseudogene & pseudogene \\\\\n";
    } else {
      ## In the case of a gene with multiple transcripts, each will have
      ## its line on the table but the first two columns will be empty
      my @acc_cols = map {
        HistoneCatalogue::mk_latex_string ($_ || "n/a") . " & " . HistoneCatalogue::mk_latex_string ($$gene{"transcripts"}{$_} || "n/a") . "\\\\\n"
      } (sort keys %{$$gene{'transcripts'}});
      print {$table} join ("      & & &", @acc_cols);
    }
  }

  say {$table} "  \\bottomrule";
  say {$table} "\\end{ctabular}";
  close ($table)
    or die "Couldn't close $path after writing: $!";
}

sub make_csv_catalogue {
  my $fpath = shift;
  my $csv = Text::CSV->new ({
    binary => 1,
    eol    => $/,
  }) or die "Cannot use Text::CSV: ". Text::CSV->error_diag ();

  open (my $fh, ">:encoding(utf8)", $fpath)
    or die "Could not open $fpath for writing: $!";

  $csv->print ($fh, ["Gene name", "Gene UID", "Transcript accession", "Protein accession"]);
  foreach my $gene (@_) {
    if ($$gene{'pseudo'}) {
      $csv->print ($fh, [$$gene{'histone'}, $$gene{'symbol'}, $$gene{'uid'}, 'n/a', 'n/a']);
    } else {
      while (my ($mrna, $prot) = each %{$$gene{"transcripts"}}) {
        $csv->print ($fh, [$$gene{'histone'}, $$gene{'symbol'}, $$gene{'uid'}, $mrna, $prot]);
      }
    }
  }
  close $fh or die "Could not close $fpath after writing: $!";
}

1; # a package must return true
