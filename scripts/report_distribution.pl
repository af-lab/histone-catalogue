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

use 5.010;
use strict;
use warnings;
use List::Util qw(max min);
use Bio::DB::EUtilities;

use MyLib; # FIXME: we should stop using this.

## report_distribution.pl --sequences path/for/sequences --cluster cluster_n --email address

# Check input options
my %opts = MyLib::parse_argv("sequences", "cluster", "email");

## Load only the desired cluster
my @data = grep {
  $_->{'cluster'} == $opts{cluster}
} MyLib::load_canonical ($opts{sequences});

##
## sort by start coordinates
##
## we will need to download the sequences to obtain the real genomic
## coordinates, since the ones in file are actually the coordinates for
## what was downloaded (and we are downloading extra 500bp upstream and
## downstream of the gene). Also, we need information on the strand,
## as the direction information is lost on the CSV file

## set up fetcher for current sequence with annotations
my $fetcher = Bio::DB::EUtilities->new(
  -eutil    => 'efetch',
  -db       => 'nucleotide',
  -retmode  => 'text',
  -rettype  => 'gb',
  -email    => $opts{email},
);

foreach my $gene (@data) {
  $fetcher->set_parameters (
    -id         => $gene->{chr_acc},
    -seq_start  => $gene->{start},
    -seq_stop   => $gene->{end},
  );
  my $seq = fetch_seq ($fetcher);

  ## update chromosome acc with the version number
  $gene->{chr_acc} .= "." . $seq->version;

  foreach my $feat ($seq->get_SeqFeatures) {
    ## what we really want is the coordinates for the CDS, because the
    ## gene may be longer. However, pseudo genes do not have a CDS
    if ($gene->{pseudo}) {
      next unless $feat->primary_tag eq "gene";
    } else {
      next unless $feat->primary_tag eq "CDS";
    }

    ## there may be multiple genes in the sequence since the UTR sometimes
    ## overlap. Make sure we get the feature for the right gene
    my $symbol = ($feat->get_tag_values("gene"))[0];
    next unless uc ($symbol) eq uc ($gene->{symbol});

    ## get the location information where the correct start and end
    ## coordinates for this feature, and adjust our information
    my $location = $feat->location;
    $gene->{strand} = $location->strand;
    $gene->{end}    = $gene->{start} + $location->end -1;
    $gene->{start} += $location->start -1;
  }
}

## sort by chromossome first, and then by start coordinates
my @sorted = sort {
  $a->{chr_acc} cmp $b->{chr_acc}
    or
  $a->{start} <=> $b->{start}
} @data;


## variables that will be used in the printed reports
my $type;     # the histone type
my $symbol;   # if the gene is already annotated, its current name
my $strand;   # string + or -
my $start;    # coordinates in the subject where the HSP starts
my $end;      # coordinates in the subject where the HSP ends
my $dist;     # distance between the beginning of a HSP and the end of the last HSP
my $chr;      # accession for the chromosome (NC_000006.12 for human chromosome 6)

format CLUSTER_TOP =
Histone  Gene          Strand      Start        End      Distance  Chromosome
type     symbol                                       to previous  accession
--------------------------------------------------------------------------------
.
format CLUSTER =
@<<<     @<<<<<<<<<<<<   @    @######### @#########   @##########  @<<<<<<<<<<<<
$type,   $symbol,    $strand,     $start,      $end,        $dist, $chr
.

$^ = 'CLUSTER_TOP';
$~ = 'CLUSTER';
$= = 80;  # set page to 80 characters wide
$- = 0;   # print header again

## set end of the last HSP to the start of the first HSP
my $last_end = $sorted[0]->{start};
foreach my $gene (@sorted) {
  $type   = $gene->{histone};
  $symbol = $gene->{symbol};
  $strand = $gene->{strand} > 0 ? "+" : "-";
  $start  = $gene->{start};
  $end    = $gene->{end};
  $chr    = $gene->{chr_acc};

  $dist   = $end - $last_end;

  write;

  ## calculate for the next iteration
  $last_end = $end;
}


## takes fetcher from Bio::DB::EUtilities, fetches the first sequence,
## and returns a Bio::Seq object without any temporary file
sub fetch_seq {
  my $fetcher = shift;
  open(my $fh, "<", \$fetcher->get_Response->content)
    or die "Could not open response content string for reading: $!";
  my $seq = Bio::SeqIO->new(
    -fh      => $fh,
    -format  => "genbank",
  )->next_seq();
  close ($fh);
  return $seq;
}
