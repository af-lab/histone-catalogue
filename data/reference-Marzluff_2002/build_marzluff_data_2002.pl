#!/usr/bin/perl
## Copyright (C) 2013 Carnë Draug <carandraug+dev@gmail.com>
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
use Text::CSV;                  # Comma-separated values manipulator
use File::Spec;                 # Perform operation on file names
use Bio::Seq;                   # Sequence object, with features
use Bio::SeqIO;                 # Handler for SeqIO formats
use FindBin;                    # Locate directory of original perl script
use Bio::DB::EUtilities;        # Retrieve entries from Entrez
use List::Util;                 # General-utility list subroutines

## This script builds the directory structure and CSV file with the data from
## our reference (Marzluff paper from 2002), as if it was downloaded with
## our tool. By doing this, our code to compare the latest data with this
## reference, can be used to compare against our own data later, or to see
## how much it changed over time.

my $email = ''; # used to download sequences from Genbank

## Because this is meant to be only used once anyway, we don't remove the
## directory if it's already there. Actually, exiting with an error means
## we are doing something funny. Remove it manually if that's really what
## you want to do.
my $proteins_dir = File::Spec->catdir($FindBin::Bin, "proteins");
mkdir ($proteins_dir) or die "Could not mkdir $proteins_dir: $!";

my %histone; # will hold all the histones

my %table1 = (
  ## from Table 1 (page 488 and 489)
  coding => {
    H2A => {
      HIST1H2AA => "AY131982",
      HIST1H2AB => "AY131983",
      HIST1H2AC => "AY131984",
      HIST1H2AD => "AY131985",
      HIST1H2AE => "AY131986",
      HIST1H2AG => "AY131987",
      HIST1H2AH => "AY131988",
      HIST1H2AI => "AY131989",
      HIST1H2AJ => "AY131990",
      HIST1H2AK => "AY131991",
      HIST1H2AL => "AY131992",
      HIST1H2AM => "AY131993",
      HIST2H2AA => "AY131971",
      HIST2H2AB => "AY131972",
      HIST2H2AC => "AY131973",
      HIST3H2A  => "AY131974",
    },
    H2B => {
      HIST1H2BA => "AF531284",
      HIST1H2BB => "AF531285",
      HIST1H2BC => "AF531286",
      HIST1H2BD => "AF531287",
      HIST1H2BE => "AF531288",
      HIST1H2BF => "AF531289",
      HIST1H2BG => "AF531290",
      HIST1H2BH => "AF531291",
      HIST1H2BI => "AF531292",
      HIST1H2BJ => "AF531293",
      HIST1H2BK => "AF531294",
      HIST1H2BL => "AF531295",
      HIST1H2BM => "AF531296",
      HIST1H2BN => "AF531297",
      HIST1H2BO => "AF531298",
      HIST2H2BE => "AY131979",
      HIST3H2BB => "AY131981",
    },
    H3 => {
      HIST1H3A  => "AF531274",
      HIST1H3B  => "AF531275",
      HIST1H3C  => "AF531276",
      HIST1H3D  => "AF531277",
      HIST1H3E  => "AF531278",
      HIST1H3F  => "AF531279",
      HIST1H3G  => "AF531280",
      HIST1H3H  => "AF531281",
      HIST1H3I  => "AF531282",
      HIST1H3J  => "AF531283",
      HIST2H3C  => "AF531307",
      HIST3H3   => "AF531308",
    },
    H4 => {
      HIST1H4A  => "AY128654",
      HIST1H4B  => "AY128655",
      HIST1H4C  => "AY128656",
      HIST1H4D  => "AY128657",
      HIST1H4E  => "AY128658",
      HIST1H4F  => "AY128659",
      HIST1H4G  => "AY128660",
      HIST1H4H  => "AY128661",
      HIST1H4I  => "AY128662",
      HIST1H4J  => "AY128663",
      HIST1H4K  => "AY128664",
      HIST1H4L  => "AY128665",
      HIST2H4   => "AF525682",
      HIST4H4   => "AY128653",
    },
  },

  ## from Table 1 (page 488 and 489)
  pseudo => {
    H2A => {},
    H2B => {
      HIST2H2BA => "AY131975",
      HIST2H2BB => "AY131976",
      HIST2H2BC => "AY131977",
      HIST2H2BD => "AY131978",
      HIST3H2BA => "AY131980",
    },
    H3 => {
      HIST2H3A  => "AF531305",
      HIST1H3B  => "AF531306",
    },
    H4 => {},
  },
);

## H2A consensus sequence from Table 2: Histone H2A protein variants (page 494)
my $H2A_consensus = "SGRGKQGGKARAKAKTRSSRAGLQFPVGRVHRLLRKGNYAERVGAGAPVYLAAVLEYLTAEILE".
                    "LAGNAARDNKKTRIIPRHLQLAIRNDEELNKLLGKVTIAQGGVLPNIQAVLLPKKTESHHKAKGK";

my @table2 = qw(HIST1H2AA HIST1H2AB HIST1H2AC HIST1H2AD HIST1H2AE HIST1H2AF HIST1H2AH
                HIST1H2AI HIST1H2AJ HIST1H2AK HIST1H2AL HIST1H2AM HIST2H2AA HIST2H2AB
                HIST2H2AC HIST3H2A);

## Fill all the H2A sequences with the consensus sequence, and then
## apply the differences to each of them.
$histone{$_} = $H2A_consensus foreach @table2;

for ("HIST1H2AA") {
  substr ($histone{$_},  13, 1) = "S";
  substr ($histone{$_},  15, 1) = "S";
  substr ($histone{$_},  29, 1) = "I";
  substr ($histone{$_},  42, 1) = "I";
  substr ($histone{$_},  70, 1) = "S";
  substr ($histone{$_},  98, 1) = "G";
  substr ($histone{$_}, 124, 6) = "HKAQSK";
}
substr ($histone{HIST1H2AB},  98, 1) = "R";
substr ($histone{HIST1H2AC},  15, 1) = "S";
substr ($histone{HIST1H2AC},  98, 1) = "R";
substr ($histone{HIST1H2AE},  98, 1) = "R";
substr ($histone{HIST1H2AH}, 123, 4) = "HKAK";
chop   ($histone{HIST1H2AH});
substr ($histone{HIST1H2AJ}, 123, 4) = "HKTK";
chop   ($histone{HIST1H2AJ});
substr ($histone{HIST2H2AA},  15, 1) = "S";
substr ($histone{HIST2H2AA},  50, 1) = "M";
for ("HIST2H2AB") {
  substr ($histone{$_},  15, 1) = "S";
  substr ($histone{$_},  50, 1) = "M";
  substr ($histone{$_},  86, 1) = "V";
  substr ($histone{$_},  98, 1) = "G";
  substr ($histone{$_}, 123, 6) = "KPGPNK";
}
substr ($histone{HIST2H2AC},  15, 1) = "S";
substr ($histone{HIST2H2AC},  50, 1) = "M";
substr ($histone{HIST2H2AC}, 127, 1) = "S";
substr ($histone{HIST3H2A},   15, 1) = "S";
substr ($histone{HIST3H2A},   50, 1) = "M";
substr ($histone{HIST3H2A},   98, 1) = "R";

## H2B consensus sequence from Table 3: Histone H2B protein variants (page 495)
my $H2B_consensus = "PEPAKSAPAPKKGSKKAVTKAQKKDGKKRKRSRKESYSVYVYKVLKQVHPDTGISSKAMGIM".
                    "NSFVNDIFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSSK";

my @table3 = qw(HIST1H2BA HIST1H2BB HIST1H2BC HIST1H2BD HIST1H2BE HIST1H2BF HIST1H2BG
                HIST1H2BH HIST1H2BI HIST1H2BJ HIST1H2BK HIST1H2BL HIST1H2BM HIST1H2BN
                HIST1H2BO HIST2H2BE HIST3H2BB);

# Fill all the H2B sequences with the consensus sequence, and then
# apply the differences to each of them.
$histone{$_} = $H2B_consensus foreach @table3;

for ("HIST1H2BA") {
  substr ($histone{$_},   5, 1) = "G";
  substr ($histone{$_},   7, 3) = "TIS";
  substr ($histone{$_},  13, 1) = "F";
  substr ($histone{$_},  18, 1) = "V";
  substr ($histone{$_},  20, 1) = "T";
  substr ($histone{$_},  24, 1) = "E";
  substr ($histone{$_},  31, 1) = "T";
  substr ($histone{$_},  38, 1) = "I";
  substr ($histone{$_},  40, 1) = "I";
  substr ($histone{$_},  59, 1) = "S";
  substr ($histone{$_},  66, 1) = "T";
  substr ($histone{$_},  74, 1) = "S";
  substr ($histone{$_},  83, 1) = "S";
  substr ($histone{$_},  89, 1) = "S";
  substr ($histone{$_}, 124, 1) = "N";
  substr ($histone{$_},   2, 2) = "VSS"; # at the end because of the addition
}
substr ($histone{HIST1H2BB},  17, 1) = "I";
substr ($histone{HIST1H2BB},  38, 1) = "I";
substr ($histone{HIST1H2BD},   3, 1) = "T";
substr ($histone{HIST1H2BH},   1, 1) = "D";
substr ($histone{HIST1H2BJ}, 123, 1) = "A";
substr ($histone{HIST1H2BK}, 123, 1) = "A";
substr ($histone{HIST1H2BL},   2, 1) = "L";
substr ($histone{HIST1H2BL},  74, 1) = "S";
substr ($histone{HIST1H2BM},   3, 1) = "V";
substr ($histone{HIST1H2BM},   8, 1) = "V";
substr ($histone{HIST1H2BM},  17, 1) = "I";
substr ($histone{HIST1H2BM},  18, 1) = "N";
substr ($histone{HIST1H2BN},   3, 1) = "S";
substr ($histone{HIST1H2BO},   1, 1) = "D";
substr ($histone{HIST1H2BO},  38, 1) = "I";
substr ($histone{HIST2H2BE},  38, 1) = "I";
for ("HIST3H2BB") {
  substr ($histone{$_},   1, 1) = "D";
  substr ($histone{$_},   3, 1) = "S";
  substr ($histone{$_},  31, 1) = "G";
  substr ($histone{$_},  38, 1) = "I";
  substr ($histone{$_},  74, 1) = "S";
  substr ($histone{$_},  93, 1) = "V";
}

## There was no actual H3 protein sequence on the paper, we had to deduce and
## track it down on the references.
##
## It is written on page 489:
##
##    ... the histone H3 genes encode the two variants previously described [15] ...
##
## The reference [15] is "Zweidler, A. (1984). Core histone variants of the mouse: primary
## structure and expression. In Histone Genes: Structure, Organization and Regulation (G. Stein,
## W.Stein W. F. Marzluff, Eds), pp.373-395. John Wiley and Sons, New York" (the correct page
## numbers are actually 339-371).
##
## Later, on page 494:
##
##    There are two variants encoded by the histone H3 genes, H3.1 and H3.2. These
##    differ in a single amino acid change at position 96, a cysteine in H3.1 genes
##    and a serine in H3.2 genes. The functional histone H3 gene (HIST2H3C) in the
##    HIST2 cluster on chromosome 1 encodes the H3.2 protein and all 11 histone H3
##    genes in the HIST1 cluster encode the H3.1 protein.
##
## However, there are only 10 histone H3 genes for the HIST1 cluster on Table 1. Also, it says
## nothing about the protein encoded by HIST3H3. That is mentioned on page 492:
##
##    HIST3H3, described previously as a solitary gene [19], encodes a novel histone
##    H3 protein, with the serine at position 96 and four additional changes (A24V,
##    V71M, A98S, A111V)
##
## This leads to a total of 3 unique H3 protein sequences: all the H3 on HIST1 encoding the same
## thing; a single coding H3 gene on HIST2 with just one amino acid difference; and the lone
## H3 gene on HIST3 with the five changes.

foreach (keys $table1{'coding'}{'H3'}) {
  ## sequence from Zweidler, A (1984) chapter on the histone genes book
  $histone{$_} = "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKST".
                 "ELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEACEAYLVGLFEDTNLCAIHAKR".
                 "VTIMPKDIQLARRIRGERA";
}
substr ($histone{HIST2H3C}, 95, 1) = "S";
for ("HIST3H3") {
  substr ($histone{$_},  23, 1) = "V";
  substr ($histone{$_},  70, 1) = "M";
  substr ($histone{$_},  95, 1) = "S";
  substr ($histone{$_},  97, 1) = "S";
  substr ($histone{$_}, 110, 1) = "V";
}

## There was also no H4 protein sequence on the paper. Instead, on the caption of Table 1 it is said
##
##    The human H4 gene described by Akasaka and coworkers [56] is HISTH4I[sic]
##
## Then, on page 489, it is written
##
##    all the histone H4 genes encode the same protein
##
## thus logically stating that all have the same sequence from reference 56 which targets "Akasaka, T.,
## et al. (1997). A recurring translocation, t(3;6)(q27;p21), in non-Hodgkin's lymphoma results in
## replacement of the 5' regulatory region of BCL6 with a novel H4 histone gene. Cancer Res. 57: 7--12."
foreach (keys $table1{'coding'}{'H4'}) {
  ## sequence from Akasaka, T et al (1997)
  $histone{$_} = "SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVL".
                 "KVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG";
}

## Fill the proteins directory with sequence files
while (my ($symbol, $seq) = each %histone) {
  ## Histones sequences are usually given without the first methionine.
  ## Put it back before saving it as a file.
  $seq = "M" . $seq;
  ## Remove newlines from the sequences because it confuses Bio::Seq
  ## We name the file based on the symbol because we don't have UIDs or
  ## accession numbers for all of them.
  Bio::SeqIO->new(
    -format => 'genbank',
    -file   => ">" . File::Spec->catdir ($proteins_dir, $symbol.".gb"),
  )->write_seq(Bio::Seq->new(
    -seq        => $seq,
    -display_id => $symbol,
  ));
}

## Create a data.csv file as if generated by out script. But this time, the
## only data we have is the gene symbol and whether they are coding or pseudo
## genes.
my $data_path = File::Spec->catdir($FindBin::Bin, "data.csv");
open (my $data, ">", $data_path) or die "Could not open $data_path for writing: $!";
my $csv = Text::CSV->new ({
  binary => 1,
  eol    => $/,
}) or die "Cannot use Text::CSV: ". Text::CSV->error_diag ();

$csv->print ($data, ['gene symbol', 'pseudo']);
$csv->print ($data, [$_, 0]) foreach keys %histone;
foreach my $type (keys $table1{'pseudo'}) {
  $csv->print ($data, [$_, 1]) foreach keys $table1{'pseudo'}{$type};
}
close ($data) or die "Couldn't close $data_path after writing: $!";

## Just for curiosity, we compared the sequences mentioned throughout the
## paper against the sequences obtained from their accession numbers.
## Does the sequences from the accession number on table 1, match what was
## actually written on the paper? The answer is no, there are 16 mistakes.
my $mistakes = 0;
foreach my $type (sort keys $table1{'coding'}) {
  my $fetcher = Bio::DB::EUtilities->new(
    -eutil   => "efetch",
    -db      => "nucleotide",
    -id      => [sort values $table1{'coding'}{$type}],
    -email   => $email,
    -retmode => "text",
    -rettype => "gb",
  );
  my $response = $fetcher->get_Response->content;
  open (my $seq_fh, "<", \$response)
    or die "Could not open sequences string for reading: $!";
  my $seqs = Bio::SeqIO->new(
    -fh     => $seq_fh,
    -format => "genbank",
  );
  while (my $seq = $seqs->next_seq()) {
    my $accession = $seq->display_id();
    my $symbol = List::Util::first {
      $table1{'coding'}{$type}{$_} eq $accession
    } keys $table1{'coding'}{$type};
    my $protein = (($seq->get_SeqFeatures("CDS"))[0]->get_tag_values("translation"))[0];
    if (! exists $histone{$symbol}) {
      say "coding gene $symbol with accession $accession has no protein sequence";
      $mistakes++;
    ## substr - because when dealing with histones, we never count the first
    ##residue and so we need to "unfix" the annotation
    } elsif (substr ($protein, 1) ne $histone{$symbol}) {
      say "different protein sequence for $symbol";
      $mistakes++;
    }
  }
}

## reverse check. Do all the proteins mentioned on this table, appear on the
## big table 1 with all genes and teir accession numbers? We already know this
## is not true. Table 1 says that HIST1H2AF is not  present in humans but
## Table 2 says it has the same sequence as the consensus. And HIST1H2AG
## appears as coding protein on Table 1 but is absent from Table 2 where its
## sequence should be.
sub find_on_table1 {
  my $type = shift;
  foreach my $gene (@_) {
    if (! grep {$gene eq $_} keys $table1{'coding'}{$type}) {
      say "gene $gene has protein sequence but no accession number";
      $mistakes++;
    }
  }
}
find_on_table1 ("H2A", @table2);
find_on_table1 ("H2B", @table3);

say "Found a total of $mistakes mistakes";
