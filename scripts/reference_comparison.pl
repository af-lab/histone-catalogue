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
use File::Spec;                 # Perform operation on file names
use Getopt::Long;               # Parse program arguments
use Bio::SeqIO;                 # Handler for SeqIO formats

use FindBin;                    # Locate directory of original perl script
use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyLib;                      # Load functions

## This script compares the new data against some reference. In our case, the reference
## is Marzluff, W.F., Gongidi, P., Woods, K.R., Jin, J., Maltais, l.J. (2002) The human
## and mouse replication-dependent histone genes. Genomics (80) 5:487--498
## doi:10.1006/geno.2002.6850

## Check input options
my %path = ("sequences" => "",
            "results"   => "");
my $email;
GetOptions(
            "sequences=s" => \$path{sequences},
            "results=s"   => \$path{results},
            'email=s'     => \$email,
          ) or die "Error processing options.";
for (keys %path) {
  die "No path for $_ specified. Use the --$_ option." unless $path{$_};
}
## we don't check the email, that's optional later on

my %reference = (
  consensus => {
    ## H2A consensus sequence from Table 2: Histone H2A protein
    ## variants (page 494)
    H2A => "SGRGKQGGKARAKAKTRSSRAGLQFPVGRVHRLLRKGNYAERVGAGAPVYLAAVLEYLTAEILELA".
           "GNAARDNKKTRIIPRHLQLAIRNDEELNKLLGKVTIAQGGVLPNIQAVLLPKKTESHHKAKGK",
    ## H2B consensus sequence from Table 3: Histone H2B protein
    ## variants (page 495)
    H2B => "PEPAKSAPAPKKGSKKAVTKAQKKDGKKRKRSRKESYSVYVYKVLKQVHPDTGISSKAMGIMNSFV".
           "NDIFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSSK",
    ## no consensus sequence for the other histones
    },

  ## from Table 1 (page 488 and 489)
  coding => [qw(HIST1H2AA HIST1H2AB HIST1H2AC HIST1H2AD HIST1H2AE HIST1H2AG HIST1H2AH
                HIST1H2AI HIST1H2AJ HIST1H2AK HIST1H2AL HIST1H2AM HIST2H2AA HIST2H2AB
                HIST2H2AC HIST3H2A
                
                HIST1H2BA HIST1H2BB HIST1H2BC HIST1H2BD HIST1H2BE HIST1H2BF HIST1H2BG
                HIST1H2BH HIST1H2BI HIST1H2BJ HIST1H2BK HIST1H2BL HIST1H2BM HIST1H2BN
                HIST1H2BO HIST2H2BE HIST3H2BB
                
                HIST1H3A  HIST1H3B  HIST1H3C  HIST1H3D  HIST1H3E  HIST1H3F  HIST1H3G
                HIST1H3H  HIST1H3I  HIST1H3J  HIST2H3C  HIST3H3
                
                HIST1H4A  HIST1H4B  HIST1H4C  HIST1H4D  HIST1H4E  HIST1H4F  HIST1H4G
                HIST1H4H  HIST1H4I  HIST1H4J  HIST1H4K  HIST1H4L  HIST2H4   HIST4H4)],

  ## from Table 1 (page 488 and 489)
  pseudo => [qw(HIST2H2BA HIST2H2BB HIST2H2BC HIST2H2BD HIST3H2BA HIST2H3A  HIST1H3B)],
);

## all H2A sequences from Table 2: Histone H2A protein variants (page 494)
foreach (@{$reference{coding}}) {
  next unless $_ =~ m/^HIST\dH2A/;
  $reference{$_} = $reference{consensus}{H2A};
}

## these are swapped but which one was correct at the time? Table 1 mentions
## only --G, table 2, only mentions --F
delete $reference{HIST1H2AG};
$reference{HIST1H2AF} = $reference{consensus}{H2A};

for ("HIST1H2AA") {
  substr ($reference{$_},  13, 1) = "S";
  substr ($reference{$_},  15, 1) = "S";
  substr ($reference{$_},  29, 1) = "I";
  substr ($reference{$_},  42, 1) = "I";
  substr ($reference{$_},  70, 1) = "S";
  substr ($reference{$_},  98, 1) = "G";
  substr ($reference{$_}, 124, 6) = "HKAQSK";
}
substr ($reference{HIST1H2AB},  98, 1) = "R";
substr ($reference{HIST1H2AC},  15, 1) = "S";
substr ($reference{HIST1H2AC},  98, 1) = "R";
substr ($reference{HIST1H2AE},  98, 1) = "R";
substr ($reference{HIST1H2AH}, 123, 4) = "HKAK";
chop   ($reference{HIST1H2AH});
substr ($reference{HIST1H2AJ}, 123, 4) = "HKTK";
chop   ($reference{HIST1H2AJ});
substr ($reference{HIST2H2AA},  15, 1) = "S";
substr ($reference{HIST2H2AA},  50, 1) = "M";
for ("HIST2H2AB") {
  substr ($reference{$_},  15, 1) = "S";
  substr ($reference{$_},  50, 1) = "M";
  substr ($reference{$_},  86, 1) = "V";
  substr ($reference{$_},  98, 1) = "G";
  substr ($reference{$_}, 123, 6) = "KPGPNK";
}
substr ($reference{HIST2H2AC},  15, 1) = "S";
substr ($reference{HIST2H2AC},  50, 1) = "M";
substr ($reference{HIST2H2AC}, 127, 1) = "S";
substr ($reference{HIST3H2A},   15, 1) = "S";
substr ($reference{HIST3H2A},   50, 1) = "M";
substr ($reference{HIST3H2A},   98, 1) = "R";

## all H2B sequences from Table 3: Histone H2B protein variants (page 495)
foreach (@{$reference{coding}}) {
  next unless $_ =~ m/^HIST\dH2B/;
  $reference{$_} = $reference{consensus}{H2B};
}

for ("HIST1H2BA") {
  substr ($reference{$_},   5, 1) = "G";
  substr ($reference{$_},   7, 3) = "TIS";
  substr ($reference{$_},  13, 1) = "F";
  substr ($reference{$_},  18, 1) = "V";
  substr ($reference{$_},  20, 1) = "T";
  substr ($reference{$_},  24, 1) = "E";
  substr ($reference{$_},  31, 1) = "T";
  substr ($reference{$_},  38, 1) = "I";
  substr ($reference{$_},  40, 1) = "I";
  substr ($reference{$_},  59, 1) = "S";
  substr ($reference{$_},  66, 1) = "T";
  substr ($reference{$_},  74, 1) = "S";
  substr ($reference{$_},  83, 1) = "S";
  substr ($reference{$_},  89, 1) = "S";
  substr ($reference{$_}, 124, 1) = "N";
  substr ($reference{$_},   2, 2) = "VSS"; # at the end because of the addition
}
substr ($reference{HIST1H2BB},  17, 1) = "I";
substr ($reference{HIST1H2BB},  38, 1) = "I";
substr ($reference{HIST1H2BD},   3, 1) = "T";
substr ($reference{HIST1H2BH},   1, 1) = "D";
substr ($reference{HIST1H2BJ}, 123, 1) = "A";
substr ($reference{HIST1H2BK}, 123, 1) = "A";
substr ($reference{HIST1H2BL},   2, 1) = "L";
substr ($reference{HIST1H2BL},  74, 1) = "S";
substr ($reference{HIST1H2BM},   3, 1) = "V";
substr ($reference{HIST1H2BM},   8, 1) = "V";
substr ($reference{HIST1H2BM},  17, 1) = "I";
substr ($reference{HIST1H2BM},  18, 1) = "N";
substr ($reference{HIST1H2BN},   3, 1) = "S";
substr ($reference{HIST1H2BO},   1, 1) = "D";
substr ($reference{HIST1H2BO},  38, 1) = "I";
substr ($reference{HIST2H2BE},  38, 1) = "I";
for ("HIST3H2BB") {
  substr ($reference{$_},   1, 1) = "D";
  substr ($reference{$_},   3, 1) = "S";
  substr ($reference{$_},  31, 1) = "G";
  substr ($reference{$_},  38, 1) = "I";
  substr ($reference{$_},  74, 1) = "S";
  substr ($reference{$_},  93, 1) = "V";
}


## There was no H3 protein sequence on the paper. Instead, it's written on page 489:
##
##    ... the histone H3 genes encode the two variants previously described [15] ...
##
## The reference target is "Zweidler, A. (1984). Core histone variants of the mouse: primary
## structure and expression. In Histone Genes: Structure, Organization and Regulation (G. Stein,
## W.Stein W. F. Marzluff, Eds), pp.373-395. John Wiley and Sons, New York". The correct page
## numbers are actually 339-371. And on page 494:
##
##    There are two variants encoded by the histone H3 genes, H3.1 and H3.2. These
##    differ in a single amino acid change at position 96, a cysteine in H3.1 genes
##    and a serine in H3.2 genes. The functional histone H3 gene (HIST2H3C) in the
##    HIST2 cluster on chromosome 1 encodes the H3.2 protein and all 11 histone H3
##    genes in the HIST1 cluster encode the H3.1 protein.
##
## However, on Table 1, there's only mention of 10 histone H3 genes for the HIST1 cluster. There
## is also mention of a third H3. Previously on page 492:
##
##    HIST3H3, described previously as a solitary gene [19], encodes a novel histone
##    H3 protein, with the serine at position 96 and four additional changes (A24V,
##    V71M, A98S, A111V)
##
## This leads to a total of 3 unique H3 protein sequences
foreach (@{$reference{coding}}) {
  next unless $_ =~ m/^HIST\dH3/;
  ## sequence from Zweidler, A (1984) chapter on the histone genes book
  $reference{$_} = "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKST".
                   "ELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEACEAYLVGLFEDTNLCAIHAKR".
                   "VTIMPKDIQLARRIRGERA";
}
foreach (@{$reference{coding}}) {
  next unless $_ =~ m/^HIST2H3/;
  substr ($reference{$_}, 95, 1) = "S";
}
substr ($reference{HIST3H3},  23, 1) = "V";
substr ($reference{HIST3H3},  70, 1) = "M";
substr ($reference{HIST3H3},  95, 1) = "S";
substr ($reference{HIST3H3},  97, 1) = "S";
substr ($reference{HIST3H3}, 110, 1) = "V";

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
foreach (@{$reference{coding}}) {
  next unless $_ =~ m/^HIST\dH4/;
  ## sequence from Akasaka, T et al (1997)
  $reference{$_} = "SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVL".
                   "KVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG";
}

## start comparing


my %changes = (
  current  => [], # list of all current genes
  removed  => [], # were removed
  added    => [], # were added
  pseudo   => [], # changed to pseudo genes
  coding   => [], # changed to coding genes
  sequence => [], # sequence has changed
);
my @data = MyLib::load_canonical ($path{sequences});
foreach my $gene (@data) {
  my $symbol = $$gene{'gene symbol'};
  push (@{$changes{current}}, $symbol);
  if (! grep ($_ eq $symbol, @{$reference{coding}}) &&
      ! grep ($_ eq $symbol, @{$reference{pseudo}})) {
    ## it's a new gene
    push (@{$changes{added}}, $symbol);
  } elsif ($$gene{pseudo} && grep ($_ eq $symbol, @{$reference{coding}})) {
    ## used to be coding gene, now is a pseudo gene
    push (@{$changes{pseudo}}, $symbol);
  } elsif (! $$gene{pseudo} && grep ($_ eq $symbol, @{$reference{pseudo}})) {
    ## used to be pseudo gene, now is a coding gene
    push (@{$changes{coding}}, $symbol);
  } elsif (exists $$gene{'protein accession'} && exists $reference{$symbol}) {
    ## compare the sequence
    my $seq = MyLib::load_protein($path{sequences}, $$gene{'protein accession'});
    if ($seq->seq ne $reference{$symbol}) {
      ## sequence has changed
      push (@{$changes{sequence}}, $symbol);
    }
  }
}

foreach my $symbol (@{$reference{pseudo}}, @{$reference{coding}}) {
  if (! grep ($_ eq $symbol, @{$changes{current}})) {
    ## gene no longer exists
    push (@{$changes{removed}}, $symbol);
  }
}

## Write down results
my $var_path = File::Spec->catdir($path{results}, "table-reference_comparison.tex");
open (my $var_fh, ">", $var_path) or die "Could not open $var_path for writing: $!";

say {$var_fh} "\\begin{tabular}{p{\\dimexpr\\textwidth-2\\tabcolsep\\relax}}";
say {$var_fh} "  \\toprule";

my %pairs = (
             "Changed sequences",              => \@{$changes{sequence}},
             "New genes",                      => \@{$changes{added}},
             "Removed genes",                  => \@{$changes{removed}},
             "Now identified as pseudo genes", => \@{$changes{pseudo}},
             "Now identified as coding genes", => \@{$changes{coding}},
             );

my $space = 1; # have we left a space yet?
foreach (keys %pairs) {
  if (scalar (@{$pairs{$_}})) {
    say {$var_fh} "  \\addlinespace" unless $space;
    say {$var_fh} "  $_: \\\\";
    say {$var_fh} join (", ", @{$pairs{$_}}) . "\\\\";
    $space = 0;
  }
}

say {$var_fh} "  \\bottomrule";
say {$var_fh} "\\end{tabular}";

close ($var_fh) or die "Couldn't close $var_path after writing: $!";


## This check is an interesting check, not part of the publication. It checks for the
## differences between the sequences talked about the paper, and the actual sequences
## from the referenced accession numbers. That's how we found the mistake between
## HIST1H2AF and HIST1H2AG
if ($email) {
  my %ref_acc = (
    ## accession numbers from table 1
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

    HIST1H2BA => "AF531284",
    HIST1H2BC => "AF531285",
    HIST1H2BD => "AF531286",
    HIST1H2BE => "AF531287",
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
  );
  require Bio::DB::EUtilities;
  my $mistakes = 0;
  foreach my $symbol (keys %ref_acc) {
    my $fetcher = Bio::DB::EUtilities->new(
                                           -eutil   => "efetch",
                                           -db      => "nucleotide",
                                           -id      => [$ref_acc{$symbol}],
                                           -retmode => "text",
                                           -email   => $email,
                                           -rettype => "gb",
                                           );
    my $response = $fetcher->get_Response->content;
    open(my $seq_fh, "<", \$response) or die "Could not open sequences string for reading: $!";
    my $parser = Bio::SeqIO->new(
                                 -fh     => $seq_fh,
                                 -format => "genbank",
                                 );
    my $seq = $parser->next_seq; # we know there's only 1 sequence
    my ($feature) = $seq->get_SeqFeatures("CDS");
    my ($protein) = $feature->get_tag_values("translation");

    ## substr - because when dealing with histones, we never count the first
    ##residue and so we need to "unfix" the annotation
    if (! exists $reference{$symbol}) {
      say "Found accession for $symbol but not a sequence for reference";
      $mistakes++;
    } elsif (substr ($protein, 1) ne $reference{$symbol}) {
      say "Found different sequences for $symbol";
      $mistakes++;
    }
  }
  say "Total number of mistakes found is: $mistakes";
}

