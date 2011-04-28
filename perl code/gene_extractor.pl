#!/usr/bin/perl -w
use 5.010;                      # Use Perl 5.10
use strict;                     # Enforce some good programming rules
use Bio::DB::EUtilities;        # Retrieve entries from Entrez
use Bio::SeqIO;                 # Work with sequence files

## extra bp to retrieve from the flanking sequence
## because of the way the code checks which is the right gene to name the file,
## it won't work well with negative values for these
my $bp5_extra = 20;
my $bp3_extra = 20;
## GenBank database to query
my $queried   = "gene";
## results limitations (stop if query returns more than this number of results)
my $limit     = 100;
## e-mail address
my $email     = 'd.pinto2@nuigalway.ie';
## directory to save sequences and log (creates it if it doesn't exist)
my $save_dir  = 'indu_histones';

## Open file for log
{
  my ($second,
      $minute,
      $hour,
      $day,
      $month,
      $year
      )         = (localtime)[0,1,2,3,4,5];
  my $time      = sprintf ("[%04d-%02d-%02d %02d:%02d:%02d]", $year+1900, $month+1, $day, $hour, $minute, $second);

  if (-e $save_dir && -d _) {
    ## all is good
  } elsif (-e $save_dir) {
    say STDERR "Directory to save sequences and log already exists as non-directory";
  } else {
    mkdir $save_dir or die "Couldn't mkdir $save_dir: $!\n";
  }

  my $log_file  = "$save_dir/${time}_gene_extractor.log";
  open (LOG, ">", $log_file) or die "Couldn't open $log_file for writing: $!\n";
  ## Flush LOG immediately
  select (LOG); $|=1; select(STDOUT);
}


my @all_ids;
my @genes = ("Histone H2A", "Histone H2B", "Histone H3", "Histone H4");

foreach my $gene (@genes) {
  my $query     = '"Hydra magnipapillata"[ORGANISM] AND '."${gene}*[GENE]";

  ## make the search
  my $searcher  = Bio::DB::EUtilities -> new (-eutil   => 'esearch',
                                              -db      => $queried,
                                              -term    => $query,
                                              -tool    => 'bioperl',
                                              -email   => $email,
                                              -retmax  => $limit
                                              );
  my @ids         = $searcher->get_ids;
  my $n_results   = $searcher->get_count;
  my $query_trans = $searcher->get_query_translation;

  say LOG "Finished search for \'$query_trans\'";
  say LOG "  - Found $n_results results";

  if ($n_results == 0) {
    next;
  } else {
    say LOG "  - Ids: ".join("; ", @ids);
    my %seen;
    foreach (@ids) {
      $seen{$_}++;
    }
    @ids = keys %seen;
    say LOG "  - Found ". ($#ids+1) ." unique ids";
    say LOG "  - Ids: ".join("; ", @ids);
  }

  push (@all_ids, @ids);
}

## Extract unique elements from list
{
  say LOG "Removing duplicate entries from total search";
  my %seen;
  foreach (@all_ids) {
    $seen{$_}++
  }
  @all_ids = keys %seen;
  say LOG "  - Found ". ($#all_ids+1) ." unique ids";
}

## Get data from the ids
my $summaries = Bio::DB::EUtilities -> new (-eutil  => 'esummary',
                                            -email  => $email,
                                            -db     => $queried,
                                            -id     => \@all_ids);

my $fetcher   = Bio::DB::EUtilities -> new (-eutil    => 'efetch',
                                            -db       => 'nucleotide',
                                            -rettype  => 'gb');

## variable to hold the fetched files which later will be opened as a file
my $file;
while (my $docsum = $summaries->next_DocSum) {

  ## Show all of docsum
#  say $docsum->to_string;

  my ($name)  = $docsum->get_contents_by_name("Name");

  ## some items in DocSum are also named ChrStart so we pick the genomic
  ## information item and get the corrdinates from it
  my ($genomic_info)  = $docsum->get_Items_by_name('GenomicInfoType');

  if (!$genomic_info) {
    say LOG "$name had an empty value of GenomicInfoType";
    next;
  }

  ## get coordinates of sequence
  ## get_contents_by_name always returns a list
  my ($chr_acc_ver)   = $genomic_info->get_contents_by_name("ChrAccVer");
  my ($chr_start)     = $genomic_info->get_contents_by_name("ChrStart");
  my ($chr_stop)      = $genomic_info->get_contents_by_name("ChrStop");
  my $strand;

  say LOG "$name is found at Acession number $chr_acc_ver with start and stop coordinates $chr_start and $chr_stop";

  if ($chr_start > $chr_stop) {
    $strand     = 2;
    $chr_start  = $chr_start +1 - (-$bp5_extra);
    $chr_stop   = $chr_stop  +1 + (-$bp5_extra);
  } else {
    $strand     = 1;
    $chr_start  = $chr_start +1 - $bp5_extra;
    $chr_stop   = $chr_stop  +1 + $bp5_extra;
  }

  $fetcher -> set_parameters (-id         => $chr_acc_ver,
                              -seq_start  => $chr_start,
                              -seq_stop   => $chr_stop,
                              -strand     => $strand);

  ## EUtilities cannot retrieve a seq object, only the content of the file.
  ## Concatenate all sequences in a variable/file for opening later.
  ## Not very memory efficient for large sequences or large number of sequences,
  ## so it might be a good idea to keep track of variable size and once it
  ## reachs a certain limit, take care of them, clear the variable, and then
  ## continue.
  $file  .= $fetcher->get_Response->content;
}

open (SEQUENCES, "<", \$file) or die "Couldn't open sequences variable for reading: $!\n";
my $seqio   = Bio::SeqIO->new(-fh     => \*SEQUENCES,
                              -format => 'genbank');

while (my $seq = $seqio->next_seq) {
  
  my @features = $seq->get_SeqFeatures;
  my @gene_ids;
#  my @mRNA_products;
  foreach my $feature (@features) {
    ## Get the value of 'gene' inside the 'gene' feature
    if      ($feature->primary_tag eq 'gene' && $feature->start == $bp5_extra+1 && $feature->end == $seq->length-$bp3_extra) {
      push (@gene_ids, $feature->get_tag_values('gene'));
    ## Get the value of 'product' inside the 'mRNA' feature
#    }elsif  ($feature->primary_tag eq 'mRNA') {
#      push (@mRNA_products, $feature->get_tag_values('product'));
    }
  }
#  my $filepath = "$save_dir/GENE_IDS: @gene_ids mRNA's_PRODUCTS: @mRNA_products.gb";
  my $filepath = "$save_dir/@gene_ids.gb";

  my $out = Bio::SeqIO->new(-file   => ">$filepath",
                            -format => "genbank");

  $out->write_seq($seq);
  $out->DESTROY;
}

## Destroying the object also closes the filehandle
$seqio->DESTROY;
close LOG;
exit;
