#!/usr/bin/perl -w
use 5.010;                      # Use Perl 5.10
use strict;                     # Enforce some good programming rules
use Bio::SeqIO;                 # Handler for SeqIO Formats
use File::Spec;                 # Portably perform operations on file names

my $out_format  = "fasta";    # format to save translated sequences
my $trans_skip  = 3;          # number of bp to skip at start of CDS before start translating

my $dirname     = $ARGV[0];

my $all_out = Bio::SeqIO->new(-file   =>">".File::Spec->catfile($dirname, 'all translations.faa'),
                              -format => $out_format);

opendir(DIR, $dirname) or die "Can't opendir $dirname: $!\n";
while (defined(my $filename = readdir(DIR))) {
  next if $filename =~ /^\.\.?$/;          # skip special files . and ..

  my $in    = Bio::SeqIO->new(-file   => File::Spec->catfile($dirname, $filename));

  while (my $seq = $in->next_seq) {
    my @features = $seq->get_SeqFeatures;
    ## loop trough features until finding CDS
    foreach my $feature (@features) {
      if ($feature->primary_tag eq 'CDS') {
        ## skip wrong CDS 
        my ($genename) = $feature->get_tag_values('gene');
        my $potential_name = "$genename.gb";
        if($potential_name ne $filename) {
          say "skipping $feature in $filename";
          next;
        }

        ## replace sequence of the sequence object with the subsequence of the CDS
        $seq->seq($seq->subseq( ($feature->start + $trans_skip), $feature->end)); 

        ## get protein sequence object
        my $protein = $seq->translate;
        ## replace display_ids and descriptions for something more suitable
        $protein->display_id($feature->get_tag_values('protein_id'));
        $protein->desc      ($feature->get_tag_values('product'));

        ## save protein objects on individual and global file
        my $out = Bio::SeqIO->new(-file   =>">".File::Spec->catfile($dirname, "protein_$filename"),
                                  -format => $out_format);
        $out    ->write_seq($protein);
        $out    ->DESTROY;
        $all_out->write_seq($protein);

      }
    }
  }
}
$all_out->DESTROY;
closedir(DIR);


