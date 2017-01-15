Marzluff 2002 as reference for histone catalogue
================================================

The sequences found here were built automatically from the data in our
reference, the published "Marzluff, W.F., Gongidi, P., Woods, K.R., Jin, J.,
Maltais, l.J. (2002) The human and mouse replication-dependent histone genes.
Genomics (80) 5:487--498 doi:10.1006/geno.2002.6850"

There were some issues while doing this since a lot of the information
required to follow many references and there are some small contradictions
between. For more details, check the comments in the source of the companion
script `build_marzluff_data_2002.pl` which has a detailed explanation of
the deducing process.

This script is able to build the `data.csv` file in this directory as well as the
individual sequence files for each protein. Simply run it without any extra
arguments.

Aside building the reference data, it also performs check of the data against
the online database. Table 1 of the paper refers accession number for all the
genes. We compared those sequences against the protein sequences used in the
text for H2A (Table 2), H2B (Table 3), H3 (reference 15 and text), and H4
(reference 56). Comments on the source code have more detail. This comparison
found the following discrepancies within the paper and its data:

* different protein sequence for Hist3h2a
* different protein sequence for Hist3h2a
* different protein sequence for Hist1h2ak
* different protein sequence for Hist2h2ab
* coding gene Hist1h2bl with accession AY158927 has no protein sequence
* different protein sequence for Hist2h2bb
* different protein sequence for Hist2h2be
* different protein sequence for Hist3h2bb
* different protein sequence for Hist2h3ca2
* different protein sequence for Hist2h3ca1

Found a total of 10 mistakes

This check is interesting, but is not part of the publication. For use in the
comparison, it was used the data from the tables and text of the paper. The
reason is that such table of differences serves two purposes:

1. offer a view that online databases are not static, that things do change
with time, and it's the reference paper that should update itself. Which one
is used is irrelevant for this.
2. for anyone that used the previous publication as basis for work, a quick
list of changes. Anyone that has done this is more likely to have used the
data given in the text than the online sequences, and even if the online
sequences were used, the differences would have already been noticed. For this
purpose a comparison against the text is more helpful.

