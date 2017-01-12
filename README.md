Histone Catalogue
=================

This project provides a catalogue of canonical core histone genes,
encoded proteins, and pseudogenes based on reference genome
annotations.  It also provides a discussion on the definition of
histone variables, isoforms, clusters, and their nomenclature.

Since curation and annotation are dynamic and evolving, the catalogue
was made a live publication so it can provide always up to date
information in an accessible format.  Inspired by the ideals of
reproducible research, this project contains all the code required to
create a new build of the catalogue from current annotations as well
as the code to automate such build.  SCons, a software build system,
is used to automate the build.


Build instructions
==================

To build the catalogue, run `scons' from the root of the project.
This will check if all required sofware is installed, search the
databases for the histone genes, download all required sequences,
analyze the sequences, generate figures, and finally compile a
catalogue in pdf format.

See `scons -h' for a complete list of targets and options.

Data
----

There is no data on this repository, all of it is download from the
Entrez databases as part of the build.  To download only the sequence
data and skip all the analysis, use the `data' target:

    scons data

Other organisms
---------------

By default, a catalogue of the human histones is generated.  Other
organisms can be specified via the `--organism' option.  This is
heavily dependent on the annotation state of the organism reference
genome and it has only been tested in human and mice.

    scons --organism='mus musculus'

Manuscript
----------

By default, only the catalogue --- a pdf with multiple tables and
figures --- is built.  There is also a manuscript which builds
something akin to a publication of the catalogue which includes a
discussion of the histone genes as well as some analysis not present
on the catalogue.

Update
------

If there is previously downloaded sequence data, a new build will
not automatically download new data.  Use the `update' target for
that.  Note that this will only update the data, if you want to
rebuild

Email
-----

The Entrez databases are searched via E-utilities requires, but does
not enforce, an email address.  This email can then be used by NCBI
staff to contact you in case you accidentally overload their servers.

> In order not to overload the E-utility servers, NCBI recommends that
> users post no more than three URL requests per second and limit
> large jobs to either weekends or between 9:00 PM and 5:00 AM Eastern
> time during weekdays. Failure to comply with this policy may result
> in an IP address being blocked from accessing NCBI.
> [...]
> The value of email will be used only to contact developers if NCBI
> observes requests that violate our policies, and we will attempt
> such contact prior to blocking access.

For more details, see section *"Usage guidelines and requirements"*,
on [A General Introduction to the E-utilities](http://www.ncbi.nlm.nih.gov/books/NBK25497/).

To set an email, use the `--email' option like so:

    scons --email=example@domain.top

Examples
--------

* Build the catalogue for human histones:

        scons

* Build the manuscript for human histones:

        scons manuscript

* Only download the sequences data for humans:

        scons data

* Only download the sequences data for mice:

        scons --organism='mus musculus' data

* Update previously downloaded data:

        scons update

* Rebuild catalogue with new data

        scons update catalogue


Dependencies
============

Several pieces of software are required to build the histone
catalogue:

* [SCons](www.scons.org) which provides the build system.
* pdflatex, bibtex, epstopdf, and several other latex packages are
  required to build the catalogue and manuscript pdf files.  Simplest
  method is to install [TeX Live](http://www.tug.org/tex-live/) which
  provides all of them in a single distribution.
* [Perl](https://www.perl.org/) as well as several perl modules.
* `bp_genbank_ref_extractor' which is used for search and download of
  sequences is part of bioperl's
  [Bio-EUtilities](https://metacpan.org/release/Bio-EUtilities)
  distribution.
* [weblogo](weblogo.threeplusone.com) to generate the sequence logos.
* [inkscape](https://inkscape.org/) for figure svg to pdf conversion.

A complete list of required perl modules and latex packages is listed
via `scons -h'.


Directory structure
===================

* data - data that is not automatically generated such as the data
  from Marzluff 2002 paper which we use as reference for comparison.
* figs - figures generated during the build.
* lib-perl5 - library for handling sequences and needed by our perl
  scripts.
* results - data after processing.  Includes aligned sequences, as
  well as LaTeX tables and variable definitions.
* results/sequences - sequences downloaded as part of the build.
* scripts - collection of scripts for data analysis.
* sections - LaTeX source for the different manuscript sections.
* site_scons - scons configuration for this project.
* t - tests for lib-perl5.
