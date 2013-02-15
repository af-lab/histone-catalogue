Human Histone Catalogue
=======================

A catalogue of all the human histone genes and proteins, with discussion
on the definition of histone variables, histone isoforms, histone clusters,
and what it means to its nomenclature.

This paper is generated automatically. The concept is that when the databases
are updated, this paper updates itself automatically. There is no actual data on
this repository (data used for the paper version is on a separate branch. See
below for details), only code that automates the downloading and analysis of
current data, the generation of figures and tables, and building of a PDF.

Directory structure
-------------------

Note that most of the directories *should* be empty before the build. Even some
of the directories that start with some files will fill up after running scons.


* figs/ - figures
* results/ - data after processing (all generated automatically so should be empty)
  * results/sequences/ - automatically downloaded sequences
* scripts/ - directory for scripts that will analyse the data
* tables/ - LaTeX tables which will be generated automatically
* histone_catalog.tex - LaTeX source for the paper
* README.md - this file, in Markdown syntax
* references.bib - BibTeX database for the references in the paper
* SConstruct - build instructions for SCons

Building instructions
---------------------

By default, SCons will download new data from the Entrez Gene database, analyse
it, and prepare PDF's of both the short report and published paper. Each of
these is a different target and SCons can be used to process only part of them.
The different target names are:

    data        - download new data (requires email)
    analysis    - analyse data
    report      - prepare PDF of short report of the analysis
    publication - prepare PDF of the published paper

Not specifying any target, calling `scons' on its own, is equivalent to
selecting all of them. That means, if all dependencies are installed, the
following command will take care of everything:

    scons --email example@domain.com

Note than in order to prevent overload of the NCBI servers, a valid email is
required to connect to the Entrez Gene database using E-utilities. Other
limitations may be in place. See "Usage guidelines and requirements", on
[A General Introduction to the E-utilities](http://www.ncbi.nlm.nih.gov/books/NBK25497/).

> In order not to overload the E-utility servers, NCBI recommends that users
> post no more than three URL requests per second and limit large jobs to either
> weekends or between 9:00 PM and 5:00 AM Eastern time during weekdays. Failure
> to comply with this policy may result in an IP address being blocked from
> accessing NCBI.
>
> ...
>
> The value of email will be used only to contact developers if NCBI observes
> requests that violate our policies, and we will attempt such contact prior to
> blocking access.

Run `scons -h` for more details on the building targets and options.

Dependencies
------------

The following dependencies need to be installed for a successful build of all
targets (see Building instructions). Individual targets will require only a
subset of these. Run `scons -h` for a complete list, divided by target.

* [SCons](www.scons.org) - the build system
* a TeX distribution with pdfTeX such as [TeX Live](http://www.tug.org/tex-live/)
* perl - at least version 5.010 and the following modules
  * Bio::EUtilities (includes the bp_genbank_ref_extractor script)
  * Bio::SeqIO
  * Bio::Tools::Run::Alignment::TCoffee
  * Email::Valid
  * Lingua::EN::Numbers
  * Text::CSV
* [weblogo](weblogo.threeplusone.com) - generate the sequence logos

The plan is to add a tag for the source version that got accepted for
publication, and a branch with the data used for it. A checkout of that branch
should only require pdfTeX and BibTeX.
