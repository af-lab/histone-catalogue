Human Histone Catalogue
=======================

A catalogue of all the human histone genes and proteins, with discussion
on the definition of histone variables, histone isoforms, histone clusters,
and what it means to its nomenclature.

This paper is generated automatically. The concept is that when the databases
are updated, this paper updates itself automatically. There is no actual data on
this repository, only the reference data which is used for comparison. All the
code required to build a recent version of the paper is available as well as
build instructions which have been set with SCons. There are targets for
automatic download and analysis of current data, generation of figures and
tables, and building of a ready-to-publish paper.

Patches are welcome
-------------------

Please send us improvements to our text and code. We see this as a live
publication, both data, analysis, and text should change with time. And if
changes belong to upstream projects, so be it. Please give it back to the
community.


Directory structure
-------------------

Note that most of the directories *should* be empty before the build. Even some
of the directories that start with some files will fill up after running scons.


* data/ - data that is not automatically generated (data from the Marzluff 2002
paper which we use as reference for comparison)
* figs/ - figures
* results/ - data after processing including LaTeX tables
  * results/sequences/ - automatically downloaded sequences
* scripts/ - directory for scripts that will analyse the data
* sections/ - directory whith the source for the paper
* histone_catalog.tex - LaTeX source for the paper
* README.md - this file, in Markdown syntax
* references.bib - BibTeX database for the references in the paper
* SConstruct - build instructions for SCons

Building instructions
---------------------

By default, SCons will download new data from the Entrez Gene database, analyse
it, and prepare a ready to publish PDF of the paper. Each of
these is a different target and SCons can be used to process only part of them.
The different target names are:

    data        - download new data (requires email)
    analysis    - analyse data
    publication - prepare PDF of the published paper

Not specifying any target, calling `scons` on its own, is equivalent to
selecting all of them. That means, if all dependencies are installed, the
following command will take care of everything:

    scons --email=example@domain.com

Note than in order to prevent overload of the NCBI servers, a valid email is
required to connect to the Entrez Gene database using E-utilities. Other
limitations may be in place. See section *"Usage guidelines and requirements"*,
on [A General Introduction to the E-utilities](http://www.ncbi.nlm.nih.gov/books/NBK25497/).

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
targets (see *Building instructions*). Individual targets will require only a
subset of these. Run `scons -h` for a complete list, divided by target.

* [SCons](www.scons.org) - the build system
* a TeX distribution with pdfTeX such as [TeX Live](http://www.tug.org/tex-live/)
* perl - at least version 5.010 and the following modules
  * Bio::EUtilities (includes the bp_genbank_ref_extractor script)
  * Bio::SeqIO
  * Bio::Tools::Run::Alignment::TCoffee
  * Email::Valid
  * Text::CSV
  * Statistics::Basic
* [weblogo](weblogo.threeplusone.com) - generate the sequence logos

The plan is to add a tag for the source version that got accepted for
publication, and a branch with the data used for it. A checkout of that branch
should only require pdfTeX and BibTeX.
