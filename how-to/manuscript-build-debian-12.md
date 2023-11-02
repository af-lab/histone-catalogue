# Building in Debian 12


## Install all build dependencies

The histone catalogue makes use of a lot of nice code written by other
people.  We've made sure that it's all available in the Debian
repositories so you only need to:

    sudo apt install \
        libbio-asn1-entrezgene-perl \
        libbio-asn1-entrezgene-perl \
        libbio-eutilities-perl \
        libbio-eutilities-perl \
        libbio-perl-perl \
        libbio-perl-perl \
        libbio-tools-phylo-paml-perl \
        libbio-tools-phylo-paml-perl \
        libbio-tools-run-alignment-clustalw-perl \
        libbio-tools-run-alignment-tcoffee-perl \
        libbio-tools-run-alignment-tcoffee-perl \
        libfile-which-perl \
        libmodule-scandeps-perl \
        libmodule-scandeps-perl \
        libmoose-perl \
        libmoosex-strictconstructor-perl \
        libnamespace-autoclean-perl \
        libstatistics-basic-perl \
        libtest-output-perl \
        libtext-csv-perl \
        python3-weblogo \
        scons \
        texlive-base \
        texlive-font-utils \
        texlive-fonts-extra \
        texlive-latex-base \
        texlive-latex-extra \
        texlive-latex-recommended \
        texlive-publishers \
        texlive-science


## Get a clone of the repository

All code is in a git repository.  First install git and then clone the
repository locally:

    sudo apt install git
    git clone https://github.com/af-lab/histone-catalogue.git
    cd histone-catalogue


## Build manuscript

    scons update manuscript

### Building manuscript with an NCBI account

If building fails during the step of downloading data, try again
during the weekend or night (relative to Eastern Time Zone).  In
addition, create and use an E-utilities API key.

To create an E-utilities API key, first register for a NCBI account at
https://www.ncbi.nlm.nih.gov/ .  You can create the API key in the
"Settings" page of your NCBI account (after signing in, click on your
NCBI username in the upper right corner of any NCBI page).  There is
an "API Key Management" area.  Click the "Create an API Key" button,
and copy the resulting key.  Then build the manuscript like so:

    scons \
        --email='example@domain.top' \
        --api_key='xxxxxxxxxxxxxxxxxxx' \
        update manuscript


## Finding the manuscript

There should appear a `manuscript.pdf` and `catalogue.pdf` document in
the `histone-catalogue` directory.
