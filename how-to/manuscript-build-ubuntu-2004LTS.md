# Manuscript builds using Ubuntu 20.04 LTS
This has been tested using Ubuntu 20.04 LTS beta in Virtualbox 6.1 running on macOS 10.13.6 High Sierra. Ubuntu 20.04 LTS is in beta release only at time of writing so be sure to update regularly.

All commands are typed into command line after  `$` in Ubuntu Terminal app.

## 1. Ensure Ubuntu is updated
`sudo add-apt-repository universe`

`sudo add-apt-repository multiverse`

`sudo apt update`

`sudo apt upgrade`

## 2. Install dependencies
`sudo apt install scons`

`sudo apt install git`

`sudo apt install inkscape`

`sudo apt install bioperl`

`sudo apt install paml`

`sudo apt install scipy`

`sudo apt --yes install libbio-eutilities-perl 
libbio-tools-phylo-paml-perl 
libbio-tools-run-alignment-tcoffee-perl`

`sudo apt --yes install texlive 
texlive-latex-extra 
texlive-fonts-extra 
texlive-font-utils 
texlive-science 
texlive-publishers 
texlive-lang-greek`

`sudo apt --yes install libmoose-perl 
libmoosex-strictconstructor-perl 
libmodule-scandeps-perl 
libtest-output-perl 
libstatistics-basic-perl 
libtext-csv-perl 
libfile-which-perl`

`sudo apt install python3-pip`

`sudo pip3 install weblogo`

## 3. Clone histone-catalogue
`git clone https://github.com/af-lab/histone-catalogue.git`

## 4. Build manuscript
`cd histone-catalogue`

`git pull`

`scons --email=youremail@youruniversity.edu manuscript`

Print or copy the manuscript.pdf file generated in the histone-catalogue directory.

### Notes
* To refresh all RefSeq data use `scons -c --email=youremail@youruniversity.edu`
  * Executing scons without manuscript will generate figures only as catalogue.pdf
  * Email address is required but is simply recorded in the log file. It is not transferred remotely or used for any other reason.
* First or refreshed builds take a while to download fresh data from NCBI so it is often convenient to start them in a virtual screen session using `screen -DR build`


