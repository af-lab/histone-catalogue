# Manuscript builds using Ubuntu 16.04 LTS
This has been tested using Ubuntu 16.04.2 LTS as a live persistent installation on a 32 GB USB stick. It is easier to use Ubuntu 17.04 but this How To is provided because 16.04 is the official long term support (LTS) version. All commands are typed into command line after  `$` in Ubuntu Terminal app.
## 1. Ensure Ubuntu is updated
`sudo add-apt-repository universe`

`sudo add-apt-repository multiverse`

`sudo apt update`

`sudo apt upgrade`

It is not possible to dist-upgrade a live USB installation because the kernal is read-only.
### Appstream cache update error
A known critical error warning occurs with the Appstream system component. For details see:

* https://askubuntu.com/questions/854168/how-i-can-fix-appstream-cache-update-completed-but-some-metadata-was-ignored-d
* https://askubuntu.com/questions/25717/how-do-i-enable-the-backports-repository

To fix this, do the following:

`sudo nano /etc.apt/source.list` then remove # character from front of lines for xenial-backports to uncomment

`sudo apt install appstream/xenial-backports`

`sudo appstreamcli refresh --force`

`sudo apt update`

A notice also occurs after apt commands about ignoring unattended upgrades that can be ignored, or the offending files could be removed.

* https://askubuntu.com/questions/829370/n-ignoring-file-50unattended-upgrades-ucf-dist-in-directory-etc-apt-apt-con

## 2. Install dependencies
`sudo apt install scons`

`sudo apt install git`

`sudo apt install inkscape`

`sudo apt install bioperl`

`sudo apt install paml`

`sudo apt --yes install texlive texlive-latex-extra texlive-fonts-extra 
texlive-science texlive-publishers texlive-lang-greek`

`sudo apt --yes install libmoose-perl libmoosex-strictconstructor-perl libmodule-scandeps-perl libtest-output-perl libstatistics-basic-perl libtext-csv-perl`

`sudo apt install python-setuptools`

`sudo easy_install weblogo`

## 3. Install libbio-eutilities-perl
The libbio-eutilities-perl package is not in the official repository so needs to be included manually.

However, the of thsi package version available claims a bioperl dependency ahead of the version available in the Ubuntu 16.04 repository. In fact the bioperl version in the 16.04 repository is sufficient so we need to edit the bioperl dependency specification in the libbio-eutilities-perl package before installing.

'wget https://launchpad.net/ubuntu/+archive/primary/+files/libbio-eutilities-perl_1.75-2_all.deb`

`ar x libbio-eutilities-perl_1.75-2_all.deb`

`tar xzf control.tar.gz`

`nano control` then edit line 'Depends: perl, bioperl (>= 1.7.1)' to read 'Depends: perl, bioperl (>= 1.6.924)'

`tar c {post,pre}{inst,rm} md5sums control | gzip -c > control.tar.gz` and ignore the warnings

`ar rcs libbio-eutilities-perl_1.75-2_all2.deb debian-binary control.tar.gz data.tar.xz`

`sudo apt install ~/libbio-eutilities-perl_1.75-2_all2.deb`

## 4. Clone histone-catalogue
`git clone https://github.com/af-lab/histone-catalogue.git`

## 5. Build manuscript
`cd histone-catalogue`

`git pull`

`scons --email=youremail@youruniversity.edu manuscript`

Print or copy the manuscript.pdf file generated in the histone-catalogue directory.

### Notes
* To refresh all RefSeq data use `scons -c --email=youremail@youruniversity.edu`
  * Executing scons without manuscript will generate figures only as catalogue.pdf
  * Email address is required but is simply recorded in the log file. It is not transferred remotely or used for any other reason.
* First or refreshed builds take a while to download fresh data from NCBI so it is often convenient to start them in a virtual screen session using `screen -DR build`


