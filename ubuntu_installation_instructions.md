# Andrew's Ubuntu Installation Instructions for SEQUIN

Run on 10/21/2024.

## Install R

Get OS information:

```console
(base) weismanal@NCATS-2260152-P:~/notebook/2024-10-18/testing_sequin_installation_on_linux$ uname -a
Linux NCATS-2260152-P 5.10.16.3-microsoft-standard-WSL2 #1 SMP Fri Apr 2 22:23:49 UTC 2021 x86_64 x86_64 x86_64 GNU/Linux
(base) weismanal@NCATS-2260152-P:~/notebook/2024-10-18/testing_sequin_installation_on_linux$ cat /etc/os-release
PRETTY_NAME="Ubuntu 22.04.5 LTS"
NAME="Ubuntu"
VERSION_ID="22.04"
VERSION="22.04.5 LTS (Jammy Jellyfish)"
VERSION_CODENAME=jammy
ID=ubuntu
ID_LIKE=debian
HOME_URL="https://www.ubuntu.com/"
SUPPORT_URL="https://help.ubuntu.com/"
BUG_REPORT_URL="https://bugs.launchpad.net/ubuntu/"
PRIVACY_POLICY_URL="https://www.ubuntu.com/legal/terms-and-policies/privacy-policy"
UBUNTU_CODENAME=jammy
```

Create new conda environment:

```console
(base) weismanal@NCATS-2260152-P:~/notebook/2024-10-18/testing_sequin_installation_on_linux$ conda create -n sequin
Channels:
 - conda-forge
 - defaults
Platform: linux-64
Collecting package metadata (repodata.json): done
Solving environment: done

## Package Plan ##

  environment location: /home/weismanal/programs/miniconda3/envs/sequin



Proceed ([y]/n)?

Preparing transaction: done
Verifying transaction: done
Executing transaction: done
#
# To activate this environment, use
#
#     $ conda activate sequin
#
# To deactivate an active environment, use
#
#     $ conda deactivate
```

Set up latest R packages database from [here](https://github.com/eddelbuettel/r2u?tab=readme-ov-file) based on instructions [here](https://cran.r-project.org/bin/linux/ubuntu/fullREADME.html):

```console
(sequin) weismanal@NCATS-2260152-P:~/notebook/2024-10-18/testing_sequin_installation_on_linux$ sudo bash add_cranapt_jammy.sh
8 packages can be upgraded. Run 'apt list --upgradable' to see them.
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
ca-certificates is already the newest version (20240203~22.04.1).
gnupg is already the newest version (2.2.27-3ubuntu2.1).
wget is already the newest version (1.21.2-2ubuntu1.1).
0 upgraded, 0 newly installed, 0 to remove and 8 not upgraded.
-----BEGIN PGP PUBLIC KEY BLOCK-----

<andrew snip>

-----END PGP PUBLIC KEY BLOCK-----
Hit:1 http://archive.ubuntu.com/ubuntu jammy InRelease
Hit:2 http://security.ubuntu.com/ubuntu jammy-security InRelease
Hit:3 http://archive.ubuntu.com/ubuntu jammy-updates InRelease
Hit:4 http://archive.ubuntu.com/ubuntu jammy-backports InRelease
Ign:5 https://r2u.stat.illinois.edu/ubuntu jammy InRelease
Get:6 https://r2u.stat.illinois.edu/ubuntu jammy Release [5713 B]
Get:7 https://r2u.stat.illinois.edu/ubuntu jammy Release.gpg [793 B]
Get:8 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 Packages [1876 kB]
Get:9 https://r2u.stat.illinois.edu/ubuntu jammy/main all Packages [6144 kB]
Fetched 8026 kB in 2s (3673 kB/s)
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
8 packages can be upgraded. Run 'apt list --upgradable' to see them.
-----BEGIN PGP PUBLIC KEY BLOCK-----

<andrew snip>

-----END PGP PUBLIC KEY BLOCK-----
Warning: apt-key is deprecated. Manage keyring files in trusted.gpg.d instead (see apt-key(8)).
Executing: /tmp/apt-key-gpghome.1h6z6u3dUR/gpg.1.sh --keyserver keyserver.ubuntu.com --recv-keys 67C2D66C4B1D4339 51716619E084DAB9
gpg: key 51716619E084DAB9: "Michael Rutter <marutter@gmail.com>" 1 new signature
gpg: key 67C2D66C4B1D4339: public key "Launchpad PPA for Dirk Eddelbuettel" imported
gpg: Total number processed: 2
gpg:               imported: 1
gpg:         new signatures: 1
8 packages can be upgraded. Run 'apt list --upgradable' to see them.
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
The following additional packages will be installed:
  libblas3 libgfortran5 liblapack3 libtcl8.6 libtk8.6 libxss1
Suggested packages:
  tcl8.6 tk8.6 elpa-ess r-doc-info | r-doc-pdf r-mathlib r-base-html
Recommended packages:
  r-recommended r-base-dev r-doc-html
The following NEW packages will be installed:
  libblas3 libgfortran5 liblapack3 libtcl8.6 libtk8.6 libxss1 r-base-core
0 upgraded, 7 newly installed, 0 to remove and 8 not upgraded.
Need to get 34.2 MB of archives.
After this operation, 65.4 MB of additional disk space will be used.
Get:1 http://archive.ubuntu.com/ubuntu jammy/main amd64 libblas3 amd64 3.10.0-2ubuntu1 [228 kB]
Get:2 http://archive.ubuntu.com/ubuntu jammy-updates/main amd64 libgfortran5 amd64 12.3.0-1ubuntu1~22.04 [879 kB]
Get:3 http://archive.ubuntu.com/ubuntu jammy/main amd64 liblapack3 amd64 3.10.0-2ubuntu1 [2504 kB]
Get:4 http://archive.ubuntu.com/ubuntu jammy/main amd64 libtcl8.6 amd64 8.6.12+dfsg-1build1 [990 kB]
Get:5 http://archive.ubuntu.com/ubuntu jammy/main amd64 libxss1 amd64 1:1.2.3-1build2 [8476 B]
Get:6 http://archive.ubuntu.com/ubuntu jammy/main amd64 libtk8.6 amd64 8.6.12-1build1 [784 kB]
Get:7 https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/ r-base-core 4.4.1-3.2204.0 [28.8 MB]
Fetched 34.2 MB in 1s (43.1 MB/s)
Selecting previously unselected package libblas3:amd64.
(Reading database ... 56677 files and directories currently installed.)
Preparing to unpack .../0-libblas3_3.10.0-2ubuntu1_amd64.deb ...
Unpacking libblas3:amd64 (3.10.0-2ubuntu1) ...
Selecting previously unselected package libgfortran5:amd64.
Preparing to unpack .../1-libgfortran5_12.3.0-1ubuntu1~22.04_amd64.deb ...
Unpacking libgfortran5:amd64 (12.3.0-1ubuntu1~22.04) ...
Selecting previously unselected package liblapack3:amd64.
Preparing to unpack .../2-liblapack3_3.10.0-2ubuntu1_amd64.deb ...
Unpacking liblapack3:amd64 (3.10.0-2ubuntu1) ...
Selecting previously unselected package libtcl8.6:amd64.
Preparing to unpack .../3-libtcl8.6_8.6.12+dfsg-1build1_amd64.deb ...
Unpacking libtcl8.6:amd64 (8.6.12+dfsg-1build1) ...
Selecting previously unselected package libxss1:amd64.
Preparing to unpack .../4-libxss1_1%3a1.2.3-1build2_amd64.deb ...
Unpacking libxss1:amd64 (1:1.2.3-1build2) ...
Selecting previously unselected package libtk8.6:amd64.
Preparing to unpack .../5-libtk8.6_8.6.12-1build1_amd64.deb ...
Unpacking libtk8.6:amd64 (8.6.12-1build1) ...
Selecting previously unselected package r-base-core.
Preparing to unpack .../6-r-base-core_4.4.1-3.2204.0_amd64.deb ...
Unpacking r-base-core (4.4.1-3.2204.0) ...
Setting up libblas3:amd64 (3.10.0-2ubuntu1) ...
update-alternatives: using /usr/lib/x86_64-linux-gnu/blas/libblas.so.3 to provide /usr/lib/x86_64-linux-gnu/libblas.so.3 (libblas.so.3-x86_64-linux-gnu) in auto mode
Setting up libtcl8.6:amd64 (8.6.12+dfsg-1build1) ...
Setting up libgfortran5:amd64 (12.3.0-1ubuntu1~22.04) ...
Setting up libxss1:amd64 (1:1.2.3-1build2) ...
Setting up liblapack3:amd64 (3.10.0-2ubuntu1) ...
update-alternatives: using /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3 to provide /usr/lib/x86_64-linux-gnu/liblapack.so.3 (liblapack.so.3-x86_64-linux-gnu) in auto mode
Setting up libtk8.6:amd64 (8.6.12-1build1) ...
Setting up r-base-core (4.4.1-3.2204.0) ...
Installing new version of config file /etc/R/Makeconf ...
Installing new version of config file /etc/R/Renviron.site ...
Installing new version of config file /etc/R/repositories ...

Creating config file /etc/R/Renviron with new version
Processing triggers for hicolor-icon-theme (0.17-2) ...
Processing triggers for libc-bin (2.35-0ubuntu3.8) ...
Processing triggers for man-db (2.10.2-1) ...
Processing triggers for install-info (6.8-4build1) ...
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
make is already the newest version (4.3-4.1build1).
make set to manually installed.
python3-dbus is already the newest version (1.2.18-3build1).
python3-dbus set to manually installed.
python3-apt is already the newest version (2.4.0ubuntu4).
python3-apt set to manually installed.
python3-gi is already the newest version (3.42.1-0ubuntu1).
python3-gi set to manually installed.
0 upgraded, 0 newly installed, 0 to remove and 8 not upgraded.
Installing package into ‘/usr/local/lib/R/site-library’
(as ‘lib’ is unspecified)
trying URL 'https://cloud.r-project.org/src/contrib/bspm_0.5.7.tar.gz'
Content type 'application/x-gzip' length 27097 bytes (26 KB)
==================================================
downloaded 26 KB

* installing *source* package ‘bspm’ ...
** package ‘bspm’ successfully unpacked and MD5 sums checked
** using staged installation
* installing /usr/share/dbus-1/system-services/org.r_project.linux1.service
* installing /etc/dbus-1/system.d/org.r_project.linux1.conf
** R
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (bspm)

The downloaded source packages are in
        ‘/tmp/Rtmp0SxCej/downloaded_packages’
```

Install suggested and recommended packages from above output:

```console
(sequin) weismanal@NCATS-2260152-P:~/notebook/2024-10-18/testing_sequin_installation_on_linux$ sudo apt install r-recommended r-base-dev r-doc-html tcl8.6 tk8.6 elpa-ess r-doc-info r-doc-pdf r-mathlib r-base-html
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
The following additional packages will be installed:
  gfortran gfortran-11 icu-devtools libblas-dev libbz2-dev libdeflate-dev libgfortran-11-dev libicu-dev libjpeg-dev libjpeg-turbo8-dev libjpeg8-dev liblapack-dev liblzma-dev libncurses-dev libpcre2-16-0 libpcre2-32-0 libpcre2-dev libpcre2-posix3 libpkgconf3 libpng-dev
  libpng-tools libreadline-dev pkgconf r-cran-boot r-cran-class r-cran-cluster r-cran-codetools r-cran-foreign r-cran-kernsmooth r-cran-lattice r-cran-mass r-cran-matrix r-cran-mgcv r-cran-nlme r-cran-nnet r-cran-rpart r-cran-spatial r-cran-survival
Suggested packages:
  xlispstat pspp jags julia gfortran-multilib gfortran-doc gfortran-11-multilib gfortran-11-doc libcoarrays-dev liblapack-doc icu-doc liblzma-doc ncurses-doc readline-doc texlive-base texlive-latex-base texlive-plain-generic texlive-fonts-recommended
  texlive-fonts-extra texlive-extra-utils texlive-latex-recommended texlive-latex-extra texinfo mozilla | www-browser r-cran-cardata r-cran-latticeextra r-cran-colorspace r-cran-sfsmisc r-cran-sasmixed tcl-tclreadline
The following NEW packages will be installed:
  elpa-ess gfortran gfortran-11 icu-devtools libblas-dev libbz2-dev libdeflate-dev libgfortran-11-dev libicu-dev libjpeg-dev libjpeg-turbo8-dev libjpeg8-dev liblapack-dev liblzma-dev libncurses-dev libpcre2-16-0 libpcre2-32-0 libpcre2-dev libpcre2-posix3 libpkgconf3
  libpng-dev libpng-tools libreadline-dev pkgconf r-base-dev r-base-html r-cran-boot r-cran-class r-cran-cluster r-cran-codetools r-cran-foreign r-cran-kernsmooth r-cran-lattice r-cran-mass r-cran-matrix r-cran-mgcv r-cran-nlme r-cran-nnet r-cran-rpart r-cran-spatial
  r-cran-survival r-doc-html r-doc-info r-doc-pdf r-mathlib r-recommended tcl8.6 tk8.6
0 upgraded, 48 newly installed, 0 to remove and 8 not upgraded.
Need to get 67.6 MB of archives.
After this operation, 177 MB of additional disk space will be used.
Do you want to continue? [Y/n]
Get:1 http://archive.ubuntu.com/ubuntu jammy/universe amd64 elpa-ess all 18.10.2-2 [1217 kB]
Get:2 https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/ r-base-dev 4.4.1-3.2204.0 [4296 B]
Get:3 https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/ r-base-html 4.4.1-3.2204.0 [95.8 kB]
Get:4 https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/ r-doc-html 4.4.1-3.2204.0 [607 kB]
Get:5 https://r2u.stat.illinois.edu/ubuntu jammy/main all r-cran-boot all 1.3-31-1.ca2204.1 [638 kB]
Get:6 https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/ r-doc-info 4.4.1-3.2204.0 [669 kB]
Get:7 https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/ r-doc-pdf 4.4.1-3.2204.0 [9990 kB]
Get:8 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-mass amd64 7.3-61-1.ca2204.1 [1114 kB]
Get:9 https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/ r-mathlib 4.4.1-3.2204.0 [2644 kB]
Get:10 https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/ r-recommended 4.4.1-3.2204.0 [2590 B]
Get:11 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-class amd64 7.3-22-1.ca2204.1 [84.7 kB]
Get:12 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-cluster amd64 2.1.6-1.ca2204.1 [554 kB]
Get:13 https://r2u.stat.illinois.edu/ubuntu jammy/main all r-cran-codetools all 0.2-20-1.ca2204.1 [89.6 kB]
Get:14 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-foreign amd64 0.8.87-1.ca2204.1 [247 kB]
Get:15 http://archive.ubuntu.com/ubuntu jammy-updates/main amd64 libgfortran-11-dev amd64 11.4.0-1ubuntu1~22.04 [842 kB]
Get:16 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-kernsmooth amd64 2.23-24-1.ca2204.1 [91.4 kB]
Get:17 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-lattice amd64 0.22-6-1.ca2204.1 [1340 kB]
Get:18 http://archive.ubuntu.com/ubuntu jammy-updates/main amd64 gfortran-11 amd64 11.4.0-1ubuntu1~22.04 [11.2 MB]
Get:19 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-matrix amd64 1.7-1-1.ca2204.1 [4185 kB]
Get:20 http://archive.ubuntu.com/ubuntu jammy/main amd64 gfortran amd64 4:11.2.0-1ubuntu1 [1182 B]
Get:21 http://archive.ubuntu.com/ubuntu jammy/main amd64 icu-devtools amd64 70.1-2 [197 kB]
Get:22 http://archive.ubuntu.com/ubuntu jammy/main amd64 libblas-dev amd64 3.10.0-2ubuntu1 [164 kB]
Get:23 http://archive.ubuntu.com/ubuntu jammy/main amd64 libbz2-dev amd64 1.0.8-5build1 [32.5 kB]
Get:24 http://archive.ubuntu.com/ubuntu jammy/main amd64 libdeflate-dev amd64 1.10-2 [59.2 kB]
Get:25 http://archive.ubuntu.com/ubuntu jammy/main amd64 libicu-dev amd64 70.1-2 [11.6 MB]
Get:26 http://archive.ubuntu.com/ubuntu jammy/main amd64 libjpeg-turbo8-dev amd64 2.1.2-0ubuntu1 [257 kB]
Get:27 http://archive.ubuntu.com/ubuntu jammy/main amd64 libjpeg8-dev amd64 8c-2ubuntu10 [1476 B]
Get:28 http://archive.ubuntu.com/ubuntu jammy/main amd64 libjpeg-dev amd64 8c-2ubuntu10 [1472 B]
Get:29 http://archive.ubuntu.com/ubuntu jammy/main amd64 liblapack-dev amd64 3.10.0-2ubuntu1 [4774 kB]
Get:30 http://archive.ubuntu.com/ubuntu jammy-updates/main amd64 libncurses-dev amd64 6.3-2ubuntu0.1 [381 kB]
Get:31 http://archive.ubuntu.com/ubuntu jammy-updates/main amd64 libpcre2-16-0 amd64 10.39-3ubuntu0.1 [203 kB]
Get:32 http://archive.ubuntu.com/ubuntu jammy-updates/main amd64 libpcre2-32-0 amd64 10.39-3ubuntu0.1 [194 kB]
Get:33 http://archive.ubuntu.com/ubuntu jammy-updates/main amd64 libpcre2-posix3 amd64 10.39-3ubuntu0.1 [6130 B]
Get:34 http://archive.ubuntu.com/ubuntu jammy-updates/main amd64 libpcre2-dev amd64 10.39-3ubuntu0.1 [730 kB]
Get:35 http://archive.ubuntu.com/ubuntu jammy/universe amd64 libpkgconf3 amd64 1.8.0-1 [30.3 kB]
Get:36 http://archive.ubuntu.com/ubuntu jammy/main amd64 libpng-dev amd64 1.6.37-3build5 [192 kB]
Get:37 http://archive.ubuntu.com/ubuntu jammy/main amd64 libpng-tools amd64 1.6.37-3build5 [28.7 kB]
Get:38 http://archive.ubuntu.com/ubuntu jammy/main amd64 libreadline-dev amd64 8.1.2-1 [166 kB]
Get:39 http://archive.ubuntu.com/ubuntu jammy/universe amd64 pkgconf amd64 1.8.0-1 [35.3 kB]
Get:40 http://archive.ubuntu.com/ubuntu jammy/main amd64 liblzma-dev amd64 5.2.5-2ubuntu1 [159 kB]
Get:41 http://archive.ubuntu.com/ubuntu jammy/main amd64 tcl8.6 amd64 8.6.12+dfsg-1build1 [15.0 kB]
Get:42 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-nlme amd64 3.1.166-1.ca2204.1 [2318 kB]
Get:43 http://archive.ubuntu.com/ubuntu jammy/main amd64 tk8.6 amd64 8.6.12-1build1 [12.8 kB]
Get:44 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-mgcv amd64 1.9-1-1.ca2204.1 [3195 kB]
Get:45 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-nnet amd64 7.3-19-1.ca2204.1 [110 kB]
Get:46 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-rpart amd64 4.1.23-1.ca2204.1 [677 kB]
Get:47 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-spatial amd64 7.3-17-1.ca2204.1 [131 kB]
Get:48 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-survival amd64 3.7-0-1.ca2204.1 [6300 kB]
Fetched 67.6 MB in 4s (16.6 MB/s)
Extracting templates from packages: 100%
Selecting previously unselected package elpa-ess.
(Reading database ... 58470 files and directories currently installed.)
Preparing to unpack .../00-elpa-ess_18.10.2-2_all.deb ...
Unpacking elpa-ess (18.10.2-2) ...
Selecting previously unselected package libgfortran-11-dev:amd64.
Preparing to unpack .../01-libgfortran-11-dev_11.4.0-1ubuntu1~22.04_amd64.deb ...
Unpacking libgfortran-11-dev:amd64 (11.4.0-1ubuntu1~22.04) ...
Selecting previously unselected package gfortran-11.
Preparing to unpack .../02-gfortran-11_11.4.0-1ubuntu1~22.04_amd64.deb ...
Unpacking gfortran-11 (11.4.0-1ubuntu1~22.04) ...
Selecting previously unselected package gfortran.
Preparing to unpack .../03-gfortran_4%3a11.2.0-1ubuntu1_amd64.deb ...
Unpacking gfortran (4:11.2.0-1ubuntu1) ...
Selecting previously unselected package icu-devtools.
Preparing to unpack .../04-icu-devtools_70.1-2_amd64.deb ...
Unpacking icu-devtools (70.1-2) ...
Selecting previously unselected package libblas-dev:amd64.
Preparing to unpack .../05-libblas-dev_3.10.0-2ubuntu1_amd64.deb ...
Unpacking libblas-dev:amd64 (3.10.0-2ubuntu1) ...
Selecting previously unselected package libbz2-dev:amd64.
Preparing to unpack .../06-libbz2-dev_1.0.8-5build1_amd64.deb ...
Unpacking libbz2-dev:amd64 (1.0.8-5build1) ...
Selecting previously unselected package libdeflate-dev:amd64.
Preparing to unpack .../07-libdeflate-dev_1.10-2_amd64.deb ...
Unpacking libdeflate-dev:amd64 (1.10-2) ...
Selecting previously unselected package libicu-dev:amd64.
Preparing to unpack .../08-libicu-dev_70.1-2_amd64.deb ...
Unpacking libicu-dev:amd64 (70.1-2) ...
Selecting previously unselected package libjpeg-turbo8-dev:amd64.
Preparing to unpack .../09-libjpeg-turbo8-dev_2.1.2-0ubuntu1_amd64.deb ...
Unpacking libjpeg-turbo8-dev:amd64 (2.1.2-0ubuntu1) ...
Selecting previously unselected package libjpeg8-dev:amd64.
Preparing to unpack .../10-libjpeg8-dev_8c-2ubuntu10_amd64.deb ...
Unpacking libjpeg8-dev:amd64 (8c-2ubuntu10) ...
Selecting previously unselected package libjpeg-dev:amd64.
Preparing to unpack .../11-libjpeg-dev_8c-2ubuntu10_amd64.deb ...
Unpacking libjpeg-dev:amd64 (8c-2ubuntu10) ...
Selecting previously unselected package liblapack-dev:amd64.
Preparing to unpack .../12-liblapack-dev_3.10.0-2ubuntu1_amd64.deb ...
Unpacking liblapack-dev:amd64 (3.10.0-2ubuntu1) ...
Selecting previously unselected package libncurses-dev:amd64.
Preparing to unpack .../13-libncurses-dev_6.3-2ubuntu0.1_amd64.deb ...
Unpacking libncurses-dev:amd64 (6.3-2ubuntu0.1) ...
Selecting previously unselected package libpcre2-16-0:amd64.
Preparing to unpack .../14-libpcre2-16-0_10.39-3ubuntu0.1_amd64.deb ...
Unpacking libpcre2-16-0:amd64 (10.39-3ubuntu0.1) ...
Selecting previously unselected package libpcre2-32-0:amd64.
Preparing to unpack .../15-libpcre2-32-0_10.39-3ubuntu0.1_amd64.deb ...
Unpacking libpcre2-32-0:amd64 (10.39-3ubuntu0.1) ...
Selecting previously unselected package libpcre2-posix3:amd64.
Preparing to unpack .../16-libpcre2-posix3_10.39-3ubuntu0.1_amd64.deb ...
Unpacking libpcre2-posix3:amd64 (10.39-3ubuntu0.1) ...
Selecting previously unselected package libpcre2-dev:amd64.
Preparing to unpack .../17-libpcre2-dev_10.39-3ubuntu0.1_amd64.deb ...
Unpacking libpcre2-dev:amd64 (10.39-3ubuntu0.1) ...
Selecting previously unselected package libpkgconf3:amd64.
Preparing to unpack .../18-libpkgconf3_1.8.0-1_amd64.deb ...
Unpacking libpkgconf3:amd64 (1.8.0-1) ...
Selecting previously unselected package libpng-dev:amd64.
Preparing to unpack .../19-libpng-dev_1.6.37-3build5_amd64.deb ...
Unpacking libpng-dev:amd64 (1.6.37-3build5) ...
Selecting previously unselected package libpng-tools.
Preparing to unpack .../20-libpng-tools_1.6.37-3build5_amd64.deb ...
Unpacking libpng-tools (1.6.37-3build5) ...
Selecting previously unselected package libreadline-dev:amd64.
Preparing to unpack .../21-libreadline-dev_8.1.2-1_amd64.deb ...
Unpacking libreadline-dev:amd64 (8.1.2-1) ...
Selecting previously unselected package pkgconf.
Preparing to unpack .../22-pkgconf_1.8.0-1_amd64.deb ...
Adding 'diversion of /usr/bin/pkg-config to /usr/bin/pkg-config.real by pkgconf'
Adding 'diversion of /usr/share/aclocal/pkg.m4 to /usr/share/aclocal/pkg.real.m4 by pkgconf'
Adding 'diversion of /usr/share/man/man1/pkg-config.1.gz to /usr/share/man/man1/pkg-config.real.1.gz by pkgconf'
Adding 'diversion of /usr/share/pkg-config-crosswrapper to /usr/share/pkg-config-crosswrapper.real by pkgconf'
Unpacking pkgconf (1.8.0-1) ...
Selecting previously unselected package liblzma-dev:amd64.
Preparing to unpack .../23-liblzma-dev_5.2.5-2ubuntu1_amd64.deb ...
Unpacking liblzma-dev:amd64 (5.2.5-2ubuntu1) ...
Selecting previously unselected package r-base-dev.
Preparing to unpack .../24-r-base-dev_4.4.1-3.2204.0_all.deb ...
Unpacking r-base-dev (4.4.1-3.2204.0) ...
Selecting previously unselected package r-base-html.
Preparing to unpack .../25-r-base-html_4.4.1-3.2204.0_all.deb ...
Unpacking r-base-html (4.4.1-3.2204.0) ...
Selecting previously unselected package r-cran-boot.
Preparing to unpack .../26-r-cran-boot_1.3-31-1.ca2204.1_all.deb ...
Unpacking r-cran-boot (1.3-31-1.ca2204.1) ...
Selecting previously unselected package r-cran-mass.
Preparing to unpack .../27-r-cran-mass_7.3-61-1.ca2204.1_amd64.deb ...
Unpacking r-cran-mass (7.3-61-1.ca2204.1) ...
Selecting previously unselected package r-cran-class.
Preparing to unpack .../28-r-cran-class_7.3-22-1.ca2204.1_amd64.deb ...
Unpacking r-cran-class (7.3-22-1.ca2204.1) ...
Selecting previously unselected package r-cran-cluster.
Preparing to unpack .../29-r-cran-cluster_2.1.6-1.ca2204.1_amd64.deb ...
Unpacking r-cran-cluster (2.1.6-1.ca2204.1) ...
Selecting previously unselected package r-cran-codetools.
Preparing to unpack .../30-r-cran-codetools_0.2-20-1.ca2204.1_all.deb ...
Unpacking r-cran-codetools (0.2-20-1.ca2204.1) ...
Selecting previously unselected package r-cran-foreign.
Preparing to unpack .../31-r-cran-foreign_0.8.87-1.ca2204.1_amd64.deb ...
Unpacking r-cran-foreign (0.8.87-1.ca2204.1) ...
Selecting previously unselected package r-cran-kernsmooth.
Preparing to unpack .../32-r-cran-kernsmooth_2.23-24-1.ca2204.1_amd64.deb ...
Unpacking r-cran-kernsmooth (2.23-24-1.ca2204.1) ...
Selecting previously unselected package r-cran-lattice.
Preparing to unpack .../33-r-cran-lattice_0.22-6-1.ca2204.1_amd64.deb ...
Unpacking r-cran-lattice (0.22-6-1.ca2204.1) ...
Selecting previously unselected package r-cran-matrix.
Preparing to unpack .../34-r-cran-matrix_1.7-1-1.ca2204.1_amd64.deb ...
Unpacking r-cran-matrix (1.7-1-1.ca2204.1) ...
Selecting previously unselected package r-cran-nlme.
Preparing to unpack .../35-r-cran-nlme_3.1.166-1.ca2204.1_amd64.deb ...
Unpacking r-cran-nlme (3.1.166-1.ca2204.1) ...
Selecting previously unselected package r-cran-mgcv.
Preparing to unpack .../36-r-cran-mgcv_1.9-1-1.ca2204.1_amd64.deb ...
Unpacking r-cran-mgcv (1.9-1-1.ca2204.1) ...
Selecting previously unselected package r-cran-nnet.
Preparing to unpack .../37-r-cran-nnet_7.3-19-1.ca2204.1_amd64.deb ...
Unpacking r-cran-nnet (7.3-19-1.ca2204.1) ...
Selecting previously unselected package r-cran-rpart.
Preparing to unpack .../38-r-cran-rpart_4.1.23-1.ca2204.1_amd64.deb ...
Unpacking r-cran-rpart (4.1.23-1.ca2204.1) ...
Selecting previously unselected package r-cran-spatial.
Preparing to unpack .../39-r-cran-spatial_7.3-17-1.ca2204.1_amd64.deb ...
Unpacking r-cran-spatial (7.3-17-1.ca2204.1) ...
Selecting previously unselected package r-cran-survival.
Preparing to unpack .../40-r-cran-survival_3.7-0-1.ca2204.1_amd64.deb ...
Unpacking r-cran-survival (3.7-0-1.ca2204.1) ...
Selecting previously unselected package r-doc-html.
Preparing to unpack .../41-r-doc-html_4.4.1-3.2204.0_all.deb ...
Unpacking r-doc-html (4.4.1-3.2204.0) ...
Selecting previously unselected package r-doc-info.
Preparing to unpack .../42-r-doc-info_4.4.1-3.2204.0_all.deb ...
Unpacking r-doc-info (4.4.1-3.2204.0) ...
Selecting previously unselected package r-doc-pdf.
Preparing to unpack .../43-r-doc-pdf_4.4.1-3.2204.0_all.deb ...
Unpacking r-doc-pdf (4.4.1-3.2204.0) ...
Selecting previously unselected package r-mathlib.
Preparing to unpack .../44-r-mathlib_4.4.1-3.2204.0_amd64.deb ...
Unpacking r-mathlib (4.4.1-3.2204.0) ...
Selecting previously unselected package r-recommended.
Preparing to unpack .../45-r-recommended_4.4.1-3.2204.0_all.deb ...
Unpacking r-recommended (4.4.1-3.2204.0) ...
Selecting previously unselected package tcl8.6.
Preparing to unpack .../46-tcl8.6_8.6.12+dfsg-1build1_amd64.deb ...
Unpacking tcl8.6 (8.6.12+dfsg-1build1) ...
Selecting previously unselected package tk8.6.
Preparing to unpack .../47-tk8.6_8.6.12-1build1_amd64.deb ...
Unpacking tk8.6 (8.6.12-1build1) ...
Setting up tk8.6 (8.6.12-1build1) ...
Setting up r-cran-codetools (0.2-20-1.ca2204.1) ...
Setting up r-cran-boot (1.3-31-1.ca2204.1) ...
Setting up libjpeg-turbo8-dev:amd64 (2.1.2-0ubuntu1) ...
Setting up tcl8.6 (8.6.12+dfsg-1build1) ...
Setting up libncurses-dev:amd64 (6.3-2ubuntu0.1) ...
Setting up r-doc-pdf (4.4.1-3.2204.0) ...
Setting up r-doc-info (4.4.1-3.2204.0) ...
Setting up r-doc-html (4.4.1-3.2204.0) ...
Setting up libpng-tools (1.6.37-3build5) ...
Setting up libgfortran-11-dev:amd64 (11.4.0-1ubuntu1~22.04) ...
Setting up r-cran-spatial (7.3-17-1.ca2204.1) ...
Setting up libpng-dev:amd64 (1.6.37-3build5) ...
Setting up elpa-ess (18.10.2-2) ...
Install emacsen-common for emacs
emacsen-common: Handling install of emacsen flavor emacs
Install elpa-ess for emacs
install/ess-18.10.2: Handling install of emacsen flavor emacs
install/ess-18.10.2: byte-compiling for emacs
Setting up r-cran-rpart (4.1.23-1.ca2204.1) ...
Setting up r-cran-mass (7.3-61-1.ca2204.1) ...
Setting up libreadline-dev:amd64 (8.1.2-1) ...
Setting up r-cran-foreign (0.8.87-1.ca2204.1) ...
Setting up libpcre2-16-0:amd64 (10.39-3ubuntu0.1) ...
Setting up libpcre2-32-0:amd64 (10.39-3ubuntu0.1) ...
Setting up libpkgconf3:amd64 (1.8.0-1) ...
Setting up r-base-html (4.4.1-3.2204.0) ...
Setting up icu-devtools (70.1-2) ...
Setting up r-cran-lattice (0.22-6-1.ca2204.1) ...
Setting up r-cran-nlme (3.1.166-1.ca2204.1) ...
Setting up gfortran-11 (11.4.0-1ubuntu1~22.04) ...
Setting up r-cran-kernsmooth (2.23-24-1.ca2204.1) ...
Setting up liblzma-dev:amd64 (5.2.5-2ubuntu1) ...
Setting up libpcre2-posix3:amd64 (10.39-3ubuntu0.1) ...
Setting up r-cran-cluster (2.1.6-1.ca2204.1) ...
Setting up r-mathlib (4.4.1-3.2204.0) ...
Setting up r-cran-nnet (7.3-19-1.ca2204.1) ...
Setting up libjpeg8-dev:amd64 (8c-2ubuntu10) ...
Setting up libdeflate-dev:amd64 (1.10-2) ...
Setting up libicu-dev:amd64 (70.1-2) ...
Setting up libblas-dev:amd64 (3.10.0-2ubuntu1) ...
update-alternatives: using /usr/lib/x86_64-linux-gnu/blas/libblas.so to provide /usr/lib/x86_64-linux-gnu/libblas.so (libblas.so-x86_64-linux-gnu) in auto mode
Setting up libbz2-dev:amd64 (1.0.8-5build1) ...
Setting up r-cran-class (7.3-22-1.ca2204.1) ...
Setting up libpcre2-dev:amd64 (10.39-3ubuntu0.1) ...
Setting up libjpeg-dev:amd64 (8c-2ubuntu10) ...
Setting up gfortran (4:11.2.0-1ubuntu1) ...
update-alternatives: using /usr/bin/gfortran to provide /usr/bin/f95 (f95) in auto mode
update-alternatives: using /usr/bin/gfortran to provide /usr/bin/f77 (f77) in auto mode
Setting up pkgconf (1.8.0-1) ...
Setting up liblapack-dev:amd64 (3.10.0-2ubuntu1) ...
update-alternatives: using /usr/lib/x86_64-linux-gnu/lapack/liblapack.so to provide /usr/lib/x86_64-linux-gnu/liblapack.so (liblapack.so-x86_64-linux-gnu) in auto mode
Setting up r-cran-matrix (1.7-1-1.ca2204.1) ...
Setting up r-cran-mgcv (1.9-1-1.ca2204.1) ...
Setting up r-base-dev (4.4.1-3.2204.0) ...
Setting up r-cran-survival (3.7-0-1.ca2204.1) ...
Setting up r-recommended (4.4.1-3.2204.0) ...
Processing triggers for man-db (2.10.2-1) ...
Processing triggers for install-info (6.8-4build1) ...
Processing triggers for libc-bin (2.35-0ubuntu3.8) ...
```

Run R to see what we get:

```console
(sequin) weismanal@NCATS-2260152-P:~/notebook/2024-10-18/testing_sequin_installation_on_linux$ R

R version 4.4.1 (2024-06-14) -- "Race for Your Life"
Copyright (C) 2024 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

Warning message:
D-Bus service not found!

- If you are in a container environment, please consider adding the
  following to your configuration to silence this warning:

  options(bspm.sudo = TRUE)

- If you are in a desktop/server environment, please remove any 'bspm'
  installation from the user library and force a new system
  installation as follows:

  $ sudo Rscript --vanilla -e 'install.packages("bspm", repos="https://cran.r-project.org")'
```

Reinstall `bspm` per the above message, though that didn't seem to fix the warning above:

```console
(sequin) weismanal@NCATS-2260152-P:~/notebook/2024-10-18/testing_sequin_installation_on_linux$ sudo Rscript --vanilla -e 'install.packages("bspm", repos="https://cran.r-project.org")'
Installing package into ‘/usr/local/lib/R/site-library’
(as ‘lib’ is unspecified)
trying URL 'https://cran.r-project.org/src/contrib/bspm_0.5.7.tar.gz'
Content type 'application/x-gzip' length 27097 bytes (26 KB)
==================================================
downloaded 26 KB

* installing *source* package ‘bspm’ ...
** package ‘bspm’ successfully unpacked and MD5 sums checked
** using staged installation
* installing /usr/share/dbus-1/system-services/org.r_project.linux1.service
* installing /etc/dbus-1/system.d/org.r_project.linux1.conf
** R
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (bspm)

The downloaded source packages are in
        ‘/tmp/Rtmpwj3rH5/downloaded_packages’
```

Still per [this website](https://cran.r-project.org/bin/linux/ubuntu/fullREADME.html), install `r-base`:

```console
(sequin) weismanal@NCATS-2260152-P:~/notebook/2024-10-18/testing_sequin_installation_on_linux$ sudo apt-get install r-base
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
The following NEW packages will be installed:
  r-base
0 upgraded, 1 newly installed, 0 to remove and 8 not upgraded.
Need to get 47.2 kB of archives.
After this operation, 66.6 kB of additional disk space will be used.
Get:1 https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/ r-base 4.4.1-3.2204.0 [47.2 kB]
Fetched 47.2 kB in 0s (122 kB/s)
Selecting previously unselected package r-base.
(Reading database ... 60416 files and directories currently installed.)
Preparing to unpack .../r-base_4.4.1-3.2204.0_all.deb ...
Unpacking r-base (4.4.1-3.2204.0) ...
Setting up r-base (4.4.1-3.2204.0) ...
```

Note that `r-base-dev` was already installed above. Now R should be fully set up and we can move on to Andrei's [SEQUIN installation instructions here](https://github.com/ncats/public_sequin/tree/main).

## Install SEQUIN

Clone the `public_sequin` repository using:

```bash
git clone https://github.com/ncats/public_sequin.git
```

DO NOT enter the cloned directory before the next steps; just say in the current directory.

Follow the [SEQUIN installation instructions](https://github.com/ncats/public_sequin/tree/main). Note that R needs to be run with `sudo` to install packages using `install_sequin.R`, which also addresses the above warning when starting R:

```console
(sequin) weismanal@NCATS-2260152-P:~/notebook/2024-10-18/testing_sequin_installation_on_linux$ sudo R

R version 4.4.1 (2024-06-14) -- "Race for Your Life"
Copyright (C) 2024 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> setwd("/home/weismanal/notebook/2024-10-18/testing_sequin_installation_on_linux/public_sequin")
> source("install_sequin.R")
Install system packages as root...
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
Fetched 0 B in 0s (0 B/s)
Install system packages as root...
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
Fetched 0 B in 0s (0 B/s)

renv: Project Environments for R

Welcome to renv! It looks like this is your first time using renv.
This is a one-time message, briefly describing some of renv's functionality.

renv will write to files within the active project folder, including:

  - A folder 'renv' in the project directory, and
  - A lockfile called 'renv.lock' in the project directory.

In particular, projects using renv will normally use a private, per-project
R library, in which new packages will be installed. This project library is
isolated from other R libraries on your system.

In addition, renv will update files within your project directory, including:

  - .gitignore
  - .Rbuildignore
  - .Rprofile

Finally, renv maintains a local cache of data on the filesystem, located at:

  - "~/.cache/R/renv"

This path can be customized: please see the documentation in `?renv::paths`.

Please read the introduction vignette with `vignette("renv")` for more information.
You can browse the package documentation online at https://rstudio.github.io/renv/.
Do you want to proceed? [y/N]: y

- "~/.cache/R/renv" has been created.
- renv activated -- please restart the R session.
'getOption("repos")' replaces Bioconductor standard repositories, see
'help("repositories", package = "BiocManager")' for details.
Replacement repositories:
    CRAN: https://mirrors.nics.utk.edu/cran
# Downloading packages -------------------------------------------------------
- Downloading AnnotationDbi from BioCsoft ...   OK [4.2 Mb in 0.21s]
- Downloading DBI from CRAN ...                 OK [1.1 Mb in 0.36s]
...
```

Eventually, the script dies with a `curl` installation error. Quit R and run:

```console
(sequin) weismanal@NCATS-2260152-P:~/notebook/2024-10-18/testing_sequin_installation_on_linux$ sudo apt install libcurl4-openssl-dev
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
Suggested packages:
  libcurl4-doc libidn11-dev libkrb5-dev libldap2-dev librtmp-dev libssh2-1-dev libssl-dev
The following NEW packages will be installed:
  libcurl4-openssl-dev
0 upgraded, 1 newly installed, 0 to remove and 8 not upgraded.
Need to get 386 kB of archives.
After this operation, 1698 kB of additional disk space will be used.
Get:1 http://archive.ubuntu.com/ubuntu jammy-updates/main amd64 libcurl4-openssl-dev amd64 7.81.0-1ubuntu1.18 [386 kB]
Fetched 386 kB in 0s (1000 kB/s)
Selecting previously unselected package libcurl4-openssl-dev:amd64.
(Reading database ... 60548 files and directories currently installed.)
Preparing to unpack .../libcurl4-openssl-dev_7.81.0-1ubuntu1.18_amd64.deb ...
Unpacking libcurl4-openssl-dev:amd64 (7.81.0-1ubuntu1.18) ...
Setting up libcurl4-openssl-dev:amd64 (7.81.0-1ubuntu1.18) ...
Processing triggers for man-db (2.10.2-1) ...
```

Doing this sort of thing for many additional failed R packages,

- `openssl` requires `sudo apt install libssl-dev`
- `xml2` requires `sudo apt install libxml2-dev`
- `Cairo` requires `sudo apt install libcairo2-dev`
- `RMariaDB` requires `sudo apt install libmariadb-dev`
- `textshaping` requires `sudo apt install libharfbuzz-dev libfribidi-dev`
- `ragg` requires `sudo apt install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev`

required quitting R and installing the above operating system packages, each time re-running:

```console
sudo R
setwd("/home/weismanal/notebook/2024-10-18/testing_sequin_installation_on_linux/public_sequin")
source("install_sequin.R")
```

After all this, it seemed to work:

```console
(sequin) weismanal@NCATS-2260152-P:~/notebook/2024-10-18/testing_sequin_installation_on_linux$ sudo R

R version 4.4.1 (2024-06-14) -- "Race for Your Life"
Copyright (C) 2024 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> setwd("/home/weismanal/notebook/2024-10-18/testing_sequin_installation_on_linux/public_sequin")
source("install_sequin.R")
Install system packages as root...
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
Fetched 0 B in 0s (0 B/s)
Install system packages as root...
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
Fetched 0 B in 0s (0 B/s)
- renv activated -- please restart the R session.
'getOption("repos")' replaces Bioconductor standard repositories, see
'help("repositories", package = "BiocManager")' for details.
Replacement repositories:
    CRAN: https://mirrors.nics.utk.edu/cran
The following package(s) will be installed:
- annotate             [1.82.0]
- AnnotationDbi        [1.66.0]
- aroma.light          [3.34.0]
- Biobase              [2.64.0]
- BiocFileCache        [2.12.0]
- BiocGenerics         [0.50.0]
- BiocIO               [1.14.0]
- BiocParallel         [1.38.0]
- BiocVersion          [3.19.1]
- biomaRt              [2.60.1]
- Biostrings           [2.72.1]
- ComplexHeatmap       [2.20.0]
- DelayedArray         [0.30.1]
- DESeq2               [1.44.0]
- EDASeq               [2.38.0]
- edgeR                [4.2.2]
- geneplotter          [1.82.0]
- GenomeInfoDb         [1.40.1]
- GenomeInfoDbData     [1.2.12]
- GenomicAlignments    [1.40.0]
- GenomicFeatures      [1.56.0]
- GenomicRanges        [1.56.2]
- GO.db                [3.19.1]
- impute               [1.78.0]
- IRanges              [2.38.1]
- KEGGREST             [1.44.1]
- limma                [3.60.6]
- MAST                 [1.30.0]
- MatrixGenerics       [1.16.0]
- org.Hs.eg.db         [3.19.1]
- org.Mm.eg.db         [3.19.1]
- preprocessCore       [1.66.0]
- Rhtslib              [3.0.0]
- Rsamtools            [2.20.0]
- rtracklayer          [1.64.0]
- RUVSeq               [1.38.0]
- S4Arrays             [1.4.1]
- S4Vectors            [0.42.1]
- ShortRead            [1.62.0]
- SingleCellExperiment [1.26.0]
- SparseArray          [1.4.8]
- SummarizedExperiment [1.34.0]
- XVector              [0.44.0]
- zlibbioc             [1.50.0]
These packages will be installed into "/usr/local/lib/R/site-library".

Do you want to proceed? [Y/n]: Y

# Installing packages --------------------------------------------------------
- Installing BiocVersion ...                    OK [copied from cache]
- Installing BiocGenerics ...                   OK [copied from cache]
- Installing S4Vectors ...                      OK [copied from cache]
- Installing IRanges ...                        OK [copied from cache]
- Installing zlibbioc ...                       OK [copied from cache]
- Installing XVector ...                        OK [copied from cache]
- Installing GenomeInfoDbData ...               OK [copied from cache]
- Installing GenomeInfoDb ...                   OK [copied from cache]
- Installing Biostrings ...                     OK [copied from cache]
- Installing KEGGREST ...                       OK [copied from cache]
- Installing Biobase ...                        OK [copied from cache]
- Installing AnnotationDbi ...                  OK [copied from cache]
- Installing BiocFileCache ...                  OK [copied from cache]
- Installing BiocIO ...                         OK [copied from cache]
- Installing BiocParallel ...                   OK [copied from cache]
- Installing ComplexHeatmap ...                 OK [copied from cache]
- Installing MatrixGenerics ...                 OK [copied from cache]
- Installing GenomicRanges ...                  OK [copied from cache]
- Installing S4Arrays ...                       OK [copied from cache]
- Installing SparseArray ...                    OK [copied from cache]
- Installing DelayedArray ...                   OK [copied from cache]
- Installing SummarizedExperiment ...           OK [copied from cache]
- Installing DESeq2 ...                         OK [copied from cache]
- Installing aroma.light ...                    OK [copied from cache]
- Installing Rhtslib ...                        OK [copied from cache]
- Installing Rsamtools ...                      OK [copied from cache]
- Installing biomaRt ...                        OK [copied from cache]
- Installing GenomicAlignments ...              OK [copied from cache]
- Installing rtracklayer ...                    OK [copied from cache]
- Installing GenomicFeatures ...                OK [copied from cache]
- Installing ShortRead ...                      OK [copied from cache]
- Installing EDASeq ...                         OK [copied from cache]
- Installing GO.db ...                          OK [copied from cache]
- Installing SingleCellExperiment ...           OK [copied from cache]
- Installing MAST ...                           OK [copied from cache]
- Installing limma ...                          OK [copied from cache]
- Installing edgeR ...                          OK [copied from cache]
- Installing RUVSeq ...                         OK [copied from cache]
- Installing annotate ...                       OK [copied from cache]
- Installing geneplotter ...                    OK [copied from cache]
- Installing impute ...                         OK [copied from cache]
- Installing org.Hs.eg.db ...                   OK [copied from cache in 0.2s]
- Installing org.Mm.eg.db ...                   OK [copied from cache in 0.17s]
- Installing preprocessCore ...                 OK [copied from cache]
Successfully installed 44 packages in 1.1 seconds.
# Downloading packages -------------------------------------------------------
- Downloading devtools from CRAN ...            OK [file is up to date]
- Downloading pkgdown from CRAN ...             OK [file is up to date]
- Downloading ragg from CRAN ...                OK [file is up to date]
- Downloading whisker from CRAN ...             OK [file is up to date]
- Downloading pkgload from CRAN ...             OK [file is up to date]
- Downloading profvis from CRAN ...             OK [file is up to date]
- Downloading rcmdcheck from CRAN ...           OK [file is up to date]
- Downloading sessioninfo from CRAN ...         OK [file is up to date]
- Downloading xopen from CRAN ...               OK [file is up to date]
- Downloading remotes from CRAN ...             OK [file is up to date]
- Downloading roxygen2 from CRAN ...            OK [file is up to date]
- Downloading rversions from CRAN ...           OK [file is up to date]
- Downloading testthat from CRAN ...            OK [file is up to date]
- Downloading praise from CRAN ...              OK [file is up to date]
- Downloading waldo from CRAN ...               OK [file is up to date]
- Downloading diffobj from CRAN ...             OK [file is up to date]
- Downloading rematch2 from CRAN ...            OK [file is up to date]
- Downloading urlchecker from CRAN ...          OK [file is up to date]
- Downloading usethis from CRAN ...             OK [file is up to date]
- Downloading gert from CRAN ...                OK [file is up to date]
- Downloading zip from CRAN ...                 OK [file is up to date]
- Downloading gh from CRAN ...                  OK [file is up to date]
- Downloading gitcreds from CRAN ...            OK [file is up to date]
- Downloading ini from CRAN ...                 OK [file is up to date]
- Downloading egg from CRAN ...                 OK [file is up to date]
- Downloading enrichR from CRAN ...             OK [file is up to date]
- Downloading fields from CRAN ...              OK [file is up to date]
- Downloading maps from CRAN ...                OK [file is up to date]
- Downloading flashClust from CRAN ...          OK [file is up to date]
- Downloading gclus from CRAN ...               OK [file is up to date]
- Downloading gfonts from CRAN ...              OK [file is up to date]
- Downloading ggdendro from CRAN ...            OK [file is up to date]
- Downloading heatmaply from CRAN ...           OK [file is up to date]
- Downloading seriation from CRAN ...           OK [file is up to date]
- Downloading qap from CRAN ...                 OK [file is up to date]
- Downloading registry from CRAN ...            OK [file is up to date]
- Downloading vegan from CRAN ...               OK [file is up to date]
- Downloading permute from CRAN ...             OK [file is up to date]
- Downloading webshot from CRAN ...             OK [file is up to date]
- Downloading kmed from CRAN ...                OK [file is up to date]
- Downloading markdown from CRAN ...            OK [file is up to date]
- Downloading mnormt from CRAN ...              OK [file is up to date]
- Downloading openxlsx from CRAN ...            OK [file is up to date]
- Downloading pheatmap from CRAN ...            OK [file is up to date]
- Downloading plotrix from CRAN ...             OK [file is up to date]
- Downloading psych from CRAN ...               OK [file is up to date]
- Downloading GPArotation from CRAN ...         OK [file is up to date]
- Downloading renv from CRAN ...                OK [file is up to date]
- Downloading reshape from CRAN ...             OK [file is up to date]
- Downloading rvest from CRAN ...               OK [file is up to date]
- Downloading selectr from CRAN ...             OK [file is up to date]
- Downloading shinyBS from CRAN ...             OK [file is up to date]
- Downloading shinyWidgets from CRAN ...        OK [file is up to date]
- Downloading shinycssloaders from CRAN ...     OK [file is up to date]
- Downloading shinyjs from CRAN ...             OK [file is up to date]
- Downloading shinythemes from CRAN ...         OK [file is up to date]
- Downloading showtext from CRAN ...            OK [file is up to date]
- Downloading sysfonts from CRAN ...            OK [file is up to date]
- Downloading showtextdb from CRAN ...          OK [file is up to date]
- Downloading spatial from CRAN ...             OK [file is up to date]
- Downloading umap from CRAN ...                OK [file is up to date]
Successfully downloaded 61 packages in 21 seconds.

The following package(s) will be installed:
- abind            [1.4-8]
- anytime          [0.3.9]
- ape              [5.8]
- askpass          [1.2.1]
- assertthat       [0.2.1]
- backports        [1.5.0]
- base64enc        [0.1-3]
- BH               [1.84.0-0]
- BiocManager      [1.30.25]
- bit              [4.5.0]
- bit64            [4.5.2]
- bitops           [1.0-9]
- blob             [1.2.4]
- boot             [1.3-31]
- brew             [1.0-10]
- brio             [1.1.5]
- broom            [1.0.7]
- broom.helpers    [1.17.0]
- bslib            [0.8.0]
- ca               [0.71.1]
- cachem           [1.1.0]
- Cairo            [1.6-2]
- callr            [3.7.6]
- cards            [0.3.0]
- caTools          [1.18.3]
- checkmate        [2.3.2]
- circlize         [0.4.16]
- class            [7.3-22]
- cli              [3.6.3]
- clipr            [0.8.0]
- clue             [0.3-65]
- cluster          [2.1.6]
- clustree         [0.5.1]
- codetools        [0.2-20]
- colorspace       [2.1-1]
- commonmark       [1.9.2]
- cowplot          [1.1.3]
- cpp11            [0.5.0]
- crayon           [1.5.3]
- credentials      [2.0.2]
- crosstalk        [1.2.1]
- crul             [1.5.0]
- curl             [5.2.3]
- data.table       [1.16.2]
- DBI              [1.2.3]
- dbplyr           [2.5.0]
- deldir           [2.0-4]
- dendextend       [1.18.1]
- desc             [1.4.3]
- devtools         [2.4.5]
- diffobj          [0.3.5]
- digest           [0.6.37]
- doParallel       [1.0.17]
- dotCall64        [1.2]
- downlit          [0.4.4]
- dplyr            [1.1.4]
- dqrng            [0.4.1]
- DT               [0.33]
- dynamicTreeCut   [1.63-1]
- egg              [0.4.5]
- ellipsis         [0.3.2]
- enrichR          [3.2]
- evaluate         [1.0.1]
- expm             [1.0-0]
- fansi            [1.0.6]
- farver           [2.1.2]
- fastcluster      [1.2.6]
- fastDummies      [1.7.4]
- fastmap          [1.2.0]
- fields           [16.3]
- filelock         [1.0.3]
- fitdistrplus     [1.2-1]
- flashClust       [1.01-2]
- FNN              [1.1.4.1]
- fontawesome      [0.5.2]
- forcats          [1.0.0]
- foreach          [1.5.2]
- foreign          [0.8-87]
- formatR          [1.14]
- Formula          [1.2-5]
- fs               [1.6.4]
- futile.logger    [1.4.3]
- futile.options   [1.0.1]
- future           [1.34.0]
- future.apply     [1.11.2]
- gclus            [1.3.2]
- generics         [0.1.3]
- gert             [2.1.4]
- GetoptLong       [1.0.5]
- gfonts           [0.2.0]
- GGally           [2.2.1]
- ggdendro         [0.2.0]
- ggforce          [0.4.2]
- ggplot2          [3.5.1]
- ggraph           [2.2.1]
- ggrepel          [0.9.6]
- ggridges         [0.5.6]
- ggstats          [0.7.0]
- gh               [1.4.1]
- gitcreds         [0.1.2]
- GlobalOptions    [0.1.2]
- globals          [0.16.3]
- glue             [1.8.0]
- goftest          [1.2-3]
- GPArotation      [2024.3-1]
- gplots           [3.2.0]
- graphlayouts     [1.2.0]
- gridExtra        [2.3]
- gtable           [0.3.5]
- gtools           [3.9.5]
- haven            [2.5.4]
- heatmaply        [1.5.0]
- here             [1.0.1]
- highr            [0.11]
- Hmisc            [5.1-3]
- hms              [1.1.3]
- htmlTable        [2.4.3]
- htmltools        [0.5.8.1]
- htmlwidgets      [1.6.4]
- httpcode         [0.3.0]
- httpuv           [1.6.15]
- httr             [1.4.7]
- httr2            [1.0.5]
- hwriter          [1.3.2.1]
- ica              [1.0-3]
- igraph           [2.1.1]
- ini              [0.3.1]
- interp           [1.1-6]
- irlba            [2.3.5.1]
- isoband          [0.2.7]
- iterators        [1.0.14]
- jpeg             [0.1-10]
- jquerylib        [0.1.4]
- jsonlite         [1.8.9]
- KernSmooth       [2.23-24]
- kmed             [0.4.2]
- knitr            [1.48]
- labeling         [0.4.3]
- labelled         [2.13.0]
- lambda.r         [1.2.4]
- later            [1.3.2]
- lattice          [0.22-6]
- latticeExtra     [0.6-30]
- lazyeval         [0.2.2]
- leiden           [0.4.3.1]
- lifecycle        [1.0.4]
- listenv          [0.9.1]
- lmtest           [0.9-40]
- locfit           [1.5-9.10]
- lubridate        [1.9.3]
- magrittr         [2.0.3]
- maps             [3.4.2]
- markdown         [1.13]
- MASS             [7.3-61]
- Matrix           [1.7-1]
- matrixStats      [1.4.1]
- MCL              [1.0]
- memoise          [2.0.1]
- mgcv             [1.9-1]
- mime             [0.12]
- miniUI           [0.1.1.1]
- mnormt           [2.1.1]
- munsell          [0.5.1]
- nlme             [3.1-166]
- nnet             [7.3-19]
- openssl          [2.2.2]
- openxlsx         [4.2.7.1]
- parallelly       [1.38.0]
- patchwork        [1.3.0]
- pbapply          [1.7-2]
- permute          [0.9-7]
- pheatmap         [1.0.12]
- pillar           [1.9.0]
- pkgbuild         [1.4.4]
- pkgconfig        [2.0.3]
- pkgdown          [2.1.1]
- pkgload          [1.4.0]
- plogr            [0.2.0]
- plotly           [4.10.4]
- plotrix          [3.8-4]
- plyr             [1.8.9]
- png              [0.1-8]
- polyclip         [1.10-7]
- praise           [1.0.0]
- prettyunits      [1.2.0]
- processx         [3.8.4]
- profvis          [0.4.0]
- progress         [1.2.3]
- progressr        [0.14.0]
- promises         [1.3.0]
- ps               [1.8.0]
- psych            [2.4.6.26]
- purrr            [1.0.2]
- qap              [0.1-2]
- R.methodsS3      [1.8.2]
- R.oo             [1.26.0]
- R.utils          [2.12.3]
- R6               [2.5.1]
- ragg             [1.3.3]
- RANN             [2.6.2]
- rappdirs         [0.3.3]
- rcmdcheck        [1.4.0]
- RColorBrewer     [1.1-3]
- Rcpp             [1.0.13]
- RcppAnnoy        [0.0.22]
- RcppArmadillo    [14.0.2-1]
- RcppEigen        [0.3.4.0.2]
- RcppHNSW         [0.6.0]
- RcppProgress     [0.4.2]
- RcppTOML         [0.2.2]
- RCurl            [1.98-1.16]
- readr            [2.1.5]
- registry         [0.5-1]
- rematch2         [2.1.2]
- remotes          [2.5.0]
- renv             [1.0.11]
- reshape          [0.8.9]
- reshape2         [1.4.4]
- restfulr         [0.0.15]
- reticulate       [1.39.0]
- rjson            [0.2.23]
- rlang            [1.1.4]
- RMariaDB         [1.3.2]
- rmarkdown        [2.28]
- ROCR             [1.0-11]
- roxygen2         [7.3.2]
- rpart            [4.1.23]
- rprojroot        [2.0.4]
- RSpectra         [0.16-2]
- RSQLite          [2.3.7]
- rstudioapi       [0.17.0]
- Rtsne            [0.17]
- rversions        [2.1.2]
- rvest            [1.0.4]
- sass             [0.4.9]
- scales           [1.3.0]
- scattermore      [1.2]
- sctransform      [0.4.1]
- selectr          [0.4-2]
- seriation        [1.5.6]
- sessioninfo      [1.2.2]
- Seurat           [5.1.0]
- SeuratObject     [5.0.2]
- shape            [1.4.6.1]
- shiny            [1.9.1]
- shinyBS          [0.61.1]
- shinycssloaders  [1.1.0]
- shinyjs          [2.1.0]
- shinythemes      [1.2.0]
- shinyWidgets     [0.8.7]
- showtext         [0.9-7]
- showtextdb       [3.0]
- sitmo            [2.0.2]
- snow             [0.4-4]
- sourcetools      [0.1.7-1]
- sp               [2.1-4]
- spam             [2.11-0]
- spatial          [7.3-17]
- spatstat.data    [3.1-2]
- spatstat.explore [3.3-2]
- spatstat.geom    [3.3-3]
- spatstat.random  [3.3-2]
- spatstat.sparse  [3.1-0]
- spatstat.univar  [3.0-1]
- spatstat.utils   [3.1-0]
- statmod          [1.5.0]
- stringi          [1.8.4]
- stringr          [1.5.1]
- survival         [3.7-0]
- sys              [3.4.3]
- sysfonts         [0.8.9]
- systemfonts      [1.1.0]
- TeachingDemos    [2.13]
- tensor           [1.5]
- testthat         [3.2.1.1]
- textshaping      [0.4.0]
- tibble           [3.2.1]
- tidygraph        [1.3.1]
- tidyr            [1.3.1]
- tidyselect       [1.2.1]
- timechange       [0.3.0]
- tinytex          [0.53]
- triebeard        [0.4.1]
- TSP              [1.2-4]
- tweenr           [2.0.3]
- tzdb             [0.4.0]
- umap             [0.2.10.0]
- urlchecker       [1.0.1]
- urltools         [1.7.3]
- usethis          [3.0.0]
- utf8             [1.2.4]
- uwot             [0.2.2]
- vctrs            [0.6.5]
- vegan            [2.6-8]
- viridis          [0.6.5]
- viridisLite      [0.4.2]
- vroom            [1.6.5]
- waldo            [0.5.3]
- webshot          [0.5.5]
- WGCNA            [1.73]
- whisker          [0.4.1]
- withr            [3.0.1]
- WriteXLS         [6.7.0]
- xfun             [0.48]
- XML              [3.99-0.17]
- xml2             [1.3.6]
- xopen            [1.0.1]
- xtable           [1.8-4]
- yaml             [2.3.10]
- zip              [2.3.1]
- zoo              [1.8-12]
These packages will be installed into "/usr/local/lib/R/site-library".

Do you want to proceed? [Y/n]: Y

# Installing packages --------------------------------------------------------
- Installing BH ...                             OK [copied from cache in 0.23s]
- Installing BiocManager ...                    OK [copied from cache]
- Installing Cairo ...                          OK [copied from cache]
- Installing DBI ...                            OK [copied from cache]
- Installing base64enc ...                      OK [copied from cache]
- Installing digest ...                         OK [copied from cache]
- Installing fastmap ...                        OK [copied from cache]
- Installing rlang ...                          OK [copied from cache]
- Installing htmltools ...                      OK [copied from cache]
- Installing jsonlite ...                       OK [copied from cache]
- Installing evaluate ...                       OK [copied from cache]
- Installing xfun ...                           OK [copied from cache]
- Installing highr ...                          OK [copied from cache]
- Installing yaml ...                           OK [copied from cache]
- Installing knitr ...                          OK [copied from cache]
- Installing cachem ...                         OK [copied from cache]
- Installing jquerylib ...                      OK [copied from cache]
- Installing cli ...                            OK [copied from cache]
- Installing glue ...                           OK [copied from cache]
- Installing lifecycle ...                      OK [copied from cache]
- Installing memoise ...                        OK [copied from cache]
- Installing mime ...                           OK [copied from cache]
- Installing fs ...                             OK [copied from cache]
- Installing R6 ...                             OK [copied from cache]
- Installing rappdirs ...                       OK [copied from cache]
- Installing sass ...                           OK [copied from cache]
- Installing bslib ...                          OK [copied from cache]
- Installing fontawesome ...                    OK [copied from cache]
- Installing tinytex ...                        OK [copied from cache]
- Installing rmarkdown ...                      OK [copied from cache]
- Installing htmlwidgets ...                    OK [copied from cache]
- Installing Rcpp ...                           OK [copied from cache]
- Installing later ...                          OK [copied from cache]
- Installing magrittr ...                       OK [copied from cache]
- Installing promises ...                       OK [copied from cache]
- Installing httpuv ...                         OK [copied from cache]
- Installing lazyeval ...                       OK [copied from cache]
- Installing crosstalk ...                      OK [copied from cache]
- Installing DT ...                             OK [copied from cache]
- Installing FNN ...                            OK [copied from cache]
- Installing Formula ...                        OK [copied from cache]
- Installing generics ...                       OK [copied from cache]
- Installing fansi ...                          OK [copied from cache]
- Installing utf8 ...                           OK [copied from cache]
- Installing vctrs ...                          OK [copied from cache]
- Installing pillar ...                         OK [copied from cache]
- Installing pkgconfig ...                      OK [copied from cache]
- Installing tibble ...                         OK [copied from cache]
- Installing withr ...                          OK [copied from cache]
- Installing tidyselect ...                     OK [copied from cache]
- Installing dplyr ...                          OK [copied from cache]
- Installing purrr ...                          OK [copied from cache]
- Installing stringi ...                        OK [copied from cache]
- Installing stringr ...                        OK [copied from cache]
- Installing cpp11 ...                          OK [copied from cache]
- Installing tidyr ...                          OK [copied from cache]
- Installing forcats ...                        OK [copied from cache]
- Installing gtable ...                         OK [copied from cache]
- Installing isoband ...                        OK [copied from cache]
- Installing MASS ...                           OK [copied from cache]
- Installing lattice ...                        OK [copied from cache]
- Installing Matrix ...                         OK [copied from cache]
- Installing nlme ...                           OK [copied from cache]
- Installing mgcv ...                           OK [copied from cache]
- Installing farver ...                         OK [copied from cache]
- Installing labeling ...                       OK [copied from cache]
- Installing colorspace ...                     OK [copied from cache]
- Installing munsell ...                        OK [copied from cache]
- Installing RColorBrewer ...                   OK [copied from cache]
- Installing viridisLite ...                    OK [copied from cache]
- Installing scales ...                         OK [copied from cache]
- Installing ggplot2 ...                        OK [copied from cache]
- Installing patchwork ...                      OK [copied from cache]
- Installing ggstats ...                        OK [copied from cache]
- Installing plyr ...                           OK [copied from cache]
- Installing crayon ...                         OK [copied from cache]
- Installing hms ...                            OK [copied from cache]
- Installing prettyunits ...                    OK [copied from cache]
- Installing progress ...                       OK [copied from cache]
- Installing GGally ...                         OK [copied from cache]
- Installing rjson ...                          OK [copied from cache]
- Installing GlobalOptions ...                  OK [copied from cache]
- Installing GetoptLong ...                     OK [copied from cache]
- Installing cluster ...                        OK [copied from cache]
- Installing rpart ...                          OK [copied from cache]
- Installing nnet ...                           OK [copied from cache]
- Installing foreign ...                        OK [copied from cache]
- Installing gridExtra ...                      OK [copied from cache]
- Installing data.table ...                     OK [copied from cache]
- Installing backports ...                      OK [copied from cache]
- Installing checkmate ...                      OK [copied from cache]
- Installing rstudioapi ...                     OK [copied from cache]
- Installing htmlTable ...                      OK [copied from cache]
- Installing viridis ...                        OK [copied from cache]
- Installing Hmisc ...                          OK [copied from cache]
- Installing KernSmooth ...                     OK [copied from cache]
- Installing expm ...                           OK [copied from cache]
- Installing MCL ...                            OK [copied from cache]
- Installing R.methodsS3 ...                    OK [copied from cache]
- Installing R.oo ...                           OK [copied from cache]
- Installing R.utils ...                        OK [copied from cache]
- Installing RANN ...                           OK [copied from cache]
- Installing bitops ...                         OK [copied from cache]
- Installing RCurl ...                          OK [copied from cache]
- Installing bit ...                            OK [copied from cache]
- Installing bit64 ...                          OK [copied from cache]
- Installing blob ...                           OK [copied from cache]
- Installing timechange ...                     OK [copied from cache]
- Installing lubridate ...                      OK [copied from cache]
- Installing plogr ...                          OK [copied from cache]
- Installing RMariaDB ...                       OK [copied from cache]
- Installing gtools ...                         OK [copied from cache]
- Installing caTools ...                        OK [copied from cache]
- Installing gplots ...                         OK [copied from cache]
- Installing ROCR ...                           OK [copied from cache]
- Installing RSQLite ...                        OK [copied from cache]
- Installing RcppEigen ...                      OK [copied from cache]
- Installing RSpectra ...                       OK [copied from cache]
- Installing RcppAnnoy ...                      OK [copied from cache]
- Installing RcppArmadillo ...                  OK [copied from cache]
- Installing RcppHNSW ...                       OK [copied from cache]
- Installing RcppProgress ...                   OK [copied from cache]
- Installing RcppTOML ...                       OK [copied from cache]
- Installing Rtsne ...                          OK [copied from cache]
- Installing cowplot ...                        OK [copied from cache]
- Installing fastDummies ...                    OK [copied from cache]
- Installing survival ...                       OK [copied from cache]
- Installing fitdistrplus ...                   OK [copied from cache]
- Installing codetools ...                      OK [copied from cache]
- Installing globals ...                        OK [copied from cache]
- Installing listenv ...                        OK [copied from cache]
- Installing parallelly ...                     OK [copied from cache]
- Installing future ...                         OK [copied from cache]
- Installing future.apply ...                   OK [copied from cache]
- Installing ggrepel ...                        OK [copied from cache]
- Installing ggridges ...                       OK [copied from cache]
- Installing curl ...                           OK [copied from cache]
- Installing sys ...                            OK [copied from cache]
- Installing askpass ...                        OK [copied from cache]
- Installing openssl ...                        OK [copied from cache]
- Installing httr ...                           OK [copied from cache]
- Installing ica ...                            OK [copied from cache]
- Installing igraph ...                         OK [copied from cache]
- Installing irlba ...                          OK [copied from cache]
- Installing rprojroot ...                      OK [copied from cache]
- Installing here ...                           OK [copied from cache]
- Installing png ...                            OK [copied from cache]
- Installing reticulate ...                     OK [copied from cache]
- Installing leiden ...                         OK [copied from cache]
- Installing zoo ...                            OK [copied from cache]
- Installing lmtest ...                         OK [copied from cache]
- Installing matrixStats ...                    OK [copied from cache]
- Installing xtable ...                         OK [copied from cache]
- Installing sourcetools ...                    OK [copied from cache]
- Installing commonmark ...                     OK [copied from cache]
- Installing shiny ...                          OK [copied from cache]
- Installing miniUI ...                         OK [copied from cache]
- Installing pbapply ...                        OK [copied from cache]
- Installing plotly ...                         OK [copied from cache]
- Installing progressr ...                      OK [copied from cache]
- Installing scattermore ...                    OK [copied from cache]
- Installing reshape2 ...                       OK [copied from cache]
- Installing sctransform ...                    OK [copied from cache]
- Installing spatstat.utils ...                 OK [copied from cache]
- Installing abind ...                          OK [copied from cache]
- Installing tensor ...                         OK [copied from cache]
- Installing spatstat.sparse ...                OK [copied from cache]
- Installing goftest ...                        OK [copied from cache]
- Installing spatstat.data ...                  OK [copied from cache]
- Installing spatstat.univar ...                OK [copied from cache]
- Installing deldir ...                         OK [copied from cache]
- Installing polyclip ...                       OK [copied from cache]
- Installing spatstat.geom ...                  OK [copied from cache]
- Installing spatstat.random ...                OK [copied from cache]
- Installing spatstat.explore ...               OK [copied from cache]
- Installing sitmo ...                          OK [copied from cache]
- Installing dqrng ...                          OK [copied from cache]
- Installing uwot ...                           OK [copied from cache]
- Installing dotCall64 ...                      OK [copied from cache]
- Installing spam ...                           OK [copied from cache]
- Installing sp ...                             OK [copied from cache]
- Installing SeuratObject ...                   OK [copied from cache]
- Installing Seurat ...                         OK [copied from cache]
- Installing iterators ...                      OK [copied from cache]
- Installing foreach ...                        OK [copied from cache]
- Installing TSP ...                            OK [copied from cache]
- Installing TeachingDemos ...                  OK [copied from cache]
- Installing doParallel ...                     OK [copied from cache]
- Installing dynamicTreeCut ...                 OK [copied from cache]
- Installing fastcluster ...                    OK [copied from cache]
- Installing WGCNA ...                          OK [copied from cache]
- Installing WriteXLS ...                       OK [copied from cache]
- Installing XML ...                            OK [copied from cache]
- Installing anytime ...                        OK [copied from cache]
- Installing ape ...                            OK [copied from cache]
- Installing assertthat ...                     OK [copied from cache]
- Installing boot ...                           OK [copied from cache]
- Installing brew ...                           OK [copied from cache]
- Installing brio ...                           OK [copied from cache]
- Installing broom ...                          OK [copied from cache]
- Installing cards ...                          OK [copied from cache]
- Installing clipr ...                          OK [copied from cache]
- Installing tzdb ...                           OK [copied from cache]
- Installing vroom ...                          OK [copied from cache]
- Installing readr ...                          OK [copied from cache]
- Installing haven ...                          OK [copied from cache]
- Installing labelled ...                       OK [copied from cache]
- Installing broom.helpers ...                  OK [copied from cache]
- Installing ca ...                             OK [copied from cache]
- Installing ps ...                             OK [copied from cache]
- Installing processx ...                       OK [copied from cache]
- Installing callr ...                          OK [copied from cache]
- Installing shape ...                          OK [copied from cache]
- Installing circlize ...                       OK [copied from cache]
- Installing class ...                          OK [copied from cache]
- Installing clue ...                           OK [copied from cache]
- Installing tidygraph ...                      OK [copied from cache]
- Installing tweenr ...                         OK [copied from cache]
- Installing systemfonts ...                    OK [copied from cache]
- Installing ggforce ...                        OK [copied from cache]
- Installing graphlayouts ...                   OK [copied from cache]
- Installing ggraph ...                         OK [copied from cache]
- Installing clustree ...                       OK [copied from cache]
- Installing credentials ...                    OK [copied from cache]
- Installing triebeard ...                      OK [copied from cache]
- Installing urltools ...                       OK [copied from cache]
- Installing httpcode ...                       OK [copied from cache]
- Installing crul ...                           OK [copied from cache]
- Installing dbplyr ...                         OK [copied from cache]
- Installing dendextend ...                     OK [copied from cache]
- Installing desc ...                           OK [copied from cache]
- Installing ellipsis ...                       OK [copied from cache]
- Installing pkgbuild ...                       OK [copied from cache]
- Installing downlit ...                        OK [copied from cache]
- Installing httr2 ...                          OK [copied from cache]
- Installing textshaping ...                    OK [copied from cache]
- Installing ragg ...                           OK [built from source and cached in 1.4m]
- Installing whisker ...                        OK [built from source and cached in 0.58s]
- Installing xml2 ...                           OK [copied from cache]
- Installing pkgdown ...                        OK [built from source and cached in 2.7s]
- Installing pkgload ...                        OK [built from source and cached in 1.4s]
- Installing profvis ...                        OK [built from source and cached in 0.93s]
- Installing sessioninfo ...                    OK [built from source and cached in 0.85s]
- Installing xopen ...                          OK [built from source and cached in 0.48s]
- Installing rcmdcheck ...                      OK [built from source and cached in 1.2s]
- Installing remotes ...                        OK [built from source and cached in 1.7s]
- Installing roxygen2 ...                       OK [built from source and cached in 6.3s]
- Installing rversions ...                      OK [built from source and cached in 0.71s]
- Installing praise ...                         OK [built from source and cached in 0.47s]
- Installing diffobj ...                        OK [built from source and cached in 4.5s]
- Installing rematch2 ...                       OK [built from source and cached in 0.99s]
- Installing waldo ...                          OK [built from source and cached in 1.1s]
- Installing testthat ...                       OK [built from source and cached in 18s]
- Installing urlchecker ...                     OK [built from source and cached in 0.62s]
- Installing zip ...                            OK [built from source and cached in 5.2s]
- Installing gert ...                           OK [built from source and cached in 2.9s]
- Installing gitcreds ...                       OK [built from source and cached in 0.61s]
- Installing ini ...                            OK [built from source and cached in 0.47s]
- Installing gh ...                             OK [built from source and cached in 0.95s]
- Installing usethis ...                        OK [built from source and cached in 2.8s]
- Installing devtools ...                       OK [built from source and cached in 2.1s]
- Installing egg ...                            OK [built from source and cached in 1.4s]
- Installing enrichR ...                        OK [built from source and cached in 1.6s]
- Installing maps ...                           OK [built from source and cached in 2.4s]
- Installing fields ...                         OK [built from source and cached in 5.1s]
- Installing filelock ...                       OK [copied from cache]
- Installing flashClust ...                     OK [built from source and cached in 0.59s]
- Installing formatR ...                        OK [copied from cache]
- Installing lambda.r ...                       OK [copied from cache]
- Installing futile.options ...                 OK [copied from cache]
- Installing futile.logger ...                  OK [copied from cache]
- Installing gclus ...                          OK [built from source and cached in 0.68s]
- Installing gfonts ...                         OK [built from source and cached in 1.2s]
- Installing ggdendro ...                       OK [built from source and cached in 1.4s]
- Installing qap ...                            OK [built from source and cached in 0.61s]
- Installing registry ...                       OK [built from source and cached in 0.54s]
- Installing permute ...                        OK [built from source and cached in 0.95s]
- Installing vegan ...                          OK [built from source and cached in 12s]
- Installing seriation ...                      OK [built from source and cached in 4.4s]
- Installing webshot ...                        OK [built from source and cached in 0.61s]
- Installing heatmaply ...                      OK [built from source and cached in 2.1s]
- Installing hwriter ...                        OK [copied from cache]
- Installing interp ...                         OK [copied from cache]
- Installing jpeg ...                           OK [copied from cache]
- Installing kmed ...                           OK [built from source and cached in 1.5s]
- Installing latticeExtra ...                   OK [copied from cache]
- Installing locfit ...                         OK [copied from cache]
- Installing markdown ...                       OK [built from source and cached in 0.89s]
- Installing mnormt ...                         OK [built from source and cached in 1.7s]
- Installing openxlsx ...                       OK [built from source and cached in 32s]
- Installing pheatmap ...                       OK [built from source and cached in 1.1s]
- Installing plotrix ...                        OK [built from source and cached in 4.7s]
- Installing GPArotation ...                    OK [built from source and cached in 0.8s]
- Installing psych ...                          OK [built from source and cached in 21s]
- Installing renv ...                           OK [built from source and cached in 6.7s]
- Installing reshape ...                        OK [built from source and cached in 1.1s]
- Installing restfulr ...                       OK [copied from cache]
- Installing selectr ...                        OK [built from source and cached in 5.3s]
- Installing rvest ...                          OK [built from source and cached in 2.8s]
- Installing shinyBS ...                        OK [built from source and cached in 1.5s]
- Installing shinyWidgets ...                   OK [built from source and cached in 6.4s]
- Installing shinycssloaders ...                OK [built from source and cached in 0.91s]
- Installing shinyjs ...                        OK [built from source and cached in 1.6s]
- Installing shinythemes ...                    OK [built from source and cached in 1.8s]
- Installing sysfonts ...                       OK [built from source and cached in 1.7s]
- Installing showtextdb ...                     OK [built from source and cached in 0.98s]
- Installing showtext ...                       OK [built from source and cached in 2.0s]
- Installing snow ...                           OK [copied from cache]
- Installing spatial ...                        OK [built from source and cached in 3.5s]
- Installing statmod ...                        OK [copied from cache]
- Installing umap ...                           OK [built from source and cached in 9.4s]
Successfully installed 311 packages in 4.9 minutes.
- GitHub authentication credentials are not available.
- Please set GITHUB_PAT, or ensure the 'gitcreds' package is installed.
- See https://usethis.r-lib.org/articles/git-credentials.html for more details.
# Downloading packages -------------------------------------------------------
- Downloading presto from GitHub ...            OK [662.3 Kb in 0.4s]
- Downloading scClustViz from GitHub ...        OK [94.4 Kb in 0.44s]
Successfully downloaded 2 packages in 1.5 seconds.

The following package(s) will be installed:
- presto     [immunogenomics/presto]
- scClustViz [BaderLab/scClustViz]
These packages will be installed into "/usr/local/lib/R/site-library".

Do you want to proceed? [Y/n]: Y

# Installing packages --------------------------------------------------------
- Installing presto ...                         OK [built from source and cached in 15s]
- Installing scClustViz ...                     OK [built from source and cached in 3.2s]
Successfully installed 2 packages in 18 seconds.
The following package(s) will be updated in the lockfile:

# Bioconductor ---------------------------------------------------------------
- GenomeInfoDbData       [1.2.11 -> 1.2.12]
- GO.db                  [3.18.0 -> 3.19.1]
- org.Hs.eg.db           [3.18.0 -> 3.19.1]
- org.Mm.eg.db           [3.18.0 -> 3.19.1]

# Bioconductor 3.19 ----------------------------------------------------------
- annotate               [repo: * -> Bioconductor 3.19; ver: 1.80.0 -> 1.82.0]
- AnnotationDbi          [repo: * -> Bioconductor 3.19; ver: 1.64.1 -> 1.66.0]
- aroma.light            [repo: * -> Bioconductor 3.19; ver: 3.32.0 -> 3.34.0]
- Biobase                [repo: * -> Bioconductor 3.19; ver: 2.62.0 -> 2.64.0]
- BiocFileCache          [repo: * -> Bioconductor 3.19; ver: 2.10.1 -> 2.12.0]
- BiocGenerics           [repo: * -> Bioconductor 3.19; ver: 0.48.1 -> 0.50.0]
- BiocIO                 [repo: * -> Bioconductor 3.19; ver: 1.12.0 -> 1.14.0]
- BiocParallel           [repo: * -> Bioconductor 3.19; ver: 1.36.0 -> 1.38.0]
- BiocVersion            [repo: * -> Bioconductor 3.19; ver: 3.18.1 -> 3.19.1]
- biomaRt                [repo: Bioconductor 3.18 -> Bioconductor 3.19; ver: 2.58.1 -> 2.60.1]
- Biostrings             [repo: Bioconductor 3.18 -> Bioconductor 3.19; ver: 2.70.2 -> 2.72.1]
- ComplexHeatmap         [repo: * -> Bioconductor 3.19; ver: 2.18.0 -> 2.20.0]
- DelayedArray           [repo: * -> Bioconductor 3.19; ver: 0.28.0 -> 0.30.1]
- DESeq2                 [repo: * -> Bioconductor 3.19; ver: 1.42.0 -> 1.44.0]
- EDASeq                 [repo: * -> Bioconductor 3.19; ver: 2.36.0 -> 2.38.0]
- edgeR                  [repo: Bioconductor 3.18 -> Bioconductor 3.19; ver: 4.0.12 -> 4.2.2]
- geneplotter            [repo: * -> Bioconductor 3.19; ver: 1.80.0 -> 1.82.0]
- GenomeInfoDb           [repo: Bioconductor 3.18 -> Bioconductor 3.19; ver: 1.38.5 -> 1.40.1]
- GenomicAlignments      [repo: Bioconductor 3.18 -> Bioconductor 3.19; ver: 1.38.2 -> 1.40.0]
- GenomicFeatures        [repo: * -> Bioconductor 3.19; ver: 1.54.1 -> 1.56.0]
- GenomicRanges          [repo: * -> Bioconductor 3.19; ver: 1.54.1 -> 1.56.2]
- impute                 [repo: * -> Bioconductor 3.19; ver: 1.76.0 -> 1.78.0]
- IRanges                [repo: * -> Bioconductor 3.19; ver: 2.36.0 -> 2.38.1]
- KEGGREST               [repo: * -> Bioconductor 3.19; ver: 1.42.0 -> 1.44.1]
- limma                  [repo: * -> Bioconductor 3.19; ver: 3.58.1 -> 3.60.6]
- MAST                   [repo: * -> Bioconductor 3.19; ver: 1.28.0 -> 1.30.0]
- MatrixGenerics         [repo: * -> Bioconductor 3.19; ver: 1.14.0 -> 1.16.0]
- preprocessCore         [repo: * -> Bioconductor 3.19; ver: 1.64.0 -> 1.66.0]
- Rhtslib                [repo: Bioconductor 3.18 -> Bioconductor 3.19; ver: 2.4.1 -> 3.0.0]
- Rsamtools              [repo: * -> Bioconductor 3.19; ver: 2.18.0 -> 2.20.0]
- rtracklayer            [repo: * -> Bioconductor 3.19; ver: 1.62.0 -> 1.64.0]
- RUVSeq                 [repo: * -> Bioconductor 3.19; ver: 1.36.0 -> 1.38.0]
- S4Arrays               [repo: * -> Bioconductor 3.19; ver: 1.2.0 -> 1.4.1]
- S4Vectors              [repo: * -> Bioconductor 3.19; ver: 0.40.2 -> 0.42.1]
- ShortRead              [repo: * -> Bioconductor 3.19; ver: 1.60.0 -> 1.62.0]
- SingleCellExperiment   [repo: * -> Bioconductor 3.19; ver: 1.24.0 -> 1.26.0]
- SparseArray            [repo: Bioconductor 3.18 -> Bioconductor 3.19; ver: 1.2.3 -> 1.4.8]
- SummarizedExperiment   [repo: * -> Bioconductor 3.19; ver: 1.32.0 -> 1.34.0]
- XVector                [repo: * -> Bioconductor 3.19; ver: 0.42.0 -> 0.44.0]
- zlibbioc               [repo: * -> Bioconductor 3.19; ver: 1.48.0 -> 1.50.0]
- pwalign                [* -> 1.0.0]
- UCSC.utils             [* -> 1.0.0]

# CRAN -----------------------------------------------------------------------
- abind                  [repo: RSPM -> CRAN; ver: 1.4-5 -> 1.4-8]
- ape                    [5.7-1 -> 5.8]
- askpass                [1.2.0 -> 1.2.1]
- backports              [1.4.1 -> 1.5.0]
- BiocManager            [1.30.22 -> 1.30.25]
- bit                    [4.0.5 -> 4.5.0]
- bit64                  [4.0.5 -> 4.5.2]
- bitops                 [repo: RSPM -> CRAN; ver: 1.0-7 -> 1.0-9]
- boot                   [1.3-28.1 -> 1.3-31]
- brio                   [1.1.4 -> 1.1.5]
- broom                  [1.0.5 -> 1.0.7]
- broom.helpers          [1.14.0 -> 1.17.0]
- bslib                  [0.6.1 -> 0.8.0]
- cachem                 [1.0.8 -> 1.1.0]
- callr                  [3.7.3 -> 3.7.6]
- caTools                [1.18.2 -> 1.18.3]
- checkmate              [2.3.1 -> 2.3.2]
- circlize               [repo: RSPM -> CRAN; ver: 0.4.15 -> 0.4.16]
- cli                    [3.6.2 -> 3.6.3]
- codetools              [0.2-19 -> 0.2-20]
- colorspace             [2.1-0 -> 2.1-1]
- commonmark             [1.9.0 -> 1.9.2]
- cpp11                  [0.4.7 -> 0.5.0]
- crayon                 [1.5.2 -> 1.5.3]
- credentials            [repo: RSPM -> CRAN; ver: 2.0.1 -> 2.0.2]
- crul                   [1.4.0 -> 1.5.0]
- curl                   [5.2.0 -> 5.2.3]
- data.table             [1.14.10 -> 1.16.2]
- DBI                    [1.2.1 -> 1.2.3]
- dbplyr                 [2.4.0 -> 2.5.0]
- deldir                 [2.0-2 -> 2.0-4]
- dendextend             [1.17.1 -> 1.18.1]
- digest                 [0.6.34 -> 0.6.37]
- dotCall64              [1.1-1 -> 1.2]
- downlit                [repo: RSPM -> CRAN; ver: 0.4.3 -> 0.4.4]
- dqrng                  [0.3.2 -> 0.4.1]
- DT                     [0.31 -> 0.33]
- evaluate               [0.23 -> 1.0.1]
- expm                   [0.999-9 -> 1.0-0]
- farver                 [2.1.1 -> 2.1.2]
- fastDummies            [1.7.3 -> 1.7.4]
- fastmap                [1.1.1 -> 1.2.0]
- fields                 [15.2 -> 16.3]
- fitdistrplus           [1.1-11 -> 1.2-1]
- FNN                    [1.1.4 -> 1.1.4.1]
- foreign                [0.8-86 -> 0.8-87]
- fs                     [1.6.3 -> 1.6.4]
- future                 [1.33.1 -> 1.34.0]
- future.apply           [1.11.1 -> 1.11.2]
- gert                   [2.0.1 -> 2.1.4]
- GGally                 [2.2.0 -> 2.2.1]
- ggdendro               [0.1.23 -> 0.2.0]
- ggforce                [0.4.1 -> 0.4.2]
- ggplot2                [3.4.4 -> 3.5.1]
- ggraph                 [2.1.0 -> 2.2.1]
- ggrepel                [0.9.5 -> 0.9.6]
- ggstats                [0.5.1 -> 0.7.0]
- gh                     [repo: RSPM -> CRAN; ver: 1.4.0 -> 1.4.1]
- globals                [0.16.2 -> 0.16.3]
- glue                   [1.7.0 -> 1.8.0]
- gplots                 [3.1.3 -> 3.2.0]
- graphlayouts           [repo: RSPM -> CRAN; ver: 1.1.0 -> 1.2.0]
- gtable                 [0.3.4 -> 0.3.5]
- highr                  [0.10 -> 0.11]
- Hmisc                  [repo: RSPM -> CRAN; ver: 5.1-1 -> 5.1-3]
- htmlTable              [2.4.2 -> 2.4.3]
- htmltools              [0.5.7 -> 0.5.8.1]
- httpuv                 [repo: RSPM -> CRAN; ver: 1.6.14 -> 1.6.15]
- httr2                  [1.0.0 -> 1.0.5]
- igraph                 [1.6.0 -> 2.1.1]
- jsonlite               [1.8.8 -> 1.8.9]
- KernSmooth             [repo: RSPM -> CRAN; ver: 2.23-22 -> 2.23-24]
- knitr                  [1.45 -> 1.48]
- labelled               [2.12.0 -> 2.13.0]
- lattice                [0.22-5 -> 0.22-6]
- listenv                [0.9.0 -> 0.9.1]
- locfit                 [1.5-9.8 -> 1.5-9.10]
- markdown               [1.12 -> 1.13]
- MASS                   [7.3-60.0.1 -> 7.3-61]
- Matrix                 [1.6-5 -> 1.7-1]
- matrixStats            [1.2.0 -> 1.4.1]
- munsell                [0.5.0 -> 0.5.1]
- nlme                   [3.1-164 -> 3.1-166]
- openssl                [2.1.1 -> 2.2.2]
- openxlsx               [4.2.5.2 -> 4.2.7.1]
- parallelly             [1.36.0 -> 1.38.0]
- patchwork              [1.2.0 -> 1.3.0]
- pkgbuild               [1.4.3 -> 1.4.4]
- pkgdown                [repo: RSPM -> CRAN; ver: 2.0.7 -> 2.1.1]
- pkgload                [repo: RSPM -> CRAN; ver: 1.3.4 -> 1.4.0]
- polyclip               [1.10-6 -> 1.10-7]
- processx               [3.8.3 -> 3.8.4]
- profvis                [repo: RSPM -> CRAN; ver: 0.3.8 -> 0.4.0]
- promises               [1.2.1 -> 1.3.0]
- ps                     [repo: RSPM -> CRAN; ver: 1.7.6 -> 1.8.0]
- psych                  [repo: RSPM -> CRAN; ver: 2.4.1 -> 2.4.6.26]
- ragg                   [1.2.7 -> 1.3.3]
- RANN                   [2.6.1 -> 2.6.2]
- Rcpp                   [1.0.12 -> 1.0.13]
- RcppArmadillo          [0.12.6.6.1 -> 14.0.2-1]
- RcppEigen              [0.3.3.9.4 -> 0.3.4.0.2]
- RcppHNSW               [0.5.0 -> 0.6.0]
- RCurl                  [1.98-1.14 -> 1.98-1.16]
- remotes                [repo: RSPM -> CRAN; ver: 2.4.2.1 -> 2.5.0]
- renv                   [1.0.3 -> 1.0.11]
- reticulate             [1.34.0 -> 1.39.0]
- rjson                  [0.2.21 -> 0.2.23]
- rlang                  [1.1.3 -> 1.1.4]
- RMariaDB               [1.3.1 -> 1.3.2]
- rmarkdown              [2.25 -> 2.28]
- roxygen2               [repo: RSPM -> CRAN; ver: 7.3.1 -> 7.3.2]
- RSpectra               [0.16-1 -> 0.16-2]
- RSQLite                [repo: RSPM -> CRAN; ver: 2.3.5 -> 2.3.7]
- rstudioapi             [0.15.0 -> 0.17.0]
- rvest                  [1.0.3 -> 1.0.4]
- sass                   [0.4.8 -> 0.4.9]
- seriation              [1.5.4 -> 1.5.6]
- Seurat                 [5.0.1 -> 5.1.0]
- SeuratObject           [5.0.1 -> 5.0.2]
- shape                  [repo: RSPM -> CRAN; ver: 1.4.6 -> 1.4.6.1]
- shiny                  [1.8.0 -> 1.9.1]
- shinycssloaders        [1.0.0 -> 1.1.0]
- shinyWidgets           [0.8.1 -> 0.8.7]
- showtext               [0.9-6 -> 0.9-7]
- sp                     [2.1-2 -> 2.1-4]
- spam                   [2.10-0 -> 2.11-0]
- spatstat.data          [repo: RSPM -> CRAN; ver: 3.0-4 -> 3.1-2]
- spatstat.explore       [3.2-5 -> 3.3-2]
- spatstat.geom          [repo: RSPM -> CRAN; ver: 3.2-8 -> 3.3-3]
- spatstat.random        [3.2-2 -> 3.3-2]
- spatstat.sparse        [3.0-3 -> 3.1-0]
- spatstat.utils         [3.0-4 -> 3.1-0]
- stringi                [1.8.3 -> 1.8.4]
- survival               [repo: RSPM -> CRAN; ver: 3.5-7 -> 3.7-0]
- sys                    [3.4.2 -> 3.4.3]
- sysfonts               [0.8.8 -> 0.8.9]
- systemfonts            [1.0.5 -> 1.1.0]
- TeachingDemos          [repo: RSPM -> CRAN; ver: 2.12 -> 2.13]
- testthat               [3.2.1 -> 3.2.1.1]
- textshaping            [0.3.7 -> 0.4.0]
- tidygraph              [1.3.0 -> 1.3.1]
- tidyselect             [1.2.0 -> 1.2.1]
- tinytex                [0.49 -> 0.53]
- tweenr                 [2.0.2 -> 2.0.3]
- usethis                [repo: RSPM -> CRAN; ver: 2.2.2 -> 3.0.0]
- uwot                   [0.1.16 -> 0.2.2]
- vegan                  [2.6-4 -> 2.6-8]
- viridis                [0.6.4 -> 0.6.5]
- waldo                  [0.5.2 -> 0.5.3]
- WGCNA                  [1.72-5 -> 1.73]
- withr                  [repo: RSPM -> CRAN; ver: 3.0.0 -> 3.0.1]
- WriteXLS               [6.5.0 -> 6.7.0]
- xfun                   [0.41 -> 0.48]
- XML                    [repo: RSPM -> CRAN; ver: 3.99-0.16.1 -> 3.99-0.17]
- xopen                  [repo: RSPM -> CRAN; ver: 1.0.0 -> 1.0.1]
- yaml                   [2.3.8 -> 2.3.10]
- bspm                   [* -> 0.5.7]
- cards                  [* -> 0.3.0]
- GPArotation            [* -> 2024.3-1]
- spatstat.univar        [* -> 3.0-1]

# GitHub ---------------------------------------------------------------------
- presto                 [immunogenomics/presto: 31dc97fe -> 7636b3d0]

The version of R recorded in the lockfile will be updated:
- R                      [4.3.2 -> 4.4.1]

- Lockfile written to "/home/weismanal/notebook/2024-10-18/testing_sequin_installation_on_linux/public_sequin/renv.lock".
Warning messages:
1: In value[[3L]](cond) : Could not find any credentials
2: In value[[3L]](cond) : Could not find any credentials
3: In value[[3L]](cond) : Could not find any credentials
```

Restart R (note, without `sudo`) and open the app via:

```console
$ R
> setwd("/home/weismanal/notebook/2024-10-18/testing_sequin_installation_on_linux/public_sequin")
> shiny::runApp(launch.browser = T)
```

to eventually get in the console (upon opening the `example_sc` dataset in the browser):

```console
Warning: The 'plotly_click' event tied a source ID of 'corplot' is not registered. In order to obtain this event data, please add `event_register(p, 'plotly_click')` to the plot (`p`) that you wish to obtain event data from.
Warning: The 'plotly_click' event tied a source ID of 'A' is not registered. In order to obtain this event data, please add `event_register(p, 'plotly_click')` to the plot (`p`) that you wish to obtain event data from.
Warning: The 'plotly_click' event tied a source ID of 'A' is not registered. In order to obtain this event data, please add `event_register(p, 'plotly_click')` to the plot (`p`) that you wish to obtain event data from.
Warning: Data is of class data.frame. Coercing to dgCMatrix.
Warning: The `slot` argument of `FetchData()` is deprecated as of SeuratObject 5.0.0.
ℹ Please use the `layer` argument instead.
This warning is displayed once every 8 hours.
Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
converting counts to integer mode
Warning in irlba(A = t(x = object), nv = npcs, ...) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Warning: The `slot` argument of `GetAssayData()` is deprecated as of SeuratObject 5.0.0.
ℹ Please use the `layer` argument instead.
This warning is displayed once every 8 hours.
Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
This message will be shown once per session
Computing nearest neighbor graph
Warning: Error in getGlobalsAndPackages: The total size of the 7 globals exported for future expression (‘FUN()’) is 685.56 MiB.. This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize'). The three largest globals are ‘FUN’ (342.85 MiB of class ‘function’), ‘index’ (342.68 MiB of class ‘S4’) and ‘query’ (28.44 KiB of class ‘numeric’)
  1: shiny::runApp [/tmp/RtmpNy79O8/renv-package-new-7c4e69fc5b25/shiny/R/runapp.R#388]
```

after which the app dies.

`Ctrl-C` back to the R prompt, update the total size for exported global variables to 2 GB via:

```R
options(future.globals.maxSize = 2000 * 1024^2)
```

and run the app again:

```R
shiny::runApp(launch.browser = T)
```

Using the app seems to work in all naively tested ways, aside from warnings:

```console
Warning: The 'plotly_click' event tied a source ID of 'corplot' is not registered. In order to obtain this event data, please add `event_register(p, 'plotly_click')` to the plot (`p`) that you wish to obtain event data from.
Warning: The 'plotly_click' event tied a source ID of 'A' is not registered. In order to obtain this event data, please add `event_register(p, 'plotly_click')` to the plot (`p`) that you wish to obtain event data from.
Warning: The 'plotly_click' event tied a source ID of 'A' is not registered. In order to obtain this event data, please add `event_register(p, 'plotly_click')` to the plot (`p`) that you wish to obtain event data from.
Counts matrix provided is not sparse; vreating v5 assay in Seurat object
Warning: Data is of class data.frame. Coercing to dgCMatrix.
converting counts to integer mode
Warning in irlba(A = t(x = object), nv = npcs, ...) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Computing nearest neighbor graph
Computing SNN

--------------------------------------
--------------------------------------
Processing cluster solution: Clust
--------------------------------------
Loading required package: cluster
-- Calculating cluster-wise gene statistics --
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=00s
-- Calculating differential expression cluster vs rest effect size --
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=00s
-- Testing differential expression cluster vs rest --
-- Testing differential expression between clusters --
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=00s
`summarise()` has grouped output by 'orig.ident'. You can override using the `.groups` argument.
Joining with `by = join_by(orig.ident, seurat_clusters)`
`summarise()` has grouped output by 'key'. You can override using the `.groups` argument.
`summarise()` has grouped output by 'key'. You can override using the `.groups` argument.
$rect
$rect$w
[1] 0.6958997

$rect$h
[1] 47.42069

$rect$left
[1] 6.8441

$rect$top
[1] 435.2166


$text
$text$x
[1] 6.94156

$text$y
[1] 411.5062


Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
ℹ Please use the `linewidth` argument instead.
This warning is displayed once every 8 hours.
Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
ℹ Please use the `linewidth` argument instead.
This warning is displayed once every 8 hours.
Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
Warning: The select input "cgGene" contains a large number of options; consider using server-side selectize for massively improved performance. See the Details section of the ?selectizeInput help topic.
Warning: The `trans` argument of `sec_axis()` is deprecated as of ggplot2 3.5.0.
ℹ Please use the `transform` argument instead.
This warning is displayed once every 8 hours.
Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
Warning: The select input "goi1" contains a large number of options; consider using server-side selectize for massively improved performance. See the Details section of the ?selectizeInput help topic.
Scale for colour is already present.
Adding another scale for colour, which will replace the existing scale.
Warning: Vectorized input to `element_text()` is not officially supported.
ℹ Results may be unexpected or may change in future versions of ggplot2.
Uploading data to Enrichr... Done.
  Querying Genome_Browser_PWMs... Done.
  Querying TRANSFAC_and_JASPAR_PWMs... Done.
  Querying Transcription_Factor_PPIs... Done.
  Querying ChEA_2013... Done.
  Querying Drug_Perturbations_from_GEO_2014... Done.
  Querying ENCODE_TF_ChIP-seq_2014... Done.
Parsing results... Done.
Allowing parallel execution with up to 10 working processes.
 Flagging genes and samples with too many missing values...
  ..step 1
..connectivity..
..matrix multiplication (system BLAS)..
..normalization..
..done.
 ..cutHeight not given, setting it to 0.997  ===>  99% of the (truncated) height range in dendro.
 ..done.
 mergeCloseModules: Merging modules whose distance is less than 0
   multiSetMEs: Calculating module MEs.
     Working on set 1 ...
     moduleEigengenes: Calculating 4 module eigengenes in given set.
   Calculating new MEs...
   multiSetMEs: Calculating module MEs.
     Working on set 1 ...
     moduleEigengenes: Calculating 4 module eigengenes in given set.
Scale for y is already present.
Adding another scale for y, which will replace the existing scale.
Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
ℹ Please use `linewidth` instead.
This warning is displayed once every 8 hours.
Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family 'noto-sans-jp' not found in PostScript font database
Warning: The select input "gene" contains a large number of options; consider using server-side selectize for massively improved performance. See the Details section of the ?selectizeInput help topic.
Warning: The `<scale>` argument of `guides()` cannot be `FALSE`. Use "none" instead as of ggplot2 3.3.4.
This warning is displayed once every 8 hours.
Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
Scale for fill is already present.
Adding another scale for fill, which will replace the existing scale.
Warning: Removed 123 rows containing missing values or values outside the scale range (`geom_point()`).
Scale for fill is already present.
Adding another scale for fill, which will replace the existing scale.
Warning: Removed 123 rows containing missing values or values outside the scale range (`geom_point()`).
Scale for fill is already present.
Adding another scale for fill, which will replace the existing scale.
Warning: Removed 1230 rows containing missing values or values outside the scale range (`geom_point()`).
Scale for colour is already present.
Adding another scale for colour, which will replace the existing scale.
Scale for colour is already present.
Adding another scale for colour, which will replace the existing scale.
Scale for colour is already present.
Adding another scale for colour, which will replace the existing scale.
Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
Also defined by ‘spam’
Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
Also defined by ‘spam’
using ntop=500 top features by variance
using ntop=500 top features by variance
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
Warning in melt.default(al.df) :
  The melt generic in data.table has been passed a data.frame and will attempt to redirect to the relevant reshape2 method; please note that reshape2 is superseded and is no longer actively developed, and this redirection is now deprecated. To continue using melt methods from reshape2 while both libraries are attached, e.g. melt.list, you can prepend the namespace, i.e. reshape2::melt(al.df). In the next version, this warning will become an error.
Using contrast as id variables
Uploading data to Enrichr... Done.
  Querying Genome_Browser_PWMs... Done.
  Querying TRANSFAC_and_JASPAR_PWMs... Done.
  Querying Transcription_Factor_PPIs... Done.
  Querying ChEA_2013... Done.
  Querying Drug_Perturbations_from_GEO_2014... Done.
  Querying ENCODE_TF_ChIP-seq_2014... Done.
Parsing results... Done.
Allowing parallel execution with up to 10 working processes.
..connectivity..
..matrix multiplication (system BLAS)..
..normalization..
..done.
 ..cutHeight not given, setting it to 0.996  ===>  99% of the (truncated) height range in dendro.
 ..done.
 mergeCloseModules: Merging modules whose distance is less than 0
   multiSetMEs: Calculating module MEs.
     Working on set 1 ...
     moduleEigengenes: Calculating 8 module eigengenes in given set.
   Calculating new MEs...
   multiSetMEs: Calculating module MEs.
     Working on set 1 ...
     moduleEigengenes: Calculating 8 module eigengenes in given set.
Scale for y is already present.
Adding another scale for y, which will replace the existing scale.
Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family 'noto-sans-jp' not found in PostScript font database
Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family 'noto-sans-jp' not found in PostScript font database
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |======================================================================| 100%
Warning in irlba(A = t(x = object), nv = npcs, ...) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
PC_ 1
Positive:  GABARAPL1, IFIT1, NTN4, ARID5B, SIPA1L2, AKR1C2, ADGRL2, FIBIN, PAMR1, ADM
           LRRK2, FBXO21, OLR1, C1RL, CPM, DCN, C1R, GBP2, PLXNC1, C1S
           ADAMTS8, MIPEP, VCAM1, TRIM22, PIK3R3, CD46, TSPAN9, APLP2, RHOBTB1, OPTN
Negative:  HMGA2, DDX21, S100A16, DKK1, CCT2, FEN1, EBNA1BP2, BCAR3, PPP1R14B, CTPS1
           TMPO, PA2G4, CLSPN, SKA3, KIF11, MRTO4, KNTC1, USP1, IPO7, NAV3
           MKI67, DTL, NASP, UTP20, GTPBP4, SRM, SERBP1, DDIAS, CDK1, KIF2C
PC_ 2
Positive:  ATP2B4, EHD1, PRDX1, ENDOD1, VCL, KLF6, THY1, A2M, FAM171A1, LIMA1
           MGST1, TALDO1, CYB5R1, MYPN, SPOCD1, TXNRD1, RRAS2, STK40, HERC4, KITLG
           PCSK7, STX12, CCND1, TMEM109, DDAH1, SCN4B, TXNIP, MYCBP, NMT2, IGFBP6
Negative:  OLFML3, CTSK, DDIT4, PLXDC2, TENM4, TSKU, PCSK9, CTSS, DCHS1, COL11A1
           DDIT3, TBX3, HSD17B7, SC5D, TMEM41B, FADS1, HYOU1, HTRA1, PGM1, CEP57
           SV2A, MKX, SYT7, TMEM132B, IDI1, HSP90B1, PTPRU, NTM, ELMOD1, KIF26B
PC_ 3
Positive:  CCDC85B, GJA3, DUSP6, AVPI1, PC, STAMBPL1, ATP1A1, OAF, PPT1, EPS8
           ZC3H12C, PHYH, ORAI1, SLC6A15, GALE, SOCS2, ST3GAL4, TNFRSF1B, SLC43A3, IGSF8
           SMAP2, PPP2R1B, SSH1, BTG1, CERS2, GLMP, EMP1, NR4A1, RUSC1, VPS51
Negative:  LOXL4, PRICKLE1, TUFT1, MSRB3, ACTA2, HHAT, KRT18, TMEM263, ITGB1, RSU1
           P4HA1, EXTL1, GBP1, DNAJB4, FRS2, ANKRD1, TRAF5, NUAK1, CNN3, ANO6
           TAGLN, EEA1, PAWR, ETV6, TSPAN2, NEK7, SGPL1, HSD17B6, GXYLT1, FGD6
PC_ 4
Positive:  EDEM3, PGM2L1, TMTC3, MALAT1, MEF2D, CNNM2, NR4A1, MMP3, ID3, KRAS
           HSPA6, MMP10, SSH1, RPE65, FRMD8, WLS, SLC38A2, FAM102B, IKBIP, HEYL
           LAMC2, APAF1, CD53, ZDHHC20, ZMYM1, GREM2, IL24, CDC42SE1, P4HA3, ZBED6
Negative:  SYNC, MVK, SH2D5, TMEM52, ITPKB, ADGRB2, HPS1, FDPS, TM7SF2, LRRN4CL
           HMGCL, GSTM2, KCNJ1, ARHGAP29, FADS2, NTPCR, SIRT3, MDK, NAV2, TUFT1
           DHCR7, HRK, BATF2, VEGFB, CDK4, DUSP5, SDSL, COMMD3, ZNF385A, EFEMP2
PC_ 5
Positive:  GUCY1A2, LCE3A, GLYATL2, C1QL4, LCE2A, OLFML2B, CD9, SSX2IP, ANKRD22, ITGA10
           ABCD2, NELL2, LDHC, LCE3D, FOXR1, SLC38A2, C2CD4D, FMN2, DEUP1, COL2A1
           KIF5A, SLC22A9, SH2D2A, MMP20, HES2, KRT72, PRPH, MTNR1B, DISP3, RGS16
Negative:  CD82, HSD11B1, FTH1, TMEM132A, FMOD, EHF, BIRC3, SAA2, PLAU, IKBKE
           CLMP, NFKB2, RGS5, DKK3, ELF3, TPCN1, CHI3L2, CCNL2, SLC11A2, DDX10
           ITIH5, CH25H, PDE4DIP, G0S2, CD69, CREG1, YIF1A, TNFSF18, ZC3H12A, IRAK3
Scale for fill is already present.
Adding another scale for fill, which will replace the existing scale.
Warning: Removed 117 rows containing missing values or values outside the scale range (`geom_point()`).
Scale for fill is already present.
Adding another scale for fill, which will replace the existing scale.
Warning: Removed 117 rows containing missing values or values outside the scale range (`geom_point()`).
Scale for fill is already present.
Adding another scale for fill, which will replace the existing scale.
Warning: Removed 117 rows containing missing values or values outside the scale range (`geom_point()`).
Scale for fill is already present.
Adding another scale for fill, which will replace the existing scale.
Warning: Removed 936 rows containing missing values or values outside the scale range (`geom_point()`).
Scale for fill is already present.
Adding another scale for fill, which will replace the existing scale.
Warning: Removed 936 rows containing missing values or values outside the scale range (`geom_point()`).
Scale for fill is already present.
Adding another scale for fill, which will replace the existing scale.
Warning: Removed 468 rows containing missing values or values outside the scale range (`geom_point()`).
Scale for colour is already present.
Adding another scale for colour, which will replace the existing scale.
Scale for colour is already present.
Adding another scale for colour, which will replace the existing scale.
Scale for colour is already present.
Adding another scale for colour, which will replace the existing scale.
Scale for colour is already present.
Adding another scale for colour, which will replace the existing scale.
```

As an FYI for future app development, note in particular this warning:

```console
Warning in melt.default(al.df) :
  The melt generic in data.table has been passed a data.frame and will attempt to redirect to the relevant reshape2 method; please note that reshape2 is superseded and is no longer actively developed, and this redirection is now deprecated. To continue using melt methods from reshape2 while both libraries are attached, e.g. melt.list, you can prepend the namespace, i.e. reshape2::melt(al.df). In the next version, this warning will become an error.
```

Even though the app appears to work, note that `error_log.txt` does get created with the content:

```text
2024-10-21 15:58:33.534482  Error in  SC-DGE-GSE_Heatmap observe : Error:

2024-10-21 15:58:54.416241  Error in  SC-DGE-GSE_Heatmap observe : Error:

2024-10-21 16:12:32.617576  Error in  SC-DGE-GSE_Heatmap observe : Error:

2024-10-21 16:13:11.816594  Error in  SC-DGE-GSE_Heatmap observe : Error:

```

**@Andrei, please investigate, I think this is something you had worked on previously. Maybe it's something we can ignore? Thanks!**

When I click on More > Session info in the app, I get:

```text
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

time zone: America/New_York
tzcode source: system (glibc)

attached base packages:
 [1] parallel  grid      tools     stats4    stats     graphics  grDevices
 [8] datasets  utils     methods   base     

other attached packages:
  [1] cluster_2.1.6               MAST_1.30.0                
  [3] SingleCellExperiment_1.26.0 patchwork_1.3.0            
  [5] Cairo_1.6-2                 gfonts_0.2.0               
  [7] showtext_0.9-7              showtextdb_3.0             
  [9] sysfonts_0.8.9              clustree_0.5.1             
 [11] ggraph_2.2.1                shinyWidgets_0.8.7         
 [13] gridExtra_2.3               dendextend_1.18.1          
 [15] cowplot_1.1.3               ggdendro_0.2.0             
 [17] plotrix_3.8-4               ggrepel_0.9.6              
 [19] rvest_1.0.4                 xml2_1.3.6                 
 [21] tidyr_1.3.1                 purrr_1.0.2                
 [23] ComplexHeatmap_2.20.0       presto_1.0.0               
 [25] stringi_1.8.4               umap_0.2.10.0              
 [27] shinyjs_2.1.0               org.Mm.eg.db_3.19.1        
 [29] org.Hs.eg.db_3.19.1         biomaRt_2.60.1             
 [31] RUVSeq_1.38.0               EDASeq_2.38.0              
 [33] ShortRead_1.62.0            GenomicAlignments_1.40.0   
 [35] Rsamtools_2.20.0            Biostrings_2.72.1          
 [37] XVector_0.44.0              BiocParallel_1.38.0        
 [39] Seurat_5.1.0                SeuratObject_5.0.2         
 [41] sp_2.1-4                    heatmaply_1.5.0            
 [43] scales_1.3.0                viridis_0.6.5              
 [45] data.table_1.16.2           scClustViz_1.3.12          
 [47] GO.db_3.19.1                enrichR_3.2                
 [49] RColorBrewer_1.1-3          gplots_3.2.0               
 [51] rmarkdown_2.28              readr_2.1.5                
 [53] ape_5.8                     kmed_0.4.2                 
 [55] MCL_1.0                     flashClust_1.01-2          
 [57] WGCNA_1.73                  fastcluster_1.2.6          
 [59] dynamicTreeCut_1.63-1       Rtsne_0.17                 
 [61] openxlsx_4.2.7.1            stringr_1.5.1              
 [63] edgeR_4.2.2                 limma_3.60.6               
 [65] DESeq2_1.44.0               SummarizedExperiment_1.34.0
 [67] MatrixGenerics_1.16.0       matrixStats_1.4.1          
 [69] GenomicRanges_1.56.2        GenomeInfoDb_1.40.1        
 [71] psych_2.4.6.26              fields_16.3                
 [73] viridisLite_0.4.2           spam_2.11-0                
 [75] digest_0.6.37               backports_1.5.0            
 [77] reshape2_1.4.4              pheatmap_1.0.12            
 [79] GGally_2.2.1                locfit_1.5-9.10            
 [81] geneplotter_1.82.0          annotate_1.82.0            
 [83] XML_3.99-0.17               AnnotationDbi_1.66.0       
 [85] IRanges_2.38.1              S4Vectors_0.42.1           
 [87] lattice_0.22-6              Biobase_2.64.0             
 [89] BiocGenerics_0.50.0         Hmisc_5.1-3                
 [91] Rcpp_1.0.13                 tibble_3.2.1               
 [93] shinythemes_1.2.0           shinycssloaders_1.1.0      
 [95] shinyBS_0.61.1              plotly_4.10.4              
 [97] ggplot2_3.5.1               gtools_3.9.5               
 [99] DT_0.33                     dplyr_1.1.4                
[101] crosstalk_1.2.1             RMariaDB_1.3.2             
[103] shiny_1.9.1                

loaded via a namespace (and not attached):
  [1] R.methodsS3_1.8.2       progress_1.2.3          nnet_7.3-19            
  [4] goftest_1.2-3           vctrs_0.6.5             spatstat.random_3.3-2  
  [7] shape_1.4.6.1           png_0.1-8               registry_0.5-1         
 [10] httpcode_0.3.0          deldir_2.0-4            parallelly_1.38.0      
 [13] MASS_7.3-61             httpuv_1.6.15           foreach_1.5.2          
 [16] withr_3.0.1             xfun_0.48               survival_3.7-0         
 [19] commonmark_1.9.2        crul_1.5.0              memoise_2.0.1          
 [22] systemfonts_1.1.0       ragg_1.3.3              GlobalOptions_0.1.2    
 [25] zoo_1.8-12              pbapply_1.7-2           R.oo_1.26.0            
 [28] Formula_1.2-5           prettyunits_1.2.0       KEGGREST_1.44.1        
 [31] promises_1.3.0          httr_1.4.7              restfulr_0.0.15        
 [34] globals_0.16.3          fitdistrplus_1.2-1      rstudioapi_0.17.0      
 [37] UCSC.utils_1.0.0        miniUI_0.1.1.1          generics_0.1.3         
 [40] base64enc_0.1-3         curl_5.2.3              zlibbioc_1.50.0        
 [43] polyclip_1.10-7         ca_0.71.1               GenomeInfoDbData_1.2.12
 [46] SparseArray_1.4.8       xtable_1.8-4            doParallel_1.0.17      
 [49] evaluate_1.0.1          S4Arrays_1.4.1          BiocFileCache_2.12.0   
 [52] preprocessCore_1.66.0   hms_1.1.3               irlba_2.3.5.1          
 [55] colorspace_2.1-1        filelock_1.0.3          ROCR_1.0-11            
 [58] reticulate_1.39.0       spatstat.data_3.1-2     magrittr_2.0.3         
 [61] lmtest_0.9-40           later_1.3.2             spatstat.geom_3.3-3    
 [64] future.apply_1.11.2     bspm_0.5.7              scattermore_1.2        
 [67] RcppAnnoy_0.0.22        pillar_1.9.0            nlme_3.1-166           
 [70] iterators_1.0.14        pwalign_1.0.0           caTools_1.18.3         
 [73] compiler_4.4.1          RSpectra_0.16-2         TSP_1.2-4              
 [76] tensor_1.5              plyr_1.8.9              crayon_1.5.3           
 [79] abind_1.4-8             BiocIO_1.14.0           graphlayouts_1.2.0     
 [82] bit_4.5.0               textshaping_0.4.0       codetools_0.2-20       
 [85] openssl_2.2.2           bslib_0.8.0             GetoptLong_1.0.5       
 [88] mime_0.12               splines_4.4.1           markdown_1.13          
 [91] circlize_0.4.16         fastDummies_1.7.4       dbplyr_2.5.0           
 [94] interp_1.1-6            knitr_1.48              blob_1.2.4             
 [97] utf8_1.2.4              clue_0.3-65             WriteXLS_6.7.0         
[100] listenv_0.9.1           checkmate_2.3.2         expm_1.0-0             
[103] Matrix_1.7-1            statmod_1.5.0           tzdb_0.4.0             
[106] tweenr_2.0.3            pkgconfig_2.0.3         cachem_1.1.0           
[109] RSQLite_2.3.7           DBI_1.2.3               impute_1.78.0          
[112] fastmap_1.2.0           ica_1.0-3               sass_0.4.9             
[115] ggstats_0.7.0           dotCall64_1.2           fontawesome_0.5.2      
[118] RANN_2.6.2              rpart_4.1.23            farver_2.1.2           
[121] tidygraph_1.3.1         yaml_2.3.10             latticeExtra_0.6-30    
[124] foreign_0.8-87          rtracklayer_1.64.0      cli_3.6.3              
[127] webshot_0.5.5           leiden_0.4.3.1          lifecycle_1.0.4        
[130] askpass_1.2.1           uwot_0.2.2              gtable_0.3.5           
[133] rjson_0.2.23            ggridges_0.5.6          progressr_0.14.0       
[136] jsonlite_1.8.9          RcppHNSW_0.6.0          seriation_1.5.6        
[139] bitops_1.0-9            bit64_4.5.2             assertthat_0.2.1       
[142] spatstat.utils_3.1-0    zip_2.3.1               jquerylib_0.1.4        
[145] spatstat.univar_3.0-1   R.utils_2.12.3          lazyeval_0.2.2         
[148] htmltools_0.5.8.1       sctransform_0.4.1       rappdirs_0.3.3         
[151] glue_1.8.0              httr2_1.0.5             RCurl_1.98-1.16        
[154] mnormt_2.1.1            jpeg_0.1-10             igraph_2.1.1           
[157] R6_2.5.1                labeling_0.4.3          GenomicFeatures_1.56.0 
[160] DelayedArray_0.30.1     tidyselect_1.2.1        htmlTable_2.4.3        
[163] maps_3.4.2              ggforce_0.4.2           future_1.34.0          
[166] munsell_0.5.1           KernSmooth_2.23-24      htmlwidgets_1.6.4      
[169] aroma.light_3.34.0      hwriter_1.3.2.1         rlang_1.1.4            
[172] spatstat.sparse_3.1-0   spatstat.explore_3.3-2  fansi_1.0.6
```

## Generate `manifest.json`

`Ctrl-C` out of the app and run:

```R
> install.packages("rsconnect")
[sudo] password for weismanal:
Install system packages as root...
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
Hit http://archive.ubuntu.com/ubuntu jammy InRelease
Get:1 http://security.ubuntu.com/ubuntu jammy-security InRelease [129 kB]
Get:2 http://archive.ubuntu.com/ubuntu jammy-updates InRelease [128 kB]
Hit http://archive.ubuntu.com/ubuntu jammy-backports InRelease
Ign https://r2u.stat.illinois.edu/ubuntu jammy InRelease
Hit https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/ InRelease
Hit https://r2u.stat.illinois.edu/ubuntu jammy Release
Fetched 257 kB in 0s (0 B/s)
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
Get:1 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-sys amd64 3.4.3-1.ca2204.1 [41.5 kB]
Get:2 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-askpass amd64 1.2.1-1.ca2204.1 [23.3 kB]
Get:3 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-base64enc amd64 0.1-3-1.ca2204.1 [26.9 kB]
Get:4 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-cli amd64 3.6.3-1.ca2204.1 [1267 kB]
Get:5 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-curl amd64 5.2.3-1.ca2204.1 [441 kB]
Get:6 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-digest amd64 0.6.37-1.ca2204.1 [208 kB]
Get:7 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-glue amd64 1.8.0-1.ca2204.1 [158 kB]
Get:8 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-jsonlite amd64 1.8.9-1.ca2204.1 [631 kB]
Get:9 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-rlang amd64 1.1.4-1.ca2204.1 [1556 kB]
Get:10 https://r2u.stat.illinois.edu/ubuntu jammy/main all r-cran-lifecycle all 1.0.4-1.ca2204.1 [113 kB]
Get:11 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-openssl amd64 2.2.2-1.ca2204.1 [649 kB]
Get:12 https://r2u.stat.illinois.edu/ubuntu jammy/main all r-cran-packrat all 0.9.2-1.ca2204.1 [662 kB]
Get:13 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-pki amd64 0.1-14-1.ca2204.1 [117 kB]
Get:14 https://r2u.stat.illinois.edu/ubuntu jammy/main all r-cran-rstudioapi all 0.17.0-1.ca2204.1 [309 kB]
Get:15 https://r2u.stat.illinois.edu/ubuntu jammy/main amd64 r-cran-yaml amd64 2.3.10-1.ca2204.1 [110 kB]
Get:16 https://r2u.stat.illinois.edu/ubuntu jammy/main all r-cran-rsconnect all 1.3.1-1.ca2204.1 [783 kB]
Fetched 7096 kB in 0s (0 B/s)
Selecting previously unselected package r-cran-sys.
(Reading database ... 62472 files and directories currently installed.)
Preparing to unpack .../00-r-cran-sys_3.4.3-1.ca2204.1_amd64.deb ...
Unpacking r-cran-sys (3.4.3-1.ca2204.1) ...
Selecting previously unselected package r-cran-askpass.
Preparing to unpack .../01-r-cran-askpass_1.2.1-1.ca2204.1_amd64.deb ...
Unpacking r-cran-askpass (1.2.1-1.ca2204.1) ...
Selecting previously unselected package r-cran-base64enc.
Preparing to unpack .../02-r-cran-base64enc_0.1-3-1.ca2204.1_amd64.deb ...
Unpacking r-cran-base64enc (0.1-3-1.ca2204.1) ...
Selecting previously unselected package r-cran-cli.
Preparing to unpack .../03-r-cran-cli_3.6.3-1.ca2204.1_amd64.deb ...
Unpacking r-cran-cli (3.6.3-1.ca2204.1) ...
Selecting previously unselected package r-cran-curl.
Preparing to unpack .../04-r-cran-curl_5.2.3-1.ca2204.1_amd64.deb ...
Unpacking r-cran-curl (5.2.3-1.ca2204.1) ...
Selecting previously unselected package r-cran-digest.
Preparing to unpack .../05-r-cran-digest_0.6.37-1.ca2204.1_amd64.deb ...
Unpacking r-cran-digest (0.6.37-1.ca2204.1) ...
Selecting previously unselected package r-cran-glue.
Preparing to unpack .../06-r-cran-glue_1.8.0-1.ca2204.1_amd64.deb ...
Unpacking r-cran-glue (1.8.0-1.ca2204.1) ...
Selecting previously unselected package r-cran-jsonlite.
Preparing to unpack .../07-r-cran-jsonlite_1.8.9-1.ca2204.1_amd64.deb ...
Unpacking r-cran-jsonlite (1.8.9-1.ca2204.1) ...
Selecting previously unselected package r-cran-rlang.
Preparing to unpack .../08-r-cran-rlang_1.1.4-1.ca2204.1_amd64.deb ...
Unpacking r-cran-rlang (1.1.4-1.ca2204.1) ...
Selecting previously unselected package r-cran-lifecycle.
Preparing to unpack .../09-r-cran-lifecycle_1.0.4-1.ca2204.1_all.deb ...
Unpacking r-cran-lifecycle (1.0.4-1.ca2204.1) ...
Selecting previously unselected package r-cran-openssl.
Preparing to unpack .../10-r-cran-openssl_2.2.2-1.ca2204.1_amd64.deb ...
Unpacking r-cran-openssl (2.2.2-1.ca2204.1) ...
Selecting previously unselected package r-cran-packrat.
Preparing to unpack .../11-r-cran-packrat_0.9.2-1.ca2204.1_all.deb ...
Unpacking r-cran-packrat (0.9.2-1.ca2204.1) ...
Selecting previously unselected package r-cran-pki.
Preparing to unpack .../12-r-cran-pki_0.1-14-1.ca2204.1_amd64.deb ...
Unpacking r-cran-pki (0.1-14-1.ca2204.1) ...
Selecting previously unselected package r-cran-rstudioapi.
Preparing to unpack .../13-r-cran-rstudioapi_0.17.0-1.ca2204.1_all.deb ...
Unpacking r-cran-rstudioapi (0.17.0-1.ca2204.1) ...
Selecting previously unselected package r-cran-yaml.
Preparing to unpack .../14-r-cran-yaml_2.3.10-1.ca2204.1_amd64.deb ...
Unpacking r-cran-yaml (2.3.10-1.ca2204.1) ...
Selecting previously unselected package r-cran-rsconnect.
Preparing to unpack .../15-r-cran-rsconnect_1.3.1-1.ca2204.1_all.deb ...
Unpacking r-cran-rsconnect (1.3.1-1.ca2204.1) ...
Setting up r-cran-curl (5.2.3-1.ca2204.1) ...
Setting up r-cran-rlang (1.1.4-1.ca2204.1) ...
Setting up r-cran-sys (3.4.3-1.ca2204.1) ...
Setting up r-cran-base64enc (0.1-3-1.ca2204.1) ...
Setting up r-cran-digest (0.6.37-1.ca2204.1) ...
Setting up r-cran-yaml (2.3.10-1.ca2204.1) ...
Setting up r-cran-packrat (0.9.2-1.ca2204.1) ...
Setting up r-cran-glue (1.8.0-1.ca2204.1) ...
Setting up r-cran-pki (0.1-14-1.ca2204.1) ...
Setting up r-cran-cli (3.6.3-1.ca2204.1) ...
Setting up r-cran-lifecycle (1.0.4-1.ca2204.1) ...
Setting up r-cran-askpass (1.2.1-1.ca2204.1) ...
Setting up r-cran-jsonlite (1.8.9-1.ca2204.1) ...
Setting up r-cran-rstudioapi (0.17.0-1.ca2204.1) ...
Setting up r-cran-openssl (2.2.2-1.ca2204.1) ...
Setting up r-cran-rsconnect (1.3.1-1.ca2204.1) ...
> library(rsconnect)

Attaching package: ‘rsconnect’

The following object is masked from ‘package:shiny’:

    serverInfo

> writeManifest(appDir = "/home/weismanal/notebook/2024-10-18/testing_sequin_installation_on_linux/public_sequin")
ℹ Capturing R dependencies from renv.lock
✔ Found 360 dependencies
```

The updated `manifest.json`, `renv.lock`, `renv/activate.R`, and `error_log.txt` from the above steps are now present in the `ubuntu_installation` branch of the [`public_sequin` repository](https://github.com/ncats/public_sequin).
