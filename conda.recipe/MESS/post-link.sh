#!/bin/bash
## Complete hackish bullshit to get required r packages
## installed for the makeMeta.r function.
##
## This is a mess right now. pika needs devtools installed and
## devtools needs git2r which needs libssl, which there's no
## straightforward way to install automatically. ARG!
##
## libssl-dev    (package on e.g. Debian and Ubuntu)
## openssl-devel (package on e.g. Fedora, CentOS and RHEL)
## openssl       (Homebrew package on OS X)
##
## There is some magic for how this is incorporated into conda
## packages that I don't understand yet bcz this isn't getting
## run on install....

echo "Running ${PKG_NAME} post-link in PREFIX=${PREFIX}"

Rcode="
checkinstall <-function(package){\n
  if (!package %in% installed.packages()) {\n
    install.packages(package, repos=\"http://cran.us.r-project.org\")
  }\n
}\n
checkinstall(\"ape\")\n
checkinstall(\"TreeSim\")\n
checkinstall(\"sads\")
"

echo -e $Rcode > /tmp/rpt.r
Rscript /tmp/rpt.r >> $PREFIX/.messages.txt 2>&1 
