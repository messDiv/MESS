#!/bin/bash
## Complete hackish bullshit to get required r packages
## installed for the makeMeta.r function.
##

echo "Running ${PKG_NAME} post-link in PREFIX=${PREFIX}"

## Fix a bug in R 3.6.0 package installs to fail like this:
## mv: cannot move <somefile> to <somedir> file exists 
## ERROR:   moving to final location failed 
export R_INSTALL_STAGED=false

Rcode="
checkinstall <-function(package){\n
  if (!package %in% installed.packages()) {\n
    install.packages(package, repos=\"http://cran.us.r-project.org\")
  }\n
}\n
checkinstall(\"ape\")\n
checkinstall(\"TreeSim\")\n
checkinstall(\"nleqslv\")\n

#checkinstall(\"phyloTop\")
"

echo -e $Rcode > /tmp/rpt.r
## Not sure how i feel about this, messages.txt gets piped to stdout
## in the event nothing bad happens so it shits all over the screen,
## somewhat annoying.
Rscript /tmp/rpt.r >> $PREFIX/.messages.txt 2>&1 
