#!/bin/bash
## Complete hackish bullshit to get required r packages
## installed for the makeMeta.r function.
##

echo "Running ${PKG_NAME} post-link in PREFIX=${PREFIX}"

Rcode="
checkinstall <-function(package){\n
  if (!package %in% installed.packages()) {\n
    install.packages(package, repos=\"http://cran.us.r-project.org\")
  }\n
}\n
checkinstall(\"ape\")\n
checkinstall(\"TreeSim\")\n
checkinstall(\"meteR\")
"

echo -e $Rcode > /tmp/rpt.r
## Not sure how i feel about this, messages.txt gets piped to stdout
## in the event nothing bad happens so it shits all over the screen,
## somewhat annoying.
Rscript /tmp/rpt.r >> $PREFIX/.messages.txt 2>&1 
