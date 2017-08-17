#!/bin/bash
cd ..
rm -rf rankSVMcompare-release
cp -r rankSVMcompare rankSVMcompare-release
##grep -v rankSVMcompareData rankSVMcompare/DESCRIPTION | grep -v Remotes > rankSVMcompare-release/DESCRIPTION
##rm rankSVMcompare-release/tests/testthat/test-rankSVMcompareData.R
PKG_TGZ=$(R CMD build rankSVMcompare-release|grep building|sed 's/.*‘//'|sed 's/’.*//')
R CMD INSTALL $PKG_TGZ
R CMD check --as-cran $PKG_TGZ
