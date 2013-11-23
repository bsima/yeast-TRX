#!/bin/bash 
#

ROOT="`pwd`"

set -e
set -x

# This script will setup the Perl invironment by installing
# the Statistics::Descriptives module, instal necessary R packages
# such as Shiny, etc., and unpack the YeastGenome data.
#
# The code will look something like this:
sudo cpan
install Statistics::Descriptives
exit
# Basically just pipe the install command to cpan, then exit

# Next we will need to unpack the YeastGenome data:
unzip $ROOT/data/YeastGenome.zip
