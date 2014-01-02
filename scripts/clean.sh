#!/bin/bash
#
# Just get rid of crap

ROOT="`pwd`"

rm -rf $ROOT/data/bayanus*;
echo "Removed data/bayanus*";

rm -rf $ROOT/data/cerevisiae*;
echo "Removed data/cerevisiae*";

rm -rf $ROOT/data/martinae*;
echo "Removed data/martinae*";

rm -rf $ROOT/data/paradoxus*;
echo "Removed data/paradoxus*";

rm -rf $ROOT/figures/*;
echo "Removed figures/*";

rm $ROOT/Rplots.pdf;
echo "Removed Rplots.pdf";

echo "All clean!";
