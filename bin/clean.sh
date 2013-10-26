#!/bin/bash
#
# Just get rid of crap

rm -rf data/bayanus*;
echo "Removed data/bayanus*";

rm -rf data/cerevisiae*;
echo "Removed data/cerevisiae*";

rm -rf data/martinae*;
echo "Removed data/martinae*";

rm -rf data/paradoxus*;
echo "Removed data/paradoxus*";

rm -rf plots/*;
echo "Removed plots/*";

rm Rplots.pdf;
echo "Removed Rplots.pdf";

echo "All clean!";
