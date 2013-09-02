Notes
=====

First we need to open the yeast genome files and read the data. The *C files (the "Crick" files, as opposed to the *W or "Watson" files) are read in a bottom-to-top, right-to-left fashion. Thus, I will need to first run a reverse compliment function on these files as a whole.

The files themselves are organized with the following criteria:

* Line 0 is a header with information about what the file contains
* Line 3 begins a 4-line block of information. Each line describes part of a different species of yeast, described by a 4-letter code (Scer, Spar, Smik, Sbay)
* Following the 4-letter code, there are 12 spaces. On column 17 begins a 59-letter sequence of genetic code. Capital letters are coding proteins.
* After this block is a line of asterisks, the meaning of which I don't yet know. Then there is a blank line.
'

With each file, I need to do two things:

1. Calculate the TRX value and store this data somewhere (probably in an R data set)
2. Calculate the delta E value and store this data somewhere

Then, with the aggregate data, plot both values on a graph (or pair of graphs).

This will allow us to visualize and correlate the differences between the two measures.

## 13/9/2013 ##


1. Make the output script only show the gene name and not the whole filepath. I think this can be done with just regular expressions and such
2. Make header line print to the output file properly
