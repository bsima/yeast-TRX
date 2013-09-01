# Now let's do something with the raw data...
# --------------------------------------------
# We need to separate each species. Then we can input the raw data from
# each into an R data frame and run the TRX calculation separately on
# each species.
#
# I think the way to do this is to generate a temporary R script for each
# respective species that creates a separate data frame of the TRX sequences
# associated with that species. Then, each R script can be run from the main
# script.pl. We could even pass each species' R script to trx.r, thus
# calculating the TRX values and drawing the necessary graphs concurrently.
