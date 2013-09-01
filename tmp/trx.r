# Now let's do something with the raw data...
# --------------------------------------------
# We need to separate each species. Then we can input the raw data from
# each into an R data frame and run the TRX calculation separately on
# each species.
#
# I think the way to do this is to generate a temporary R script for each
# respective species that passess the $species variable to the script as
# the key for the data. The genetic data should all be stored together.
#
# First I will just pass it to R and calculate the average, then figure
# out the TRX calculation later.

	
