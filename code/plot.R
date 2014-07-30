library(plyr)
library(ggplot2)

# Some ideas for sensible coding:
#
# * Prefix my variable names with something, such as `my` or `var`. This doesn't actually scope them, it just signifies that I set them.
# * Setup the root directory var again.
# * Open this script in RStudio for testing.
# * I think this will become the `server.ui` logic of the Shiny app. It might be useful to check that out, and finish writing
#   the script with the Shiny app in mind.

# Get the root of the repo
myPath = getwd()
myRoot = sub('code$', '', myPath)

# mySpeciesNamesDirs = paste(myRoot, '/data', sep = '')
# mySpeciesNames     = dir(path = mySpeciesNamesDirs, pattern = '([A-Za-z]+)\b(?<!Genome)(?<!tmp)', recursive = FALSE)
# mySpeciesNames     = list.dirs(path = mySpeciesNamesDirs, recursive = FALSE)
# mySpeciesNames     = sub('(.+YeastGenome.)|(.+data.")', '', mySpeciesNames)

mySpeciesNames = c('bayanus', 'cerevisiae', 'martinae', 'paradoxus')

myGeneNameDirs = paste(myRoot, '/data/YeastGenome-tmp', sep = '')
myGeneNames    = list.files(path = myGeneNameDirs, pattern = '(e.{1,})[[:punct:]]aln')
myGeneNames    = sub('[[:punct:]]aln', '', myGeneNames) # This gets rid of the filename .aln ending

myWeightedCorrelation = function( x, y, w = rep(1,length(x)) ) {
    stopifnot( length(x) == dim(y)[2] )
    w = w / sum(w)
    # Center x and y, using the weighted means
    x  = x - sum(x * w)
    ty = t(y - colSums(t(y) * w))
    # Compute the variance
    vy = sum(w * x * x)
    vy = colSums(w * ty * ty)
    # Compute the covariance
    vxy = colsums(ty * x * w)
    # Compute the correlation
    vxy / sqrt(vx * vy)
}

myMakePlot = function(species) {

    mySpeciesFolder = paste(myRoot, '/figures/', species, sep = '')
    dir.create(mySpeciesFolder)

    for (myGene in myGeneNames) {

        mySmoothData = paste(myRoot, '/data/', species, '/', myGene, '/smooth.csv', sep = '')

        myData = read.table(mySmoothData, sep = ',', header = TRUE)

        # plot = qplot(data = df, x = 1:length(df$position), y = trxMean, facets =~ gene, color = gene) + geom_line()
        #genePosition = df$position
        #trxMean = df$trxMean
        #energyMean = df$energyMean
        myCorrelation = cor(x = myData$trxMean, y = abs(myData$energyMean), use = "na.or.complete")

        myWeightedData = myWeightedCorrelation( myData$postion, myData$  )

        myPlot = ggplot() +
                   geom_point(data = myData, aes(x = position, y = trxMean), color = 'blue') +
                   geom_point(data = myData, aes(x = position, y = abs(energyMean)), color = 'red') +
                   geom_line(data = myData, aes(x = position, y = trxMean)) +
                   geom_line(data = myData, aes(x = position, y = abs(energyMean))) +
                   geom_line(data = myWeightedData, aes(x = position, y = IDK)) +
                   annotate("text", label = paste('r=', myCorrelation, sep=''),x=1000,y=5,size=5, color = 'red')
        myPlot # Generates the plot

        myPlotPath = paste('./plots/', species, '/', gene, '.png', sep = '')
        ggsave(plotPath)
    }
}

for (x in mySpeciesNames) {
    myMakePlot(x)
}

#q()
