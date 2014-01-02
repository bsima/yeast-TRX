library(plyr)
library(ggplot2)

# Some ideas for sensible coding:
#
# * Prefix my variable names with something, such as `my` or `var`. This doesn't actually scope them, it just signifies that I set them.
# * Setup the root directory var again.
# * Open this script in RStudio for testing.
# * I think this will become the `server.ui` logic of the Shiny app. It might be useful to check that out, and finish writing
#   the script with the Shiny app in mind.

myPath = getcwd()
myRoot = sub('code$', '', myPath)


myGeneNameDirs = paste(myRoot, '/data/YeastGenome-tmp', sep = '')
myGeneNames    = list.files(path = myGeneNameDirs, pattern = '(e.{1,})[[:punct:]]aln')
myGeneNames    = sub('[[:punct:]]aln', '', myGeneNames) # This gets rid of the filename .aln ending

myMakePlot = function(species) {

    mySpeciesFolder = paste(myRoot, '/plots/', species, sep = '')
    dir.create(mySpeciesFolder)

    for (myGene in myGeneNames) {

        mySmoothData = paste(myRoot, '/data/', species, '/', gene, '/smooth.csv', sep = '')

        myData = read.table(mySmoothData, sep = ',', header = TRUE)

        # plot = qplot(data = df, x = 1:length(df$position), y = trxMean, facets =~ gene, color = gene) + geom_line()
        #genePosition = df$position
        #trxMean = df$trxMean
        #energyMean = df$energyMean
        myCorrelation = cor(x = myData$trxMean, y = abs(myData$energyMean), use = "na.or.complete")

        WC( myData$postion, myData$  )

        myPlot = ggplot() +
                   geom_point(data = myData, aes(x = position, y = trxMean), color = 'blue') +
                   geom_point(data = myData, aes(x = position, y = abs(energyMean)), color = 'red') +
                   geom_line(data = myData, aes(x = position, y = trxMean)) +
                   geom_line(data = myData, aes(x = position, y = abs(energyMean))) +
                   annotate("text", label = paste('r=', correlation, sep=''),x=1000,y=5,size=5, color = 'red')
        myPlot # Generates the plot

        myPlotPath = paste('./plots/', species, '/', gene, '.png', sep = '')
        ggsave(plotPath)
    }
}

mySpeciesNames = c("bayanus","cerevisiae","martinae","paradoxus")
#geneNames    = c("eYAL001C","eYCL031C","eYDR181C","eYER111C","eYGR100W")

for (x in mySpeciesNames) {
    myMakePlot(x)
}

#q()
