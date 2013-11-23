library(plyr)
library(ggplot2)

geneNameDirs = './data/YeastGenome-tmp' 
geneNames    = list.files(path = geneNameDirs, pattern = '(e.{1,})[[:punct:]]aln')
geneNames    = sub('[[:punct:]]aln', '', geneNames) # This gets rid of the filename .aln ending

plotMe = function(species) {

    speciesFolder = paste('./plots/', species, sep = '')
    dir.create(speciesFolder)

    for (gene in geneNames) {

        # geneFolder = paste('./plots/', species, '/', gene, sep = '')
        # dir.create(geneFolder)

        smoothData = paste('./data/', species, '/', gene, '/smooth.csv', sep = '') 
    
        df = read.table(smoothData, sep = ',', header = TRUE)

        # plot = qplot(data = df, x = 1:length(df$position), y = trxMean, facets =~ gene, color = gene) + geom_line()
        #genePosition = df$position
        #trxMean = df$trxMean
        #energyMean = df$energyMean
        correlation = cor(x = df$trxMean, y = abs(df$energyMean), use = "na.or.complete")
        plot = ggplot() +
                   geom_point(data = df, aes(x = position, y = trxMean), color = 'blue') +
                   geom_point(data = df, aes(x = position, y = abs(energyMean)), color = 'red') +
                   geom_line(data = df, aes(x = position, y = trxMean)) +
                   geom_line(data = df, aes(x = position, y = abs(energyMean))) +
                   annotate("text", label = paste('r=', correlation, sep=''),x=1000,y=5,size=5, color = 'red')
        plot

        plotPath = paste('./plots/', species, '/', gene, '.png', sep = '')
        ggsave(plotPath)
    }
}

speciesNames = c("bayanus","cerevisiae","martinae","paradoxus")
#geneNames    = c("eYAL001C","eYCL031C","eYDR181C","eYER111C","eYGR100W")

for (x in speciesNames) {    
    plotMe(x)
}

#q()
