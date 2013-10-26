library(ggplot2)

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
        correlation = cor(x = df$trxMean, y = df$energyMean, method = 'pearson')
        plot = ggplot() +
                   geom_point(data = df, aes(x = position, y = trxMean), color = 'blue') +
                   geom_point(data = df, aes(x = position, y = energyMean), color = 'red') +
                   geom_line(data = df, aes(x = position, y = trxMean)) +
                   geom_line(data = df, aes(x = position, y = energyMean)) +
                   guide_legend( title = legend ) +
                   stat_function( fun = correlation )
                   #geom_line(cor(df$trxMean, df$energyMean, method = "pearson"))

        plot

        plotPath = paste('./plots/', species, '/', gene, '.png', sep = '')
        ggsave(plotPath)
    }
}

speciesNames = c("bayanus","cerevisiae","martinae","paradoxus")
geneNames    = c("eYAL001C","eYCL031C","eYDR181C","eYER111C","eYGR100W")

for (x in speciesNames) {
    plotMe(x)
}

q()
