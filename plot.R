library(ggplot2)

plotMe = function(species) {

    path = paste('./data/', species, '-TRXscore.csv', sep = '')
    class(path)
    
    df = read.table(path, sep = ',', header = TRUE)
    class(df)

    plot = qplot(data = df, x = position, y = trx.score, facets =~ gene, color = gene)
    class(plot)

    plot

    plotPath = paste('./plots/', species, '.png', sep = '')
    ggsave(plotPath)

}

speciesNames = c("bayanus","cerevisiae","martinae","paradoxus")

for (x in speciesNames) {
    plotMe(x)
}

q()
