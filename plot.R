library(ggplot2)

plotMe = function(species) {

    path = paste('./data/', species, '/calc.csv', sep = '') 
    
    df = read.table(path, sep = ',', header = TRUE)

    plot = qplot(data = df, x = 1:length(df$position), y = trx.score, facets =~ gene, color = gene) + geom_line()

    plot

    plotPath = paste('./plots/', species, '.png', sep = '')
    ggsave(plotPath)

}

speciesNames = c("bayanus","cerevisiae","martinae","paradoxus")

for (x in speciesNames) {
    plotMe(x)
}

q()
