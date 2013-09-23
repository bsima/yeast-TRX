library(ggplot2)

makePlot <- function(species) {
    x = "read.table('./data/"
    y = "-TRXscore.csv',sep=',',header=TRUE)" 
    species.TRX = paste(x, species, y, sep="") 
    
    # species.TRX = paste("read.table('./data/","species",,sep="")
    # species.TRX <- read.table("./data/"species"-TRXscore.csv",sep=",",header=TRUE)
    species.Plot <- qplot(data = species.TRX, x=position, y=trx.score, facets=~gene, color=gene)
    species.Plot
    
    a = "'plot/"
    b = "Plot.png'"
    plotPath = paste(a, species, b, sep="")
    ggsave(plotPath)
}

# setwd("~/Code/Babbit")
# bayanusTRX <- read.table("./data/bayanus-TRXscore.csv",sep=",",header=TRUE)
# cerevisiaeTRX <- read.table("./data/cerevisiae-TRXscore.csv",sep=",",header=TRUE)
# martinaeTRX <- read.table("./data/martinae-TRXscore.csv",sep=",",header=TRUE)
# paradoxusTRX <- read.table("./data/paradoxus-TRXscore.csv",sep=",",header=TRUE)

speciesNames = c("bayanus","cerevisiae","martinae","paradoxus")

for (species in speciesNames) {
    makePlot(species)
}


# bayanusPlot<-qplot(data = bayanusTRX, x=position, y=trx.score, facets=~gene, color=gene)
# bayanusPlot
# ggsave("plot/bayanusPlot.png")
# 
# cerevisiaePlot<-qplot(data = cerevisiaeTRX, x=position, y=trx.score, facets=~gene, color=gene)
# cerevisiaePlot
# ggsave("plot/cerevisiaePlot.png")
# 
# martinaePlot<-qplot(data = martinaeTRX, x=position, y=trx.score, facets=~gene, color=gene)
# martinaePlot
# ggsave("plot/martinaePlot.png")
# 
# paradoxusPlot<-qplot(data = paradoxusTRX, x=position, y=trx.score, facets=~gene, color=gene)
# paradoxusPlot
# ggsave("plot/paradoxusPlot.png")

q()

# rm(bayanusTRX)
# rm(cerevisiaeTRX)
# rm(martinaeTRX)
# rm(paradoxusTRX)
# rm(bayanusPlot)
