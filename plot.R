# setwd("~/Code/Babbit")
bayanusTRX <- read.table("bayanus-TRXscore.csv",sep=",",header=TRUE)
cerevisiaeTRX <- read.table("cerevisiae-TRXscore.csv",sep=",",header=TRUE)
martinaeTRX <- read.table("martinae-TRXscore.csv",sep=",",header=TRUE)
paradoxusTRX <- read.table("paradoxus-TRXscore.csv",sep=",",header=TRUE)

library(ggplot2)
bayanusPlot<-qplot(data = bayanusTRX, x=position, y=trx.score, facets=~gene, color=gene)
bayanusPlot
ggsave("bayanusPlot.png")

cerevisiaePlot<-qplot(data = cerevisiaeTRX, x=position, y=trx.score, facets=~gene, color=gene)
cerevisiaePlot
ggsave("cerevisiaePlot.png")

martinaePlot<-qplot(data = martinaeTRX, x=position, y=trx.score, facets=~gene, color=gene)
martinaePlot
ggsave("martinaePlot.png")

paradoxusPlot<-qplot(data = paradoxusTRX, x=position, y=trx.score, facets=~gene, color=gene)
paradoxusPlot
ggsave("paradoxusPlot.png")



# rm(bayanusTRX)
# rm(cerevisiaeTRX)
# rm(martinaeTRX)
# rm(paradoxusTRX)
# rm(bayanusPlot)
