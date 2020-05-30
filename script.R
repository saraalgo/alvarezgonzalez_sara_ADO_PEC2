
library(tweeDEseqCountData)
library(tweeDEseq)
library(edgeR)


1.

datos_f <- normalizeCounts(datos)
maPlot(datos_f[,1], datos_f[,2],
       pch=19, cex=.5, ylim=c(-8,8), 
       allCol="darkgrey", lowess=TRUE)
grid(col="black")
title("normalización TMM")