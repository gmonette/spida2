# 2019_12_22
# Stacked Bar Plot with Colors and Legend
# Turn this into a method for tables to do simple 'mosaic plots'


counts <- table(mtcars$vs, mtcars$gear)
barplot(counts, main="Car Distribution by Gears and VS",
  xlab="Number of Gears", col=c("darkblue","red"),
  width = c(.2,.3,3),
  space = .01,
  legend = rownames(counts))

library(latticeExtra)
barchart(counts)
