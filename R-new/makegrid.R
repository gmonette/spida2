#' 
#' 
#' make a grid from a fit based on getD
#' 
library(spida2)
fit <- lm(mathach ~ ses * Sex, hs)
getD(fit) %>% head
# get makegrid from UTFA
