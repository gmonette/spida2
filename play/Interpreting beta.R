#' What does beta mean?
#' 
betaq <- function(cmat) {
  solve(cbind(1,cmat))
}

zf <- factor(sample(LETTERS[1:4],10,T))
levels(zf)
library(spida2)
contr.poly(4) %>% betaq %>%  svd
contr.treatment(4) %>% betaq %>% svd
contr.sum(4) %>% betaq %>% svd
contr.helmert(4) %>% betaq %>% svd
contr.nhelmert(4) %>% betaq %>% svd

# Variance patterns from a diagonal G matrix
# 
vg <- function(cmat, G = diag(1:nrow(cmat))) {
  V <- cbind(1,cmat) %*% G %*% t(cbind(1,cmat))
  V
}
contr.treatment(3) %>% vg
contr.treatment(3) %>% vg(diag(3))
contr.treatment(3) %>% vg(diag(3:1))
contr.treatment(3) %>% vg(diag(c(1,0,0)))
contr.treatment(3) %>% vg(diag(c(1,0.1,0.1)))
contr.treatment(3) %>% vg(diag(c(1,0.01,0.01)))
