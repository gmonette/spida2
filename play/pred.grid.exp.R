library(spida2)
library(car)
library(lattice)
library(latticeExtra)

hs
pred.grid <- function(...) {
  nams <- as.character(as.list(substitute(list(...)))[-1L])
  x <- list(...)
  names(x)[names(x) == ''] <- nams[names(x) == '']
  x <- lapply(x, unique)
  do.call(expand.grid, x)
}

fit4<-lm(mathach ~ Sector*ses*Sex, hs)
summary(fit4)


hsm <- within(
  hs,
  {
    id1 <- factor(paste(Sector, ses, Sex))
    id <- reorder(id1, ses + I(Sector == 'Public')*1000)
  }
)

up(hs,~school)

pred <- with(hsm, pred.grid(id, ses = seq(-3,3,.01)))
pred <- merge(pred, up(hs,~school), by = 'id', all.x = TRUE)
dim(pred)
head(pred)


pred$fit4 <- predict(fit4, newdata = pred)


xyplot(mathach ~ Sector*ses*Sex, hsm , groups = Sex,
       # skip = rep(c(F,T,F), c(21,3,30)),
       layout = c(7,6),
       between = list(y =c(0,0,.3,0,0,0)),
       par.strip.text = list(cex=.6),
       auto.key = list(space = 'right', lines = T)) +
  xyplot(fit4 ~ ses | id1, pred, groups = Sex, type = 'l')
