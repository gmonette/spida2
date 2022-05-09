library(spida2)
library(nlme)
fit <- lme(mathach ~ ses , hs, random = ~ 1 + ses| school )
summary(fit)
fit <- lme(mathach ~ ses , hs, random = list( school = ~ 1 + ses ))
fits <- lme(mathach ~ ses , hs, random = list( school = pdSymm(form= ~ 1 )))
library(car)
compareCoefs(fit,fits)
summary(fits)
summary(fit)
hs %>% 
  split(~ school+Sex) -> zz
str(zz)
?split
lme.formula
as.call((mathach ~ ses * Sex)[-2])
fits <- lme(mathach ~ ses , hs, random = list( school = pdDiag(form= ~ 1 + ses )))
summary(fits)
fits <- lme(mathach ~ ses , hs, random = list( school = pdInd(form= ~ 1 , cov = 1)))
fits <- lme(mathach ~ ses , hs, random = list( school = pdSymm(form= ~ 1 )))
debug(pdInd)
debug(pdConstruct.pdInd)
debug(pdInd)
debug(pdConstruct.pdInd)
debug(pdSymm)
debug(nlme:::pdConstruct.pdSymm)

?pdMat
nlme:::pdConstruct.pdSymm


