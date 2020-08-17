

library(car)
library(spida2)

fit <- lm(prestige ~ income*type*education, Prestige)
summary(fit)

?lht
?Anova
Anova(fit)
Anova(fit, test = 'Chisq')
getD
car:::linearHypothesis.lm
car:::linearHypothesis.default
