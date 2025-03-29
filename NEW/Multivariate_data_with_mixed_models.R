# y ~ Var/(Treat * Time) - 1,
# 
# random = ~ Var - 1 | id,
# 
# weights = varIdent(form = ~ 1 | Var),
# 
# Sometimes 1 is added, subtracted, or follows a tilde. I am guessing that these mean very different things, but I am not sure what these meanings are.
# 
# 3) Can you explain the structure of the M Matrix a bit more...for example, how do you know what values should be 1, 0 and, -1 respectively? I have some guesses at the general idea...but some things confuse me. For example, why are the first two values equal to 1 for "Var 1 Treat B at baseline". If the first column is Treatment A at baseline, shouldn't it be 0, just as the second column was 0 when you were estimating Treatment A at baseline in the first row?
# 
# 4) What does "DinD" stand for? I googled around but was not able to find anything on this. 
# 
# 
# On Fri, Apr 2, 2021 at 4:02 PM Georges A Monette <georges@yorku.ca> wrote:
# Hi Andrew,
# 
# Here's some sample code for a 'doubly' multivariate analysis using lme. It's using the equivalent of repeated measures anova. Doing repeated measures manova would require coding a new corClass in lme. 
# 
# All the best,
# Georges


library(spida2)
library(nlme)
Nid <- 50
dt <- expand.grid(
  id = 1:Nid,
  time = 1:3,
  var = 1:2)
set.seed(123)
dt <-
  within(
    dt,
    {
      treat <- 1 + (id %% 2)
      Treat <- factor(paste0(' ', c('A','B')[treat]))
      Var <- factor(paste0(' ', var))
      y <- rnorm(id) + (1 + 1*(var-1))*rnorm(Nid)[id] + 1 * (var - 1)
      Time <- factor(paste0(' ',time))
    }
  )
xqplot(dt)
head(dt)

fit <- lme(
  y ~ Var/(Treat * Time) - 1,
  data = dt,
  random = ~ Var - 1 | id,
  correlation = corSymm(form = ~ var | id/time),
  weights = varIdent(form = ~ 1 | Var),
  control = list(msMaxIter = 1000, msMaxEval = 1000, msVerbose = TRUE))
summary(fit)
L1 <- rbind(
  'Var 1 Treat A at baseline'                            = c(1,0,  0, 0,  0,0,  0,0,  0,0, 0,0),                                
  'Var 1 Treat B at baseline'                            = c(1,0,  1, 0,  0,0,  0,0,  0,0, 0,0),                                
  'Var 1 Treat A at time 2'                              = c(1,0,  0, 0,  1,0,  0,0,  0,0, 0,0),                                
  'Var 1 Treat B at time 2'                              = c(1,0,  1, 0,  1,0,  0,0,  1,0, 0,0),                                
  'Var 1 Treat A at time 3'                              = c(1,0,  0, 0,  0,0,  1,0,  0,0, 0,0),                                
  'Var 1 Treat B at time 3'                              = c(1,0,  1, 0,  0,0,  1,0,  0,0, 1,0),                                
  'Var 1 Treatment difference at baseline'               = c(0,0,  1, 0,  0,0,  0,0,  0,0, 0,0),                                
  'Var 1 Treatment effect (B-A) at time 2 (DinD)'        = c(0,0,  0, 0,  0,0,  0,0,  1,0, 0,0),
  'Var 1 Treatment effect (B-A) at time 3 (DinD)'        = c(0,0,  0, 0,  0,0,  0,0,  0,0, 1,0),
  'Var 1 Treatment effect (B-A) in change from 2 to 3'   = c(0,0,  0, 0,  0,0,  0,0, -1,0, 1,0))
wald(fit, L1)  
L2 <- cbind(0,L1)[,-(ncol(L1)+1)]
rownames(L2) <- sub('Var 1','Var 2', rownames(L1))
wald(fit, L2)
Lall <- rbind(L1,L2)
#'
#' Multivariate tests
#'
for(i in 1:10 ) {
  cat('\n\n')
  print(wald(fit, Lall[c(i,i+10),]))
}
y ~ Var/(Treat * Time) - 1,

random = ~ Var - 1 | id,

weights = varIdent(form = ~ 1 | Var),

Sometimes 1 is added, subtracted, or follows a tilde. I am guessing that these mean very different things, but I am not sure what these meanings are.

3) Can you explain the structure of the M Matrix a bit more...for example, how do you know what values should be 1, 0 and, -1 respectively? I have some guesses at the general idea...but some things confuse me. For example, why are the first two values equal to 1 for "Var 1 Treat B at baseline". If the first column is Treatment A at baseline, shouldn't it be 0, just as the second column was 0 when you were estimating Treatment A at baseline in the first row?

4) What does "DinD" stand for? I googled around but was not able to find anything on this. 


On Fri, Apr 2, 2021 at 4:02 PM Georges A Monette <georges@yorku.ca> wrote:
Hi Andrew,

Here's some sample code for a 'doubly' multivariate analysis using lme. 
It's using the equivalent of repeated measures anova. 
Doing repeated measures manova would require coding a new corClass in lme. 

All the best,
Georges


library(spida2)
library(nlme)
Nid <- 50
dt <- expand.grid(
  id = 1:Nid,
  time = 1:3,
  var = 1:2)
set.seed(123)
dt <-
  within(
    dt,
    {
      treat <- 1 + (id %% 2)
      Treat <- factor(paste0(' ', c('A','B')[treat]))
      Var <- factor(paste0(' ', var))
      y <- rnorm(id) + (1 + 1*(var-1))*rnorm(Nid)[id] + 1 * (var - 1)
      Time <- factor(paste0(' ',time))
    }
  )
xqplot(dt)
head(dt)

fit <- lme(
  y ~ Var/(Treat * Time) - 1,
  data = dt,
  random = ~ Var - 1 | id,
  correlation = corSymm(form = ~ var | id/time),
  weights = varIdent(form = ~ 1 | Var),
  control = list(msMaxIter = 1000, msMaxEval = 1000, msVerbose = TRUE))
summary(fit)
L1 <- rbind(
  'Var 1 Treat A at baseline'                            = c(1,0,  0, 0,  0,0,  0,0,  0,0, 0,0),                                
  'Var 1 Treat B at baseline'                            = c(1,0,  1, 0,  0,0,  0,0,  0,0, 0,0),                                
  'Var 1 Treat A at time 2'                              = c(1,0,  0, 0,  1,0,  0,0,  0,0, 0,0),                                
  'Var 1 Treat B at time 2'                              = c(1,0,  1, 0,  1,0,  0,0,  1,0, 0,0),                                
  'Var 1 Treat A at time 3'                              = c(1,0,  0, 0,  0,0,  1,0,  0,0, 0,0),                                
  'Var 1 Treat B at time 3'                              = c(1,0,  1, 0,  0,0,  1,0,  0,0, 1,0),                                
  'Var 1 Treatment difference at baseline'               = c(0,0,  1, 0,  0,0,  0,0,  0,0, 0,0),                                
  'Var 1 Treatment effect (B-A) at time 2 (DinD)'        = c(0,0,  0, 0,  0,0,  0,0,  1,0, 0,0),
  'Var 1 Treatment effect (B-A) at time 3 (DinD)'        = c(0,0,  0, 0,  0,0,  0,0,  0,0, 1,0),
  'Var 1 Treatment effect (B-A) in change from 2 to 3'   = c(0,0,  0, 0,  0,0,  0,0, -1,0, 1,0))
wald(fit, L1)  
L2 <- cbind(0,L1)[,-(ncol(L1)+1)]
rownames(L2) <- sub('Var 1','Var 2', rownames(L1))
wald(fit, L2)
Lall <- rbind(L1,L2)
#'
#' Multivariate tests
#'
for(i in 1:10 ) {
  cat('\n\n')
  print(wald(fit, Lall[c(i,i+10),]))
}
