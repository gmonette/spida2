head(hs)
# hs$csize <- with(hs, capply(school, school, length))
sp <- function(x) gsp(x, 0, 1, 0)
fit <- lme( mathach ~ sp(ses) + cvar(ses, school) + Sector*csize + Sector*sp(ses), hs,
            random = ~ 1 | school)
summary(fit)
fit2 <- lme( mathach ~ sp(dvar(ses,school)) + cvar(ses, school) + 
               csize + Sector*sp(dvar(ses,school)), hs,
            random = ~ 1 | school)
summary(fit)
summary(fit2)
AIC(fit,fit2)
summary(lm(csize ~ Sector, up(hs, ~ school)))

pred <- with(hs, pred.grid(csize = c(40,50,60),ses.m = c(-1,0,1), ses.d = c(-1,0,1), Sector))
xqplot(pred)
pred$school <- with(pred, paste(ses.m, Sector, csize))
pred$ses <- with(pred, ses.m + ses.d)
pred$fit <- predict(fit, newdata = pred, level = 0)
pred$fit2 <- predict(fit2, newdata = pred, level = 0)
pred
(xyplot(fit ~ ses | Sector * csize, pred, groups = school, type = 'b' )+
  layer(panel.grid(h=-1,v=-1))) %>% 
  useOuterStrips
xyplot(fit ~ ses | Sector * csize, pred, groups = school, type = 'b' )+
  layer(panel.grid(h=-1,v=-1))
%>% 
  useOuterStrips
xyplot(fit2 ~ ses | Sector * csize, pred, groups = school, type = 'b' ) %>% 
  useOuterStrips

---
title: "R Script Template"
author: "Georges Monette"
date: "2024-11-16"
output: pdf_document
---

```{r setup, include=FALSE}

library(pacman)

pacman::p_load('rstatix',
               'readxl',
               'kableExtra',
               'dplyr',
               'tidyr',             #### added GM for pivot_wider
               'rstudioapi',        #### added GM to set working directory
               'lattice',
               'latticeExtra',
               'latex2exp',
               'ggdist',
               'ggridges',
               'ggbeeswarm',
               'ggthemes',
               'EnvStats',
               'lavaan',
               'semPlot',
               'lme4',
               'nnet',
               'foreign',
               'spida2',
               'flextable',
               'nlme')

# setwd_here()
options(knitr.table.format = "latex")

knitr::opts_chunk$set(echo = FALSE)
knitr::opts_chunk$set(echo = TRUE, comment = '    ')
knitr::opts_chunk$set(cache = TRUE)


# 
# Functions:
# 

layer <- latticeExtra::layer
vif <- car:::vif.default                    # GM: may need to install car package
head.tbl_df <- function(x, ...) head(as.data.frame(x), ...)
h <- function(x,...) head(x,...)
xq <- function(x, name = paste0(substitute(x),'-xq.pdf'), ...) {
  pdf(name, title = name)
  xqplot(x,...)
  dev.off()
}

#    
summary.gls <- function (object, verbose = FALSE, ...) 
{
    # modified summary.gls method to suppress printing of
    # beta correlation matrix
    # vif is more informative for large models with
    #    multiparameter factors and effects
    fixSig <- attr(object[["modelStruct"]], "fixedSigma")
    fixSig <- !is.null(fixSig) && fixSig
    stdBeta <- sqrt(diag(as.matrix(object$varBeta)))
    corBeta <- t(object$varBeta/stdBeta)/stdBeta
    beta <- coef(object)
    dims <- object$dims
    dimnames(corBeta) <- list(names(beta), names(beta))
    object$corBeta <- matrix('Printing Supressed',1,1)
    tTable <- data.frame(beta, stdBeta, beta/stdBeta, beta)
    dimnames(tTable) <- list(names(beta), c("Value", "Std.Error", 
        "t-value", "p-value"))
    tTable[, "p-value"] <- 2 * pt(-abs(tTable[, "t-value"]), 
        dims$N - dims$p)
    object$tTable <- as.matrix(tTable)
    resd <- resid(object, type = "pearson")
    if (length(resd) > 5) {
        resd <- quantile(resd, na.rm = TRUE)
        names(resd) <- c("Min", "Q1", "Med", "Q3", "Max")
    }
    object$residuals <- resd
    aux <- logLik(object)
    structure(c(object, list(BIC = BIC(aux), AIC = AIC(aux))), 
        verbose = verbose, class = c("summary.gls", class(object)))
}

```

\tableofcontents

# Notes

- Removed 

# Data Input

```{r data read}

setwd(dirname(getActiveDocumentContext()$path))
                 
dir <- "C:/Users/sharm/OneDrive/Documents/Productivity/Research/Open/UHN/CANTBI/"           # for Bhanu's PC
dir <- '../data/'                                                                           # for Georges's PC

gose <- read_excel(paste0(dir,"cantbi-gose-dec72023.xlsx"))
gose$mult <- with(gose, capply(`Injury Date`,`Pat Id`, function(x) length(unique(x))))
tab_(gose, ~ mult)
tab_(gose, ~ `Pat Id`) %>% c %>% tab



ct <- read_excel(paste0(dir,"CT_findings.xlsx"))
tab_(ct, ~ `Pat Id`) %>% c %>% tab    
ct$mult <- with(ct, ave(`Pat Id`,`Pat Id`, FUN = length))
subset(ct, mult > 1) %>% h
# 
# Patient 3342 has multiple assessments in 'ct' with a date for
# Updates but there is no ax date.
# 
subset(gose, `Pat Id` == 3342) %>% h

# One Pat Id '3342' has two Ax in ct.  
# The variable 'Injury Date' can be used, with 'Pat Id', 
# to create a 'stint' id in both the 'gose' and the 'dem'
# data frames but not in the 'ct' data frame.
# 
# 
names(ct) %and% names(gose)
names(ct) %and% names(dem)

# names(gose) %>% intersect(names(ct))

dem <- read_excel(paste0(dir,"cantbi-dem-dec72023.xlsx"))
tab_(dem, ~ `Pat Id`) %>%c %>%  tab # No Pat Id is duplicated

# Problem: GONE
# There's a Pat Id with multiple stints in gose which can
# be identified using Injury Date 
# but the Pat Id is also duplicated in ct without an
# obvious variable to disambiguate the stints.
# 
# 

```

## Data exploration and derived variables

- Is there a variable that measures elapsed time since injury?
  - Concern: 'Gose Date' is available for all occasions but
    'Injury Date' is missing for a notable proportion of
    Adult Versions.  Could it be imputed by backdating from
    the date of the first 'Gose Date'? Answer: Problem because
    there's a large variability in the elapsed time from the
    first 'Gose Date' and 'Injury Date' ranging from almost
    immediate to over 400 days.
    __However__ these dates are clumpy perhaps corresponding to
    a categorical variable on time category which could then be used
    to better impute injury date if that variable is available
    for cases with missing injury date.
    
    ... on further analysis: Once we've eliminated all the individuals
    who are dead throughout the study, or those who die during the study,
    there are only 4 individuals excluded due to incomplete data to
    calculate 'etime' (Elapsed time since injury)
    
```{r imputing-injury-date}
gose <- within(
 gose,
  {
    ID     <- `Pat Id`
    id     <- paste0(ID,":", `Injury Date`) # one subject has two stints ??
    Days_since_injury <- as.numeric(`Gose Date`) - as.numeric(`Injury Date`)
    Months_since_injury <- Days_since_injury * (12/365)
    etime  <- Months_since_injury
    doa    <- capply(`Gose Score`, id, function(x) max(x) == 1)
    died   <- capply(`Gose Score`, id, function(x) min(x) == 1)
    Version <- tr(`Gose Version`,
                  c('adult', 'ped','pwee'),
                  c('Adult','Pediatric','PWEE'))
    gose <- ifelse(
      Version != 'Adult', 
      9 - as.numeric(`Gose Score`), 
      as.numeric(`Gose Score`)
      )
    
    # adding a jitter per id

    gose_j <- gose + capply(id, id, function(i) .15 * rnorm(1))
    
  }
)
up(gose, ~ id) %>% tab_(~ID) %>% c %>% tab  # OK 
gose$id
# xq(gose)
```
### Longitudinal data set


```{r longitudinal}
dl <- left_join(gose, ct, by = 'Pat Id', multiple = 'last', suffix = c(".gose",".ct"))
subset(dl, select = c('Pat Id', grepv('*\\.gose$|*\\.ct$', names(dl)))) %>% h
dl <- left_join(dl, dem, by = 'Pat Id', suffix = c('.dl','.dem'))
xq(dl)
subset(dl, `Injury Date.x` != `Injury Date.y`)
names(dl)

dl <- within(
  dl,
  {
    age <- as.numeric(`Age Y`) + 
      ifelse(
        `Age M` == '-', 
        .5 / 12,
        as.numeric(`Age M`)/12
      ) +
      ifelse(
        `Age D` == '-',
        15/30 * 1/12, 
        as.numeric(`Age D`)/365)
    Sev <- `Gcs 24.x`
    Severity <- tr(Sev, 
                   c('Mild','Moderate','Severe'),
                   c('Mild', 'msTBI', 'msTBI')
    )
    time <- as.numeric(tr(`Gose Timepoint`,
               c('days7to10','month_3','month_6','month_12'),
               c(1,2,3,4))
    )
    etime[etime <0] <- NA
    gose[gose<2.5] <- NA
    etime2 <- etime
  }
)

# data checks and corrections
dl <- as.data.frame(dl)

xyplot(etime ~ time, dl)
subset(dl, time == 1 & etime > 5)
id1 <- subset(dl, time == 1 & etime > 5)$id
id2 <- subset(dl, time == 4 & etime < 5)$id
subset(dl, id == id1, c(`Gose Timepoint`, time, etime, `Gose Date`, `Injury Date.y`))
subset(dl, id == id2, c(`Gose Timepoint`, time, etime, `Gose Date`, `Injury Date.y`))
dl[,c('etime','time')]

# post-facto test:
stopifnot(with(dl, all(`Gcs 24.x`==`Gcs 24.y`)))

```

## Longitudinal Analysis

```{r longitudinal2}

xyplot(gose_j ~ etime|Version*Severity, dl , groups = id, type = 'l') %>% 
  useOuterStrips

sp <- function(x) gsp(x, 10, c(0,1),0)

# mod1 <- nlme(gose ~ theta1 + theta2*exp(-theta3*etime) + alpha *  sp(etime2), data=dl, subset = Severity == "Mild",
#              fixed=list(
#         theta1 ~ 1 ,
#         theta2 ~ 1 ,
#         theta3 ~ 1,
#         alpha ~ 1),
#     random=list(id = list(theta1 ~ 1, theta2 ~ 1, alpha ~ 1)),
#     start=list(fixed=c(8, - 3, .2,  -.1)),
#     na.action = na.omit,
#     control=nlmeControl(msMaxIter=200,maxIter = 100, msVerbose = TRUE,
#                         returnObject =TRUE))
# 
# mod1
# coef(mod1)
# pred <- data.frame(etime = seq(0,36,.1), etime2 = 0)
# pred$y <- predict(mod1, level = 0,  newdata= pred)
# 
# xyplot(y ~ etime, pred, ylim = c(0,10))

# mod 2
# 
dl$Sev_mild <- with(dl, 1*(Severity == 'Mild'))
mod2 <- nlme(
  gose ~ theta1 + theta2*exp(-theta3*etime) , 
  data=dl, 
  fixed=list(
        theta1 ~ 1 + Sev_mild ,
        theta2 ~ 1 + Sev_mild,
        theta3 ~ 1+ Sev_mild) ,
    random=list(id = list(theta1 ~ 1, theta2 ~ 1)),
    start=list(fixed=c(8,0, - 3,0,.2,0)),
    na.action = na.omit,
    control=nlmeControl(msMaxIter=200,maxIter = 100, msVerbose = TRUE,
                        returnObject =TRUE))

mod2
coef(mod2)
pred <- expand.grid(etime = seq(0,36,.1), Sev_mild = c(0,1))
pred$y <- predict(mod2, level = 0,  newdata= pred)

xyplot(y ~ etime, pred, groups = Sev_mild, ylim = c(0,10), type = 'l',
       auto.key = T)


# 
# mod 3
# 
# Simple recovery spline
# 
sp <- function(x) gsp(x, 8, 1, 0)
xyplot(gose ~ etime|Severity*Sex*Version, dl, groups = id, type = 'b')
dl$Sev_mild <- with(dl, 1*(Severity == 'Mild'))
dlnopwee <- droplevels(subset(dl, Version != 'PWEE'))
dlnopwee %>% tab(~ Severity+Sex+Version)
mod2 <- lme(
  gose ~ (sp(etime) + Severity + Sex + Version)^3, 
  data=dlnopwee, 
  random = list(ID = ~ sp(etime)),
    na.action = na.omit,
    control=lmeControl(maxIter = 100, msVerbose = TRUE,
                       msMaxIter = 500, msMaxEval = 20000,
                        returnObject =TRUE))

summary(mod2)
Anova(mod2)
wald(mod2, 'Sex')
coef(mod2)

## Drop Sex

mod3 <- lme(
  gose ~ (sp(etime) + Severity + Version)^3, 
  data=dlnopwee, 
  random = list(ID = ~ sp(etime)),
    na.action = na.omit,
    control=lmeControl(maxIter = 100, msVerbose = TRUE,
                       msMaxIter = 500, msMaxEval = 20000,
                        returnObject =TRUE))

Anova(mod3)

pred3 <- with(dlnopwee, pred.grid(etime = seq(0,18,.1), Severity, Version ))
pred3$y <- predict(mod3, level = 0,  newdata= pred3)

xyplot(y ~ etime | Severity, pred3, groups = Version, ylim = c(0,10), type = 'l',
       auto.key = T, lwd = 4) +
  xyplot(gose ~ etime | Severity, dlnopwee, groups = id, type = 'l')

predI <- with(dl, pred.grid(etime = seq(0,36,.1), Severity, Sex, ID))
predI$y <- predict(mod2, level = 1,  newdata= predI)
tab(up(dl, ~ ID), ~ Severity + Sex)
xyplot(y ~ etime | Severity*Sex, predI, groups = ID, ylim = c(0,10), type = 'l',
       auto.key = T)


```

##   Using baseline as covariate

```{r baseline-as-covariate-models}

dl$baseline <- with(dl, capply(dl, ~id, with, c(gose[time == '1'],NA)[1]))
dl$mintime <- with(dl, capply(dl, ~id, with, c(min(time))))
zork
sp <- function(x) gsp(x, 8, 1, 0)
xyplot(gose ~ etime|Severity*Sex*Version, dl, groups = id, type = 'b')
dl$Sev_mild <- with(dl, 1*(Severity == 'Mild'))
dlnopwee <- droplevels(subset(dl, Version != 'PWEE'))
dlnopwee %>% tab(~ Severity+Sex+Version)
mod2 <- lme(
  gose ~ (sp(etime) + Severity + Sex + Version)^3, 
  data=dlnopwee, 
  random = list(ID = ~ sp(etime)),
    na.action = na.omit,
    control=lmeControl(maxIter = 100, msVerbose = TRUE,
                       msMaxIter = 500, msMaxEval = 20000,
                        returnObject =TRUE))

summary(mod2)
Anova(mod2)
wald(mod2, 'Sex')
coef(mod2)

## Drop Sex

mod3 <- lme(
  gose ~ (sp(etime) + Severity + Version)^3, 
  data=dlnopwee, 
  random = list(ID = ~ sp(etime)),
    na.action = na.omit,
    control=lmeControl(maxIter = 100, msVerbose = TRUE,
                       msMaxIter = 500, msMaxEval = 20000,
                        returnObject =TRUE))

Anova(mod3)

pred3 <- with(dlnopwee, pred.grid(etime = seq(0,18,.1), Severity, Version ))
pred3$y <- predict(mod3, level = 0,  newdata= pred3)

xyplot(y ~ etime | Severity, pred3, groups = Version, ylim = c(0,10), type = 'l',
       auto.key = T, lwd = 4) +
  xyplot(gose ~ etime | Severity, dlnopwee, groups = id, type = 'l')

predI <- with(dl, pred.grid(etime = seq(0,36,.1), Severity, Sex, ID))
predI$y <- predict(mod2, level = 1,  newdata= predI)
tab(up(dl, ~ ID), ~ Severity + Sex)
xyplot(y ~ etime | Severity*Sex, predI, groups = ID, ylim = c(0,10), type = 'l',
       auto.key = T)










```


```{r}
options(max.print=100000)
# subset(gose, select = c('id','Injury Date','Gose Date', 'etime','Gose Timepoint','Gose Score')) %>% as.data.frame %>% sortdf(~id)
#
# Q: Are missing etime's worth salvaging?    
# A: Only 4 subjects. Probably not worth imputing.
# 
subset(gose, select = c('id','Injury Date','Gose Date', 'etime','Gose Timepoint','Gose Score','doa','died')) %>% 
  subset(`Gose Score` != '-') %>% 
  as.data.frame %>% 
  subset(doa == FALSE & died == FALSE) %>% 
  sortdf(~id) 

tab(~is.na(etime) + I(`Gose Score`=='-'), gose)
tab(gose$`Injury Date`=='-')
tab(is.na(etime)+is.na(`Injury Date`)+is.na(`Gose Date`), gose)
head(gose)
goseID <- up(gose, ~ ID)
head(goseID)
xyplot(minetime ~  |factor(`Gose Version`), gose)
tab(~is.na(`Gose Timepoint`), gose)
tab(~`Gose Timepoint`, gose)

tab(~ is.na(etime)+`Gose Version`, gose)
tab(~ is.na(`Gose Date`)+`Gose Version`, gose) # complete
tab(~ I(`Injury Date`=='-')+`Gose Version`, gose) # complete
tab(~ is.na(etime), gose) # complete
gose[, c('etime','Injury Date','Gose Date')] %>% as.data.frame
```

### Visualizing GOSE vs elapsed time since injury

```{r imputing-injury-date2}


tab(gose, ~`Gose Version`)

xyplot(gose_j ~ etime|`Gose Version`, sortdf(gose, ~etime), groups = id, type = 'b')
xyplot(gose_j ~ etime|`Gose Version`, sortdf(gose, ~etime), groups = id, type = 'b')
xyplot(gose_j ~ etime|`Gose Version`, sortdf(gose, ~etime), groups = id, type = 'b')

```



```{r data clean}

# selecting variables of interest
gose_small <- gose %>%
  select(2, 4, 5, 8, 26) %>%
  filter(`Gose Version` %in% c("ped", "adult")) %>%  
  rename("Timepoint" = "Gose Timepoint",
         "Version" = "Gose Version",
         "GOSE" = "Gose Score",
         "Severity" = "Gcs 24",
         "ID" = "Pat Id")
ct <- ct %>%
  #select(2, 38) %>%
  rename("ID" = "Pat Id")
dem <- dem %>%
  #select(2:9) %>%
  rename("ID" = "Pat Id",
         "Severity" = "Gcs 24")
# 
# Elapsed time from injury to Gose Date
# 
gose <- left_join(gose, ct, by = 'ID', multiple = 'first')

ct %>% within({
  multiplicity <- capply(ID, ID, length)
}) %>% subset(multiplicity > 1) %>% Head
gose %>% subset(ID == '3342') %>% Head
# ID 3342 occurs twice in 'gose' data frame and
# if ct data frame
# 

# creating a dataset for just the pwee (i.e., a sub-population not required for the main analyses)
gose_pwee <- gose %>%
  select(2, 4, 5, 8, 26) %>%
  filter(`Gose Version` %in% "pwee") %>%
  rename("Timepoint" = "Gose Timepoint",
         "Version" = "Gose Version",
         "GOSE" = "Gose Score",
         "Severity" = "Gcs 24",
         "ID" = "Pat Id")
gose_pwee$Timepoint <- recode_factor(
  gose_pwee$Timepoint,
  "days7to10" = "Days 7 to 10",
  "month_3" = "Month 3",
  "month_6" = "Month 6",
  "month_12" = "Month 12"
)
gose_pwee$Timepoint <- gose_pwee$Timepoint %>%
  factor(levels = c("Days 7 to 10", "Month 3", "Month 6", "Month 12"))

# the pediatric and adult GOSE data are reverse coded
# creating variable GOSEv2 to get all data in the same orientation
gose_small <- gose_small %>%
  mutate(GOSEv2 = GOSE) %>%
  mutate(GOSEv2 = case_when(
    (GOSE == 1 & Version == "ped") ~ 8,
    (GOSE == 2 & Version == "ped") ~ 7,
    (GOSE == 3 & Version == "ped") ~ 6,
    (GOSE == 4 & Version == "ped") ~ 5,
    (GOSE == 5 & Version == "ped") ~ 4,
    (GOSE == 6 & Version == "ped") ~ 3,
    (GOSE == 7 & Version == "ped") ~ 2,
    (GOSE == 8 & Version == "ped") ~ 1,
    (GOSE == 1 & Version == "adult") ~ 1,
    (GOSE == 2 & Version == "adult") ~ 2,
    (GOSE == 3 & Version == "adult") ~ 3,
    (GOSE == 4 & Version == "adult") ~ 4,
    (GOSE == 5 & Version == "adult") ~ 5,
    (GOSE == 6 & Version == "adult") ~ 6,
    (GOSE == 7 & Version == "adult") ~ 7,
    (GOSE == 8 & Version == "adult") ~ 8
  ))

# removing '-' and replacing them with blanks so R can read them as NAs
gose_small[3] <- lapply(gose_small[3], function(col) as.numeric( gsub("-$|\\,", "", col) ) )

# giving timepoints numeric values for ease of use
gose_small <- gose_small %>%
  mutate(Time = case_when(
    Timepoint == "days7to10" ~ 1,
    Timepoint == "month_3" ~ 2,
    Timepoint == "month_6" ~ 3,
    Timepoint == "month_12" ~ 4
  )) %>%
  mutate(Time = as.factor(Time))

# creating the msTBI category
gose_small <- gose_small %>%
  mutate(Sev = case_when(
    Severity == "Mild" ~ "Mild",
    Severity == "Moderate" ~ "msTBI",
    Severity == "Severe" ~ "msTBI"
  ))

# recoding the Version labels
gose_small$Version <- recode(
  gose_small$Version,
  "ped" = "Pediatric",
  "adult" = "Adult"
)

# recoding the Timepoint labels
gose_small$Timepoint <- recode_factor(
  gose_small$Timepoint,
  "days7to10" = "Days 7 to 10",
  "month_3" = "Month 3",
  "month_6" = "Month 6",
  "month_12" = "Month 12"
)
gose_small$Timepoint <- gose_small$Timepoint %>%
  factor(levels = c("Days 7 to 10", "Month 3", "Month 6", "Month 12"))

# merging dem data to gose
gose_small <- full_join(
  gose_small, dem, by = c("ID")
)

# merging ct data now
gose_small <- full_join(
  gose_small, ct, by = c("ID")
)

# removing '-' and replacing them with blanks
gose_small[11:13] <- lapply(gose_small[11:13], function(col) as.numeric( gsub("-$|\\,", 0, col) ) )

# replacing NA with 9, to code as "CT not done"
gose_small$CT <- replace(gose_small$CT, is.na(gose_small$CT), 9)
gose_small$CT_factor <- as.factor(gose_small$CT)

# creating an age variable
gose_small <- gose_small %>%
  mutate(Age = `Age Y` + (`Age M`/12) + (`Age D`/365)) %>%
  select(-c("Severity.y")) %>%
  rename("Severity" = "Severity.x")

subset(gose_small, Age < 1) %>% as.data.frame
subset(gose_small, Age < 1) %>% as.data.frame
# 
# Remove subjects < 1.0 year of age
# 
gose_small <- subset(gose_small, Age > 1.0)


# identifying the number of datapoints per subject
gose_timepoints <- gose_small %>%
  group_by(ID) %>%
  summarise(Datapoints = length(Time))

# dropping all GOSEv2=1 data for analysis, as GOSEv2 = 1 = deceased
# NOTE the group wants to still report GOSE = 1 for descriptives
gose_clean <- gose_small %>%
  subset(GOSEv2 %in% c(2:8, NA)) 

# merging with the gose_timepoints dataset
gose_clean <- left_join(
  gose_clean, gose_timepoints, by = c("ID")
)



```

# Data visualization

## Data summary

```{r counts}

counts <- gose_small %>%
  group_by(Timepoint, Version, Sev) %>%
  summarise(n = n())
flextable(counts)

```

```{r lineplot}

overview <- gose_small %>%
  filter(Time != 'NA') %>%
  group_by(ID) %>%
  ggplot(aes(x = Time, y = GOSEv2, group = ID)) +
  geom_line(aes(color = ID)) +
  facet_grid(rows = vars(Version), cols = vars(Sev)) +
  theme_hc()
overview + theme(legend.position = 'none')

```



## Descriptive tables

Currently, these tables include GOSE = 1. If this needs to be changed again, the tables can be re-produced.

```{r descriptive table function}

descriptive_table <- function(data, grouping_vars, outcome) {
  gose_tab <- data %>%
    group_by(all_of(across(grouping_vars))) %>%
    get_summary_stats(outcome)
  gose_tab %>%
    kbl() %>%
    kable_styling(latex_options = c("scale_down", "hold_position"))
}

```

**Table 1a**. GOSE by TIMEPOINT.

```{r gose tab1a, warning=FALSE}

descriptive_table(gose_small, c('Timepoint'), 'GOSEv2')

```

**Table 1b**. GOSE by TIMEPOINT and AGE GROUP.

```{r gose tab1b, warning=FALSE}

descriptive_table(gose_small, c('Timepoint', 'Version'), 'GOSEv2')

```

**Table 1c**. GOSE by TIMEPOINT and SEVERITY.  

```{r gose tab1c, warning=FALSE}

descriptive_table(gose_small, c('Timepoint', 'Severity'), 'GOSEv2')

```

\newpage

**Table 1d**. PWEE by TIMEPOINT.  

```{r gose tab1d, echo=FALSE, warning=FALSE}

gose_pwee[3] <- lapply(gose_pwee[3], function(col) as.numeric( gsub("-$|\\,", "", col) ) )

descriptive_table(gose_pwee, c('Timepoint'), 'GOSE')

```

\newpage

## Boxplots

```{r boxplot function}

descriptive_table <- function(data, grouping_vars, outcome) {
  gose_tab <- data %>%
    group_by(all_of(across(grouping_vars))) %>%
    get_summary_stats(outcome)
  gose_tab %>%
    kbl() %>%
    kable_styling(latex_options = c("scale_down", "hold_position"))
}

boxplot_density <- function(data, grouping_vars, title, facet_var) {
  gose_box <- data %>%
  subset(Version %in% grouping_vars) %>%
  ggplot(aes(Time, GOSEv2, color = Time)) +
  ggtitle(title) +
  stat_n_text() +
  labs(x = "Time", y = "GOSE") +
  guides(color = "none") +
  coord_cartesian(ylim = c(0,8)) +
  facet_wrap(facet_var) +
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA
  ) +
  geom_point(
    size = 1.5,
    alpha = .2,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  theme_hc()
gose_box
}

```

Currently, these figures include GOSE = 1. If this needs to be changed again, the figures can be re-produced.

**Figure 1a** Box-and-whisker, density, and raincloud plot of GOSE (ADULT vs. PEDIATRIC)

```{r gose fig1a, echo=FALSE, warning=FALSE}

boxplot_density(gose_small, c('Adult', 'Pediatric'), 'GOSE by TIME & AGE', ~Version)

```

\newpage

**Figure 1b** Box-and-whisker, density, and raincloud plot of ADULT GOSE (by SEVERITY)

```{r gose fig1b, echo=FALSE, warning=FALSE}

boxplot_density(gose_small, c('Adult'), 'ADULT: GOSE by TIME & SEVERITY', ~Sev)

```

\newpage

**Figure 1c** Box-and-whisker, density, and raincloud plot of PEDIATRIC GOSE (by SEVERITY)  

```{r gose fig1c, echo=FALSE, warning=FALSE}

boxplot_density(gose_small, c('Pediatric'), 'PEDIATRIC: GOSE by TIME & SEVERITY', ~Sev)

```

\newpage

## Histograms of GOSE at 12-months

```{r histogram function}

hist_plot <- function(data, grouping_vars, sev_vars, title) {
  gose_hist <- data %>%
    subset(Version %in% grouping_vars) %>%
    subset(Sev %in% sev_vars) %>%
    subset(Time == 4) %>%
    ggplot(aes(GOSEv2)) +
    ggtitle(title) +
    geom_histogram(bins = 8, color = "#000000", fill = "#0099F8") + 
    stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1.5) +
    coord_cartesian(ylim = c(0,150)) +
    theme_hc()
  gose_hist
}

```


**Fig 2a** GOSE scores at 12-months (all patients)

```{r gose his fig2a, echo=FALSE, warning=FALSE}

hist_plot(gose_small, c('Adult', 'Pediatric'), c('Mild', 'msTBI'), '12-mo GOSE SCORES: ALL participants')

```

\newpage

**Fig 2b** GOSE scores at 12-months (all MILD TBI patients)

```{r gose his fig2b, echo=FALSE, warning=FALSE}

hist_plot(gose_small, c('Adult', 'Pediatric'), c('Mild'), '12-mo GOSE SCORES: MTBI participants')

```

\newpage

**Fig 2c** GOSE scores at 12-months (all msTBI patients)

```{r gose his fig2c, echo=FALSE, warning=FALSE}

hist_plot(gose_small, c('Adult', 'Pediatric'), c('msTBI'), '12-mo GOSE SCORES: MS TBI participants')

```

\newpage  

**Fig 2d** GOSE scores at 12-months (PEDIATRIC only)

```{r gose his fig2d, echo=FALSE, warning=FALSE}

hist_plot(gose_small, c('Pediatric'), c('Mild', 'msTBI'), '12-mo GOSE SCORES: PEDIATRIC participants')

```

\newpage

**Fig 2e** GOSE scores at 12-months (all PEDIATRIC MILD TBI patients)

```{r gose his fig2e, echo=FALSE, warning=FALSE}

hist_plot(gose_small, c('Pediatric'), c('Mild'), '12-mo GOSE SCORES: PEDIATRIC MTBI participants')

```

\newpage

**Fig 2f** GOSE scores at 12-months (all PEDIATRIC msTBI patients)

```{r gose his fig2f, echo=FALSE, warning=FALSE}

hist_plot(gose_small, c('Pediatric'), c('msTBI'), '12-mo GOSE SCORES: PEDIATRIC MS TBI participants')

```

\newpage

**Fig 2g** GOSE scores at 12-months (ADULT only)

```{r gose his fig2g, echo=FALSE, warning=FALSE}

hist_plot(gose_small, c('Adult'), c('Mild', 'msTBI'), '12-mo GOSE SCORES: ADULT participants')

```

\newpage

**Fig 2h** GOSE scores at 12-months (All ADULT MILD TBI participants)

```{r gose his fig2h, echo=FALSE, warning=FALSE}

hist_plot(gose_small, c('Adult'), c('Mild'), '12-mo GOSE SCORES: ADULT MTBI participants')

```

\newpage

**Fig 2i** GOSE scores at 12-months (All ADULT msTBI participants)

```{r gose his fig2i, echo=FALSE, warning=FALSE}

hist_plot(gose_small, c('Adult'), c('msTBI'), '12-mo GOSE SCORES: ADULT MS TBI participants')

```

\newpage


# Models using linear splines


- The following models allows for different variances at different times which
  is suggested by the appearance of a plot against time.

- 'correlation' and 'weights' have the effect of creating a unconstrained matrix for the marginal variance
  of within-subject responses.
  
- 'corSymm' allows for an arbitrary/free correlation matrix (form = ~time) makes sure time is matched up across ids
  
- 'weights' is also set to have an arbitrary structure and
  ultimately allows for an arbitrary variance-covariance matrix
  
- This is analogous to multi-variate repeated measures
  
- 'na.omit' omits the instance not the case

```{r spline based model function}

spm <- function(x)  gsp(x, c(2,3), 1, 0) 
spms <- function(x)  gsp(x, c(3), 1, 0)

```

```{r gls}
gose_clean_gls <- gose_clean
gose_clean_gls$nTimeID <- with(gose_clean_gls, capply(ID, paste(Time, ID), length))
gose_clean_gls <- gose_clean_gls %>%
  subset(!is.na(Time) & !is.na(Sev) & nTimeID == 1 & GOSEv2 > 1) %>%
  mutate(time = as.numeric(Time)) %>%
  mutate(time1 = time - 1)

gose_clean_mstbi <- gose_clean_gls %>%
  subset(Sev == 'msTBI' & time > 1) %>%
  mutate(time1 = time - 1)

gose_clean_mild <- gose_clean_gls %>%
  subset(Sev == 'Mild')

# mild TBI
# Q4GM - should the varIdent time variable also be `time` not `Time`?    
#   
#   A: It could be. The argument for corSymm needs to be an integer
#      with lowest value 1. The 'Time' variable is a factor with
#      levels '1', '2', etc. and I'm in the habit of using a factor
#      to denote distinct levels not numerically related, which is
#      the intention in producing arbitrary variances for each time point.
#      However, in this case, gls would interpret a numeric variable
#      as if it were categorical and 'time' would work equally well.
#      
fit_mild <- gls(GOSEv2 ~ spm(time) * Version, data = gose_clean_mild, 
              correlation = corSymm(form = ~ time | ID), 
              weights = varIdent(form = ~ 1 | Time), 
              na.action=na.omit)
summary(fit_mild)

# checking for possible interaction with Version - does not seem to be
wald(fit_mild, ':')
```

GM note:

Although there is no evidence of interactions, it might, nevertheless, not be reasonable
to adopt a model that assumes parallel effects among adult and pediatric participants.
We would be assuming parallel effects due to a lack of evidence to the contrary, but this
is reasonable only if it makes sense to believe that effects could be parallel in these
two populations. That is, failure to reject the null does not justify assuming
the null, unless the null is plausible to begin with.

I'll show plots for the non-additive models as well as the additive models below.

They are quite similar but the non-additive plots show that the trajectories are
relatively parallel not because the model constained them to be parallel but
because the data show relatively parallel trajectories in the two groups. 



## Additive model with inference on change per time point 

```{r gls2}
# additive model
fit_mild_add <- gls(GOSEv2 ~ spm(time) + Version, data = gose_clean_mild,
                    correlation = corSymm(form = ~ time | ID), 
                    weights = varIdent(form = ~ 1 | Time), 
                    na.action=na.omit)
summary(fit_mild_add)
# evidence of effects over time
wald(fit_mild_add, 'spm')

# estimates of slopes and changes in slopes
L <- cbind(0, sc(spm, x = c(1,2,3,2,3), D=c(1,1,1,1,1), type = c(1,1,1,2,2) ), 0)
rownames(L) <- c('Change from Week1toMonth3','Change from Month3toMonth6', 'Change from Month6toMonth12', 'Change in change at Month3', 'Change in change at Month6')
  

L_levels <-
  cbind(1, sc(spm, x = c(1,2,3,4,1,2,3,4), D = 0), rep(0:1, each = 4))
rownames(L_levels) <- paste0(rep(c('Adult','Pediatric'), each = 4), ' ', c('Week 1','Month 3','Month 6','Month 12'))
```

### Estimated levels -- additive model - mild

```{r gls2p}
wald(fit_mild_add, L_levels)
#'
#' Estimated levels -- additive model - mild - graph
#'
ww <- waldf(fit_mild_add, L_levels) 
ww$months <- c(0,3,6,12)
ww$Version <- rep(c('Adult','Pediatric'), each = 4)
xyplot(coef ~ months, ww, type = 'l', groups = Version,
       main = 'Mild to Severe',
       ylab = TeX('GOSEv2 ($\\pm$SE)'),
       xlim = c(0,12),
       auto.key = list(reverse.rows = T),
       fit = ww$coef,
       upper = ww$coef + ww$se,
       lower = ww$coef - ww$se,
       subscripts = T) +
  latticeExtra::glayer(panel.fit(...))

L_diffs <-
  rbind(
    'Adult M3 - W1' =     c(-1, 1,  0, 0,    0, 0, 0, 0),
    'Adult M6 - W1' =     c(-1, 0,  1, 0,    0, 0, 0, 0),
    'Adult M12 - W1' =    c(-1, 0, -0, 1,    0, 0, 0, 0),
    'Adult M6 - M3' =     c(0, -1,  1, 0,    0, 0, 0, 0),
    'Adult M12 - M3' =    c(0, -1,  0, 1,    0, 0, 0, 0),
    'Adult M12 - M6' =    c(0, 0,  -1, 1,    0, 0, 0, 0),
    'Pediatric M3 - W1' = c(0,0,0,0,  -1, 1, 0, 0),
    'Pediatric M6 - W1' = c(0,0,0,0,  -1, 0, 1, 0),
    'Pediatric M12 - W1' = c(0,0,0,0,  -1, 0, 0, 1),
    'Pediatric M6 - M3' = c(0,0,0,0,  0, -1, 1, 0),
    'Pediatric M12 - M3' = c(0,0,0,0,  0, -1, 0, 1),
    'Pediatric M12 - M6' = c(0,0,0,0,  0, 0, -1, 1),
    'Pediatric - Adult at W1' = c(-1, 0,0,0, 1, 0, 0, 0),
    'Pediatric - Adult at M3' = c(0, -1, 0,0,0, 1, 0, 0),
    'Pediatric - Adult at M6' = c(0, 0, -1, 0,0,0, 1, 0),
    'Pediatric - Adult at M12' = c(0, 0, 0, -1, 0,0,0, 1)
  
  ) %*% L_levels

#
# Levels:
#
wald(fit_mild_add, L_levels)
#'
#' Differences
#'
wald(fit_mild_add, L_diffs)
#'
#'
```

## Additive model with inference on change per **month**

Same model but reporting inferences on rate of change per month instead of change
from one time point to another.


```{r gls3}
#'
# additive model
head(gose_clean_mild)
gose_clean_mild <-
  within(
    gose_clean_mild,
    {
      month <- tr(time, 1:4, c(0, 3, 6, 12))
    }
  )

spmmon <- function(month) gsp(month, c(3,6), 1, 0) 
fit_mild_add_mo <- gls(GOSEv2 ~ spmmon(month) + Version, data = gose_clean_mild,
                    correlation = corSymm(form = ~ time | ID), 
                    weights = varIdent(form = ~ 1 | Time), 
                    na.action=na.omit)
summary(fit_mild_add_mo)
# evidence of effects over time
wald(fit_mild_add_mo, 'spmmon')

# estimates of slopes and changes in slopes
L <- cbind(0, sc(spmmon, x = c(0,3,6,3,6), D=c(1,1,1,1,1), type = c(1,1,1,2,2) ), 0)
rownames(L) <- c('rate/month W1 to M3','rate/month M3 to M6', 'rate/month M6 to M12', 
                 'change in rate/month at M3', 'change in rate/month at M6')


wald(fit_mild_add_mo, L)

L_levels <-
  cbind(1, sc(spmmon, x = c(0,3,6,12,0,3,6,12), D = 0), rep(0:1, each = 4))
rownames(L_levels) <- paste0(rep(c('Adult','Pediatric'), each = 4), ' ', c('Week 1','Month 3','Month 6','Month 12'))
#'
#' Estimated levels -- additive model - mild
#'
wald(fit_mild_add_mo, L_levels)
```

### Estimated levels -- additive model - mild - graph

```{r gls3p}
ww <- waldf(fit_mild_add_mo, L_levels) 
ww$months <- c(0,3,6,12)
ww$Version <- rep(c('Adult','Pediatric'), each = 4)
xyplot(coef ~ months, ww, type = 'l', groups = Version,
       main = 'Mild to Severe',
       ylab = TeX('GOSEv2 ($\\pm$SE)'),
       xlim = c(0,12),
       auto.key = list(reverse.rows = T),
       fit = ww$coef,
       upper = ww$coef + ww$se,
       lower = ww$coef - ww$se,
       subscripts = T) +
  latticeExtra::glayer(panel.fit(...))

L_diffs <-
  rbind(
    'Adult M3 - W1' =     c(-1, 1,  0, 0,    0, 0, 0, 0),
    'Adult M6 - W1' =     c(-1, 0,  1, 0,    0, 0, 0, 0),
    'Adult M12 - W1' =    c(-1, 0, -0, 1,    0, 0, 0, 0),
    'Adult M6 - M3' =     c(0, -1,  1, 0,    0, 0, 0, 0),
    'Adult M12 - M3' =    c(0, -1,  0, 1,    0, 0, 0, 0),
    'Adult M12 - M6' =    c(0, 0,  -1, 1,    0, 0, 0, 0),
    'Pediatric M3 - W1' = c(0,0,0,0,  -1, 1, 0, 0),
    'Pediatric M6 - W1' = c(0,0,0,0,  -1, 0, 1, 0),
    'Pediatric M12 - W1' = c(0,0,0,0,  -1, 0, 0, 1),
    'Pediatric M6 - M3' = c(0,0,0,0,  0, -1, 1, 0),
    'Pediatric M12 - M3' = c(0,0,0,0,  0, -1, 0, 1),
    'Pediatric M12 - M6' = c(0,0,0,0,  0, 0, -1, 1),
    'Pediatric - Adult at W1' = c(-1, 0,0,0, 1, 0, 0, 0),
    'Pediatric - Adult at M3' = c(0, -1, 0,0,0, 1, 0, 0),
    'Pediatric - Adult at M6' = c(0, 0, -1, 0,0,0, 1, 0),
    'Pediatric - Adult at M12' = c(0, 0, 0, -1, 0,0,0, 1)
    
  ) %*% L_levels

#'
#' Levels:
#'
wald(fit_mild_add_mo, L_levels)
#'
#' Differences
#'
wald(fit_mild_add_mo, L_diffs)
```

## Interaction model with inference on change per **month**


```{r gls4}
#'
#'
# interaction model
gose_clean_mild %>% as.data.frame %>% names %>% grepv('CT', .)
gose_clean_mild %>% as.data.frame %>% .[ ,c("CT GROUP"  ,"CT" ,       "CT_factor")] %>% head(30)
gose_clean_mild <-
  within(
    gose_clean_mild,
    {
      month <- tr(time, 1:4, c(0, 3, 6, 12))
    }
  )

spmmon <- function(month) gsp(month, c(3,6), 1, 0) 
fit_mild_int_mo <- gls(GOSEv2 ~ spmmon(month) * Version, data = gose_clean_mild,
                       correlation = corSymm(form = ~ time | ID), 
                       weights = varIdent(form = ~ 1 | Time), 
                       na.action=na.omit)
summary(fit_mild_int_mo)
# evidence of effects over time
wald(fit_mild_int_mo, 'spmmon')

# estimates of slopes and changes in slopes
L <- rbind(
  cbind(0, sc(spmmon, x = c(0,3,6,3,6), D=c(1,1,1,1,1), type = c(1,1,1,2,2) ), 0, 0*sc(spmmon, x = c(0,3,6,3,6), D=c(1,1,1,1,1), type = c(1,1,1,2,2) )),
  cbind(0, sc(spmmon, x = c(0,3,6,3,6), D=c(1,1,1,1,1), type = c(1,1,1,2,2) ), 1,sc(spmmon, x = c(0,3,6,3,6), D=c(1,1,1,1,1), type = c(1,1,1,2,2) ))
) 
  
ratenames <- c('rate/month W1 to M3','rate/month M3 to M6', 'rate/month M6 to M12', 
                 'change in rate/month at M3', 'change in rate/month at M6')
rownames(L) <- paste(rep(c('Adult','Pediatric'), each = length(ratenames)), ratenames)

wald(fit_mild_int_mo, L)

L_levels <-
  cbind(1, sc(spmmon, x = c(0,3,6,12,0,3,6,12), D = 0), rep(0:1, each = 4))
L_levels <- cbind(L_levels, L_levels[,5] * L_levels[,2:4])


rownames(L_levels) <- paste0(rep(c('Adult','Pediatric'), each = 4), ' ', c('Week 1','Month 3','Month 6','Month 12'))
#'
#' Estimated levels -- additive model - mild
#'
wald(fit_mild_int_mo, L_levels)
#'
#' Estimated levels -- additive model - mild - graph
#'
ww <- waldf(fit_mild_int_mo, L_levels) 
ww$months <- c(0,3,6,12)
ww$Version <- rep(c('Adult','Pediatric'), each = 4)
xyplot(coef ~ months, ww, type = 'l', groups = Version,
       main = 'Mild to Severe',
       ylab = TeX('GOSEv2 ($\\pm$SE)'),
       xlim = c(0,12),
       auto.key = list(reverse.rows = T),
       fit = ww$coef,
       upper = ww$coef + ww$se,
       lower = ww$coef - ww$se,
       subscripts = T) +
  latticeExtra::glayer(panel.fit(...))

L_diffs <-
  rbind(
    'Adult M3 - W1' =     c(-1, 1,  0, 0,    0, 0, 0, 0),
    'Adult M6 - W1' =     c(-1, 0,  1, 0,    0, 0, 0, 0),
    'Adult M12 - W1' =    c(-1, 0, -0, 1,    0, 0, 0, 0),
    'Adult M6 - M3' =     c(0, -1,  1, 0,    0, 0, 0, 0),
    'Adult M12 - M3' =    c(0, -1,  0, 1,    0, 0, 0, 0),
    'Adult M12 - M6' =    c(0, 0,  -1, 1,    0, 0, 0, 0),
    'Pediatric M3 - W1' = c(0,0,0,0,  -1, 1, 0, 0),
    'Pediatric M6 - W1' = c(0,0,0,0,  -1, 0, 1, 0),
    'Pediatric M12 - W1' = c(0,0,0,0,  -1, 0, 0, 1),
    'Pediatric M6 - M3' = c(0,0,0,0,  0, -1, 1, 0),
    'Pediatric M12 - M3' = c(0,0,0,0,  0, -1, 0, 1),
    'Pediatric M12 - M6' = c(0,0,0,0,  0, 0, -1, 1),
    'Pediatric - Adult at W1' = c(-1, 0,0,0, 1, 0, 0, 0),
    'Pediatric - Adult at M3' = c(0, -1, 0,0,0, 1, 0, 0),
    'Pediatric - Adult at M6' = c(0, 0, -1, 0,0,0, 1, 0),
    'Pediatric - Adult at M12' = c(0, 0, 0, -1, 0,0,0, 1)
    
  ) %*% L_levels

#'
#' Levels:
#'
wald(fit_mild_int_mo, L_levels)
#'
#' Differences
#'
wald(fit_mild_int_mo, L_diffs)



  
pred_mild <- with(gose_clean_mild, pred.grid(time, Version))
pred_mild <- waldf(fit_mild_add, pred = pred_mild)
# 
# Q4GM what exactly is this reporting? Predicted values at each time, by Version?
# 
# A: right!
pred_mild

pred_mild_plot <- ggplot(
  pred_mild, aes(x = time, y = pred_mild$coef, fill = Version)) +
  geom_line() +
  geom_ribbon(aes(x = time, 
                  y = pred_mild$coef, 
                  ymax = (pred_mild$coef + pred_mild$se),
                  ymin = (pred_mild$coef - pred_mild$se),
                  alpha = 0.25)) +
  facet_grid(cols = vars(Version)) + 
  xlab('Time') + 
  ylab('GOSE') +
  ggtitle('GOSE MILD predicted trajectories - Side-by-side') +
  theme_hc()
pred_mild_plot + theme(legend.position = 'none')

pred_mild_plot2 <- ggplot(
  pred_mild, aes(x = time, y = pred_mild$coef, fill = Version)) +
  geom_line() +
  geom_ribbon(aes(x = time, 
                  y = pred_mild$coef, 
                  ymax = (pred_mild$coef + pred_mild$se),
                  ymin = (pred_mild$coef - pred_mild$se),
                  alpha = 0.10)) +
  #facet_grid(cols = vars(Version)) + 
  xlab('Time') + 
  ylab('GOSE') +
  ggtitle('GOSE MILD predicted trajectories - Overlaid') +
  theme_hc()
pred_mild_plot2
```

## Controlling for Age

Since ages differ considerably in the two groups (Adult vs Pediatric)
we centre age within groups to reduce collinearity and so estimates
of the effect of 'Version' are comparisons of children of average
age with adults of average age.

The following is only a promising start since there appears to
be possible quadratic effects on rates of recovery which need
to be explored more carefully.


```{r mild-int-age}
gose_clean_mild <-
  within(
    gose_clean_mild,
    {
      # treating 0 M as missing although there are no 12's for months
      age <- `Age Y` + ifelse(`Age M` == 0, .5, `Age M`/12)  
      age_mean <- capply(age, Version, mean, na.action = na.omit)
      age_cen <- age - age_mean
      age_cen2 <- age_cen^2
    }
  )

poly2 <- function(agec) gsp(agec, 0, 2, 2)
sp3 <- function(agec) gsp(agec, 0, 3, 2)    # cubic spline

fit_mild_int_mo_age2 <- update(fit_mild_int_mo, . ~ (. + Version) * poly2(age_cen))
fit_mild_int_mo_agesp3 <- update(fit_mild_int_mo, . ~ (. + Version) * sp3(age_cen))
fit_mild_int_mo_age2 %>% summary
fit_mild_int_mo_agesp3 %>% summary
anova(update(fit_mild_int_mo_age2, method = 'ML'),
      update(fit_mild_int_mo_agesp3, method = 'ML'))
xyplot(GOSEv2 ~ age_cen | Version,gose_clean_mild )
xyplot(jitter(as.numeric(factor(Version))) ~ age_cen, gose_clean_mild)
tab(gose_clean_mild$Version)
gose_clean_mild$age_cen
```

#### Test of quadratic term

```{r mild-int-age2}
#
# 
wald(fit_mild_int_mo_age2, 'cen.*D2') # illustrating convenience of grepping a term to get overall F test

# in Adults only:

ww <- wald(fit_mild_int_mo_age2, 'poly2.age_cen.D2')
print(ww)

# in Pediatric only:

L <- ww[[1]]$L  # L matrix 
Lped <- L[1:4,] + L[5:8,]

wald(fit_mild_int_mo_age2, Lped)

# Plot

xyplot(GOSEv2 ~ age_cen | Version,gose_clean_mild )


unique(gose_clean_mild$month)
pred <- with(gose_clean_mild, pred.grid(Version, month = c(0,3,6,12), age_cen = seq(-20,40,5)))
pred <- subset(pred, !(abs(age_cen) >10 & Version == 'Pediatric'))                                        
ww %>% head
ww <- waldf(fit_mild_int_mo_age2, pred = pred)
ww$time <- with(ww, age_cen +  month/12)
xyplot(coef ~ time | Version, ww, groups = age_cen, type = 'l', 
       layout = c(1,2))
ww$time
# FIX ABOVE
# 
gose_clean_mild$time
summary(fit_mild_int_mo_age2)
ga(summary.gls)

ga(summary.gls)

summary(fit_mild_int_mo_age2)
vif(fit_mild_int_mo_age2)

wald(fit_mild_int_mo_age2, ':.*:')

```

There are significant 3-way interaction in both the Adult
and in the Pediatric groups

## HERE ############

```{r mild-int-age2-contd}

tab(gose_clean_mild, ~ Version)


# 
# 
# 
                                        
```




Note that age has a stronger effect in the Pediatric group

```{r ct-status}

```



## MS TBI


```{r msTBI-gls}
# msTBI
fit_ms <- gls(GOSEv2 ~ spms(time) * Version, data = gose_clean_mstbi,
              correlation = corSymm(form = ~ time1 | ID), 
              weights = varIdent(form = ~ 1 | time), 
              na.action=na.omit)
summary(fit_ms)

# evidence of interaction with Version?
wald(fit_ms, ':')  # No

# additive model
fit_ms_add <- gls(GOSEv2 ~ spms(time) + Version, data = gose_clean_mstbi,
                  correlation = corSymm(form = ~ time1 | ID), 
                  weights = varIdent(form = ~ 1 | time), 
                  na.action=na.omit)
summary(fit_ms_add)
# Q4GM - should this be spms?
# 
# A: When the second argument to 'wald' is a character, it
# is used as a regular expression to grep into the
# names of coefficients of the model. So 'spms'
# would produce the same output in this case.
wald(fit_ms_add, 'spm')
  
# Estimates of slopes and changes in slopes
L <- cbind(0,sc(spms, x = c(2,3,3), D=c(1,1,1), type = c(1,1,2) ), 0)
rownames(L) <- c('Change from Month3toMonth6','Change from Month6toMonth12',  
                   'Change in slope at Month6')
wald(fit_ms_add, L)

# predictions for ms
pred_ms <- with(gose_clean_mstbi, pred.grid(time, Version))
pred_ms <- waldf(fit_ms_add, pred = pred_ms)
pred_ms

pred_ms_plot <- ggplot(
  pred_ms, aes(x = time, y = pred_ms$coef, fill = Version)) +
  geom_line() +
  geom_ribbon(aes(x = time, 
                  y = pred_ms$coef, 
                  ymax = (pred_ms$coef + pred_ms$se),
                  ymin = (pred_ms$coef - pred_ms$se),
                  alpha = 0.25)) +
  facet_grid(cols = vars(Version)) + 
  xlab('Time') + 
  ylab('GOSE') +
  ggtitle('GOSE MSTBI predicted trajectories - Side-by-side') +
  theme_hc()
pred_ms_plot + theme(legend.position = 'none')

pred_ms_plot2 <- ggplot(
  pred_ms, aes(x = time, y = pred_ms$coef, fill = Version)) +
  geom_line() +
  geom_ribbon(aes(x = time, 
                  y = pred_ms$coef, 
                  ymax = (pred_ms$coef + pred_ms$se),
                  ymin = (pred_ms$coef - pred_ms$se),
                  alpha = 0.25)) +
  #facet_grid(cols = vars(Version)) + 
  xlab('Time') + 
  ylab('GOSE') +
  ggtitle('GOSE MSTBI predicted trajectories - Overlaid') +
  theme_hc()
pred_ms_plot2 + theme(legend.position = 'none')  

# merged plot for all severities
pred_mild$Sev <- 'Mild'
pred_ms$Sev <- "msTBI"
pred_mild$L <- NULL
pred_ms$L <- NULL
pred_all <- full_join(pred_mild, pred_ms)
pred_all$Sev <- factor(pred_all$Sev)

pred_all_plot <- ggplot(
  pred_all, aes(x = time, y = pred_all$coef, fill = Version)) +
  geom_line() +
  geom_ribbon(aes(x = time, 
                  y = pred_all$coef, 
                  ymax = (pred_all$coef + pred_all$se),
                  ymin = (pred_all$coef - pred_all$se),
                  alpha = 0.25)) +
  facet_grid(cols = vars(Sev)) + 
  xlab('Time') + 
  ylab('GOSE') +
  ggtitle('GOSE predicted trajectories (all severity) - Side-by-side') +
  theme_hc()
pred_all_plot #+ theme(legend.position = 'none')  

```

# GOSE trajectory analysis{.tabset}

## Mild TBI trajectory analysis

The bootstrapping method in the `traj` package identified **1 trajectory group** among all mild cases. Below the one-group trajectory is plotted. The plots for 2, 3, and 4 groups (not plotted) show near overlap of 95% CI intervals about each trajectory (suggesting they are the same).

```{r gose traj analysis mild, echo=FALSE, warning=FALSE}

# creating a wide dataset
gose_clean_wide <- gose_clean %>%
  group_by(ID) %>%
  pivot_wider(id_cols = ID, names_from = Time, values_from = c(GOSEv2, Sev, Age, Sex, CT)) %>%
  rename("T1" = "GOSEv2_1", 
         "T2" = "GOSEv2_2",
         "T3" = "GOSEv2_3",
         "T4" = "GOSEv2_4")

gose_clean_wide <- gose_clean_wide %>%
  mutate(Sev = case_when(
          (Sev_1 == "Mild" | Sev_2 == "Mild" | Sev_3 == "Mild" | Sev_4 == "Mild") ~ "Mild",
          (Sev_1 == "msTBI" | Sev_2 == "msTBI" | Sev_3 == "msTBI" | Sev_4 == "msTBI") ~ "msTBI"),
         Sex = case_when(
           (Sex_1 == "Female" | Sex_2 == "Female" | Sex_3 == "Female" | Sex_4 == "Female") ~ "Female",
           (Sex_1 == "Male" | Sex_2 == "Male" | Sex_3 == "Male" | Sex_4 == "Male") ~ "Male"))

gose_clean_wide <- gose_clean_wide %>% 
  rowwise() %>% 
  mutate(Age = mean(c(Age_1, Age_2, Age_3, Age_4), na.rm = TRUE)) %>%
  mutate(CT = mean(c(CT_1, CT_2, CT_3, CT_4), na.rm = TRUE))

gose_clean_wide <- gose_clean_wide[, c(1:5, 27:30)]
gose_clean_wide[2:5] <- lapply(gose_clean_wide[2:5], function(col) as.numeric( gsub("NULL", "", col) ) )

# trajectory analysis
gose_traj_mild <- gose_clean_wide %>%
  select(1, 2:6) %>%
  subset(Sev %in% c("Mild")) %>%
  rename("GOSE_1" = "T1",
         "GOSE_2" = "T2",
         "GOSE_3" = "T3",
         "GOSE_4" = "T4") %>%
  mutate(T1 = 0,
         T2 = 3,
         T3 = 6,
         T4 = 12) 

#gose_traj_colorder <- c("ID", "GOSE_1", "GOSE_2", "GOSE_3", "GOSE_4", "T1", "T2", "T3", "T4")
#gose_traj <- gose_traj[, gose_traj_colorder]

gose_traj_mild <- gose_traj_mild[complete.cases(gose_traj_mild), ]
gose_traj_mild <- as.matrix(gose_traj_mild)
gose_traj_mild <- apply(gose_traj_mild, 2, as.numeric)

gose_traj_mild_data <- gose_traj_mild[, c(1, 2:5)]
gose_traj_mild_time <- gose_traj_mild[, c(1, 7:10)]

gose_traj_mild_list <- list(
  gt_mild_data = gose_traj_mild_data,
  gt_mild_time = gose_traj_mild_time)

gose_mild_s1 = traj::Step1Measures(gose_traj_mild_list$gt_mild_data, gose_traj_mild_list$gt_mild_time, ID = TRUE)
gose_mild_s2 = traj::Step2Selection(gose_mild_s1)
gose_mild_s3 = traj::Step3Clusters(gose_mild_s2, nclusters = NULL, K.max = 3)
plot(gose_mild_s3)

gose_traj_mild_clust <- as.data.frame(gose_mild_s3$partition)
gose_traj_mild_clust_group <- right_join(gose_clean, gose_traj_mild_clust, by = "ID")
gose_traj_mild_clust_group$Cluster <- as.factor(gose_traj_mild_clust_group$Cluster)

gose_traj_mild_fig <- gose_traj_mild_clust_group %>%
  #filter(Timepoint %in% c("Month 3", "Month 6", "Month 12")) %>%
  ggplot(aes(Timepoint, GOSEv2, color = Cluster, group = ID)) +
  geom_point(size = 0.005, alpha = 0.25) +
  geom_line(size= 0.005, position=position_jitter(w=0.1, h=0), alpha = 0.25) +
  geom_smooth(size = 2.5, se = TRUE, method = "loess", aes(group = Cluster)) +
  ggtitle("MILD trajectory analysis - 4 timepoint") +
  theme_hc()
gose_traj_mild_fig

```

\newpage

## Pediatric mild TBI trajectory analysis

The bootstrapping method in the `traj` package identified **3 trajectory groups** among all pediatric mild cases. Roughly, the groups are... "3-month improvers (n = 22)", "6-month improvers (n=10)", "high-level maintainers (n=30)".

```{r gose traj analysis ped mild, echo=FALSE, warning=FALSE}

# this follows the same steps as the `GOSE trajectory analysis` code chunk
gose_traj_mild_ped <- gose_clean_wide %>%
  select(1, 2:6, 8) %>%
  subset(Sev %in% c("Mild")) %>%
  subset(Age < 18) %>%
  rename("GOSE_1" = "T1",
         "GOSE_2" = "T2",
         "GOSE_3" = "T3",
         "GOSE_4" = "T4") %>%
  mutate(T1 = 0,
         T2 = 3,
         T3 = 6,
         T4 = 12) 

gose_traj_colorder <- c("ID", "GOSE_1", "GOSE_2", "GOSE_3", "GOSE_4", "T1", "T2", "T3", "T4")
gose_traj_mild_ped <- gose_traj_mild_ped[, gose_traj_colorder]

gose_traj_mild_ped <- gose_traj_mild_ped[complete.cases(gose_traj_mild_ped), ]
gose_traj_mild_ped <- as.matrix(gose_traj_mild_ped)
gose_traj_mild_ped <- apply(gose_traj_mild_ped, 2, as.numeric)

gose_traj_mild_ped_data <- gose_traj_mild_ped[, c(1, 2:5)]
gose_traj_mild_ped_time <- gose_traj_mild_ped[, c(1, 6:9)]

gose_traj_mild_ped_list <- list(
  gt_mild_ped_data = gose_traj_mild_ped_data,
  gt_mild_ped_time = gose_traj_mild_ped_time)

gose_mild_ped_s1 = traj::Step1Measures(gose_traj_mild_ped_list$gt_mild_ped_data, gose_traj_mild_ped_list$gt_mild_ped_time, ID = TRUE)
gose_mild_ped_s2 = traj::Step2Selection(gose_mild_ped_s1)
gose_mild_ped_s3 = traj::Step3Clusters(gose_mild_ped_s2, nclusters = NULL, K.max = 3)
plot(gose_mild_ped_s3)

gose_traj_mild_ped_clust <- as.data.frame(gose_mild_ped_s3$partition)
gose_traj_mild_ped_clust_group <- right_join(gose_clean, gose_traj_mild_ped_clust, by = "ID")
gose_traj_mild_ped_clust_group$Cluster <- as.factor(gose_traj_mild_ped_clust_group$Cluster)

gose_traj_mild_ped_fig <- gose_traj_mild_ped_clust_group %>%
  #filter(Timepoint %in% c("Month 3", "Month 6", "Month 12")) %>%
  ggplot(aes(Timepoint, GOSEv2, color = Cluster, group = ID)) +
  geom_point(size = 0.005, alpha = 0.25) +
  geom_line(size= 0.005, position=position_jitter(w=0.1, h=0), alpha = 0.25) +
  geom_smooth(size = 2.5, se = TRUE, method = "loess", aes(group = Cluster)) +
  ggtitle("PEDIATRIC MILD trajectory analysis - 4 timepoint") +
  theme_hc()
gose_traj_mild_ped_fig

```

\newpage

## Pediatric mild TBI trajectory group model

The multinomial logistic regression from the `nnet` package is used when the dependent variable is an unordered, categorical variable. Older **Age** was identified as a significant predictor of membership in **Group 2 ("6-month improvers")**, relative to the reference category of Group 1.

```{r ped mild tbi traj model}

# identifying predictors of group membership, based on the trajectory analyses above
gose_traj_mild_ped_clust_group$Cluster <- relevel(
  gose_traj_mild_ped_clust_group$Cluster, ref = "1"
)

# running the multinomial logistic regression
ped_traj_model <- multinom(Cluster ~ Age + Sex + CT, data = gose_traj_mild_ped_clust_group)
summary(ped_traj_model)
exp(coef(ped_traj_model))

# producing z and p scores
ped_traj_model_z <- summary(ped_traj_model)$coefficients / summary(ped_traj_model)$standard.errors
ped_traj_model_z
ped_traj_model_p <- (1 - pnorm(abs(ped_traj_model_z), 0, 1)) * 2
ped_traj_model_p

#gose_traj_mild_ped_clust_group %>%
#  group_by(Cluster) %>%
#  summarize(mean_age = mean(Age), mean_gose = mean(GOSEv2))

#gose_traj_mild_ped_clust_group %>%
#  group_by(Cluster, Timepoint) %>%
#  summarize(gose = mean(GOSEv2))

```

\newpage

## Adult mild TBI trajectory analysis

The bootstrapping method in the `traj` package identified **3 trajectory groups** among all adult mild cases. Roughly, the groups are... "3-month improvers (n=35)", "6-month improvers (n=12)", "high-level maintainers (n=41)".

```{r gose traj analysis adult mild, echo=FALSE, warning=FALSE}

# this follows the same steps as the `GOSE trajectory analysis` code chunk
gose_traj_mild_adult <- gose_clean_wide %>%
  select(1, 2:6, 8) %>%
  subset(Sev %in% c("Mild")) %>%
  subset(Age >= 18) %>%
  rename("GOSE_1" = "T1",
         "GOSE_2" = "T2",
         "GOSE_3" = "T3",
         "GOSE_4" = "T4") %>%
  mutate(T1 = 0,
         T2 = 3,
         T3 = 6,
         T4 = 12) 

gose_traj_colorder <- c("ID", "GOSE_1", "GOSE_2", "GOSE_3", "GOSE_4", "T1", "T2", "T3", "T4")
gose_traj_mild_adult <- gose_traj_mild_adult[, gose_traj_colorder]

gose_traj_mild_adult <- gose_traj_mild_adult[complete.cases(gose_traj_mild_adult), ]
gose_traj_mild_adult <- as.matrix(gose_traj_mild_adult)
gose_traj_mild_adult <- apply(gose_traj_mild_adult, 2, as.numeric)

gose_traj_mild_adult_data <- gose_traj_mild_adult[, c(1, 2:5)]
gose_traj_mild_adult_time <- gose_traj_mild_adult[, c(1, 6:9)]

gose_traj_mild_adult_list <- list(
  gt_mild_adult_data = gose_traj_mild_adult_data,
  gt_mild_adult_time = gose_traj_mild_adult_time)

gose_mild_adult_s1 = traj::Step1Measures(gose_traj_mild_adult_list$gt_mild_adult_data, gose_traj_mild_adult_list$gt_mild_adult_time, ID = TRUE)
gose_mild_adult_s2 = traj::Step2Selection(gose_mild_adult_s1)
gose_mild_adult_s3 = traj::Step3Clusters(gose_mild_adult_s2, nclusters = NULL)
# plot(gose_mild_adult_s3)                                                                  #  mod GM

gose_traj_mild_adult_clust <- as.data.frame(gose_mild_adult_s3$partition)
gose_traj_mild_adult_clust_group <- right_join(gose_clean, gose_traj_mild_adult_clust, by = "ID")
gose_traj_mild_adult_clust_group$Cluster <- as.factor(gose_traj_mild_adult_clust_group$Cluster)

gose_traj_mild_adult_fig <- gose_traj_mild_adult_clust_group %>%
  #filter(Timepoint %in% c("Month 3", "Month 6", "Month 12")) %>%
  ggplot(aes(Timepoint, GOSEv2, color = Cluster, group = ID)) +
  geom_point(size = 0.005, alpha = 0.25) +
  geom_line(size= 0.005, position=position_jitter(w=0.1, h=0), alpha = 0.5) +
  geom_smooth(size = 2.5, se = TRUE, method = "loess", aes(group = Cluster)) +
  ggtitle("ADULT MILD trajectory analysis - 4 timepoint") +
  theme_hc()
gose_traj_mild_adult_fig

```

\newpage

## Adult mild TBI trajectory group model

The multinomial logistic regression from the `nnet` package is used when the dependent variable is an unordered, categorical variable. **Male** sex trends towards being a predictor of membership in **Group 2 ("6-month improvers")**, relative to the reference category of Group 1/Female.

```{r adult mild tbi traj model}

# identifying predictors of group membership, based on the trajectory analyses above
gose_traj_mild_adult_clust_group$Cluster <- relevel(
  gose_traj_mild_adult_clust_group$Cluster, ref = "1"
)

# running the multinomial logistic regression
adult_traj_model <- multinom(Cluster ~ Age + Sex + CT, data = gose_traj_mild_adult_clust_group)
summary(adult_traj_model)
exp(coef(adult_traj_model))

# producing z and p scores
adult_traj_model_z <- summary(adult_traj_model)$coefficients / summary(adult_traj_model)$standard.errors
adult_traj_model_z
adult_traj_model_p <- (1 - pnorm(abs(adult_traj_model_z), 0, 1)) * 2
adult_traj_model_p

```

\newpage

## msTBI trajectory analysis

The bootstrapping method in the `traj` package identified **1 trajectory group** among all msTBI cases. Below the one-group trajectory is plotted. The plots for 2, 3, and 4 groups (not plotted) show near overlap of 95% CI intervals about each trajectory (suggesting they are the same).

```{r gose traj analysis sev, echo=FALSE, warning=FALSE}

# this follows the same steps as the `GOSE trajectory analysis` code chunk
gose_traj_sev <- gose_clean_wide %>%
  select(1, 2:6, 8) %>%
  subset(Sev %in% c("msTBI")) %>%
  rename(#"GOSE_1" = "T1",
         "GOSE_2" = "T2",
         "GOSE_3" = "T3",
         "GOSE_4" = "T4") %>%
  mutate(#T1 = 0,
         T2 = 3,
         T3 = 6,
         T4 = 12) 

gose_traj_colorder <- c("ID", "GOSE_2", "GOSE_3", "GOSE_4", "T2", "T3", "T4")
gose_traj_sev <- gose_traj_sev[, gose_traj_colorder]

gose_traj_sev <- gose_traj_sev[complete.cases(gose_traj_sev), ]
gose_traj_sev <- as.matrix(gose_traj_sev)
gose_traj_sev <- apply(gose_traj_sev, 2, as.numeric)

gose_traj_sev_data <- gose_traj_sev[, c(1, 2:4)]
gose_traj_sev_time <- gose_traj_sev[, c(1, 5:7)]

gose_traj_sev_list <- list(
  gt_sev_data = gose_traj_sev_data,
  gt_sev_time = gose_traj_sev_time)

gose_sev_s1 = traj::Step1Measures(gose_traj_sev_list$gt_sev_data, gose_traj_sev_list$gt_sev_time, ID = TRUE)
gose_sev_s2 = traj::Step2Selection(gose_sev_s1)
gose_sev_s3 = traj::Step3Clusters(gose_sev_s2, nclusters = NULL)
# plot(gose_sev_s3)                                                                        ##### mod GM

gose_traj_sev_clust <- as.data.frame(gose_sev_s3$partition)
gose_traj_sev_clust_group <- right_join(gose_clean, gose_traj_sev_clust, by = "ID")
gose_traj_sev_clust_group$Cluster <- as.factor(gose_traj_sev_clust_group$Cluster)

gose_traj_sev_fig <- gose_traj_sev_clust_group %>%
  filter(Timepoint %in% c("Month 3", "Month 6", "Month 12")) %>%
  ggplot(aes(Timepoint, GOSEv2, color = Cluster, group = ID)) +
  geom_point(size = 0.005, alpha = 0.25) +
  geom_line(size= 0.005, position=position_jitter(w=0.1, h=0), alpha = 0.25) +
  geom_smooth(size = 2.5, se = TRUE, method = "loess", aes(group = Cluster)) +
  ggtitle("msTBI trajectory analysis - 3 timepoint") +
  theme_hc()
gose_traj_sev_fig

```

\newpage

# Logistic models

## All patients - 7+ point cutoff

Using GOSE >= 7 as a cutoff for "favorable" status at 12-months, and using all available data, the logistic regression identified Severity and Age as significant predictors. The log odd coefficients were exponentiated and expressed as odds ratios. The odds of an **msTBI** patient (relative to mild TBI, holding Age and Sex and CT status constant) having a favorable outcome are **0.15 (0.08, 0.27)**. For every unit increase in **Age** (holding Severity and Sex constant), the odds of achieving a favorable outcome are **0.99 (0.98, 0.99)**. No Main effects of sex were observed. 

```{r logit all patients, echo=FALSE, warning=FALSE}

# creating the cutpoints for 'favorable' outcome
gose_logis <- gose_clean %>%
  mutate(cut4 = case_when(
    (GOSEv2 >= 4) ~ 1,
    (GOSEv2 < 4) ~ 0
  ))
gose_logis <- gose_logis %>%
  mutate(cut7 = case_when(
    (GOSEv2 >= 7) ~ 1,
    (GOSEv2 < 7) ~ 0
  ))
gose_logis <- gose_logis %>%
  mutate(cut8 = case_when(
    (GOSEv2 >= 8) ~ 1,
    (GOSEv2 < 8) ~ 0
  ))
gose_logis$cut4 <- as.factor(gose_logis$cut4)
gose_logis$cut7 <- as.factor(gose_logis$cut7)
gose_logis$cut8 <- as.factor(gose_logis$cut8)

# running the logistic regression
gose_logit_all <- glm(cut7 ~ Sev + Age + Sex + CT,
                      data = subset(gose_logis, Time == 4),
                      family = "binomial")
summary(gose_logit_all)
confint(gose_logit_all)
confint.default(gose_logit_all)
exp(coef(gose_logit_all))
exp(cbind(OR = coef(gose_logit_all), confint(gose_logit_all)))

```

\newpage

## All patients - Variable cutoff

Using GOSE >= 7 and GOSE >= 4 as a cutoff for "favorable" status at 12-months for *mild TBI* and **msTBI** patients, respectively. the logistic regression identified Age as a significant predictor. The log odd coefficients were exponentiated and expressed as odds ratios. For every unit increase in **Age** (holding Severity and Sex and CT status constant), the odds of achieving a favorable outcome are **0.98 (0.97, 0.99)**. No Main effects of sex were observed. 

```{r logit all patients var, echo=FALSE, warning=FALSE}

# creating the cutpoints for 'favorable' outcome
gose_logis_mild <- gose_clean %>%
  filter(Sev %in% "Mild") %>%
  mutate(cut = case_when(
    (GOSEv2 >= 7) ~ 1,
    (GOSEv2 < 7) ~ 0
  ))
gose_logis_sev <- gose_clean %>%
  filter(Sev %in% "msTBI") %>%
  mutate(cut = case_when(
    (GOSEv2 >= 4) ~ 1,
    (GOSEv2 < 4) ~ 0
  ))
gose_logis_var <- rbind(gose_logis_mild, gose_logis_sev)
gose_logis_var$cut <- as.factor(gose_logis_var$cut)

# running the logistic regressions
gose_logit_var <- glm(cut ~ Sev + Age + Sex + CT,
                      data = subset(gose_logis_var, Time == 4),
                      family = "binomial")
summary(gose_logit_var)
confint(gose_logit_var)
confint.default(gose_logit_var)
exp(coef(gose_logit_var))
exp(cbind(OR = coef(gose_logit_var), confint(gose_logit_var)))

```

\newpage

## Mild TBI - 7pt

Using GOSE >= 7 as a cutoff for "favorable" status at 12-months, and **only mild TBI data**, the logistic regression identified Age as a significant predictor. The log odd coefficients were exponentiated and expressed as odds ratios. For every unit increase in **Age** (holding Sex constant), the odds of achieving a favorable outcome are **0.98 (0.97, 0.99)**. No Main effects of sex were observed. 

```{r logit mild tbi 7, echo=FALSE, warning=FALSE}

# running the logistic regression for cut7
gose_logit_mild <- glm(cut7 ~ Age + Sex + CT,
                      data = subset(gose_logis, Time == 4 & Sev %in% "Mild"),
                      family = "binomial")
summary(gose_logit_mild)
confint(gose_logit_mild)
confint.default(gose_logit_mild)
exp(coef(gose_logit_mild))
exp(cbind(OR = coef(gose_logit_mild), confint(gose_logit_mild)))

```

## Mild TBI - 8pt

Using GOSE = 8 as a cutoff for "favorable" status at 12-months, and **only mild TBI data**, the logistic regression identified Age as a significant predictor. The log odd coefficients were exponentiated and expressed as odds ratios. For every unit increase in **Age** (holding Sex and CT status constant), the odds of achieving a favorable outcome are **0.98 (0.97, 0.99)**. No Main effects of sex were observed. 

```{r logit mild tbi, echo=FALSE, warning=FALSE}

# running the logistic regression for cut8
gose_logit_mild <- glm(cut8 ~ Age + Sex + CT,
                      data = subset(gose_logis, Time == 4 & Sev %in% "Mild"),
                      family = "binomial")
summary(gose_logit_mild)
confint(gose_logit_mild)
confint.default(gose_logit_mild)
exp(coef(gose_logit_mild))
exp(cbind(OR = coef(gose_logit_mild), confint(gose_logit_mild)))

```

\newpage

## msTBI - 4pt

Using GOSE >= 4 as a cutoff for "favorable" status at 12-months, and **only msTBI data**, the logistic regression identified Age as a significant predictor. The log odd coefficients were exponentiated and expressed as odds ratios. For every unit increase in **Age** (holding Sex and CT status constant), the odds of achieving a favorable outcome are **0.97 (0.94, 0.99)**. No Main effects of Sex were observed.

```{r logit msTBI 4pt, echo=FALSE, warning=FALSE}

# running the logistic regression for cut4
gose_logit_mstbi <- glm(cut4 ~ Age + Sex + CT,
                      data = subset(gose_logis, Time == 4 & Sev %in% "msTBI"),
                      family = "binomial")
summary(gose_logit_mstbi)
confint(gose_logit_mstbi)
confint.default(gose_logit_mstbi)
exp(coef(gose_logit_mstbi))
exp(cbind(OR = coef(gose_logit_mstbi), confint(gose_logit_mstbi)))

```

## msTBI - 7pt

Using GOSE >= 7 as a cutoff for "favorable" status at 12-months, and **only msTBI data**, the logistic regression did not identify significant predictors.

```{r logit msTBI 7pt, echo=FALSE, warning=FALSE}

# running the logistic regression for cut7
gose_logit_mstbi <- glm(cut7 ~ Age + Sex + CT,
                      data = subset(gose_logis, Time == 4 & Sev %in% "msTBI"),
                      family = "binomial")
summary(gose_logit_mstbi)
confint(gose_logit_mstbi)
confint.default(gose_logit_mstbi)
exp(coef(gose_logit_mstbi))
exp(cbind(OR = coef(gose_logit_mstbi), confint(gose_logit_mstbi)))

```

\newpage

## Pediatric patients - Variable cutoff

```{r logit pediatric var, echo=FALSE, warning=FALSE, run=FALSE, include=FALSE}

# running the logistic regression for cut
gose_logit_ped <- glm(cut ~ Sev + Sex + CT,
                      data = subset(gose_logis_var, Time == 4 & Age < 18),
                      family = "binomial")
summary(gose_logit_ped)
confint(gose_logit_ped)
confint.default(gose_logit_ped)
exp(coef(gose_logit_ped))
exp(cbind(OR = coef(gose_logit_ped), confint(gose_logit_ped)))

```

\newpage

## Adult patients - Variable cutoff

```{r logit adult var, echo=FALSE, warning=FALSE, inclue=FALSE, run=FALSE}

# running the logistic regression for cut
gose_logit_adult <- glm(cut ~ Sev + Age + Sex + CT,
                      data = subset(gose_logis_var, Time == 4 & Age >= 18),
                      family = "binomial")
summary(gose_logit_adult)
confint(gose_logit_adult)
confint.default(gose_logit_adult)
exp(coef(gose_logit_adult))
exp(cbind(OR = coef(gose_logit_adult), confint(gose_logit_adult)))

```

\newpage

# Asymptotic recovery models


## Mild TBI - Conditional recovery models

```{r asymp-explore}

dim(gose)
dim(up(gose, ~ id))
gose <- as.data.frame(gose)
head(gose)
gose$Time
names(gose)
gose$gose

gose$Time_f <- factor(gose$`Gose Timepoint`)
tab(gose, ~ Time_f + `Gose Timepoint`)
gose$Time <- tr(as.numeric(gose$Time_f), from = c(1,2,3,4), to = c(1,4,2,3))

xyplot(etime ~ Time, gose)

id1 <- gose[ with(gose, which(etime > 5 & Time ==1)), ]$id
subset(gose, id == id1)
```

