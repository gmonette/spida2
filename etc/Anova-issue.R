library("car")
Moore <- with(Moore, Moore[fcategory != "low" | partner.status != "low",], )
xtabs(~ fcategory + partner.status, data=Moore)
mod <- lm(conformity ~ fcategory*partner.status, data=Moore)
Anova(mod)
car:::Anova.default(mod) # correct!

Data <- structure(list(Practice = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 1L), .Label = c("conventional", "organic"
), class = "factor"), crop_last = structure(c(3L, 3L, 3L, 7L, 7L, 7L, 
6L, 6L, 6L, 6L, 6L, 6L, 2L, 2L, 2L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 
7L, 7L, 7L, 7L, 6L, 6L, 6L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 2L, 2L, 
2L, 6L, 6L, 6L, 1L, 1L, 1L, 7L, 7L, 7L, 7L, 5L, 2L, 2L, 2L, 7L, 7L, 
7L, 7L, 7L, 7L, 2L, 2L, 8L, 8L, 8L, 8L, 2L, 2L, 7L, 7L, 3L, 3L, 2L, 
2L, 4L, 4L, 2L), .Label = c("av", "bet", "chicoree", "fourrage", "Lin", "ma", "po", "pois", "ww"), class = "factor"),
DM = c(2.34, 3.5775, 2.7, 10.9125, 9.3825, 9.9, 6.84, 7.83,
8.0325, 4.725, 2.7, 5.4675, 9.945, 11.295, 12.465, 2.66,
2.57, 2.36, 2.27, 1.93, 2.09, 2.21, 3.6833, 5.28, 5.6133,
5.7533, 4.78, 2.39, 1.8767, 2.6, 1.5933, 1.8767, 2.14, 1.8067,
1.8567, 1.9967, 2.3333, 2.22, 6.1967, 6.52, 4.97, 2.1667,
2.4, 2.1567, 2.7967, 2.99, 2.5667, 6.4567, 8.4567, 5.5667,
6.14, 5.1667, 3.6867, 4.13, 2.8467, 4.2267, 5.4833, 4.24,
4.2533, 4.46, 5.5167, 4.74, 3.7333, 6.1367, 7.2733, 4.92,
4.71, 2.2067, 2.4067, 9.2867, 9.18, 3.6767, 5.5233, 7.0333,
5.8233, 5.4733, 2.2867, 2.9633)), .Names = c("Practice", 
"crop_last", "DM"), row.names = c("1", "2", "3", "4", "5", "6", "7", 
"8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", 
"20", "21", "22", "28", "29", "30", "31", "32", "33", "34", "35", 
 "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", 
 "53", "54", "55", "56", "57", "58", "62", "63", "64", "65", "66", 
 "67", "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", 
 "84", "85", "86", "87", "88", "89", "90", "91", "92", "93", "94", 
 "95", "96", "97", "98"), class = "data.frame")

Data

xtabs(~ Practice + crop_last, data=Data) 

mod <- lm(DM ~ Practice + crop_last, data=Data)
summary(mod)
drop1(mod) 
Anova(mod) # correct
car:::Anova.default(mod) # wrong
car:::Anova.lm(mod)

vcov(mod)
library(spida2)
wald(mod)
wald(mod, -1)
wald(mod,'crop')
wald(mod, 'Practice')
