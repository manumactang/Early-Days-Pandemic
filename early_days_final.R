

library(ggplot2)
library(lme4)
library(moderndive)
library(plyr)
library(dplyr)
library(lubridate)
library(psych)
library(lmerTest)
library(ggeffects)
library(sjPlot) # table functions
library(sjmisc)



rm(list=ls())


###Datensatz aus Dropbox laden
load("C:/Users/Manuel/Nextcloud/Jochen-Manuel/covid/analysis/all_all2.Rda")

#load("D:/Dropbox/Forschung/covid/analysis/all_all2.Rda")

##set off scientific notation
options(scipen = 999)

#### only cases with test data 
df <- df %>%  #### <---- comment out if all cases are wanted
  filter(is.na(tests_per_day) == FALSE)

df <- df[c("deaths" ,"date" , "cases" , "tests_per_day" , "gdp_pc" , 
               "gwth_5a" , "hbeds1000",  "h_ex5a" ,
               "age_15_64" , "age_65"  , "netmig"  , "touri" ,
               "enrol_snd" , "pop", "pop_dns", "code", "dsi")]

#rescale
df$pop <- df$pop/1000000
df$touri <- df$touri/1000000
df$netmig <- df$netmig/100000
df$land <- df$land/1000000
df$gdp_pc <- df$gdp_pc/10000
df$pop_dns <- df$pop_dns/1000
df$cases <- df$cases/1000
df$tests_per_day <- df$tests_per_day/1000


df$emp <- df$emp/100
df$age_15_64 <- df$age_15_64/100
df$age_65 <- df$age_65/100
df$imp_5a <- df$imp_5a/100
df$enrol_snd <- df$enrol_snd/100
df$h_ex5a <- df$h_ex5a/100



#####

#get last day of infection
df1 <- df %>%  
  group_by(code) %>% 
  filter(cases != 0 ) %>% 
  arrange(desc(date)) %>% 
  filter(row_number()==1) %>% 
  mutate(last_infect = date) %>% 
  ungroup


df1 <- df1[c("code", "last_infect")]  

df <- merge(df, df1, by="code")

df$dsi_last <- as.numeric(df$date - df$last_infect)

summary(df$dsi_last)
######################### +4 kann es nicht geben
summary(df$last_infect)
#############################################################################
last_infect_ww <- "2020-04-02"
last_infect_ww <- ymd(last_infect_ww)
df$dsi_lastr <- as.numeric(df$date - last_infect_ww)
summary(df$dsi_lastr)
#############################################################################



###days since first reported infection world wide 
first_infect_ww <- "2019-12-31"
first_infect_ww <- ymd(first_infect_ww)
df$dsfi <- as.numeric(df$date - first_infect_ww)


df <- na.omit(df)

#####grand mean centering IVs
df$gdp_pc <- df$gdp_pc-mean(df$gdp_pc, na.rm=TRUE)
df$h_ex5a <- df$h_ex5a-mean(df$h_ex5a, na.rm=TRUE)
df$hbeds1000 <- df$hbeds1000-mean(df$hbeds1000, na.rm=TRUE)
df$pop_dns <- df$pop_dns-mean(df$pop_dns, na.rm=TRUE)
df$touri <- df$touri-mean(df$touri, na.rm=TRUE)
df$age_15_64 <- df$age_15_64-mean(df$age_15_64, na.rm=TRUE)
df$age_65 <- df$age_65-mean(df$age_65, na.rm=TRUE)
df$pop <- df$pop-mean(df$pop, na.rm=TRUE)
df$gwth_5a <- df$gwth_5a-mean(df$gwth_5a, na.rm=TRUE)
df$netmig <- df$netmig-mean(df$netmig, na.rm=TRUE)
df$tests_per_day <- df$tests_per_day - mean(df$tests_per_day, na.rm=TRUE)
df$emp <- df$emp - mean(df$emp, na.rm=TRUE)
df$enrol_snd <- df$enrol_snd - mean(df$enrol_snd, na.rm=TRUE)



#class(df$dsfi)
###group mean centering cases - IV
agg_case <- aggregate(df$cases, list(df$code), FUN=mean, na.rm=TRUE)
names(agg_case) <- c("code", "agg_cases")
df <- merge(df, agg_case, by="code")


df$case_gmc <- df$agg_cases - mean(df$cases, na.rm=TRUE )


#######################################

##code still character
df$code <- as.factor(df$code)

unique(df$code)

df$dsi_lastlog <- (-1 ) * log10((df$dsi_lastr*-1)+1)
summary(df$dsi_lastlog)


#A unconditional model
n1 <- lmer(deaths ~ 1 + (1 | code), data = df)
summary(n1)
icc <- 1208   /(1208   + 6631)
icc
#ICC 0.15

#B unconditional growth model
n2 <- lmer(deaths ~ dsi_lastlog + (1 | code) + (0 + dsi_lastlog | code), data = df)
summary(n2)


df$dsi_lastlog <- (-1 ) * log10((df$dsi_lastr*-1)+1)
summary(df$dsi_lastlog)


##THIS IS THE FINAL ONE#######
########################################################################################################## ! FINAL FINAL MODEL

####demography#######

n4etestlog_1 <- lmer(deaths ~ dsi_lastlog + case_gmc + tests_per_day + gdp_pc +  
                       gwth_5a + 
                       hbeds1000 + h_ex5a +
                       age_15_64 + age_65  + 
                       netmig + 
                       touri +
                       enrol_snd + 
                       pop + 
                       pop_dns + 
                       dsi_lastlog:pop +
                       dsi_lastlog:pop_dns +
                       dsi_lastlog:enrol_snd + 
                       dsi_lastlog:case_gmc +
                       dsi_lastlog:age_65  + 
                       dsi_lastlog:age_15_64 +
                       (1 | code) + (0 + dsi_lastlog | code), data = df)
summary(n4etestlog_1)

dem <- as.data.frame(coef(summary(n4etestlog_1)))
#round
dem <- round(dem, 3)
#flag number of stars according to sig level
dem$star <- NA 
dem$star[dem$`Pr(>|t|)` > 0.01 & dem$`Pr(>|t|)` <= 0.05] <- 1
dem$star[dem$`Pr(>|t|)` > 0.001 & dem$`Pr(>|t|)` <= 0.01] <- 2
dem$star[dem$`Pr(>|t|)` <= 0.001] <- 3

##Fits
AIC(n4etestlog_1)
BIC(n4etestlog_1)
logLik(n4etestlog_1)



### Economic#####

n4etestlog_2 <- lmer(deaths ~ dsi_lastlog + case_gmc + tests_per_day + gdp_pc +  
                     gwth_5a + 
                     hbeds1000 + h_ex5a +
                     age_15_64 + age_65  + 
                     netmig + 
                     touri +
                     enrol_snd +
                     pop + 
                     pop_dns + 
                       dsi_lastlog:pop +
                       dsi_lastlog:pop_dns +
                       dsi_lastlog:enrol_snd + 
                       dsi_lastlog:age_65  + 
                       dsi_lastlog:age_15_64 + 
                       dsi_lastlog:case_gmc +
                     dsi_lastlog:gdp_pc + 
                     dsi_lastlog:gwth_5a + 
                    
                     (1 | code) + (0 + dsi_lastlog | code), data = df)
summary(n4etestlog_2)

econ <- as.data.frame(coef(summary(n4etestlog_2)))
econ <- round(econ, 3)
econ$star <- NA 
econ$star[econ$`Pr(>|t|)` > 0.01 & econ$`Pr(>|t|)` <= 0.05] <- 1
econ$star[econ$`Pr(>|t|)` > 0.001 & econ$`Pr(>|t|)` <= 0.01] <- 2
a$star[econ$`Pr(>|t|)` <= 0.001] <- 3


AIC(n4etestlog_2)
BIC(n4etestlog_2)
logLik(n4etestlog_2)


####health#######

n4etestlog_3 <- lmer(deaths ~ dsi_lastlog + case_gmc + tests_per_day + gdp_pc +  
                       gwth_5a + 
                       hbeds1000 + h_ex5a + 
                       age_15_64 + age_65  + 
                       netmig + 
                       touri +
                       enrol_snd +
                       pop + 
                       pop_dns + 
                       dsi_lastlog:pop +
                       dsi_lastlog:pop_dns +
                       dsi_lastlog:enrol_snd + 
                       dsi_lastlog:age_65  + 
                       dsi_lastlog:age_15_64 +
                       dsi_lastlog:case_gmc +
                       dsi_lastlog:hbeds1000 + 
                       dsi_lastlog:h_ex5a +
                       
                       (1 | code) + (0 + dsi_lastlog | code), data = df)
summary(n4etestlog_3)

hth <- as.data.frame(coef(summary(n4etestlog_3)))
hth <- round(hth, 3)
hth$star <- NA 
hth$star[hth$`Pr(>|t|)` > 0.01 & hth$`Pr(>|t|)` <= 0.05] <- 1
hth$star[hth$`Pr(>|t|)` > 0.001 & hth$`Pr(>|t|)` <= 0.01] <- 2
hth$star[hth$`Pr(>|t|)` <= 0.001] <- 3


AIC(n4etestlog_3)
BIC(n4etestlog_3)
logLik(n4etestlog_3)



#####globalization########
  
  n4etestlog_4 <- lmer(deaths ~ dsi_lastlog + case_gmc + tests_per_day + gdp_pc +  
                         gwth_5a + 
                         hbeds1000 + h_ex5a +
                         age_15_64 + age_65  + 
                         netmig + 
                         touri +
                         enrol_snd +
                         pop + 
                         pop_dns +  
                         dsi_lastlog:pop +
                         dsi_lastlog:pop_dns +
                         dsi_lastlog:enrol_snd + 
                         dsi_lastlog:age_65  + 
                         dsi_lastlog:age_15_64 +
                         dsi_lastlog:case_gmc +
                         dsi_lastlog:netmig +
                         dsi_lastlog:touri +
                         
                         (1 | code) + (0 + dsi_lastlog | code), data = df)
  summary(n4etestlog_4)
  
  
  glob <- as.data.frame(coef(summary(n4etestlog_4)))
  glob <- round(glob, 3)
  glob$star <- NA 
  glob$star[glob$`Pr(>|t|)` > 0.01 & glob$`Pr(>|t|)` <= 0.05] <- 1
  glob$star[glob$`Pr(>|t|)` > 0.001 & glob$`Pr(>|t|)` <= 0.01] <- 2
  glob$star[ glob$`Pr(>|t|)` <= 0.001] <- 3
  
  
  AIC(n4etestlog_4)
  BIC(n4etestlog_4)
  logLik(n4etestlog_4)



#####final alltogether#######



n4etestlog <- lmer(deaths ~ dsi_lastlog + case_gmc + tests_per_day + gdp_pc +  
                     gwth_5a + 
                     hbeds1000 + h_ex5a + 
                     age_15_64 + age_65  + 
                     netmig + 
                     touri +
                     enrol_snd +
                     pop + 
                     pop_dns +  
                     dsi_lastlog:gdp_pc + 
                     dsi_lastlog:gwth_5a + 
                     dsi_lastlog:touri + 
                     dsi_lastlog:hbeds1000 + 
                     dsi_lastlog:enrol_snd + 
                     dsi_lastlog:h_ex5a +
                     dsi_lastlog:case_gmc +
                     dsi_lastlog:age_65  + 
                     dsi_lastlog:age_15_64  + 
                     dsi_lastlog:netmig +
                     dsi_lastlog:pop +
                     dsi_lastlog:pop_dns +
                     (1 | code) + (0 + dsi_lastlog | code), data = df)
summary(n4etestlog)


a <- as.data.frame(coef(summary(n4etestlog)))
a <- round(a, 3)
a$star <- NA 
a$star[a$`Pr(>|t|)` > 0.01 & a$`Pr(>|t|)` <= 0.05] <- 1
a$star[a$`Pr(>|t|)` > 0.001 & a$`Pr(>|t|)` <= 0.01] <- 2
a$star[ a$`Pr(>|t|)` <= 0.001] <- 3


AIC(n4etestlog)
BIC(n4etestlog)
logLik(n4etestlog)


anova(n4etestlog, n2)


##### correlation matrix
vars <- df[c("dsi_last","case_gmc" , "tests_per_day" , "gdp_pc" ,  
             "gwth_5a" ,
             "hbeds1000" , "h_ex5a"  , 
             "age_15_64" , "age_65"  , "netmig" , "touri" ,
             "pop" ,
             "enrol_snd" ,
             "pop_dns")]

vars <- vars[vars$dsi_last==0,]
vars$dsi_last <- NULL
###
cor(vars)
##########################

##### Interaction Effect Plots##########################


#### Hospital Beds
mydf <- ggpredict(n4etestlog, terms = c("dsi_lastlog", "hbeds1000"))


ggplot(mydf, aes(x, predicted, shape = group, color=group)) + 
  geom_pointrange(aes(ymin=conf.low, ymax=conf.high, color=group),
                  position=position_dodge(0.2)) +
  geom_line(position=position_dodge(0.2)) +
  labs(title= "Predicted Values of Death by COVID-19 Infection",
       subtitle="Interaction Effect with hospital beds per 1000 inhabitants 
incl. 95% Confidence Interval", x = "LogDays",
       y= "Predicted Deaths") + 
  scale_color_discrete(name= "Hospital beds",
                       labels= c("- 1 SD", "Mean", "+ 1 SD")) + 
  scale_shape_discrete(name= "Hospital beds",
                       labels= c("- 1 SD", "Mean", "+ 1 SD")) + 
  
  theme_bw(base_size=12)



#### Economic Growth
mydf1 <- ggpredict(n4etestlog, terms = c("dsi_lastlog", "gwth_5a"))

ggplot(mydf1, aes(x, predicted, shape = group, color=group)) + 
  geom_pointrange(aes(ymin=conf.low, ymax=conf.high, color=group),
                  position=position_dodge(0.2)) + 
  geom_line(position=position_dodge(0.2)) +
  geom_point(aes(shape=group, color=group), position=position_dodge(0.2)) +
  labs(title= "Predicted Values of Death by COVID-19 Infection",
       subtitle="Interaction Effect with 5-year average economic growth 
incl. 95% Confidence Interval",
       x = "LogDays",
       y= "Predicted Deaths") + 
  scale_color_discrete(name= "Economic growth 
(as % of GDP)",
                       labels= c("- 1 SD", "Mean", "+ 1 SD")) + 
  scale_shape_discrete(name= "Economic growth 
(as % of GDP)",
                       labels= c("- 1 SD", "Mean", "+ 1 SD")) + 
  theme_bw(base_size=12)


#### Globalization
mydf2 <- ggpredict(n4etestlog, terms = c("dsi_lastlog", "touri"))

ggplot(mydf2, aes(x, predicted, shape = group, color=group)) + 
  geom_pointrange(aes(ymin=conf.low, ymax=conf.high, color=group),
                  position=position_dodge(0.2)) + 
  geom_line(position=position_dodge(0.2)) +
  geom_point(aes(shape=group, color=group), position=position_dodge(0.2)) +
  labs(title= "Predicted Values of Death by COVID-19 Infection",
       subtitle="Interaction Effect with total amount of annual tourists
incl. 95% Confidence Interval",
       x = "LogDays",
       y= "Predicted Deaths") + 
  scale_color_discrete(name= "Annual tourists",
                       labels= c("- 1 SD", "Mean", "+ 1 SD")) + 
  scale_shape_discrete(name= "Annual tourists",
                       labels= c("- 1 SD", "Mean", "+ 1 SD")) +
  theme_bw(base_size=12)


#### Population Age
mydf3 <- ggpredict(n4etestlog, terms = c("dsi_lastlog", "age_65"))

ggplot(mydf3, aes(x, predicted, shape = group, color=group)) + 
  geom_pointrange(aes(ymin=conf.low, ymax=conf.high, color=group),
                  position=position_dodge(0.2)) + 
  geom_line(position=position_dodge(0.2)) +
  geom_point(aes(shape=group, color=group), position=position_dodge(0.2)) +
  labs(title= "Predicted Values of Death by COVID-19 Infection",
       subtitle="Interaction Effect with % of people over the age of 65
incl. 95% Confidence Interval",
       x = "LogDays",
       y= "Predicted Deaths") + 
  scale_color_discrete(name= "People over 
the age of 65",
                       labels= c("- 1 SD", "Mean", "+ 1 SD")) + 
  scale_shape_discrete(name= "People over 
the age of 65",
                       labels= c("- 1 SD", "Mean", "+ 1 SD")) +
  theme_bw(base_size=12)


#### Health expenditure
mydf4 <- ggpredict(n4etestlog, terms = c("dsi_lastlog", "h_ex5a"))

ggplot(mydf4, aes(x, predicted, shape = group, color=group)) + 
  geom_pointrange(aes(ymin=conf.low, ymax=conf.high, color=group),
                  position=position_dodge(0.2)) + 
  geom_line(position=position_dodge(0.2)) +
  geom_point(aes(shape=group, color=group), position=position_dodge(0.2)) +
  labs(title= "Predicted Values of Death by COVID-19 Infection",
       subtitle="Interaction Effect with health expenditure as % of GDP
incl. 95% Confidence Interval",
       x = "LogDays",
       y= "Predicted Deaths") + 
  scale_color_discrete(name= "Health Expenditure",
                       labels= c("- 1 SD", "Mean", "+ 1 SD")) + 
  scale_shape_discrete(name= "Health Expenditure",
                       labels= c("- 1 SD", "Mean", "+ 1 SD")) +
  theme_bw(base_size=12)
















