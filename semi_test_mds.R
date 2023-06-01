# install.packages('C:/Users/dmanevski/Dropbox (MF Uni LJ)/Damjan Manevski/Research/relsurv/relsurv_2.2-9.tar.gz', repos=NULL, type='source')
# install.packages('C:/Users/dame_/Dropbox (MF Uni LJ)/Damjan Manevski/Research/relsurv/relsurv_2.2-9.tar.gz', repos=NULL, type='source')

install.packages('C:/Users/dmanevski/Documents/GitHub/mstate_0.3.2.tar.gz', repos=NULL, type='source')


# Packages:
library(survival)
# library(relsurv)
library(dplyr)
library(mstate)

# Load anonymised data:
j <- read.table("C:/Users/dmanevski/Dropbox (MF Uni LJ)/Damjan Manevski/Research/2019 Liesbeth/MDS_Damjan/anonymised_mds.txt",header = TRUE)
# j <- read.table("C:/Users/dame_/Dropbox (MF Uni LJ)/Damjan Manevski/Research/2019 Liesbeth/MDS_Damjan/anonymised_mds.txt",header = TRUE)
load("C:/Users/dmanevski/Dropbox (MF Uni LJ)/Damjan Manevski/Research/2019 Liesbeth/MDS_Damjan/joinpoptab.RData")

colnames(j)[8] <- 'year'
# j <- j[j$ID %in% 1:300,]
tmat <- mstate::transMat(list(c(2,3), c(4), c(), c()),
                         c('Alive relapse-free','Relapse','NRM','DaR'))

df0 <- mstate::msprep(time=c(NA, 'rel_time', 'death_time', 'death_time'),
                                  status=c(NA, 'rel_status', 'death_status', 'death_status'),
                                  trans=tmat,
                                  data=j,
                                  keep=c('age', 'sex', 'year', 'country'))
# df0 %>% slice(1:20)


# expanded covariates
df <- expand.covs(df0,c("age","sex"))

df$Tstart <- df$Tstart*365.241
df$Tstop <- df$Tstop*365.241
df$time <- df$time*365.241
df$year <- as.Date(df$year)
df$country <- ifelse(df$country=='Netherlands', 'Netherlands, The', 'Denmark')

cause_vec <- ifelse(df$status==0, 2, sample(c(0,1,2), size=nrow(df %>% filter(status==1)), replace = T))

############## #

# Cox model with different covariate
cx <- coxph(Surv(Tstart,Tstop,status)~age.1+age.2+sexmale.2+age.3+strata(trans),
            data=df,method="breslow")
summary(cx)

cx2 <- coxph.relsurv(Surv(Tstart,Tstop,status)~age.1+age.2+sexmale.2+age.3+strata(trans),
                     data=df, split.transitions = 2:3, ratetable = joinpoptab,
                     rmap = list(age=age*365.241), init=0,na.action = 'na.omit', bwin=-1,
                     centered = FALSE, cause=cause_vec)
cx2
print(cx2)
summary(cx2)
print(summary(cx2))

cx2 <- coxph.relsurv(Surv(Tstart,Tstop,status)~age.2+sexmale.2+age.3+strata(trans),
                     data=df, split.transitions = 2:3, ratetable = joinpoptab,
                     rmap = list(age=age*365.241))
cx2
print(cx2)
summary(cx2)
print(summary(cx2))

############# #

newdata <- data.frame(trans=1:3,age.1=c(65,0,0),age.2=c(0,65,0),sexmale.2=c(0,1,0), age.3=c(0,0,65),strata=1:3, 
                      age=65, sex='male', year=as.Date('2010-01-01'), country='Denmark')

cx <- coxph(Surv(Tstart,Tstop,status)~age.1+age.2+sexmale.2+age.3+strata(trans),
                     data=df)
mod <- msfit(cx,newdata,trans=tmat)

cx2 <- coxph.relsurv(Surv(Tstart,Tstop,status)~age.1+age.2+sexmale.2+age.3+strata(trans),
                     data=df, split.transitions = 2:3, ratetable = joinpoptab,
                     rmap = list(age=age*365.241))
mod_rs <- msfit.coxph.relsurv(cx2,newdata = newdata, trans = tmat)

library(ggplot2)
ggplot(mod$Haz %>% mutate(trans=factor(trans)))+
  geom_step(aes(time, Haz, color=trans, group=trans))

ggplot(mod_rs$Haz %>% mutate(trans=factor(trans)))+
  geom_step(aes(time, Haz, color=trans, group=trans))


