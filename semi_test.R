library(mstate)
library(dplyr)

# data in wide format, for transition 1 this is dataset E1 of
# Therneau & Grambsch (2000)
tg <- data.frame(illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
                 dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
                 x1=c(1,1,1,0,0,0),x2=c(6:1))


tmat <- transMat(list(c(2,3), c(4),c(), c()), c('Healthy', 'Relapse', 'NRM', 'DaR'))
# data in long format using msprep
tglong0 <- msprep(time=c(NA,"illt","dt","dt"),status=c(NA,"ills","ds","ds"),
                  data=tg,keep=c("x1","x2"),trans=tmat)
# expanded covariates
tglong <- expand.covs(tglong0,c("x1","x2"))

tglong$Tstart <- tglong$Tstart*365.241
tglong$Tstop <- tglong$Tstop*365.241
tglong$time <- tglong$time*365.241
tglong$age <- 65*365.241
tglong$sex <- 'male'
tglong$year <- as.Date('2010-01-01')
tglong <- rbind(tglong,tglong,tglong,tglong,tglong)


# Cox model with different covariate
cx <- coxph(Surv(Tstart,Tstop,status)~x1.1+x2.2+strata(trans),
            data=tglong,method="breslow")
cx
summary(cx)

cx2 <- coxph.relsurv(Surv(Tstart,Tstop,status)~x1.1+x2.2+x2.3+strata(trans),
                     data=tglong %>%
                       mutate(Tstop=Tstop+runif(nrow(.)), 
                              x1.2=x1.2+runif(nrow(.)), 
                              x1.1=x1.1+runif(nrow(.))), 
                     split.transitions = c(2:3),
                     rmap = list(age=age))
cx2
summary(cx2)

# new data, to check whether results are the same for transition 1 as
# those in appendix E.1 of Therneau & Grambsch (2000)
newdata <- data.frame(trans=1:3,x1.1=c(0,0,0),x2.2=c(0,1,0),x2.3=c(0,0,1),strata=1:3, 
                      age=65*365.241, sex='male', year=as.Date('2010-01-01'))
mod <- msfit(cx,newdata,trans=tmat)
mod

mod_rs <- msfit.coxph.relsurv(cx2,newdata = newdata, trans = tmat)









pt <- probtrans(mod,predt=0)



# Non-parametric:
tmat
msprep
coxph()
msfit()
msfit.relsurv(split.transitions)
probtrans

# Semi-parametric:
tmat
msprep
coxph.relsurv # not coxph
msfit.relsurv
probtrans

cx.relsurv <- coxph.relsurv(Surv(Tstart,Tstop,status)~x1.1+x2.2+strata(trans),
                            data=tglong, # method="breslow", - leave this in ... in coxph
                            # relsurv steps:
                            split.transitions=2, ratetable, time.format, rmap, variance,
                            # rsadd arguments: use EM by default: na.action by default, 
                            # init, bwin, centered, cause. 
                            
)
summary(cx.relsurv)
# Output: coxph object with a corresponding structure like coxph.
# Focus on what msfit needs: coefficients, baseline hazards, strata,
# SEs (maybe not necessary), , etc.

# coxph object. Using coxph for non-split and rsadd for split hazard

# Option 1:
mod <- msfit.relsurv(
  cx.relsurv, # not mod
  data, # maybe we do not need it
  newdata, # additional to data
  trans
  # here do bootstrap
)
pt <- probtrans(mod, predt = 0)















# 


library(relsurv)
fit <- rsadd(Surv(time,cens)~sex+ratetable(age=age*365.241),
             ratetable=slopop,data=rdata, method = 'EM')
summary(fit)
