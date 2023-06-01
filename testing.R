library(mstate)

# data in wide format, for transition 1 this is dataset E1 of
# Therneau & Grambsch (2000)
tg <- data.frame(illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
                 dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
                 x1=c(1,1,1,0,0,0),x2=c(6:1))


tmat <- trans.illdeath()
# data in long format using msprep
tglong <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
                 data=tg,keep=c("x1","x2"),trans=tmat)
# expanded covariates
tglong <- expand.covs(tglong,c("x1","x2"))
# Cox model with different covariate
cx <- coxph(Surv(Tstart,Tstop,status)~x1.1+x2.2+strata(trans),
            data=tglong,method="breslow")
summary(cx)
# new data, to check whether results are the same for transition 1 as
# those in appendix E.1 of Therneau & Grambsch (2000)
newdata <- data.frame(trans=1:3,x1.1=c(0,0,0),x2.2=c(0,1,0),strata=1:3)
mod <- msfit(cx,newdata,trans=tmat)

# Option 1:
mod <- msfit.relsurv(
  cx, # not mod
  data,
  newdata, # additional to data
  trans,
  # relsurv arguments:
  split.transitions, 
  ratetable, 
  rmap, 
  time.format, 
  variance
)

# Option 2:
cx.relsurv <- coxph.relsurv(
  cx, 
  # relsurv arguments:
  data,
  split.transitions, 
  ratetable, 
  rmap, 
  time.format, 
  variance)

msfit.relsurv(
  cx.relsurv, 
  newdata,
  trans
)

# Output:
[1] Haz
[2] varHaz
[3] trans

# 


library(relsurv)
fit <- rsadd(Surv(time,cens)~sex+ratetable(age=age*365.241),
             ratetable=slopop,data=rdata,, method = 'EM')
summary(fit)