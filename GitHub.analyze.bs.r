# This code is provided for illustrative purposes only and comes
# with no warranty or guarantee.

# Uses bootstrapping to estimate 95% CIs for the risk differences and
# relative risks.

library(survival)
library(cmprsk)
library(riskRegression)
library(prodlim)
library(rms)

N.boot <- 2000

################################################################################
# Read in EFFECT data.
################################################################################

zlist <- list(age=0,prevMI=0,prevCHF=0,HeartRate=0,SBP=0,
   creatinine=0,enzymes=0,stemi=0,pci=0,treat=0,event.type=0,T.event=0)

effect.df <- data.frame(scan("cohort.dat",zlist))
# NOTE: Due to privacy regulations, the data cannot be disseminated.
# Do not contact the author requesting these data.

remove(zlist)

################################################################################
# Create RCS terms for use in regression models.
################################################################################

age.mat <- rcspline.eval(effect.df$age,nk=5,inclx=F,norm=0)
HeartRate.mat <- rcspline.eval(effect.df$HeartRate,nk=5,inclx=F,norm=0)
SBP.mat <- rcspline.eval(effect.df$SBP,nk=5,inclx=F,norm=0)
creatinine.mat <- rcspline.eval(effect.df$creatinine,nk=5,inclx=F,norm=0)

colnames(age.mat) <- c("age2","age3","age4")
colnames(HeartRate.mat) <- c("HeartRate2","HeartRate3","HeartRate4")
colnames(SBP.mat) <- c("SBP2","SBP3","SBP4")
colnames(creatinine.mat) <- c("creatinine2","creatinine3","creatinine4")

age.mat.df <- data.frame(age.mat)
HeartRate.mat.df <- data.frame(HeartRate.mat)
SBP.mat.df <- data.frame(SBP.mat)
creatinine.mat.df <- data.frame(creatinine.mat)

effect.df <- cbind(effect.df,age.mat.df,HeartRate.mat.df,SBP.mat.df,
  creatinine.mat.df)

remove(age.mat,HeartRate.mat,SBP.mat,creatinine.mat)
remove(age.mat.df,HeartRate.mat.df,SBP.mat.df,creatinine.mat.df)

################################################################################
# Weighted Aalen-Johansen estimator.
################################################################################

psm <- glm(treat ~ age + prevMI + prevCHF + HeartRate + SBP + creatinine +
   enzymes + stemi + pci +
  age2 + age3 + age4 + HeartRate2 + HeartRate3 + HeartRate4 +
  SBP2 + SBP3 + SBP4 + creatinine2 + creatinine3 + creatinine4,
  family="binomial",data=effect.df)

ps <- psm$fitted

Z <- effect.df$treat

effect.df$iptw <- (Z/ps) + ((1-Z)/(1-ps))

remove(psm,ps,Z)

cif.aj <- prodlim(data=effect.df,Hist(T.event,event.type) ~ treat,
  type="risk",caseweights=effect.df$iptw)

risk.aj <- predict(cif.aj,times=365*(1:5),cause='1',
  newdata=data.frame(treat=0:1))

risk0.aj <- risk.aj$'treat=0'
risk1.aj <- risk.aj$'treat=1'

rd.aj <- risk1.aj - risk0.aj
rr.aj <- risk1.aj/risk0.aj

remove(cif.aj,risk.aj)

cat(round(rd.aj,3),file="boot.rd.aj.out",fill=T,append=F)
cat(round(rr.aj,3),file="boot.rr.aj.out",fill=T,append=F)

################################################################################
# Estimate the relative risks at years 1 to 5 using the IPTW-IPCW estimator.
################################################################################

effect.df$treat.factor <- factor(effect.df$treat)

treat_mod <- glm(treat.factor ~ age + prevMI + prevCHF + HeartRate + SBP + 
  creatinine + enzymes + stemi + pci +
  age2 + age3 + age4 + HeartRate2 + HeartRate3 + HeartRate4 +
  SBP2 + SBP3 + SBP4 + creatinine2 + creatinine3 + creatinine4,
  family="binomial",data=effect.df)

cox_mod <- CSC(Hist(T.event,event.type) ~ age + prevMI + prevCHF + HeartRate + 
  SBP + creatinine + enzymes + stemi + pci + treat.factor +
  age2 + age3 + age4 + HeartRate2 + HeartRate3 + HeartRate4 +
  SBP2 + SBP3 + SBP4 + creatinine2 + creatinine3 + creatinine4,
  data=effect.df)

censor_mod <- coxph(Surv(T.event,event.type==0) ~ 1,data=effect.df)

ate.iptw <- ate(event=cox_mod,treatment=treat_mod,censor=censor_mod,
  data=effect.df,estimator="IPTW",times=365*(1:5),cause=1,se=F)

rr.iptw <- ate.iptw$diffRisk$estimate.B/ate.iptw$diffRisk$estimate.A

remove(treat_mod,cox_mod,censor_mod,ate.iptw)

cat(round(rr.iptw,3),file="boot.rr.iptw.out",fill=T,append=F)

################################################################################
# Bootstrap
################################################################################

rd.aj.bs <- NULL
rr.aj.bs <- NULL
rr.iptw.bs <- NULL

effect.df <- subset(effect.df,select=c(age,prevMI,prevCHF,HeartRate,SBP,
  creatinine,enzymes,stemi,pci,treat,treat.factor,event.type,T.event))
# Keep only original variables. We will define RCS components in the
# bootstrap sample.

for (iter in 1:N.boot){
  cat(iter,file="boot.iter",fill=T,append=F)
  set.seed(iter)

  sample.id <- sample(1:nrow(effect.df),size=nrow(effect.df),replace=T)

  boot.df <- effect.df[sample.id,]

  remove(sample.id)

  ##############################################################################
  # Create RCS terms for use in regression models.
  ##############################################################################

  age.mat <- rcspline.eval(boot.df$age,nk=5,inclx=F,norm=0)
  HeartRate.mat <- rcspline.eval(boot.df$HeartRate,nk=5,inclx=F,norm=0)
  SBP.mat <- rcspline.eval(boot.df$SBP,nk=5,inclx=F,norm=0)
  creatinine.mat <- rcspline.eval(boot.df$creatinine,nk=5,inclx=F,norm=0)

  colnames(age.mat) <- c("age2","age3","age4")
  colnames(HeartRate.mat) <- c("HeartRate2","HeartRate3","HeartRate4")
  colnames(SBP.mat) <- c("SBP2","SBP3","SBP4")
  colnames(creatinine.mat) <- c("creatinine2","creatinine3","creatinine4")

  age.mat.df <- data.frame(age.mat)
  HeartRate.mat.df <- data.frame(HeartRate.mat)
  SBP.mat.df <- data.frame(SBP.mat)
  creatinine.mat.df <- data.frame(creatinine.mat)

  boot.df <- cbind(boot.df,age.mat.df,HeartRate.mat.df,SBP.mat.df,
    creatinine.mat.df)

  remove(age.mat,HeartRate.mat,SBP.mat,creatinine.mat)
  remove(age.mat.df,HeartRate.mat.df,SBP.mat.df,creatinine.mat.df)

  ##############################################################################
  # Weighted Aalen-Johansen estimator.
  ##############################################################################

  psm <- glm(treat ~ age + prevMI + prevCHF + HeartRate + SBP + creatinine +
    enzymes + stemi + pci +
    age2 + age3 + age4 + HeartRate2 + HeartRate3 + HeartRate4 +
    SBP2 + SBP3 + SBP4 + creatinine2 + creatinine3 + creatinine4,
    family="binomial",data=boot.df)

  ps <- psm$fitted

  Z <- boot.df$treat

  boot.df$iptw <- (Z/ps) + ((1-Z)/(1-ps))

  remove(psm,ps,Z)

  cif.aj <- prodlim(data=boot.df,Hist(T.event,event.type) ~ treat,
    type="risk",caseweights=boot.df$iptw)

  risk.aj <- predict(cif.aj,times=365*(1:5),cause='1',
    newdata=data.frame(treat=0:1))

  risk0.aj <- risk.aj$'treat=0'
  risk1.aj <- risk.aj$'treat=1'

  rd.aj <- risk1.aj - risk0.aj
  rr.aj <- risk1.aj/risk0.aj

  rd.aj.bs <- rbind(rd.aj.bs,rd.aj)
  rr.aj.bs <- rbind(rr.aj.bs,rr.aj)

  remove(cif.aj,risk.aj,risk0.aj,risk1.aj,rd.aj,rr.aj)

  ##############################################################################
  # IPTW-IPCW estimator.
  ##############################################################################

  treat_mod_bs <- glm(treat.factor ~ age + prevMI + prevCHF + HeartRate + SBP +
    creatinine + enzymes + stemi + pci +
    age2 + age3 + age4 + HeartRate2 + HeartRate3 + HeartRate4 +
    SBP2 + SBP3 + SBP4 + creatinine2 + creatinine3 + creatinine4,
     family="binomial",data=boot.df)

  cox_mod_bs <- CSC(Hist(T.event,event.type) ~ age + prevMI + prevCHF + 
    HeartRate + SBP + creatinine + enzymes + stemi + pci + treat.factor +
    age2 + age3 + age4 + HeartRate2 + HeartRate3 + HeartRate4 +
    SBP2 + SBP3 + SBP4 + creatinine2 + creatinine3 + creatinine4,
    data=boot.df)

  censor_mod_bs <- coxph(Surv(T.event,event.type==0) ~ 1,data=boot.df)

  ate.iptw <- ate(event=cox_mod_bs,treatment=treat_mod_bs,censor=censor_mod_bs,
    data=boot.df,estimator="IPTW",times=365*(1:5),cause=1,se=F)

  rr.iptw.bs <- rbind(rr.iptw.bs,
    ate.iptw$diffRisk$estimate.B/ate.iptw$diffRisk$estimate.A)

  remove(treat_mod_bs,cox_mod_bs,censor_mod_bs,ate.iptw)
}

cat(round(apply(rd.aj.bs,MARGIN=2,FUN=quantile,probs=0.025),3),
  file="boot.rd.aj.out",fill=T,append=T)
cat(round(apply(rd.aj.bs,MARGIN=2,FUN=quantile,probs=0.975),3),
  file="boot.rd.aj.out",fill=T,append=T)

cat(round(apply(rr.aj.bs,MARGIN=2,FUN=quantile,probs=0.025),3),
  file="boot.rr.aj.out",fill=T,append=T)
cat(round(apply(rr.aj.bs,MARGIN=2,FUN=quantile,probs=0.975),3),
  file="boot.rr.aj.out",fill=T,append=T)

cat(round(apply(rr.iptw.bs,MARGIN=2,FUN=quantile,probs=0.025),3),
  file="boot.rr.iptw.out",fill=T,append=T)
cat(round(apply(rr.iptw.bs,MARGIN=2,FUN=quantile,probs=0.975),3),
  file="boot.rr.iptw.out",fill=T,append=T)
