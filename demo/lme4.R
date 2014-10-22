traceHat.sp <-function(fit) {
  fit1 <- fit
  U1 <- cBind( t(getME(fit1,"A")),as(getME(fit1,"X"),"dgCMatrix"))
  n1 <- nrow(getME(fit1,"A"))
  D1 <- cBind(Diagonal(n1), Matrix(0, n1,ncol(getME(fit1,"X"))))  
  M1 <-rBind(U1,D1)
  rho <- sum(diag(solve(crossprod(M1)) %*% crossprod(U1)))
  p <- length(fixef(fit1))
  N <- nrow(getME(fit1,"X"))
  K.corr <- N*(N-p-1)*(rho+1)/((N-p)*(N-p-2)) + N*(p+1)/((N-p)*(N-p-2))
  #REML
  #K.corr <- (N-p-1)*(rho+1)/(N-p-2) + (p+1)/(N-p-2)
  return(K.corr)
  #return(rho)
}

calAIC.lmer <- function(fit) {
  fit1 <- fit
  logL <- sum(dnorm(getME(fit1,"y"),predict(fit1),attr(VarCorr(fit1),"sc"),log=T))
  df <- traceHat.sp(fit1)
  cAIC <- -2*logL+2*df
  return( c(logL,df,cAIC) )
}



library(lme4)
library(Matrix)
library(cAIC4)
library(mlmRev)

fm1 <- lmer(Yield ~ 1|Batch, Dyestuff ,REML=F)
system.time(print(cAIC(fm1)))
system.time(print(calAIC.lmer(fm1)))

fm2 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,REML=F)
system.time(print(cAIC(fm2)))
system.time(print(cAIC(fm2,analytic = FALSE)))
system.time(print(calAIC.lmer(fm2)))

fm3  <- lmer(y ~ dept*service + (1|s) ,InstEval,REML=F)
system.time(print(cAIC(fm3)))
system.time(print(calAIC.lmer(fm3)))
            
fm4 <- lmer(score ~  (1|lea), Chem97,REML=F)
system.time(print(cAIC(fm4)))
system.time(print(calAIC.lmer(fm4)))

system.time (fm8 <- lmer(y ~ dept*service + (1|s) + (1|d),InstEval,REML=F))
system.time(print(cAIC(fm8)))
system.time(print(calAIC.lmer(fm8)))

fm9 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
system.time(print(cAIC(fm9)))
