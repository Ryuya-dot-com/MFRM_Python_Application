#!/usr/bin/env Rscript
# Offline fixed-SD PCM comparison; version and model are checked before fitting.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L, args[2] %in% c('fixed','refit'))
job <- jsonlite::read_json(args[1],simplifyVector=TRUE)
stopifnot(identical(job$schema,'mfrm_pcm_native_input_v1'), !job$scientific_inference_ready,
          as.character(packageVersion('TAM'))==job$TAM_version,
          !file.exists('fit.rds'), !file.exists('result.json'))
fixed <- args[2]=='fixed'; n <- length(job$pid)
resp <- matrix(t(job$responses),ncol=2,byrow=TRUE,dimnames=list(NULL,c('C1','C2')))
facets <- data.frame(rater=rep(rep(c('R1','R2'),each=2),n),task=rep(c('T1','T2'),2*n))
xsi <- job$xsi[c(1,2,3,4,5,7,6,8)]
warnings <- character()
fit <- withCallingHandlers(TAM::tam.mml.mfr(resp=resp,facets=facets,pid=rep(job$pid,each=4),
  formulaA=~item+rater+task+item:step,constraint='cases',delete.red.items=FALSE,
  xsi.inits=cbind(1:8,xsi),xsi.fixed=if(fixed)cbind(1:8,xsi) else NULL,
  beta.inits=matrix(c(1,1,0),nrow=1),beta.fixed=matrix(c(1,1,0),nrow=1),
  est.variance=TRUE,variance.inits=matrix(job$sigma^2,1,1),variance.fixed=matrix(c(1,1,job$sigma^2),nrow=1),
  control=list(nodes=seq(-job$bound,job$bound,length.out=job$nodes),snodes=0,
    maxiter=if(fixed)3L else 2000L,conv=1e-10,convD=1e-12,convM=1e-10,Msteps=20L,progress=FALSE),verbose=FALSE),
  warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
saveRDS(fit,'fit.rds')
jsonlite::write_json(list(runtime=list(R=R.version.string,TAM=as.character(packageVersion('TAM'))),
  mode=args[2],A=fit$A,B=fit$B,xsi_names=rownames(fit$xsi),item_names=dimnames(fit$A)[[1]],
  xsi=fit$xsi$xsi,beta=fit$beta,variance=fit$variance,resp=fit$resp,pid=fit$pid,person=fit$person,
  deviance=fit$deviance,deviance_history=fit$deviance.history,iterations=fit$iter,reached_cap=fit$iter>=(if(fixed)3L else 2000L),
  control=fit$control,warnings=warnings,scientific_inference_ready=FALSE),
  'result.json',digits=NA,auto_unbox=TRUE,na='null',null='null',pretty=TRUE)
