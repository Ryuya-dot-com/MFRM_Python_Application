#!/usr/bin/env Rscript
# Native unit-weight PCM fits with retained regression means and missing responses.
args <- commandArgs(TRUE); stopifnot(length(args)==3L)
spec <- jsonlite::read_json(args[1],simplifyVector=TRUE); folder<-args[2]; mode<-args[3]
stopifnot(!spec$scientific_inference_ready,as.character(packageVersion('TAM'))==spec$tam_version)
write_new<-function(x,path){stopifnot(!file.exists(path));jsonlite::write_json(x,path,digits=NA,auto_unbox=TRUE,pretty=TRUE,na='null',null='null')}
input<-jsonlite::read_json(file.path(folder,'input.json'),simplifyVector=TRUE)
if(mode %in% c('fixed','refit')) {
  namespace<-asNamespace('TAM'); hooks<-c('tam_mml_progress_em0','tam_mml_progress_em')
  originals<-lapply(hooks,get,envir=namespace)
  trace(hooks[1],where=namespace,print=FALSE,tracer=quote({
    frames<-Filter(function(f) exists('tam_fct',f,inherits=FALSE)&&identical(f$tam_fct,'tam.mml.mfr'),sys.frames())
    stopifnot(length(frames)==1L);e<-frames[[1]]
    .GlobalEnv$.before[[iter]]<-list(iter=iter,xsi=e$xsi,beta=e$beta,variance=e$variance)
  }))
  trace(hooks[2],where=namespace,print=FALSE,tracer=quote({
    frames<-Filter(function(f) exists('tam_fct',f,inherits=FALSE)&&identical(f$tam_fct,'tam.mml.mfr'),sys.frames())
    stopifnot(length(frames)==1L);e<-frames[[1]]
    .GlobalEnv$.after[[iter]]<-list(iter=iter,xsi=e$xsi,beta=e$beta,variance=e$variance,deviance=deviance,
      xsi_change=xsi_change,beta_change=beta_change,variance_change=variance_change,deviance_change=deviance_change)
  }))
  response<-matrix(input$responses,ncol=2,byrow=TRUE,dimnames=list(NULL,c('C1','C2')))
  facets<-data.frame(rater=rep(rep(c('R1','R2'),each=2),32),task=rep(c('T1','T2'),64))
  for(bound in spec$bounds) for(kind in mode) {
    key<-paste0('tam_b',bound,'_',kind); path<-file.path(folder,key)
    start<-if(kind=='fixed') input$continuous_point else spec$start
    # App rater sign is positive; native rater/task effects sum to zero.
    xsi<-c(start[4:5]+.2-(.35+start[2])/2,(start[2]-.35)/2,start[3]-.2,start[c(6,8,7,9)])
    .GlobalEnv$.before<-list();.GlobalEnv$.after<-list(); warnings<-character()
    fixed<-kind=='fixed'
    fit<-withCallingHandlers(TAM::tam.mml.mfr(resp=response,facets=facets,pid=rep(sprintf('P%03d',0:31),each=4),
      Y=matrix(rep(input$x,each=4),ncol=1,dimnames=list(NULL,'x')),formulaA=~item+rater+task+item:step,
      constraint='cases',delete.red.items=FALSE,est.variance=TRUE,
      xsi.inits=cbind(1:8,xsi),xsi.fixed=if(fixed)cbind(1:8,xsi) else NULL,
      beta.inits=rbind(c(1,1,0),c(2,1,start[1])),beta.fixed=if(fixed)rbind(c(1,1,0),c(2,1,start[1])) else matrix(c(1,1,0),nrow=1),
      variance.inits=matrix(exp(2*start[10]),1,1),variance.fixed=if(fixed)matrix(c(1,1,exp(2*start[10])),nrow=1) else NULL,
      control=list(nodes=seq(-bound,bound,by=spec$spacing),snodes=0,maxiter=if(fixed)3L else spec$maxiter,
        conv=1e-10,convD=1e-12,convM=1e-10,Msteps=20L,progress=FALSE),verbose=FALSE),
      warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    stopifnot(length(.before)==fit$iter,length(.after)==fit$iter,!file.exists(paste0(path,'.rds')))
    stopifnot(max(abs(fit$Y[,2]-input$x))<1e-15)
    saveRDS(fit,paste0(path,'.rds'))
    write_new(list(kind=kind,bound=bound,iter=fit$iter,reached_cap=fit$iter>=(if(fixed)3L else spec$maxiter),warnings=warnings,
      A=fit$A,B=fit$B,AXsi=fit$AXsi,xsi_names=rownames(fit$xsi),item_names=dimnames(fit$A)[[1]],xsi=fit$xsi$xsi,
      beta=fit$beta,variance=fit$variance,Y=fit$Y,resp=fit$resp,pid=fit$pid,person=fit$person,
      initial=.before[[1]],before=.before,after=.after,deviance=fit$deviance,history=fit$deviance.history,
      control=fit$control,runtime=list(R=R.version.string,TAM=as.character(packageVersion('TAM')))),paste0(path,'.json'))
    cat(key,'iterations',fit$iter,'SD',sqrt(fit$variance[1]),'\n')
  }
  for(name in hooks)untrace(name,where=namespace)
  stopifnot(all(vapply(seq_along(hooks),function(i)identical(get(hooks[i],namespace),originals[[i]]),logical(1))))
} else {
  stopifnot(mode=='grid'); cases<-jsonlite::read_json(file.path(folder,'r_cases.json'))
  idx<-input$indices; persons<-idx$person+1L; r<-idx$facets$Rater+1L;t<-idx$facets$Task+1L;c<-idx$facets$Criterion+1L
  lse<-function(v){m<-max(v);m+log(sum(exp(v-m)))}
  answer<-lapply(cases,function(case){
    par<-unlist(case$coordinates);theta<-seq(-case$bound,case$bound,by=spec$spacing);sigma<-exp(par[10])
    raters<-c(.35,par[2]);tasks<-c(par[3],.4-par[3]);criteria<-par[4:5]
    cumulative<-rbind(c(0,par[6],par[6]+par[7],0),c(0,par[8],par[8]+par[9],0))
    values<-lapply(1:32,function(p){
      rows<-which(persons==p)
      ll<-vapply(theta,function(v){
        eta<-v+raters[r[rows]]-tasks[t[rows]]-criteria[c[rows]]
        logits<-outer(eta,0:3)-cumulative[c[rows],,drop=FALSE]
        sum(logits[cbind(seq_along(rows),idx$score_k[rows]+1L)]-apply(logits,1,lse))
      },numeric(1))
      logw<-dnorm(theta,input$x[p]*par[1],sigma,log=TRUE)+log(spec$spacing)
      mass<-lse(ll+logw);post<-exp(ll+logw-mass);eap<-sum(post*theta)
      c(nll=-mass,log_prior_mass=lse(logw),eap=eap,sd=sqrt(sum(post*(theta-eap)^2)))
    }); v<-do.call(rbind,values)
    list(nll=sum(v[,'nll']),normalized_nll=sum(v[,'nll']+v[,'log_prior_mass']),eap=v[,'eap'],sd=v[,'sd'],log_prior_mass=v[,'log_prior_mass'])
  })
  write_new(answer,file.path(folder,'r_grid.json'))
}
