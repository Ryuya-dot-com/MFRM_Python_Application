#!/usr/bin/env Rscript
# Fixed-calibration case-weight calls and an independent scalar PCM calculation.
args <- commandArgs(TRUE); stopifnot(length(args)==1L)
out <- args[1]; spec <- jsonlite::read_json(file.path(out,'protocol.json'),simplifyVector=TRUE)
stopifnot(!spec$scientific_inference_ready,as.character(packageVersion('TAM'))==spec$tam_version)
write_new <- function(x,path) {
  stopifnot(!file.exists(path)); jsonlite::write_json(x,path,digits=NA,auto_unbox=TRUE,pretty=TRUE,na='null',null='null')
}
ns <- asNamespace('TAM')
for(n in c('tam.mml.mfr','tam_mml_mfr_proc_multiple_person_ids','tam_mml_compute_deviance','tam_calc_posterior'))
  writeLines(deparse(get(n,ns)),file.path(out,paste0(n,'_installed.R')))
lse <- function(v) {m<-max(v);m+log(sum(exp(v-m)))}
for(cell in spec$cells) {
  folder <- file.path(out,cell); input <- jsonlite::read_json(file.path(folder,'input.json'),simplifyVector=TRUE)
  par <- input$continuous_point; idx <- input$indices
  r <- idx$facets$Rater+1L; t <- idx$facets$Task+1L; c <- idx$facets$Criterion+1L
  theta <- seq(-spec$bound,spec$bound,by=spec$spacing); sigma <- exp(par[10])
  raters<-c(.35,par[2]);tasks<-c(par[3],.4-par[3]);criteria<-par[4:5]
  cumulative<-rbind(c(0,par[6],par[6]+par[7],0),c(0,par[8],par[8]+par[9],0))
  scalar_grid <- function(response_weights) {
    values <- lapply(1:32,function(p) {
      rows<-which(idx$person+1L==p)
      ll<-vapply(theta,function(v) {
        eta<-v+raters[r[rows]]-tasks[t[rows]]-criteria[c[rows]]
        logits<-outer(eta,0:3)-cumulative[c[rows],,drop=FALSE]
        sum(response_weights[rows]*(logits[cbind(seq_along(rows),idx$score_k[rows]+1L)]-apply(logits,1,lse)))
      },numeric(1))
      logw<-dnorm(theta,input$x[p]*par[1],sigma,log=TRUE)+log(spec$spacing)
      mass<-lse(ll+logw);post<-exp(ll+logw-mass);eap<-sum(post*theta)
      c(person_nll=-mass,log_prior_mass=lse(logw),eap=eap,sd=sqrt(sum(post*(theta-eap)^2)))
    });v<-do.call(rbind,values)
    as.list(as.data.frame(v))
  }
  reference <- list(unit=scalar_grid(rep(1,length(idx$person))))
  response<-matrix(input$responses,ncol=2,byrow=TRUE,dimnames=list(NULL,c('C1','C2')))
  facets<-data.frame(rater=rep(rep(c('R1','R2'),each=2),32),task=rep(c('T1','T2'),64))
  old<-readRDS(file.path(spec$prior,cell,'tam_b12_fixed.rds'))
  xsi<-c(par[4:5]+.2-(.35+par[2])/2,(par[2]-.35)/2,par[3]-.2,par[c(6,8,7,9)])
  for(scheme in names(spec$case_weights)) {
    w<-unlist(spec$case_weights[[scheme]]); stopifnot(length(w)==32,all(is.finite(w)),all(w>0))
    if(scheme!='unit')reference[[scheme]]<-scalar_grid(w[idx$person+1L])
    warnings<-character()
    fit<-withCallingHandlers(TAM::tam.mml.mfr(resp=response,facets=facets,pid=rep(sprintf('P%03d',0:31),each=4),
      Y=matrix(rep(input$x,each=4),ncol=1,dimnames=list(NULL,'x')),pweights=rep(w,each=4),
      formulaA=~item+rater+task+item:step,constraint='cases',delete.red.items=FALSE,est.variance=TRUE,
      xsi.inits=cbind(1:8,xsi),xsi.fixed=cbind(1:8,xsi),
      beta.inits=rbind(c(1,1,0),c(2,1,par[1])),beta.fixed=rbind(c(1,1,0),c(2,1,par[1])),
      variance.inits=matrix(exp(2*par[10]),1,1),variance.fixed=matrix(c(1,1,exp(2*par[10])),nrow=1),
      control=list(nodes=theta,snodes=0,maxiter=3L,conv=1e-10,convD=1e-12,convM=1e-10,Msteps=20L,progress=FALSE),verbose=FALSE),
      warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    checks<-c(design=identical(fit$A,old$A)&&identical(fit$B,old$B),responses=identical(fit$resp,old$resp),
      persons=identical(fit$pid,old$pid),means=max(abs(fit$Y[,2]-input$x))<1e-15,
      case_weights=max(abs(fit$pweights-32*w/sum(w)))<1e-14,
      xsi=max(abs(fit$xsi$xsi-xsi))<1e-12,beta=max(abs(fit$beta-c(0,par[1])))<1e-12,
      variance=abs(fit$variance[1]-exp(2*par[10]))<2e-10)
    path<-file.path(folder,scheme,'tam');stopifnot(!file.exists(paste0(path,'.rds')))
    saveRDS(fit,paste0(path,'.rds'))
    write_new(list(checks=as.list(checks),deviance=fit$deviance,pweights=fit$pweights,person=fit$person,
      xsi=fit$xsi$xsi,beta=fit$beta,variance=fit$variance,iter=fit$iter,reached_cap=fit$iter>=3L,warnings=warnings,
      score_change_to_previous=list(eap=max(abs(fit$person$EAP-old$person$EAP)),sd=max(abs(fit$person$SD.EAP-old$person$SD.EAP))),
      runtime=list(R=R.version.string,TAM=as.character(packageVersion('TAM')))),paste0(path,'.json'))
    stopifnot(all(checks));cat(cell,scheme,'TAM fixed controls pass\n')
  }
  write_new(reference,file.path(folder,'scalar_R.json'))
}
