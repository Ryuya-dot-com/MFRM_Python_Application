#!/usr/bin/env Rscript
# Diagnose saved-posterior scoring and the fixed-mean variance update, without refits.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
root <- args[1]
spec <- jsonlite::read_json(file.path(root,"input.json"), simplifyVector=TRUE)
python <- jsonlite::read_json(file.path(root,"python.json"))
stopifnot(!spec$scientific_inference_ready, !spec$qualification_eligible,
          !file.exists(file.path(root,"update_audit.json")))
sources <- c("tam.mml.mfr", "tam_mml_mstep_regression", "tam_mml_person_posterior")
source_text <- lapply(sources, function(n) deparse(get(n,asNamespace("TAM"))))
names(source_text) <- sources
results <- list()
for (id in spec$tam$id) {
  fit <- readRDS(file.path(root,paste0(id,".rds")))
  theta <- as.numeric(fit$theta)
  y <- fit$resp
  n <- nrow(y)
  stopifnot(identical(dim(y),c(24L,20L)), all(fit$pweights == 1),
            fit$beta[1] == 0, all(fit$Y == 1), all(fit$B[1,,1] == 0:4))
  ll <- matrix(0,n,length(theta))
  for (i in seq_len(ncol(y))) {
    logits <- outer(0:4,theta) + fit$AXsi[i,]
    shifted <- sweep(logits,2,apply(logits,2,max),"-")
    logp <- sweep(shifted,2,log(colSums(exp(shifted))),"-")
    ll <- ll + logp[y[,i]+1,,drop=FALSE]
  }
  joint <- sweep(ll,2,dnorm(theta,sd=sqrt(fit$variance[1]),log=TRUE),"+")
  centered <- exp(joint-apply(joint,1,max))
  fresh <- centered/rowSums(centered)
  eap <- as.numeric(fresh%*%theta)
  sd <- sqrt(rowSums(fresh*sweep(matrix(theta,n,length(theta),byrow=TRUE),1,eap,"-")^2))
  stored_eap <- as.numeric(fit$hwt%*%theta)
  stored_sd <- sqrt(rowSums(fit$hwt*sweep(matrix(theta,n,length(theta),byrow=TRUE),1,stored_eap,"-")^2))
  uncentered_second <- mean(as.numeric(fit$hwt%*%(theta^2)))
  posterior_mean <- mean(stored_eap)
  native_update <- TAM:::tam_mml_mstep_regression(resp=y,hwt=fit$hwt,resp.ind=fit$resp.ind,
    pweights=fit$pweights,pweightsM=fit$pweights,Y=fit$Y,theta=fit$theta,
    theta2=fit$theta^2,YYinv=solve(crossprod(fit$Y)),ndim=1,nstud=n,
    beta.fixed=fit$beta.fixed,variance=fit$variance,Variance.fixed=fit$variance.fixed,
    group=fit$group,G=1,nomiss=TRUE,iter=2001,min.variance=.001,beta=fit$beta)
  predicted_update <- uncentered_second-posterior_mean^2+1e-10
  py <- python[[id]]$finite
  differences <- c(stored_eap=max(abs(stored_eap-fit$person$EAP)),
    stored_sd=max(abs(stored_sd-fit$person$SD.EAP)),
    refreshed_python_eap=max(abs(eap-unlist(py$eap))),
    refreshed_python_sd=max(abs(sd-unlist(py$sd))),
    variance_update=abs(native_update$variance[1]-predicted_update))
  stopifnot(max(differences) < 1e-11)
  results[[id]] <- list(differences=as.list(differences),refreshed_eap=eap,refreshed_sd=sd,
    native_refreshed_eap_difference=max(abs(eap-fit$person$EAP)),
    stored_posterior_mean=posterior_mean,stored_second_moment=uncentered_second,
    native_mean_fixed=fit$beta[1],native_variance=fit$variance[1],
    native_variance_update=native_update$variance[1],
    mean_zero_expected_complete_loglik_variance_update=uncentered_second,
    centered_update_difference=uncentered_second-native_update$variance[1],
    stored_log_sigma_nll_score=n*(1-uncentered_second/fit$variance[1]),
    centered_fixed_point_score=n*(1e-10-posterior_mean^2)/fit$variance[1],
    fresh_log_sigma_nll_score=n*(1-mean(as.numeric(fresh%*%(theta^2)))/fit$variance[1]))
}
jsonlite::write_json(list(classification="OBSERVED_DEVELOPMENT_ONLY",scientific_inference_ready=FALSE,
  qualification_eligible=FALSE,arithmetic_checks_pass=TRUE,results=results,
  installed_TAM_sources=source_text,runtime=list(R=R.version.string,TAM=as.character(packageVersion("TAM"))),
  source_sha256=digest::digest(file="validation/mml_pcm_tam_update_audit.R",algo="sha256"),
  original_summary_sha256=digest::digest(file=file.path(root,"summary.json"),algo="sha256"),
  scope="Stored versus refreshed posterior; isolated native variance M-step replay, no estimator modification"),
  file.path(root,"update_audit.json"),digits=NA,pretty=TRUE,auto_unbox=TRUE)
cat("Stored/fresh posterior and fixed-mean variance update checks: 12/12\n")
