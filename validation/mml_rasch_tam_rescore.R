#!/usr/bin/env Rscript
# Recompute TAM posterior summaries at every retained endpoint, without refitting.
# Original failed checks and saved objects remain unchanged.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L, !dir.exists(args[3]), as.character(packageVersion("TAM")) == "4.3.25")
native <- args[1]; reference <- args[2]; out <- args[3]
hash <- function(p) digest::digest(file=p, algo="sha256")
verify <- function(folder) {
  m <- jsonlite::read_json(file.path(folder,"manifest.json"))
  stopifnot(all(vapply(names(m),function(n) identical(hash(file.path(folder,n)),m[[n]]),logical(1))))
}
verify(native); verify(reference)
dir.create(out,recursive=TRUE)
write_new <- function(value,name) {
  stopifnot(!file.exists(file.path(out,name)))
  jsonlite::write_json(value,file.path(out,name),digits=NA,auto_unbox=TRUE,pretty=TRUE,na="null",null="null")
}
write_new(list(classification="POST_RESULT_SCORING_DIAGNOSTIC",scientific_inference_ready=FALSE,
  no_refitting=TRUE,original_checks_unchanged=TRUE,original_moment_tolerance=1e-8,
  scope="Fresh TAM fixed-grid posterior versus independent continuous reference at all twelve saved points",
  native_manifest_sha256=hash(file.path(native,"manifest.json")),
  reference_manifest_sha256=hash(file.path(reference,"manifest.json"))),"protocol.json")
file.copy("validation/mml_rasch_tam_rescore.R",file.path(out,"source.R"))
ns <- asNamespace("TAM")
source_names <- c("tam.mml","tam_mml_calc_prob","tam_mml_calc_prob_R","tam_calc_posterior","tam_mml_person_posterior")
write_new(setNames(lapply(source_names,function(n) deparse(get(n,ns))),source_names),"installed_sources.json")
write_new(list(R=R.version.string,TAM=as.character(packageVersion("TAM"))),"runtime.json")
original <- jsonlite::read_json(file.path(reference,"summary.json"),simplifyVector=FALSE)
records <- lapply(original$cases,function(record) {
  id <- record$id
  fit <- readRDS(file.path(native,paste0(id,".rds")))
  oracle <- jsonlite::read_json(file.path(reference,paste0(id,".json")),simplifyVector=FALSE)$continuous[[2]]
  n <- nrow(fit$resp); items <- ncol(fit$resp); nk <- dim(fit$A)[2]; q <- nrow(fit$theta)
  fresh <- get("tam_mml_calc_prob",ns)(iIndex=seq_len(items),A=fit$A,AXsi=matrix(0,items,nk),
    B=fit$B,xsi=fit$xsi$xsi,theta=fit$theta,nnodes=q,maxK=nk,recalc=TRUE,use_rcpp=FALSE)
  prior <- dnorm(fit$theta[,1],mean=fit$beta[1],sd=sqrt(fit$variance[1])); prior <- prior/sum(prior)
  posterior <- get("tam_calc_posterior",ns)(rprobs=fresh$rprobs,gwt=matrix(prior,n,q,byrow=TRUE),
    resp=fit$resp,nitems=items,resp.ind.list=fit$resp.ind.list,resp.ind=fit$resp.ind,
    normalization=TRUE,avoid.zerosum=FALSE)$hwt
  stopifnot(all(is.finite(posterior)),max(abs(rowSums(posterior)-1)) < 1e-12)
  summaries <- function(hwt) {
    get("tam_mml_person_posterior",ns)(pid=fit$pid,nstud=n,pweights=fit$pweights,
      resp=fit$resp,resp.ind=fit$resp.ind,snodes=0,hwtE=hwt,hwt=hwt,ndim=1,theta=fit$theta)$person
  }
  saved <- summaries(fit$hwt); recomputed <- summaries(posterior)
  stopifnot(identical(as.character(recomputed$pid),as.character(fit$person$pid)))
  error <- list(cached_eap_reconstruction=max(abs(saved$EAP-fit$person$EAP)),
    cached_sd_reconstruction=max(abs(saved$SD.EAP-fit$person$SD.EAP)),
    fresh_vs_reference_eap=max(abs(recomputed$EAP-unlist(oracle$eap))),
    fresh_vs_reference_sd=max(abs(recomputed$SD.EAP-unlist(oracle$sd))),
    cached_vs_fresh_eap=max(abs(saved$EAP-recomputed$EAP)),
    cached_vs_fresh_sd=max(abs(saved$SD.EAP-recomputed$SD.EAP)),
    cached_vs_fresh_probabilities=max(abs(fit$rprobs-fresh$rprobs)),
    cached_vs_fresh_posterior=max(abs(fit$hwt-posterior)))
  checks <- list(cached_reconstruction=max(error$cached_eap_reconstruction,error$cached_sd_reconstruction) < 1e-12,
    fresh_moments=max(error$fresh_vs_reference_eap,error$fresh_vs_reference_sd) < 1e-8)
  result <- list(id=id,original_checks=record$checks,errors=error,checks=checks,
    fresh_eap=recomputed$EAP,fresh_sd=recomputed$SD.EAP,cached_eap=saved$EAP,cached_sd=saved$SD.EAP)
  write_new(result,paste0(id,".json"))
  cat(id,"fresh EAP error",error$fresh_vs_reference_eap,"cached/fresh EAP difference",error$cached_vs_fresh_eap,"\n")
  result[c("id","original_checks","errors","checks")]
})
verify(native); verify(reference)
passed <- length(records)==12L && all(vapply(records,function(r) all(unlist(r$checks)),logical(1)))
write_new(list(classification="POST_RESULT_SCORING_DIAGNOSTIC",scientific_inference_ready=FALSE,
  original_all_checks_pass=original$all_checks_pass,completed_points=length(records),
  fresh_scoring_checks_pass=passed,cases=records),"summary.json")
files <- list.files(out,full.names=TRUE)
write_new(setNames(lapply(files,hash),basename(files)),"manifest.json")
if (!passed) quit(status=1L)
