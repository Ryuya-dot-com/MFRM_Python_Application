#!/usr/bin/env Rscript
# Matched-model native TAM refits of all six retained development datasets.
# Explicit native-coordinate A/B matrices avoid software-dependent naming constraints.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(2L,3L), !dir.exists(args[2]))
parent <- args[1]; out <- args[2]
case_pattern <- if (length(args) == 3L) args[3] else "_Q31$"
hash <- function(p) digest::digest(file = p, algo = "sha256")
manifest <- jsonlite::read_json(file.path(parent, "manifest.json"))
stopifnot(all(vapply(names(manifest), function(n) identical(hash(file.path(parent,n)), manifest[[n]]), logical(1))))
stopifnot(as.character(packageVersion("TAM")) == "4.3.25")
dir.create(out, recursive = TRUE)
write_new <- function(value, name) {
  path <- file.path(out, name); stopifnot(!file.exists(path))
  jsonlite::write_json(value, path, digits = NA, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
}
protocol <- list(classification = "OBSERVED_DEVELOPMENT_ONLY", scientific_inference_ready = FALSE,
  parent_manifest_sha256 = hash(file.path(parent, "manifest.json")), TAM = "4.3.25",
  grids = list(list(id = "b12_h10", bound = 12, spacing = .1), list(id = "b16_h05", bound = 16, spacing = .05)),
  start = "all native structural coordinates zero; population mean zero, SD one",
  method = "tam.mml with explicit A/B; item.elim=FALSE preserves declared categories",
  maxiter = 2000L, conv = 1e-10, convD = 1e-12, convM = 1e-12, Msteps = 20L,
  limits = list(same_point_grid_nll_per_person = 1e-8, point_difference = 1e-5,
    continuous_nll_per_person = 1e-8, continuous_score_supnorm = 1e-4,
    continuous_joint_se_relative = .005),
  se_scope = "Native TAM SE saved but not equated with full joint inverse information",
  failure_policy = "Retain both planned grids for every dataset; no replacement, retuning or best-fit selection")
write_new(protocol, "protocol.json")
file.copy("validation/mml_rasch_tam_refit.R", file.path(out, "source.R"))
write_new(list(R = R.version.string, TAM = as.character(packageVersion("TAM")), session = capture.output(sessionInfo())), "runtime.json")
source_names <- c("tam.mml", "tam_mml_mstep_regression", "tam_calc_posterior", "tam_mml_se_quick",
  "tam_mml_calc_prob", "tam_mml_person_posterior")
write_new(setNames(lapply(source_names, function(n) deparse(get(n, asNamespace("TAM")))), source_names), "installed_sources.json")
cases <- jsonlite::read_json(file.path(parent, "reference_input.json"), simplifyVector = FALSE)$cases
cases <- cases[grepl(case_pattern, vapply(cases, `[[`, character(1), "id"))]
stopifnot(length(cases) > 0L, !anyDuplicated(vapply(cases, `[[`, character(1), "id")))
expected_fits <- 2L*length(cases)
records <- list(); failures <- list()
for (case in cases) {
  id <- sub("_Q31$", "", case$id); cfg <- case$config
  stopifnot(cfg$model %in% c("RSM", "PCM"))
  data <- read.csv(file.path(parent, paste0(id, ".csv")), stringsAsFactors = FALSE)
  file.copy(file.path(parent, paste0(id, ".csv")), file.path(out, paste0(id, ".csv")))
  facet_names <- unlist(cfg$facet_names)
  items <- unique(data[facet_names]); items <- items[do.call(order, items), , drop = FALSE]
  item_key <- function(d) do.call(paste, c(d[facet_names], sep = "|"))
  item_index <- match(item_key(data), item_key(items))
  persons <- sort(unique(data$Person)); person_index <- match(data$Person, persons)
  stopifnot(!anyDuplicated(paste(person_index, item_index)))
  response <- matrix(NA_real_, length(persons), nrow(items), dimnames = list(persons, paste0("cell", seq_len(nrow(items)))))
  response[cbind(person_index,item_index)] <- data$Score-cfg$rating_min
  nk <- cfg$n_cat
  counts <- vapply(facet_names, function(f) length(cfg$facet_levels[[f]])-as.integer(f != cfg$noncenter_facet), integer(1))
  step_levels <- if (cfg$model == "RSM") "shared" else unlist(cfg$facet_levels[[cfg$step_facet]])
  dimension <- sum(counts)+length(step_levels)*(nk-2)
  A <- array(0, c(nrow(items), nk, dimension), dimnames = list(colnames(response), as.character(0:(nk-1)), paste0("native",seq_len(dimension))))
  B <- array(rep(rep(0:(nk-1), each=nrow(items)), 1), c(nrow(items), nk, 1))
  cursor <- 0L
  for (f in facet_names) {
    levels <- unlist(cfg$facet_levels[[f]])
    for (j in seq_len(counts[f])) {
      effect <- as.numeric(items[[f]] == levels[j])
      if (f != cfg$noncenter_facet) effect <- effect-as.numeric(items[[f]] == tail(levels,1))
      A[,,cursor+j] <- outer(effect*cfg$facet_signs[[f]], 0:(nk-1))
    }
    cursor <- cursor+counts[f]
  }
  for (s in step_levels) {
    keep <- if (cfg$model == "RSM") rep(TRUE, nrow(items)) else items[[cfg$step_facet]] == s
    for (j in seq_len(nk-2)) A[,,cursor+j] <- -outer(as.numeric(keep), as.numeric(0:(nk-1) >= j & 0:(nk-1) < nk-1))
    cursor <- cursor+nk-2
  }
  stopifnot(cursor == dimension, max(response, na.rm=TRUE) == nk-1)
  write_new(list(config=cfg, items=items, response=response, persons=persons, A=A, B=B), paste0(id,"_design.json"))
  for (grid in protocol$grids) {
    key <- paste(id, grid$id, sep="_"); warnings <- character()
    result <- tryCatch({
      fit <- withCallingHandlers(TAM::tam.mml(resp=response, pid=persons, A=A, B=B,
        item.elim=FALSE, constraint="cases", est.variance=TRUE,
        xsi.inits=cbind(seq_len(dimension), rep(0,dimension)),
        beta.fixed=matrix(c(1,1,0), nrow=1), variance.inits=matrix(1,1,1),
        control=list(nodes=seq(-grid$bound,grid$bound,by=grid$spacing), snodes=0L,
          maxiter=protocol$maxiter, conv=protocol$conv, convD=protocol$convD,
          convM=protocol$convM, Msteps=protocol$Msteps, progress=FALSE), verbose=FALSE),
        warning=function(w) { warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning") })
      saveRDS(fit, file.path(out,paste0(key,".rds")))
      stopifnot(identical(dim(fit$A),dim(A)), max(abs(fit$A-A)) < 1e-14,
        identical(dim(fit$B),dim(B)), max(abs(fit$B-B)) < 1e-14,
        max(abs(fit$beta)) < 1e-14, all(fit$pweights == 1), length(fit$xsi$xsi) == dimension,
        fit$variance[1] > 0)
      # Retain cached summaries, and additionally score the returned parameters.
      # An EM object's final E-step can precede its last parameter update.
      ns <- asNamespace("TAM"); nn <- nrow(response); ni <- ncol(response); nq <- nrow(fit$theta)
      fresh <- get("tam_mml_calc_prob",ns)(iIndex=seq_len(ni),A=fit$A,AXsi=matrix(0,ni,nk),
        B=fit$B,xsi=fit$xsi$xsi,theta=fit$theta,nnodes=nq,maxK=nk,recalc=TRUE,use_rcpp=FALSE)
      prior <- dnorm(fit$theta[,1],mean=0,sd=sqrt(fit$variance[1])); prior <- prior/sum(prior)
      posterior <- get("tam_calc_posterior",ns)(rprobs=fresh$rprobs,gwt=matrix(prior,nn,nq,byrow=TRUE),
        resp=fit$resp,nitems=ni,resp.ind.list=fit$resp.ind.list,resp.ind=fit$resp.ind,
        normalization=TRUE,avoid.zerosum=FALSE)$hwt
      stopifnot(all(is.finite(posterior)),max(abs(rowSums(posterior)-1)) < 1e-12)
      fresh_person <- get("tam_mml_person_posterior",ns)(pid=fit$pid,nstud=nn,pweights=fit$pweights,
        resp=fit$resp,resp.ind=fit$resp.ind,snodes=0,hwtE=posterior,hwt=posterior,ndim=1,theta=fit$theta)$person
      list(id=key, point=c(fit$xsi$xsi, log(sqrt(fit$variance[1]))), grid=grid,
        xsi=fit$xsi, variance=fit$variance, beta=fit$beta, AXsi=fit$AXsi,
        theta=fit$theta, weights=fit$pi.k, deviance=fit$deviance, iter=fit$iter,
        terminated_below_cap=fit$iter < protocol$maxiter, warnings=warnings,
        person=fit$person, fresh_person=fresh_person, pid=fit$pid, deviance_history=fit$deviance.history)
    }, error=function(e) list(id=key, error=conditionMessage(e), warnings=warnings))
    write_new(result, paste0(key,".json"))
    if (!is.null(result$error)) failures[[key]] <- result else records[[key]] <- result
    cat(key, if (is.null(result$error)) paste("iterations",result$iter,"sigma",exp(tail(result$point,1))) else result$error, "\n")
    flush.console()
  }
}
write_new(list(classification=protocol$classification, scientific_inference_ready=FALSE,
  expected_fits=expected_fits, completed_fits=length(records), failures=failures,
  all_terminated_below_cap=length(records)==expected_fits && all(vapply(records, `[[`, logical(1), "terminated_below_cap"))), "summary.json")
files <- list.files(out,full.names=TRUE)
write_new(setNames(lapply(files,hash),basename(files)), "manifest.json")
if (length(records) != expected_fits) quit(status=1L)
