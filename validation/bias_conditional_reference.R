# Independent scalar conditional likelihood reference: no Python/app imports.
args <- commandArgs(trailingOnly=TRUE)
library(jsonlite)
library(numDeriv)
input <- fromJSON(args[1], simplifyVector=FALSE)
family_data <- if(is.null(input$family_file)) NULL else read.csv(file.path(dirname(args[1]),input$family_file))
check_cell <- function(cell) {
  # Legacy development requests may supply an exclusion flag. Prospective
  # family requests derive exclusion from the complete observed family data.
  if (is.null(family_data)) {
    if (isTRUE(cell$excluded_extreme)) return(list(id=cell$id,status='excluded_extreme_level'))
    dat <- read.csv(file.path(dirname(args[1]), cell$file))
  } else {
    dat <- family_data[family_data$Rater==cell$rater & family_data$Task==cell$task,,drop=FALSE]
  }
  step_names <- grep('^step[0-9]+$',names(dat),value=TRUE)
  step_names <- step_names[order(as.integer(sub('step','',step_names)))]; nk <- length(step_names)+1L
  stopifnot(nk>=2L,identical(step_names,paste0('step',seq_len(nk-1L))))
  if (!nrow(dat) && !is.null(family_data)) return(list(id=cell$id,status='unobserved_declared_cell'))
  if (!is.null(family_data) && isTRUE(input$omit_extreme)) {
    extreme <- function(x) length(x)>0L && (all(x==0)||all(x==nk-1L))
    if (extreme(family_data$score[family_data$Rater==cell$rater]) ||
        extreme(family_data$score[family_data$Task==cell$task]))
      return(list(id=cell$id,status='excluded_extreme_level'))
  }
  y <- dat$score; a <- dat$slope; w <- dat$weight; k <- seq_len(nk)-1L
  if (!length(y) || !all(y %in% k) || !all(is.finite(a)&a>0) ||
      !all(is.finite(w)&w>=0) || !any(w>0) ||
      length(cell$bound)!=1L || !is.finite(cell$bound) || cell$bound<=0)
    return(list(id=cell$id,status='invalid_input'))
  steps <- as.matrix(dat[, step_names,drop=FALSE])
  cum <- matrix(0,nrow(dat),nk)
  for (j in seq_len(nk-1L)) cum[,j+1L] <- cum[,j]+steps[,j]
  logits0 <- a * (outer(dat$eta, k)-cum)
  state <- function(b) {
    logits <- logits0 + outer(a*b,k)
    maxima <- apply(logits,1,max)
    ex <- exp(logits-maxima); prob <- ex/rowSums(ex)
    nll <- sum(w*(maxima+log(rowSums(ex))-logits[cbind(seq_along(y),y+1)]))
    mean <- drop(prob %*% k)
    variance <- rowSums(prob*(matrix(k,nrow(dat),nk,byrow=TRUE)-mean)^2)
    list(nll=nll, score=sum(w*a*(mean-y)), info=sum(w*a^2*variance), probability=prob)
  }
  answer <- list(id=cell$id, status='interior')
  if (!identical(input$include_points,FALSE)) answer$points <- lapply(c(-.4,0,.6), function(b) {
    val <- state(b)
    val$gradient_fd <- as.numeric(grad(function(x) state(x)$nll,b))
    # Absolute-coordinate differences avoid numDeriv's tiny-coordinate steps.
    # Retained original stencil at 1e-3/r=6 amplified cancellation. This
    # post-result reference recheck uses three absolute Richardson stencils.
    val$information_fd <- vapply(c(.02,.04,.08), function(h)
      as.numeric(hessian(function(u) state(b+u)$nll,0,
        method.args=list(eps=h,d=h,r=4))), numeric(1))
    val
  })
  active <- w>0; bound <- cell$bound
  if (all(y[active]==0) || all(y[active]==nk-1L)) {
    answer$status <- 'no_finite_mle'; return(answer)
  }
  if (state(-bound)$score>=0 || state(bound)$score<=0) {
    answer$status <- 'search_range_limited'; return(answer)
  }
  estimate <- uniroot(function(x)state(x)$score,c(-bound,bound),tol=1e-12)$root
  val <- state(estimate)
  answer$estimate <- estimate; answer$se <- 1/sqrt(val$info)
  null <- state(0); answer$lr_raw_nll <- 2*(null$nll-val$nll)
  lr <- answer$lr_raw_nll
  shift <- a*estimate
  if (max(abs(shift))*(nk-1L)<=.5) {
    mean0 <- drop(null$probability %*% k)
    centered <- shift*(matrix(k,nrow(dat),nk,byrow=TRUE)-mean0)
    cumulant <- log1p(rowSums(null$probability*expm1(centered)))
    lr <- 2*sum(w*(shift*(y-mean0)-cumulant))
  }
  answer$lr <- max(0,lr)
  if (all(w[active]==1)) {
    answer$p <- pchisq(answer$lr,1,lower.tail=FALSE)
    cutoff <- qchisq(.95,1)
    endpoint <- function(bound) {
      f <- function(x) 2*(state(x)$nll-val$nll)-cutoff
      if (f(bound)<=0) return(list(value=bound, limited=TRUE))
      list(value=uniroot(f,sort(c(bound,estimate)),tol=1e-10)$root, limited=FALSE)
    }
    lower <- endpoint(-bound); upper <- endpoint(bound)
    answer$ci <- c(lower$value,upper$value)
    answer$ci_status <- paste(unique(ifelse(c(lower$limited,upper$limited),'limited by search range','ok')),collapse='; ')
  }
  answer
}
case_results <- lapply(input$cases,function(cell) tryCatch(check_cell(cell),error=function(e)
  list(id=cell$id,status='reference_execution_failed',error=conditionMessage(e))))
derived <- NULL; observed_ids <- NULL
if (!is.null(family_data)) {
  present <- vapply(input$cases,function(cell)
    any(family_data$Rater==cell$rater & family_data$Task==cell$task),logical(1))
  observed_ids <- as.list(vapply(input$cases[present],function(cell)cell$id,character(1)))
  p <- vapply(case_results[present],function(x)if(is.null(x$p)) NA_real_ else x$p,numeric(1))
  derived <- list(list(holm=as.list(p.adjust(p,'holm',n=length(p))),BH=as.list(p.adjust(p,'BH',n=length(p)))))
}
result <- list(classification='INDEPENDENT_DEVELOPMENT_REFERENCE', scientific_inference_ready=FALSE,
  R=R.version.string, numDeriv=as.character(packageVersion('numDeriv')),
  cases=case_results, family_observed_ids=observed_ids, derived_adjustments=derived,
  adjustments=lapply(input$families,function(p){
    p <- vapply(p,function(x)if(is.null(x)) NA_real_ else x,numeric(1))
    list(holm=as.list(p.adjust(p,'holm',n=length(p))),BH=as.list(p.adjust(p,'BH',n=length(p))))
  }))
write_json(result,args[2],auto_unbox=TRUE,digits=17,pretty=TRUE,na='null')
