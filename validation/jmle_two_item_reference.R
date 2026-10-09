# Direct joint row likelihood and full numerical Hessian. Independent of the
# Python closed-form estimator. All counts are deterministic formula controls.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
cases <- list(c(1,1),c(3,1),c(1,3),c(40,24))
rows <- lapply(cases,function(counts) {
  y <- rbind(c(1,0),c(0,1))
  nll <- function(x) {
    eta <- outer(x[1:2],c(0,0),'+') + matrix(rep(c(x[3]/2,-x[3]/2),each=2),2,2)
    loss <- pmax(eta,0)+log1p(exp(-abs(eta)))-y*eta
    sum(counts*rowSums(loss))
  }
  grad <- function(x) {
    eta <- outer(x[1:2],c(0,0),'+') + matrix(rep(c(x[3]/2,-x[3]/2),each=2),2,2)
    residual <- plogis(eta)-y
    c(counts*rowSums(residual),sum(counts*(residual[,1]-residual[,2]))/2)
  }
  fit <- optim(c(.3,-.2,0),nll,grad,method='BFGS',hessian=TRUE,
               control=list(maxit=10000,reltol=1e-13))
  original <- list(point=unname(fit$par),nll=fit$value,gradient=unname(grad(fit$par)))
  point <- fit$par
  # Fixed numerical refinement, using only the row-likelihood derivatives.
  # No closed-form estimate, DGP truth, or altered acceptance tolerance enters.
  for (iteration in 1:3) point <- point-solve(optimHess(point,nll,grad),grad(point))
  info <- optimHess(point,nll,grad)
  list(n10=counts[1],n01=counts[2],point=unname(point),nll=nll(point),
       gradient=unname(grad(point)),convergence=fit$convergence,original_BFGS=original,
       fixed_Newton_refinements=3,
       information=unname(info),joint_se=unname(sqrt(diag(solve(info)))),
       information_min_eigenvalue=min(eigen(info,symmetric=TRUE,only.values=TRUE)$values))
})
jsonlite::write_json(list(classification='DETERMINISTIC_JMLE_ROW_LIKELIHOOD_REFERENCE',
                         cases=rows,R_version=R.version.string),args[1],auto_unbox=TRUE,digits=17,pretty=TRUE)
