#' Summary method for coxph.relsurv object
#' 
#' Produces a summary of a fitted coxph.relsurv model.
#' @param object The result of a coxph.relsurv fit
#' @param conf.int Level for computation of the confidence intervals. If set to FALSE no confidence intervals are printed
#' @param scale Vector of scale factors for the coefficients, defaults to 1. The printed coefficients, se, and confidence intervals will be associated with one scale unit
#' @param ... For future methods
#' @return An object of class summary.coxph.relsurv
#' 
#' @author Damjan Manevski \email{damjan.manevski@@mf.uni-lj.si}
#' @seealso \code{\link{coxph.relsurv}}
#' 
`summary.coxph.relsurv` <- function (object, conf.int = 0.95, scale = 1, ...) 
{
  cox <- object
  beta <- cox$coefficients * scale
  if (is.null(cox$coefficients)) {
    return(object)
  }
  nabeta <- !(is.na(beta))
  beta2 <- beta[nabeta]
  if (is.null(beta) | is.null(cox$var)) 
    stop("Input is not valid")
  se <- sqrt(diag(cox$var)) * scale
  # if (!is.null(cox$naive.var)) 
  #   nse <- sqrt(diag(cox$naive.var))
  rval <- list(call = cox$call, fail = cox$fail, na.action = cox$na.action, 
               n = cox$n, loglik = cox$loglik)
  if (!is.null(cox$nevent)) 
    rval$nevent <- cox$nevent
  # if (is.null(cox$naive.var)) {
  tmp <- cbind(beta, exp(beta), se, beta/se, stats::pchisq((beta/se)^2, 
                                                    1, lower.tail = FALSE))
  dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", 
                                       "se(coef)", "z", "Pr(>|z|)"))
  # }
  # else {
  #   tmp <- cbind(beta, exp(beta), nse, se, beta/se, stats::pchisq((beta/se)^2, 
  #                                                          1, lower.tail = FALSE))
  #   dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", 
  #                                        "se(coef)", "robust se", "z", "Pr(>|z|)"))
  # }
  rval$coefficients <- tmp
  
  
  # Relsurv part:
  st_ch <- as.character(cox$split.transitions)
  remember_names <- c()
  
  for(st in st_ch){
    coef <- cox$coefficients_relsurv[[st]]
    se <- sqrt(diag(cox$var_relsurv[[st]]))
    
    tmp_j <- cbind(coef, exp(coef), se, 
                   coef/se, stats::pchisq((coef/se)^2, 
                                   1, lower.tail = FALSE))
    if(!exists('tmp_rs')){
      tmp_rs <- tmp_j
    } else{
      tmp_rs <- rbind(tmp_rs, tmp_j)
    }
    remember_names <- c(remember_names, names(coef))
  }
  
  dimnames(tmp_rs) <- list(remember_names, c("coef", "exp(coef)", 
                                             "se(coef)", "z", "Pr(>|z|)"))
  rval$coefficients_relsurv <- tmp_rs
  
  
  if (conf.int) {
    z <- qnorm((1 + conf.int)/2, 0, 1)
    tmp <- cbind(exp(beta), exp(-beta), exp(beta - z * se), 
                 exp(beta + z * se))
    dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)", 
                                         paste("lower .", round(100 * conf.int, 2), sep = ""), 
                                         paste("upper .", round(100 * conf.int, 2), sep = "")))
    rval$conf.int <- tmp
    
    tmp_rs2 <- cbind(exp(tmp_rs[,1]), exp(-tmp_rs[,1]), exp(tmp_rs[,1] - z * tmp_rs[,3]), 
                     exp(tmp_rs[,1] + z * tmp_rs[,3]))
    colnames(tmp_rs2) <- c("exp(coef)", "exp(-coef)", 
                           paste("lower .", round(100 * conf.int, 2), sep = ""), 
                           paste("upper .", round(100 * conf.int, 2), sep = ""))
    rownames(tmp_rs2) <- remember_names
    rval$conf.int.rs <- tmp_rs2
    
  }
  df <- length(beta2)
  logtest <- -2 * (cox$loglik[1] - cox$loglik[2])
  rval$logtest <- c(test = logtest, df = df, pvalue = stats::pchisq(logtest, 
                                                             df, lower.tail = FALSE))
  # rval$sctest <- c(test = cox$score, df = df, pvalue = stats::pchisq(cox$score, 
  #                                                             df, lower.tail = FALSE))
  # rval$rsq <- c(rsq = 1 - exp(-logtest/cox$n), maxrsq = 1 - 
  #                 exp(2 * cox$loglik[1]/cox$n))
  # rval$waldtest <- c(test = as.vector(round(cox$wald.test, 
  #                                           2)), df = df, pvalue = stats::pchisq(as.vector(cox$wald.test), 
  #                                                                         df, lower.tail = FALSE))
  # if (!is.null(cox$rscore)) 
  #   rval$robscore <- c(test = cox$rscore, df = df, pvalue = stats::pchisq(cox$rscore, 
  #                                                                  df, lower.tail = FALSE))
  # rval$used.robust <- !is.null(cox$naive.var)
  # if (!is.null(cox$concordance)) {
  #   rval$concordance <- cox$concordance[6:7]
  #   names(rval$concordance) <- c("C", "se(C)")
  # }
  # if (inherits(cox, "coxphms")) {
  #   rval$cmap <- cox$cmap
  #   rval$states <- cox$states
  # }
  
  # rval$coefficients_relsurv <- cox$coefficients_relsurv
  # rval$var_relsurv <- cox$var_relsurv
  
  class(rval) <- "summary.coxph.relsurv"
  rval
}

#' Print method for summary.coxph.relsurv objects
#' 
#' Produces a printed summary of a fitted coxph.relsurv model
#' @param x The result of a call to summary.coxph.relsurv
#' @param digits significant digits to print
#' @param expand If the summary is for a multi-state coxph fit, print the results in an expanded format
#' @param ... For future methods
#' 
#' @author Damjan Manevski \email{damjan.manevski@@mf.uni-lj.si}
#' @seealso \code{\link{coxph.relsurv}}, \code{\link{summary.coxph.relsurv}}
#' 
`print.summary.coxph.relsurv` <- function (x, digits = max(getOption("digits") - 3, 3), #signif.stars = getOption("show.signif.stars"), 
                                         expand = FALSE, ...) 
{
  if (!is.null(x$call)) {
    cat("Call:\n")
    dput(x$call)
    cat("\n")
  }
  # if (!is.null(x$fail)) {
  #   cat(" Coxreg failed.", x$fail, "\n")
  #   return()
  # }
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  omit <- x$na.action
  cat("  n=", x$n)
  if (!is.null(x$nevent)) 
    cat(", number of events=", x$nevent, "\n")
  else cat("\n")
  if (length(omit)) 
    cat("   (", stats::naprint(omit), ")\n", sep = "")
  # if (nrow(x$coef) == 0) {
  #   cat("   Null model\n")
  #   return()
  # }
  
  if (!is.null(x$coefficients)) {
    cat("\n")
    cat("Non-split transitions:\n \n")
    
    stats::printCoefmat(x$coefficients, digits = digits, signif.stars = FALSE, 
                 ...)
  }
  if (!is.null(x$conf.int)) {
    cat("\n")
    print(x$conf.int)
  }
  
  cat("\n")
  cat("Excess transitions:\n \n")
  
  stats::printCoefmat(x$coefficients_relsurv, digits = digits, signif.stars = FALSE, 
               ...)
  
  if (!is.null(x$conf.int.rs)) {
    cat("\n")
    print(x$conf.int.rs)
  }
  
  # cat("\n")
  # cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
  
  invisible()
}
