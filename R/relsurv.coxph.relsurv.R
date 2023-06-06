#' Fit Extended Proportional Hazards Regression Model with Relative Survival
#' 
#' An extension of the survival::coxph function where one can split death-related transitions
#' and model the corresponding excess hazards using relative survival.
#' @param formula A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the Surv function
#' @param data The data used for fitting the model
#' @param na.action A missing-data filter function. This is applied to the model.frame after any subset argument has been used. Default is options()\$na.action
#' @param split.transitions An integer vector containing the numbered transitions that should be split. Use same numbering as in the given transition matrix
#' @param ratetable The population mortality table. A table of event rates, organized as a ratetable object, see for example relsurv::slopop. Default is slopop
#' @param time.format Define the time format which is used in the data. Possible options: c('days', 'years', 'months'). Default is 'days'
#' @param rmap An optional list to be used if the variables in the data are not organized (and named) in the same way as in the ratetable object
#' @param init Vector of initial values of the iteration used when estimating the effects for the excess hazards. Default initial value is zero for all variables
#' @param bwin Controls the bandwidth used for smoothing in the EM algorithm when estimating the effects for the excess hazards. The follow-up time is divided into quartiles and bwin specifies a factor by which the maximum between events time length on each interval is multiplied. The default bwin=-1 lets the function find an appropriate value. If bwin=0, no smoothing is applied
#' @param centered If TRUE, all variables are centered before fitting the EM algorithmg for the excess hazards and the baseline excess hazard is calculated accordingly. Default is FALSE
#' @param cause An optional vector of the same length as the number of cases if cause of death is (partially) given. Use 0 for population deaths, 1 for disease-specific deaths, 2 (default) for unknown. 
#' @param ... Other arguments that will be passed to coxph
#' @return Returns a coxph.relsurv object that contains the coefficient estimates for the extended model
#' 
#' @author Damjan Manevski \email{damjan.manevski@@mf.uni-lj.si}
#' @seealso \code{\link{coxph}}, \code{\link{msfit.relsurv}}, \code{\link[relsurv]{rsadd}}
#' 
#' @export
`coxph.relsurv` <- function(formula, data, na.action,
                            split.transitions, ratetable = relsurv::slopop, 
                            time.format = "days", rmap, 
                            init, bwin, centered, cause,
                            ...
){

  # TO DO:
  # variance?

  Call <- match.call()
  
  # Check split.transitions argument value:
  trans_df <- as.data.frame(unique(data[,c('from', 'to', 'trans')]))
  
  if(missing(split.transitions)){
    stop('The split.transitions argument is empty.')
  }
  else{
    if(inherits(split.transitions, c('numeric', 'integer'))){
      if(!all(split.transitions %in% trans_df$trans)) stop("Invalid transitions used inside argument split.transitions.")
    }
    else stop("Argument split.transitions expects values of class numeric/integer.")
    
    # Check intermediate states:
    
    to_vals <- unique(trans_df$to)[which(unique(trans_df$to) %in% unique(trans_df$from))]
    intermediate_transitions <- trans_df$trans[trans_df$to %in% to_vals]
    
    if(any(split.transitions %in% intermediate_transitions)) stop("You've listed an intermediate transition for which the hazard would have to be split in excess and population hazard. Population mortality tables for intermediate events haven't been implemented in this function. Please include only transitions that go to death states in the split.transitions argument.")
  }
  
  if(missing(centered)){
    centered <- FALSE
  }
  
  if(!missing(rmap)){
    # if(substitution){
      rmap <- substitute(rmap)
    # }
  }
  
  if(!missing(init)){
    init_arg <- init
  }
  
  if(!missing(cause)){
    cause_arg <- cause
  }
  
  # Define time-related objects:
  Year <- 365.241
  Month <- Year/12
  
  time.format.orig <- time.format
  
  if(time.format == "days"){
    if(max(data$time) < 30) warning("Your max time in the data is less than 30 days. If time is not stored in days, please use argument time.format. \n")
  }
  else if(time.format == "years"){
    if(max(data$time) > 100) warning("Your max time in the data is more than 100 years. If time is not stored in years, please use argument time.format. \n")
    
    data$Tstart <- data$Tstart*Year
    data$Tstop <- data$Tstop*Year
    data$time <- data$time*Year

    time.format <- "days" # Fix argument
  }
  else if(time.format == "months"){
    if(max(data$time) > 600) warning("Your max time in the data is more than 600 months. If time is not stored in months, please use argument time.format. \n")
    
    data$Tstart <- data$Tstart*Month
    data$Tstop <- data$Tstop*Month
    data$time <- data$time*Month

    time.format <- "days" # Fix argument
  }
  else{
    stop("Argument time.format should take values in c('days', 'years', 'months').")
  }
  
  ##### #
  # Prepare new formula for coxph:
  av <- deparse(formula[[3]])
  kovs <- strsplit(av, '\\+')[[1]]
  
  # Helper function for matching transitions:
  find_trans <- function(x){
    x <- gsub(" ", "", x, fixed = TRUE)
    substr(x, nchar(x)-1, nchar(x))
  }
  
  # Prepare relsurv formulas:
  relsurv_formula <- list() 
  relsurv_no_of_covs <- list()
  for(st in split.transitions){
    relsurv_kovs_tmp <- kovs[grep(paste0('.', st), sapply(kovs, find_trans), fixed=TRUE)]
    if(length(relsurv_kovs_tmp)==0) stop('In split.transitions you have supplied a transition for which there are no covariates in the formula.')
    relsurv_formula <- append(relsurv_formula, 
                              as.formula(paste0(deparse(formula[[2]]), '~', paste0(relsurv_kovs_tmp, collapse = '+'))))
    relsurv_no_of_covs <- append(relsurv_no_of_covs,
                                 length(relsurv_kovs_tmp))
  }
  names(relsurv_formula) <- split.transitions
  names(relsurv_no_of_covs) <- split.transitions
  

  for(st in split.transitions){
    kovs <- kovs[-grep(paste0('.', st), kovs, fixed=TRUE)]
  }
  
  kovs <- paste0(kovs, collapse='+')
  coxph_formula <- as.formula(paste0(deparse(formula[[2]]), '~', kovs))
  
  cx <- survival::coxph(formula=coxph_formula, data=data,
                        na.action=na.action, ...)
  coxph.object <- cx

  ##### #
  # relsurv part:
  relsurv.coefficients <- list()
  relsurv.var <- list()
  relsurv.mods <- list()

  ie <- 1
  for(st in split.transitions){
    if(!missing(init)){
      init <- rep(init_arg, relsurv_no_of_covs[[as.character(st)]])
    }
    
    if(!missing(cause)){
      cause <- cause_arg[data$trans==st]
    }

    mod <- relsurv::rsadd(formula = relsurv_formula[[as.character(st)]],
                   data = data[data$trans==st,],
                   ratetable = ratetable, na.action=na.action,
                   method = 'EM', init = init, bwin = bwin,
                   centered = centered, cause = cause,
                   rmap = rmap)

    relsurv.coefficients <- append(relsurv.coefficients, list(mod$coefficients))
    relsurv.var <- append(relsurv.var, list(mod$var))
    
    relsurv.mods[[ie]] <- mod
    ie <- ie+1
  }
  names(relsurv.coefficients) <- split.transitions
  names(relsurv.var) <- split.transitions
  names(relsurv.mods) <- split.transitions
  
  cx$relsurv.coefficients <- relsurv.coefficients
  cx$relsurv.var <- relsurv.var
  cx$split.transitions <- split.transitions
  
  cx2 <- cx
  cx2$call <- Call
  
  # Remove coxph objects you do not need:
  cx2$loglik <- NULL
  cx2$score <- NULL
  cx2$iter <- NULL
  cx2$means <- NULL
  cx2$first <- NULL
  cx2$info <- NULL
  cx2$method <- NULL
  cx2$assign <- NULL
  cx2$wald.test <- NULL
  cx2$concordance <- NULL
  cx2$y <- NULL
  cx2$timefix <- NULL
  cx2$formula <- NULL
  cx2$xlevels <- NULL
  # Maybe use them:
  cx2$linear.predictors <- NULL
  cx2$residuals <- NULL
  cx2$terms <- NULL
  
  cx2$coxph.object <- coxph.object
  cx2$relsurv.mods <- relsurv.mods
  
  class(cx2) <- 'coxph.relsurv'
  
  return(cx2)

  # Liesbeth:
  # coxph - should be use breslow? What about remaining arguments from coxph
  
  # Opombe:
  # Premisli output iz coxph v primeru, ko so vse covariate splittane
  
}


#' Print method for coxph.relsurv object
#' 
#' Print method for coxph.relsurv object
#' @param x Object of class coxph.relsurv to be printed
#' @param digits Number of digits
#' @param signif.stars Should significance stars be printed
#' @param ... Further arguments to print
#' 
#' @author Damjan Manevski \email{damjan.manevski@@mf.uni-lj.si}
#' @seealso \code{\link{coxph.relsurv}}
#' 
#' @export 
`print.coxph.relsurv` <- function (x, digits = max(1L, getOption("digits") - 3L), signif.stars = FALSE, 
                                 ...){
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }

  savedig <- options(digits = digits)
  on.exit(options(savedig))
  
  cat("Non-split transitions:\n \n")
  
  if(!is.null(x$coefficients)){
    coef <- x$coefficients
    se <- sqrt(diag(x$var))
    
    if (is.null(x$naive.var)) {
      tmp <- cbind(coef, exp(coef), se, coef/se, stats::pchisq((coef/se)^2, 
                                                               1, lower.tail = FALSE))
      dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
                                           "se(coef)", "z", "p"))
    }
    else {
      nse <- sqrt(diag(x$naive.var))
      tmp <- cbind(coef, exp(coef), nse, se, coef/se, stats::pchisq((coef/se)^2, 
                                                                    1, lower.tail = FALSE))
      dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
                                           "se(coef)", "robust se", "z", "p"))
    }
    
    stats::printCoefmat(tmp, digits = digits, P.values = TRUE,
                        has.Pvalue = TRUE, signif.stars = signif.stars, ...)
  }

  cat("\nExcess transitions:\n \n")
  
  st_ch <- as.character(x$split.transitions)
  remember_names <- c()
  
  for(st in st_ch){
    coef <- x$relsurv.coefficients[[st]]
    se <- sqrt(diag(x$relsurv.var[[st]]))
    
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
                                       "se(coef)", "z", "p"))
  
  stats::printCoefmat(tmp_rs, digits = digits, P.values = TRUE,
                      has.Pvalue = TRUE, signif.stars = signif.stars, ...)
  
  cat('\n')
  
  omit <- x$na.action
  cat("n=", x$n)
  if (!is.null(x$nevent)) 
    cat(", number of events=", x$nevent, "\n")
  else cat("\n")
  if (length(omit)) 
    cat("   (", stats::naprint(omit), ")\n", sep = "")
  invisible(x)
}
