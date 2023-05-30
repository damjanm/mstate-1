`coxph.relsurv` <- function(formula, data, na.action,
                            split.transitions, ratetable = relsurv::slopop, 
                            time.format = "days", rmap, 
                            init, bwin, centered, cause,
                            ...
){
  # ... other arguments will be passed to coxph
  
  # cause: A vector of the same length as the number of rows in data. 0 for population deaths, 1 for disease specific deaths, 2 (default) for unknown.

  # TO DO:
  # variance?

  Call <- match.call()
  
  if(missing(split.transitions)){
    stop('The split.transitions argument is empty.')
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
  
  ##### #
  # relsurv part:
  relsurv_coef <- list()
  relsurv_var <- list()

  for(st in split.transitions){
    if(!missing(init)){
      init <- rep(init_arg, relsurv_no_of_covs[[as.character(st)]])
    }
    
    if(!missing(cause)){
      cause <- na.omit(cause_arg[data$trans==st])
    }
    
    mod <- relsurv::rsadd(formula = relsurv_formula[[as.character(st)]],
                   data = subset(data, trans==st),
                   ratetable = ratetable, na.action=na.action,
                   method = 'EM', init = init, bwin = bwin,
                   centered = centered, cause = cause,
                   rmap = rmap)
    
 #    mod <- relsurv::rsadd(formula = relsurv_formula[[as.character(st)]],
 #                          data = subset(data, trans==st) %>%
 # mutate(Tstop=Tstop+runif(nrow(.)), x1.2=x1.2+runif(nrow(.)), x1.1=x1.1+runif(nrow(.))),
 #                          ratetable = ratetable, na.action=na.action,
 #                          method = 'EM', init = init, bwin = bwin,
 #                          centered = centered, cause = cause,
 #                          rmap = rmap)
    
    relsurv_coef <- append(relsurv_coef, list(mod$coefficients))
    relsurv_var <- append(relsurv_var, list(mod$var))
  }
  names(relsurv_coef) <- split.transitions
  names(relsurv_var) <- split.transitions
  
  cx$relsurv_coef <- relsurv_coef
  cx$relsurv_var <- relsurv_var
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
  
  class(cx2) <- 'coxph.relsurv'
  
  return(cx2)

  # Liesbeth:
  # coxph - should be use breslow? What about remaining arguments from coxph
  
  # Opombe:
  # Premisli output iz coxph v primeru, ko so vse covariate splittane
  
}


print.coxph.relsurv <- function (x, digits = max(1L, getOption("digits") - 3L), signif.stars = FALSE, 
                                 ...){
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }

  savedig <- options(digits = digits)
  on.exit(options(savedig))
  
  cat("Non-split transitions:\n \n")
  
  coef <- x$coefficients
  se <- sqrt(diag(x$var))
  
  if (is.null(coef) | is.null(se)) 
    stop("Input is not valid")
  if (is.null(x$naive.var)) {
    tmp <- cbind(coef, exp(coef), se, coef/se, pchisq((coef/se)^2, 
                                                      1, lower.tail = FALSE))
    dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
                                         "se(coef)", "z", "p"))
  }
  else {
    nse <- sqrt(diag(x$naive.var))
    tmp <- cbind(coef, exp(coef), nse, se, coef/se, pchisq((coef/se)^2, 
                                                           1, lower.tail = FALSE))
    dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
                                         "se(coef)", "robust se", "z", "p"))
  }
  
  stats::printCoefmat(tmp, digits = digits, P.values = TRUE,
                      has.Pvalue = TRUE, signif.stars = signif.stars, ...)
  
  cat("\nExcess transitions:\n \n")
  
  st_ch <- as.character(x$split.transitions)
  remember_names <- c()
  
  for(st in st_ch){
    coef <- x$relsurv_coef[[st]]
    se <- sqrt(diag(x$relsurv_var[[st]]))
    
    tmp_j <- cbind(coef, exp(coef), se, 
                 coef/se, pchisq((coef/se)^2, 
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
    cat("   (", naprint(omit), ")\n", sep = "")
  invisible(x)
}
