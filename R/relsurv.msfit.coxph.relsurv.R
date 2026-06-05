#' Compute subject-specific transition hazards with (co-)variances from coxph.relsurv 
#' 
#' An extension of the msfit function for coxph.relsurv
#' @param object A coxph.relsurv object
#' @param newdata A data frame with the same variable names as those that appear in the \code{coxph.relsurv} formula
#' @param variance Submitted to msfit for the overall hazards
#' @param vartype Submitted to msfit for the overall hazards
#' @param trans Transition matrix describing the states and transitions in the
#' multi-state model. See \code{trans} in \code{\link{msprep}} for more
#' detailed information
#' @return An object of class \code{"msfit"}, which is a list containing
#' \item{Haz }{A data frame with \code{time}, \code{Haz}, \code{trans},
#' containing the estimated subject-specific hazards for each of the
#' transitions in the multi-state model} \item{varHaz }{A data frame with
#' \code{time}, \code{Haz}, \code{trans1}, \code{trans2} containing the
#' variances (\code{trans1}=\code{trans2}) and covariances
#' (\code{trans1}<\code{trans2}) of the estimated hazards. This element is only
#' returned when \code{variance}=\code{TRUE}} \item{trans}{The extended transition
#' matrix}
#' 
#' @author Damjan Manevski \email{damjan.manevski@@mf.uni-lj.si}
#' @seealso \code{\link{coxph.relsurv}}, \code{\link{msfit.relsurv}}, \code{\link[relsurv]{rsadd}}
#' 
#' @export
`msfit.coxph.relsurv` <- function(object, 
                                newdata, 
                                variance = FALSE,
                                vartype = c("aalen", "greenwood"),
                                trans){

  trans_new <- modify_transMat(trans, object$split.transitions)
  
  # Run variance always (otherwise it doesn't work for larger data sets due to a problem in msfit):
  variance_tmp <- TRUE
  # When fixed, change variance_tmp->variance
  
  msf <- msfit(object$coxph.object, newdata = newdata, variance = variance_tmp, vartype = vartype, trans = trans)
  
  all_trans <- unique(as.numeric(msf$trans))
  all_trans <- all_trans[!is.na(all_trans)]
  all_trans <- sort(all_trans)
  wh <- (all_trans %in% object$split.transitions)
  trans_non_split <- all_trans[!wh]
  trans_split <- all_trans[wh]
  
  Haz_non_split <- subset(msf$Haz, trans %in% trans_non_split)
  if(nrow(Haz_non_split) > 0){
    # Prepare column for new transition numbers
    colnames(Haz_non_split)[3] <- 'trans_old'
    Haz_non_split$trans <- 0
  }
  
  Haz_split <- msf$Haz[0,]
  
  
  
  # Get all transitions:
  transitions <- which(!is.na(trans), arr.ind = TRUE)
  link_trans <- list()
  
  #################### #
  #  Prepare hazards
  #################### #
  
  # We go through all original transitions and match them
  # to the ones in the new transMat (trans_new).
  # We then obtain the new hazards.
  
  for(i in 1:nrow(transitions)){
    # The transition in trans:
    trans_1 <- trans[transitions[i,1], transitions[i,2]]
    
    
    # We deal differently based on the type of transition
    # (whether we have to split the transition or not):
    if(!(trans_1 %in% object$split.transitions)){
      # The adequate transition in trans_new:
      trans_2 <- trans_new[rownames(trans)[transitions[i,1]],
                           colnames(trans)[transitions[i,2]]]
      # Save the linkage:
      link_trans[[ trans_1 ]] <- trans_2

      whi <- (Haz_non_split$trans_old==trans_1)
      Haz_non_split$trans[whi] <- rep(trans_2, sum(whi))
    }
    else{
      # The adequate transition in trans_new:
      trans_2 <- c(trans_new[rownames(trans)[transitions[i,1]],
                             paste0(colnames(trans)[transitions[i,2]], ".p")],
                   trans_new[rownames(trans)[transitions[i,1]],
                             paste0(colnames(trans)[transitions[i,2]], ".e")])
      # Save the linkage:
      link_trans[[ trans_1 ]] <- trans_2
      
      
      # Prepare pop. and excess hazard objects:
      df_p <- df_e <- subset(msf$Haz, trans==trans_1)
      
      # Take the subset we need:
      # df_subset <- data[(data$from == transitions[i,1]) &
      #                     (data$to == transitions[i,2]),]
      df_subset <- df_p[1,]  #DELETE THIS

      # Find first time, when a jump happens:
      wh_jump <- which.max(df_e$Haz>0)
      is_jump <- any(df_e$Haz>0)
      if(!is_jump){
        # if(!(link_trans_ind == TRUE & substitution == FALSE)){ # If it's not called when bootstrapping
          stop(paste0("There are no events occurring in transition ", trans_1, ". Please remove it from the split.transitions argument."))
        # }
      }

      if(nrow(df_subset)==0 | nrow(df_p)==0){
        if(nrow(df_p)>0){
          df_p$trans <- trans_2[1]
          df_e$trans <- trans_2[2]
        }
      } else{
        # Calculate hazards:
        
        newdata_tmp <- newdata[newdata$trans==trans_1,]
        mod_tmp <- object$relsurv.mods[[as.character(trans_1)]]
        mod_names <- names(mod_tmp$coefficients)
        
        # check_covs <- mod_names[!(mod_names %in% colnames(newdata_tmp))]
        # if(length(check_covs)>0) stop(paste0('Please define covariate(s) ', check_covs, ' in the newdata argument.'))
        
        # Here is Owen's temporary solution
        missing_covs <- setdiff(mod_names, colnames(newdata_tmp))
        # Detect common basis expansion patterns (splines, interactions)
        spline_like <- grepl(":", missing_covs) | grepl("ns\\(", missing_covs) | grepl("poly\\(", missing_covs)
        # Warn only about real missing covariates
        missing_real <- missing_covs[!spline_like]
        # Check:
        if(length(missing_real)>0) stop(paste0('Please define covariate(s) ', missing_real, ' in the newdata argument.'))

        predict_tmp <- relsurv::predict.rsadd(mod_tmp, newdata=newdata_tmp)
        
        ####### #
        # Ta del tukaj - moras ga popraviti, zaenkrat vleces case. Namesto to, hoces te find_times dati v predict.rsadd
        find_times <- df_p$time[!(df_p$time %in% predict_tmp$time)]
        df_tmp <- data.frame(time=find_times, Haz.e=NA, Haz.p=NA)

        predict_tmp2 <- rbind(predict_tmp, df_tmp)
        
        predict_tmp2 <- predict_tmp2[order(predict_tmp2$time),]
        predict_tmp2$Haz.e <- NAfix(predict_tmp2$Haz.e, 0)
        predict_tmp2$Haz.p <- NAfix(predict_tmp2$Haz.p, 0)
        ####### #
        
        # # Calculate hazards at the wanted times
        # wh <- which(Hazs$time %in% (df_p$time))
        # wh_l <- c(NA, wh[1:(length(wh)-1)])+1
        # wh_l[1] <- 1
        # 
        # haz.pop <- sapply(1:length(wh),
        #                   function(x) sum(Hazs$haz.pop[wh_l[x]:wh[x]]))
        # haz.excess <- sapply(1:length(wh),
        #                      function(x) sum(Hazs$haz.excess[wh_l[x]:wh[x]]))
        # 
        # Checks:
        # plot(1:length(haz.pop), (df_p$Haz[1:length(haz.pop)] - cumsum(haz.pop + haz.excess)), type="l")
        # plot(1:length(haz.pop), cumsum(haz.pop), type="l")
        # plot(1:length(haz.pop), cumsum(haz.excess), type="l")

        # Population hazards:
        df_p$trans <- trans_2[1]
        df_p$Haz <- predict_tmp2$Haz.p

        # Excess hazards:
        df_e$trans <- trans_2[2]
        df_e$Haz <- predict_tmp2$Haz.e
        # Check: plot(1:length(haz.pop), (Haz[Haz$trans == trans_1, "Haz"] - df_e$Haz - df_p$Haz), type="l")

        # Make the times fully equal as in the msfit object:
        old_times <- msf$Haz$time[msf$Haz$trans == trans_1]
        df_p$time <- old_times
        df_e$time <- old_times

        # Check if cum. excess haz<0 after first event time:
        # if(wh_jump<nrow(df_e) & substitution & any(df_e$Haz[(wh_jump+1):nrow(df_e)]<0)){
        #   warning(paste0("The relative survival assumption that the observed hazard can be split in population and excess components might not hold for transition ", trans_1, ". Please check the estimated hazards. Consider removing transition ", trans_1, " from the split.transitions argument. \n"))
        # }
      }

      # Save:
      Haz_split <- rbind(Haz_split, df_p, df_e)
    }
  }
  
  # Remove old column:
  if(nrow(Haz_non_split) >0){
    Haz_non_split$trans_old <- NULL
  }

  Haz_new <- rbind(Haz_non_split, Haz_split)
  ordering <- order(Haz_new[,"trans"])
  Haz_new <- Haz_new[ordering,,drop=FALSE]
  rownames(Haz_new) <- 1:nrow(Haz_new)
  
  
  

  
  # sf0 <- summary(survfit(object))
  
  
  
  # Pripraviti Haz: torej rabis time, Haz, trans.
  # potem dobiti napovedi na podlagi modelov (rsadd hmmm)
  # pol se varianca? zaenkrat prazna I guess
  # Pogledaj se argumente, k jih je treba dodati (msfit/msfit.relsurv), recimo variance=FALSE
  
  #################### #
  # 4. Return objects:
  #################### #
  
  # Perform bootstrap if needed:
  if(!is.null(object$boot.objects)){
    all_times <- subset(Haz_new, trans==1)$time
    no_trans <- max(trans_new, na.rm = TRUE)
    
    var_obj <- matrix(NA, length(object$boot.objects), length(all_times)*no_trans)
    
    Haz.boot <- vector('list', length=length(object$boot.objects))
    
    for(ie in 1:length(object$boot.objects)){
      tmp0 <- msfit.coxph.relsurv(object$boot.objects[[ie]],
                          newdata, variance=FALSE, trans=trans)
      Haz.boot[[ie]] <- tmp0
      tmp <- tmp0$Haz
      
      these_times <- subset(tmp, trans==1)$time
      
      diff_times <- all_times[!(all_times %in% these_times)]
      
      full_Haz <- tmp
      # If needed, add times:
      if(length(diff_times) > 0){
        full_Haz <- rbind(tmp,
                          data.frame(time=rep(diff_times, no_trans),
                                     Haz=NA,
                                     trans=rep(1:no_trans, each=length(diff_times))
                          ))
        full_Haz <- full_Haz[order(full_Haz$trans, full_Haz$time), ]
        rownames(full_Haz) <- NULL
        full_Haz$Haz <- NAfix(full_Haz$Haz, 0)
        
        full_Haz <- full_Haz[full_Haz$time %in% all_times,]
      }
      
      var_obj[ie,] <- full_Haz$Haz
      
    }
    
    colVars <- function (x, na.rm = FALSE) {
      f <- function(v, na.rm = na.rm) {
        if (is.numeric(v) || is.logical(v) || is.complex(v)) 
          stats::var(v, na.rm = na.rm)
        else NA
      }
      return(unlist(lapply(x, f, na.rm = na.rm)))
    }
    
    varHaz_new <- data.frame(time=rep(all_times, no_trans),
                             varHaz=colVars(data.frame(var_obj)),
                             trans1=rep(1:no_trans, each=length(all_times)),
                             trans2=rep(1:no_trans, each=length(all_times))
    )
  }
  
  # Save the new values:
  if(variance){
    if(!exists('varHaz_new')){
      stop('Bootstrap option has to be run in coxph.relsurv.')
    }
    
    res <- list(Haz=Haz_new,varHaz=varHaz_new,trans=trans_new, Haz.boot=Haz.boot)
  } else{
    res <- list(Haz=Haz_new,trans=trans_new)
  }
  
  # If link_trans and bootstrap replications have to be added:
  # if(link_trans_ind) res$link_trans <- link_trans
  # if(bootstrap_ind) res$Haz.boot <- haz_boot
  
  class(res) <- "msfit"
  return(res)
}
