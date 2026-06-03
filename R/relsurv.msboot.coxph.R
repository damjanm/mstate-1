#' Bootstrap function in multi-state models
#' 
#' A generic nonparametric bootstrapping function for multi-state models with coxph.relsurv. Based on
#' \code{msboot} function. 
#' 
#' The function \code{msboot.coxph.relsurv} samples randomly with replacement subjects from
#' the original dataset \code{data}. The individuals are identified with
#' \code{id}, and bootstrap datasets are produced by concatenating all selected
#' rows.
#' 
#' @param theta A function of \code{data} and perhaps other arguments,
#' returning the value of the statistic to be bootstrapped; the output of theta
#' should be a scalar or numeric vector
#' @param data An object of class 'msdata', such as output from
#' \code{\link{msprep}}
#' @param B The number of bootstrap replications; the default is taken to be
#' quite small (5) since bootstrapping can be time-consuming
#' @param id Character string indicating which column identifies the subjects
#' to be resampled
#' @param verbose The level of output; default 0 = no output, 1 = print the
#' replication
#' @param ... Any further arguments to the function \code{theta}
#' @return Matrix of dimension (length of output of theta) x B, with b'th
#' column being the value of theta for the b'th bootstrap dataset
#' @author Marta Fiocco, Hein Putter <H.Putter@@lumc.nl>
#' @references Fiocco M, Putter H, van Houwelingen HC (2008). Reduced-rank
#' proportional hazards regression and simulation-based prediction for
#' multi-state models. \emph{Statistics in Medicine} \bold{27}, 4340--4358.
#' @keywords datagen
#' @examples
#' 
#' @export msboot.coxph.relsurv
`msboot.coxph.relsurv` <- function(theta,data,B=5,id="id",verbose=0,...)
{
    if (!inherits(data, "msdata"))
        stop("'data' must be a 'msdata' object")
    trans <- attr(data, "trans")
    ids <- unique(data[[id]])
    n <- length(ids)
    th <- theta(data,...) # actually only used to get the length
    
    res <- vector("list", length = B) # matrix(NA,length(th),B)
    coef_mat <- matrix(NA, length(th$coefficients), B)
    rs_coef_mat <- vector("list", B)
    
    for (b in 1:B) {
        if (verbose>0) {
            cat("\nBootstrap replication",b,"\n")
            flush.console()
        }
        bootdata <- NULL
        bids <- sample(ids,replace=TRUE)
        bidxs <- unlist(sapply(bids, function(x) which(x==data[[id]])))
        bootdata <- data[bidxs,]
        if (verbose>0) {
            print(date())
            print(events(bootdata))
            cat("applying theta ...")
        }
        thstar <- theta(bootdata,...)
        res[[b]] <- thstar
        
        coef_mat[,b] <- thstar$coefficients
        
        names(thstar$relsurv.coefficients) <- NULL
        rs_coef_mat[[b]] <- thstar$relsurv.coefficients
    }
    if (verbose) cat("\n")
    
    colVars <- function (x, na.rm = FALSE) {
      f <- function(v, na.rm = na.rm) {
        if (is.numeric(v) || is.logical(v) || is.complex(v)) 
          stats::var(v, na.rm = na.rm)
        else NA
      }
      return(unlist(lapply(x, f, na.rm = na.rm)))
    }
    
    coef_var <- colVars(data.frame(t(coef_mat)))
    names(coef_var) <- names(thstar$coefficients)
    
    rs_coef_mat <- unlist(rs_coef_mat)
    rs_coef_mat <- data.frame(variable=names(rs_coef_mat), value=rs_coef_mat)
    rs_coef_mat <- aggregate(value ~ variable, data = rs_coef_mat, FUN = var)
    
    return(list(coef_var, rs_coef_mat, res))
}
