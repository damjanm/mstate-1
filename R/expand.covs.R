#' Expand covariates in competing risks dataset in stacked format
#' 
#' Given a competing risks dataset in stacked format, and one or more
#' covariates, this function adds type-specific covariates to the dataset. The
#' original dataset with the type-specific covariates appended is returned.
#' 
#' Type-specific covariates can be used to analyse separate effects on all
#' event types in a single analysis based on a stacked data set (Putter, Fiocco
#' & Geskus (2007) and Geskus (2016)). It is only unambiguously defined for
#' numeric covariates or for explicit codings. Rows that contain the data for
#' that specific event type have the value copied from the original covariate
#' in case it is numeric. In all other rows it has the value zero. If the
#' covariate is a factor, it will be expanded on the design matrix given by
#' \code{\link[stats:model.matrix]{model.matrix}}. For standard "treatment
#' contrasts" this means that dummy variables are created.  If the covariate is
#' a factor, the column name combines the name of the covariate with the
#' specific event type. If \code{longnames}=\code{TRUE}, both parts are
#' intersected by the specific labels in the coding. Missing values in the
#' basic covariates are allowed and result in missing values in the expanded
#' covariates.
#' 
#' @aliases expand.covs expand.covs.default
#' 
#' @return An data frame object of the same class as the data argument,
#' containing the design matrix for the type-specific covariates, either on its
#' own (\code{append}=\code{FALSE}) or appended to the data
#' (\code{append}=\code{TRUE}).
#' @author Ronald Geskus and Hein Putter \email{H.Putter@@lumc.nl}
#' @seealso \code{\link{expand.covs.msdata}}.
#' @references Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics:
#' Competing risks and multi-state models. \emph{Statistics in Medicine}
#' \bold{26}, 2389--2430.
#' 
#' Geskus, Ronald B. (2016). \emph{Data Analysis with Competing Risks and
#' Intermediate States.} CRC Press, Boca Raton.
#' @keywords datagen
#' @examples
#' 
#' # small data set in stacked format
#' tg <- data.frame(time=c(5,5,1,1,9,9),status=c(1,0,2,2,0,1),failcode=rep(c("I","II"),3),
#'         x1=c(1,1,2,2,2,2),x2=c(3,3,2,2,1,1))
#' tg$x1 <- factor(tg$x1,labels=c("male","female"))
#' # expanded covariates
#' expand.covs(tg,covs=c("x1","x2"))
#' expand.covs(tg,covs=c("x1","x2"),longnames=TRUE)
#' expand.covs(tg,covs=c("x1","x2"),append=FALSE)
#' 
#' @inheritParams expand.covs.msdata
#' 
#' @export 
expand.covs <- function(data, ...) UseMethod("expand.covs")

#' @inheritParams expand.covs.msdata
#' @export
expand.covs.default <-
    function (data, covs, append = TRUE, longnames = FALSE, event.types="failcode", ...)
{
    comp.risks <- unique(data[[event.types]])
    if(length(comp.risks)==1)
        stop("Function does not create type-specific covariates with only one event type")
    trans <- trans.comprisk(K=length(comp.risks), names=comp.risks)
    trans2 <- to.trans2(trans)
    K <- nrow(trans2)
    if (is.character(covs))
        form1 <- as.formula(paste("~ ", paste(covs, collapse = " + ")))
    else form1 <- as.formula(covs)
    mm4 <- NULL
    for (j in 1:length(covs)) {
        wh <- which(!is.na(data[[covs[j]]]))
        form1 <- as.formula(paste("~ ", covs[j]))
        mm <- model.matrix(form1, data = data)
        mm <- data.frame(mm)
        mm <- mm[, -1, drop = FALSE]
        if (!longnames) {
            nc <- ncol(mm)
            if (nc == 1)
                names(mm) <- covs[j]
            else names(mm) <- paste(covs[j], 1:nc, sep = "")
        }
        nms <- names(mm)
        ms <- data.frame(trans = data[[event.types]])
        ms$trans <- factor(ms$trans, levels=comp.risks)
        ms <- cbind(ms[wh, , drop = FALSE], mm)
        mm2 <- model.matrix(as.formula(paste("~ (", paste(nms,
                                                          collapse = " + "), "):trans")), data = ms)[, -1]
        mm3 <- matrix(NA, nrow(data), ncol(mm2))
        mm3[wh, ] <- mm2
        mm3 <- data.frame(mm3)
        nms <- as.vector(t(outer(nms, comp.risks, "paste", sep = ".")))
        names(mm3) <- nms
        if (j == 1)
            mm4 <- mm3
        else mm4 <- cbind(mm4, mm3)
    }
    if (!append)
        return(mm4)
    else {
        if (!all(is.na(match(names(data), nms))))
            warning("One or more names of appended data already in data!")
        mm4 <- cbind(data, mm4)
    }
    class(mm4) <- class(data)
    return(mm4)
}





#' Expand covariates in multi-state dataset in long format
#' 
#' Given a multi-state dataset in long format, and one or more covariates, this
#' function adds transition-specific covariates, expanding the original
#' covariate(s), to the dataset. The original dataset with the
#' transition-specific covariates appended is returned.
#' 
#' For a given basic covariate \code{Z}, the transition-specific covariate for
#' transition \code{s} is called \code{Z.s}. The concept of transition-specific
#' covariates in the context of multi-state models was introduced by Andersen,
#' Hansen & Keiding (1991), see also Putter, Fiocco & Geskus (2007). It is only
#' unambiguously defined for numeric covariates or for explicit codings. Then
#' it will take the value 0 for all rows in the long format dataframe for which
#' \code{trans} does not equal \code{s}. For the rows for which \code{trans}
#' equals \code{s}, the original value of \code{Z} is copied. In
#' \code{expand.covs}, when a given covariate is a factor, it will be expanded
#' on the design matrix given by
#' \code{\link[stats:model.matrix]{model.matrix}}. Missing values in the basic
#' covariates are allowed and result in missing values in the expanded
#' covariates.
#' 
#' @param data An object of class \code{"msdata"}, such as output by
#' \code{\link{msprep}}
#' @param covs A character vector containing the names of the covariates in
#' \code{data} to be expanded
#' @param append Logical value indicating whether or not the design matrix for
#' the expanded covariates should be appended to the data (default=\code{TRUE})
#' @param longnames Logical value indicating whether or not the labels are to
#' be used for the names of the expanded covariates that are categorical
#' (default=\code{TRUE}); in case of \code{FALSE} numbers from 1 up to the
#' number of contrasts are used
#' @param \dots Further arguments to be passed to or from other methods. They
#' are ignored in this function.
#' 
#' @return An object of class 'msdata', containing the design matrix for the
#' transition- specific covariates, either on its own
#' (\code{append}=\code{FALSE}) or appended to the data
#' (\code{append}=\code{TRUE}).
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @references Andersen PK, Hansen LS, Keiding N (1991). Non- and
#' semi-parametric estimation of transition probabilities from censored
#' observation of a non-homogeneous Markov process. \emph{Scandinavian Journal
#' of Statistics} \bold{18}, 153--167.
#' 
#' Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics: Competing
#' risks and multi-state models. \emph{Statistics in Medicine} \bold{26},
#' 2389--2430.
#' @keywords datagen
#' @examples
#' 
#' # transition matrix for illness-death model
#' tmat <- trans.illdeath()
#' # small data set in wide format
#' tg <- data.frame(illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
#'         dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
#'         x1=c(1,1,1,2,2,2),x2=c(6:1))
#' tg$x1 <- factor(tg$x1,labels=c("male","female"))
#' # data in long format using msprep
#' tglong <- msprep(time=c(NA,"illt","dt"),
#'         status=c(NA,"ills","ds"),data=tg,
#'         keep=c("x1","x2"),trans=tmat)
#' # expanded covariates
#' expand.covs(tglong,c("x1","x2"),append=FALSE)
#' expand.covs(tglong,"x1")
#' expand.covs(tglong,"x1",longnames=FALSE)
#' 
#' @method expand.covs msdata
#' @export
expand.covs.msdata <- function(data, covs, append=TRUE, longnames=TRUE, ...)
{
    if (!inherits(data, "msdata"))
        stop("'data' must be an 'msdata' object")
    trans <- attr(data, "trans")
    data <- as.data.frame(data)
    trans2 <- to.trans2(trans)
    K <- nrow(trans2)
    if (is.character(covs)) form1 <- as.formula(paste("~ ",paste(covs,collapse=" + ")))
    else form1 <- as.formula(covs)
    # going to apply model.matrix, but NA's are not allowed, so have to deal with that
    mm4 <- NULL
    for (j in 1:length(covs)) {
        wh <- which(!is.na(data[[covs[j]]]))
        form1 <- as.formula(paste("~ ",covs[j]))
        mm <- model.matrix(form1,data=data)
        mm <- data.frame(mm)
        mm <- mm[,-1,drop=FALSE]
        if (!longnames) {
            nc <- ncol(mm)
            if (nc==1) names(mm) <- covs[j]
            else names(mm) <- paste(covs[j],1:nc,sep="")
        }
        nms <- names(mm)
        ms <- data.frame(trans=data[["trans"]])
        ms$trans <- factor(ms$trans)
        ms <- cbind(ms[wh,,drop=FALSE],mm)
        mm2 <- model.matrix(as.formula(paste("~ (",paste(nms,collapse=" + "),"):trans")),data=ms)[,-1]
        mm3 <- matrix(NA,nrow(data),ncol(mm2))
        mm3[wh,] <- mm2
        mm3 <- data.frame(mm3)
        nms <- as.vector(t(outer(nms,1:K,"paste",sep=".")))
        names(mm3) <- nms
        if (j==1) mm4 <- mm3 else mm4 <- cbind(mm4,mm3)
    }
    if (!append) return(mm4)
    else {
        if (!all(is.na(match(names(data),nms))))
            warning("One or more names of appended data already in data!")
        mm4 <- cbind(data,mm4)
    }
    attr(mm4, "trans") <- trans
    class(mm4) <- c("msdata", "data.frame")
    return(mm4)
}
