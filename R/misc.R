`NAfix` <- function(x, subst=-Inf) {
### Written by Christian Hoffmann; propagate last known non-NA value
### Input:
###     x: numeric vector
###     subst: scalar inidicating which value should replace NA
###         if x starts with a series of NA's
### Output:
###     (numeric) vector, with NA's replaced by last known non-NA value,
###         or 'subst'
    spec <- max(x[!is.na(x)])+1
    x <- c(spec,x)
    while (any(is.na(x))) x[is.na(x)] <- x[(1:length(x))[is.na(x)]-1]
    x[x==spec] <- subst
    x <- x[-1]
    x
}

`my.rbind` <- function(mat1, mat2)
{
    ### rbind "extended" to bind matrices with differing number of columns
    m <- max(ncol(mat1), ncol(mat2))
    if (ncol(mat1)<m) mat1 <- cbind(mat1, matrix(NA, nrow(mat1), m-ncol(mat1)))
    if (ncol(mat2)<m) mat2 <- cbind(mat2, matrix(NA, nrow(mat2), m-ncol(mat2)))
    return(rbind(mat1, mat2))
}

`my.cbind` <- function(mat1, mat2)
{
    ### cbind "extended" to bind matrices with differing number of rows
    m <- max(nrow(mat1), nrow(mat2))
    if (nrow(mat1)<m) mat1 <- rbind(mat1, matrix(NA, m-nrow(mat1), ncol(mat1)))
    if (nrow(mat2)<m) mat2 <- rbind(mat2, matrix(NA, m-nrow(mat2), ncol(mat2)))
    return(cbind(mat1, mat2))
}

`my.cbind2` <- function(x, mat)
{
    ### cbind "extended" to bind scalar x with matrix
    ### elements of x are replicated along the rows of mat
    return(cbind(rep(x, nrow(mat)), mat))
}
