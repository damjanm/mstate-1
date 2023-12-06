`msfit.aalen` <- function(object, newdata, variance=TRUE, vartype=c("aalen","greenwood"), trans){
  
  times <- object$coefficients[[1]][,1]
  lcoef <- length(object$coefficients)
  
  
  Haz <- data.frame(time=rep(times, lcoef), Haz=NA, trans=rep(1:lcoef, each=length(times)))
  
  zac <- 1
  kon <- length(times)
  
  tr_i <- 1
  for(i in object$coefficients){

    tmp_df <- i[,2:ncol(i)]
    newdata_tmp <- newdata[newdata$trans==tr_i,]
    
    if(ncol(tmp_df)>=2){
      find_names <- colnames(tmp_df)[2:ncol(tmp_df)]
      
      if(ncol(tmp_df)>=3){
        h_tmp <- rep(0, nrow(tmp_df))
        
        for(ie in 1:length(find_names)){
          h_tmp <- h_tmp+newdata_tmp[, find_names[ie]]*tmp_df[, 1+ie]
        } 
        
        Haz_tmp <- h_tmp + tmp_df[,1]
      } else{
        Haz_tmp <- newdata_tmp[, find_names]*tmp_df[,2:ncol(tmp_df)] + tmp_df[,1]
      }

    } else{
      Haz_tmp <- tmp_df[,1]
    }

    
    # Haz_tmp <- rowSums(i[,2:ncol(i)])
    
    Haz$Haz[zac:kon] <- Haz_tmp
    
    zac <- kon+1
    kon <- kon+length(times)
    
    tr_i <- tr_i+1
  }
  
  out <- list(Haz=Haz, varHaz=NA, trans=trans)
  class(out) <- 'msfit'
  
  return(out)
}