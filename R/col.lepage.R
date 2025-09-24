col.lepage <- function(x, y) {

  d <- dim(x)[2]
  statistic <- matrix(nrow = d, ncol = 6)
  colnames(statistic) <- c("L0", "L1", "L2", "L3", "L4", "L5") 
  rownames(statistic) <- colnames(x)
  pvalue <- statistic
  
  for ( i in 1:d ) {
    mod <- LePage::lepage(x[, i], y[, i])
    statistic[i, ] <- mod[1, ]
    pvalue[i, ] <- mod[2, ]     
  }
  
  list(statistic = statistic, pvalue = pvalue)
}


