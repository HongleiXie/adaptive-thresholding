
sim <- function(p = p, n = n, delta = 2, thresholding = 'hard', model = 'm1') {
  
  # compute theta
  get.theta <- function(i = i, j = j) {  
    X_bar_i <- mean(data[, i]);
    X_bar_j <- mean(data[, j]);
    
    sigma <- cov(data)[i,j];
    
    out <- sum(
      unlist(
        lapply(
          X = 1:n, # X = 1:n instead of 1:p
          FUN = function(z) { 
            ((data[z, i] - X_bar_i)*(data[z, j] - X_bar_j) - sigma)^2;
          }
        )
      )
    );
    
    
    return(out/n); # return out/n instead of out/p
    
  }
  
  # compute the estimated covariance matrix
  get.est.cov <- function(delta = delta, thresholding = thresholding) {
    
    datagrid <- expand.grid(i = 1:p, j = 1:p);
    res <- apply(datagrid, 
                 MARGIN = 1, 
                 FUN = function(z){
                   get.theta(i = z["i"], j = z["j"])
                 }
    );
    
    theta <- matrix(res, nrow = p, ncol = p);
    lambda <- delta*(sqrt(theta*log(p)/n));
    
    
    # hard thresholding
    if('hard' == thresholding) {
      cov_est <- sample_cov*(abs(sample_cov) > lambda)*1;
    }
    
    # adaptive lasso thresholding
    if('al' == thresholding) {
      d <- 1- (abs(lambda/sample_cov))^4;
      d <- d*(d > 0);
      cov_est <- sample_cov*d;
    }
    
    
    return(cov_est);
    
  }
  
  if('m1' == model) {cov_true <- m1.data(p = p, n = n); }else{ cov_true <- m2.data(p = p, n = n); }
  
  data <- rmvnorm(n = n, mean = rep(0, p), sigma = cov_true);
  sample_cov <- cov(data);
  
  cov_hard <- get.est.cov(delta = delta, thresholding = thresholding);
  
  return(
    c(
      norm( 
        x = cov_hard - cov_true,
        type = '2' 
      ),
      norm( 
        x = cov_hard - cov_true,
        type = 'F'
      )
    )
  );
  
  
}
