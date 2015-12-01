
m1.data <- function(n = n, p = p) {
  
  datagrid <- expand.grid(i = 1:(p/2), j = 1:(p/2));
  
  get.A1 <- function(i = i, j = j) {
  sigma_ij <- 1-abs(i-j)/10;
  sigma_ij <- sigma_ij*(sigma_ij > 0);
  return(sigma_ij);
  }

  A1 <- mapply(
    get.A1, 
    datagrid$i, 
    datagrid$j
    );

  A1 <- matrix(A1, nrow = p/2, ncol = p/2);
  A2 <- diag(x = 4, nrow = p/2, ncol = p/2);

  cov_true <- matrix(0, nrow = p, ncol = p);
  cov_true[1:(p/2), 1:(p/2)] <- A1;
  cov_true[((p/2)+1): p, ((p/2)+1): p] <- A2;
  
  return(cov_true);
  
}

m2.data <- function(n = n, p = p) {
  
  datagrid <- expand.grid(i = 1:(p/2), j = 1:(p/2));
  
  get.B <- function(i = i, j = j) {
    b_ij <- runif(n = 1, min = 0.3, max = 0.8)*(rbinom(n = 1, size = 1, prob = 0.2));
    return(b_ij);
  }
  
  B <- mapply(
    get.B, 
    datagrid$i, 
    datagrid$j
    );
  
  B <- matrix(B, nrow = p/2, ncol = p/2);
  
  #################################################################################################
  ### QUESTION: what if the eigenvalues are not real numbers? #####################################
  #################################################################################################
  
  epi <- max(-eigen(B)$values, 0) + 0.01;
  A1 <- B + diag(epi, p/2);
  A2 <- diag(x = 4, nrow = p/2, ncol = p/2);
  
  cov_true <- matrix(0, nrow = p, ncol = p);
  cov_true[1:(p/2), 1:(p/2)] <- A1;
  cov_true[((p/2)+1): p, ((p/2)+1): p] <- A2;
  
  return(cov_true);
  
}
