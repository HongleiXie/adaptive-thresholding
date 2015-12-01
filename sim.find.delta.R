
sim.find.delta <- function(delta = delta, 
                           fold.id = fold.id, 
                           p = p, 
                           n = n, 
                           thresholding = 'hard', 
                           model = 'm1',
                           measure = 'f'
                           ) {
  
  # split data into 5 folds
  #set.seed(12345);
  sample <- sample(1:n, size = n, replace = FALSE);
  fold <- rep(1:5, length = n);
  
  train_id <- list();
  test_id <- list();
  
  for(i in 1:5) {
    train_id[[i]] <- sample[fold != i]; 
    test_id[[i]] <- sample[fold == i]; 
  }
  
  if('m1' == model) {cov_true <- m1.data(p = p, n = n); }else{ cov_true <- m2.data(p = p, n = n); }
  
  data <- rmvnorm(n = n, mean = rep(0, p), sigma = cov_true);
  sample_cov <- cov(data);
  
  train <- train_id[[fold.id]];
  data_train <- data[train, ];
  
  # compute theta
  get.theta <- function(i = i, j = j) {
    
    X_bar_i <- mean(data_train[, i]);
    X_bar_j <- mean(data_train[, j]);
    sigma <- cov(data_train)[i,j];
    
    out <- sum(
      unlist(
        lapply(
          X = 1:length(train),
          FUN = function(z) { 
            jj <- train[z]; 
            ((data[jj, i] - X_bar_i)*(data[jj, j] - X_bar_j) - sigma)^2;
          }
        )
      )
    );
    
    
    return(out/length(train));
    
  }
  
  # find the best delta
  datagrid <- expand.grid(i = 1:p, j = 1:p);
  res <- mcmapply(get.theta, 
                  datagrid$i, 
                  datagrid$j,
                  mc.cores = detectCores() -1 
                  );
  
  theta <- matrix(res, nrow = p, ncol = p);
  lambda <- delta*sqrt(theta*log(p)/length(train));
  
  # sample variance from the training samples
  sample_sigma_train <- cov(data_train);
  
  # sample variance from the testing samples
  
  sample_sigma_test <- cov(data[test_id[[fold.id]], ]);
  
  # hard thresholding
  if('hard' == thresholding) {
    if('f' == measure) {
    norm_loss <- norm( 
      x = sample_sigma_train*(abs(sample_sigma_train) > lambda)*1 - sample_sigma_test,
      type = 'F' #Frobenius norm
      );
    }else{
      norm_loss <- norm( 
      x = sample_sigma_train*(abs(sample_sigma_train) > lambda)*1 - sample_sigma_test,
      type = '2'
      );
    }
    
  }
  
  # adaptive lasso thresholding
  if('al' == thresholding) {
    d <- 1- (abs(lambda/sample_sigma_train))^4;
    d <- d*(d > 0);
    
    if('f' == measure) {
      norm_loss <- norm( 
        x = sample_sigma_train*(abs(sample_sigma_train) > lambda)*1 - sample_sigma_test,
        type = 'F'
        );
    }else{
      norm_loss <- norm( 
        x = sample_sigma_train*(abs(sample_sigma_train) > lambda)*1 - sample_sigma_test,
        type = '2'
        );
    }

  }
  
  return(norm_loss);
  
}

