# Author: Honglei Xie
# This is the script for reproducing results in 
# "Adaptive Thresholding for Sparse Covariance Matrix Estimation", T.Cai, W.Liu (2011)
# Section 5.2: Correlation Analysis on Real Data
###################################################################################################

# load SRBC dataset containing 64 samples, 2308 gene expression values
# 4 types of tumors: 23 EWS, 8 BL-NHL, 12 NB, 21 RMS
data <- read.table(file = "~/Dropbox/RA/train.csv", 
                   header = TRUE
                   );
data <- as.data.frame(t(data));
data$group <- c(rep("EWS", 23), rep("BL-NHL", 8), rep("NB", 12), rep("RMS", 21) );
n_class <- c(23, 8, 12, 21);

# as done by Rothman, Levina, and Zhu (2009), we first rank the genes based on F-statistic
###################################################################################################
# FUNCTION: get.F.stat
# PARAMETERS: k is the number of tumor classes; 
#             n is the number of tissue samples; 
#             gene_id is the index of gene in original dataset
# OUTPUTS:    return the F-statistic described in both two papers
###################################################################################################
get.F.stat <- function(k = 4, n = 64, gene_id = gene_id) {
  gene_values <- as.numeric(data[, gene_id]);   
  x_bar <- mean(gene_values);
  mean_class <- aggregate(data[, gene_id], list(data$group), mean)$x;
  var_class <- aggregate(data[, gene_id], list(data$group), var)$x;
  numerator <- lapply( 
                      X = 1:4,
                      FUN = function(j){
                              n_class[j]*(mean_class[j]-x_bar)^2
                            }
                      );
  
  numerator <- (1/(k-1))*sum(unlist(numerator));
  
  denominator <- lapply( 
    X = 1:4,
    FUN = function(j){
      (n_class[j]-1)*var_class[j]
    }
  );
  
  denominator <- (1/(n-k))*sum(unlist(denominator));
  
  F_stat <- numerator/denominator;
  return(F_stat);
  }
# compute the F-statistic for each gene
F_stat <- as.numeric(lapply(
                        X = 1:2308, # the last variable is group label
                        FUN = function(x) {
                          get.F.stat(gene_id = x)
                          } 
                        )
                     );
# rank genes by F-statistic
gene_selected <- c(
                    order(F_stat, decreasing = TRUE)[1:40], # top 40 genes
                    order(F_stat, decreasing = TRUE)[2149:2308] # bottom 160 genes
                    );

# use the filtered dataset
data_filtered <- data[, gene_selected];
colnames(data_filtered) <- paste0("gene_", gene_selected);

# save filtered dataset containing 200 genes
write.table(data_filtered, file = "~/Dropbox/RA/data_filtered.txt");

# obtian adaptive thresholding covariance matrix estimates ########################################
# sample covariance matrix
n_class <- c(23, 8, 12, 21);

sample_cov <- (22*cov(data_filtered[1:23, ]) + 
                7*cov(data_filtered[24:31, ]) +
                11*cov(data_filtered[32:43, ]) + 
                 20*cov(data_filtered[44:64, ])
                  )/(64-4);


# split samples into equal 5 folds and let the proportions of the four types of tumors in each group
# be nearly equal

split.data <- function(tumor = tumor, train.or.test = 'train') {
  SEED <- 12345;
  set.seed(SEED);
  
  if('EWS' == tumor) {
    sample <- sample(1:23, size = 23);
    fold.id <- rep(1:5, length = 23);
  } else{
    if('BL' == tumor) {
      sample <- sample(24:31, size = 8);
      fold.id <- rep(1:5, length = 8);
    } else{
      if('NB' == tumor) {
        sample <- sample(32:43, size = 12);
        fold.id <- rep(1:5, length = 12);
        } else{
          sample <- sample(44:64, size = 21);
          fold.id <- rep(1:5, length = 21);
        
      }
      
    }
    
  } 
  
  train_id <- list();
  test_id <- list();
  
  for(i in 1:5) {
    train_id[[i]] <- sample[fold.id != i]; 
    test_id[[i]] <- sample[fold.id == i]; 
  }
  
  if('train' == train.or.test) { return(train_id) } else { return(test_id) };
  
}

train_id <- mapply(c, split.data('EWS'), split.data('BL'), split.data('NB'), split.data('RMS'), SIMPLIFY = FALSE);
test_id <- mapply(c, split.data('EWS', 'test'), split.data('BL', 'test'), split.data('NB', 'test'), split.data('RMS', 'test'), SIMPLIFY = FALSE);

#check tumor group distribution
num_BL <- rep(NA, 5);
num_EWS <- rep(NA, 5);
num_NB <- rep(NA, 5);
num_RMS <- rep(NA, 5);

get.tumor.freq <- function(pattern = pattern, fold.id = fold.id) {
  for(i in 1:5) {
    out <- length(
      grep(
        pattern = pattern,
        x = rownames(data_filtered)[train_id[[fold.id]]],
        ignore.case = TRUE
      )
    );
  }
  return(out);
  
}

tumor_gp <- c("EWS", "BL", "NB", "RMS");
freq <- matrix(NA, nrow = 5, ncol = 4);

for(i in 1:4) {
  freq[, i] <- unlist(
    lapply( 
      X = 1:5,
      FUN = function(x) { get.tumor.freq(pattern = tumor_gp[i], fold.id = x) }
      )
    );

  }

### START #########################################################################################
data_filtered <- as.matrix(data_filtered);

# compute theta hat matrix
###################################################################################################
# FUNCTION: get.theta
# PARAMETERS: 
#           i: row index; 
#           j: colum index;
#           fold.id: index of fold id if CV is used   
#           use.all: do you use all data or just training data in CV setting?
#
# OUTPUTS:    return the theta_ij which is defined in the paper
###################################################################################################
n_1 <- freq[, 1];
n_2 <- freq[, 2];
n_3 <- freq[, 3];
n_4 <- freq[, 4];

get.theta <- function(i = i, j = j, fold.id = fold.id, use.all = FALSE) {
  if(FALSE == use.all) {
    train <- train_id[[fold.id]];
    n1 <- n_1[fold.id];
    n2 <- n_2[fold.id];
    n3 <- n_3[fold.id];
    n4 <- n_4[fold.id];
    } else{ 
      train <- 1:64;
      n1 <- n_class[1];
      n2 <- n_class[2];
      n3 <- n_class[3];
      n4 <- n_class[4];
      };
  
  data <- data_filtered[train, ];
  X_bar_i <- mean(data[, i]);
  X_bar_j <- mean(data[, j]);

  sigma <- ((n1-1)*cov(data[1:n1, ]) + 
            (n2-1)*cov(data[(n1+1):(n1+n2), ])+
            (n3-1)*cov(data[(n1+n2+1):(n1+n2+n3), ]) + 
            (n4-1)*cov(data[(n1+n2+n3+1):(n1+n2+n3+n4), ])
            )/(n1+n2+n3+n4-4);
  
  sigma <- sigma[i ,j];

  out <- sum(
            unlist(
              lapply(
                    X = 1:length(train),
                    FUN = function(z) { 
                            jj <- train[z]; 
                            ((data_filtered[jj, i] - X_bar_i)*(data_filtered[jj, j] - X_bar_j) - sigma)^2;
                            }
                   )
              )
          );

  
  return(out/length(train));
  
}

###################################################################################################
# FUNCTION: function.of.delta
# PARAMETERS: 
#              delta: the tuning parameter
#              fold.id: index of fold id if CV is used
#              thresholding: thresholding function. Hard thresholding ('hard') or adaptive lasso ('al')
#
# OUTPUTS:    return the corresponding loss measured by Frobenius norm 
#             given certain training data under CV setting. 
###################################################################################################

function.of.delta <- function(delta = delta, fold.id = fold.id, thresholding = thresholding) {
  datagrid <- expand.grid(i = 1:200, j = 1:200, fold.id = fold.id, use.all = FALSE);
  library(parallel);
  res <- mcmapply(get.theta, 
                  datagrid$i, 
                  datagrid$j, 
                  datagrid$fold.id, 
                  datagrid$use.all,
                  mc.cores = detectCores() -1 
                  );

  theta <- matrix(res, nrow = 200, ncol = 200);
  
  lambda <- delta*sqrt(theta*log(200)/length(train_id[[fold.id]]));
  
  n1 <- n_1[fold.id];
  n2 <- n_2[fold.id];
  n3 <- n_3[fold.id];
  n4 <- n_4[fold.id];
  data <- data_filtered[train_id[[fold.id]], ];
  
  # sample variance from the training samples
  sample_sigma_train <- 
            ((n1-1)*cov(data[1:n1, ]) + 
            (n2-1)*cov(data[(n1+1):(n1+n2), ])+
            (n3-1)*cov(data[(n1+n2+1):(n1+n2+n3), ]) + 
            (n4-1)*cov(data[(n1+n2+n3+1):(n1+n2+n3+n4), ])
            )/(n1+n2+n3+n4-4);
  
  # sample variance from the testing samples
  n1 <- 23-n1;
  n2 <- 8-n2;
  n3 <- 12-n3;
  n4 <- 21-n4;
  data <- data_filtered[test_id[[fold.id]], ];
  
  if(1 == n2) { s2 <- 0; } else {s2 <- (n2-1)*cov(data[(n1+1):(n1+n2), ]); }
  
  sample_sigma_test <-  
                ( (n1-1)*cov(data[1:n1, ]) + 
                s2 +
                (n3-1)*cov(data[(n1+n2+1):(n1+n2+n3), ]) + 
                (n4-1)*cov(data[(n1+n2+n3+1):(n1+n2+n3+n4), ])
                 )/(n1+n2+n3+n4-4);
  
  
  # hard thresholding
  if('hard' == thresholding) {
  norm <- norm( 
            x = sample_sigma_train*(abs(sample_sigma_train) > lambda)*1 - sample_sigma_test,
            type = 'F' #Frobenius norm
            );
  norm <- norm^2;
  
  }
  # adaptive lasso thresholding
  if('al' == thresholding) {
    d <- 1- (abs(lambda/sample_sigma_train))^4;
    d <- d*(d > 0);
    norm <- norm( 
      x = sample_sigma_train*d - sample_sigma_test,
      type = 'F'
      );
    norm <- norm^2;
    
  }
  
  return(norm);
  
}

# N is the fixed integer, refer paper for details
N <- 5;
deltagrid <- expand.grid(fold.id = 1:5, delta = seq(0, 4, 1/N), thresholding = c("hard","al"));
find_delta <- deltagrid;
library(parallel);

tos <- Sys.time();
find_delta$loss <- mcmapply(
                      function.of.delta, 
                      deltagrid$delta, 
                      deltagrid$fold.id, 
                      deltagrid$thresholding,
                      mc.cores = detectCores() -1
                      );
toe <- Sys.time();
print(toe-tos);

# save results
write.table(find_delta, file = "~/Dropbox/RA/find_delta.txt");

# find the delta minimizing the loss
# al thresholding 
find_delta_al <- aggregate( 
                          loss ~ delta, 
                          data = find_delta["al" == find_delta$thresholding, ],
                          FUN = mean 
                          );
                 
delta_opt_al <- find_delta_al[which.min(find_delta_al$loss), 'delta'];

# hard thresholding
find_delta_hard <- aggregate( 
  loss ~ delta, 
  data = find_delta["hard" == find_delta$thresholding, ],
  FUN = mean 
  );

delta_opt_hard <- find_delta_hard[which.min(find_delta_hard$loss), 'delta'];

### use the optimal delta to compute the estimated covariance matrix ##############################
###################################################################################################
# FUNCTION: get.est.cov
# PARAMETERS: 
#       delta: the tuning parameter; default is 2
#       thresholding: thresholding function. Hard thresholding ('hard') or adaptive lasso ('al')
#                     default is 'hard'.
#             
# OUTPUTS:    return the estimated covariance matrix described in this paper.
###################################################################################################
get.est.cov <- function(delta = delta, thresholding = thresholding) {
  
  datagrid <- expand.grid(i = 1:200, j = 1:200);
  res <- apply(datagrid, 
               MARGIN = 1, 
               FUN = function(z){
                 get.theta(i = z["i"], j = z["j"], fold.id = NULL, use.all = TRUE)
               }
          );
  
  theta <- matrix(res, nrow = 200, ncol = 200);
  lambda <- delta*(sqrt(theta*log(200)/64));

  
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

cov_hard <- get.est.cov(delta = delta_opt_hard, thresholding = 'hard');
cov_al <- get.est.cov(delta = delta_opt_al, thresholding = 'al');

library(Matrix);
1- nnzero(x = cov_hard)/(200*200);
# 0.6169
1- nnzero(x = cov_al)/(200*200);
# 0.6169
