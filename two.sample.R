# Author: Honglei Xie
# This is the script for reproducing results in 
# "Two-sample Covariance Matrix Testing and Support Recovery in High-Dimensional and Sparse Settings", 
# T.Cai, W.Liu and Y.Xia (2014)
# Section 5.2: Real data analysis
###################################################################################################

#library(datamicroarray);
#data('singh', package = 'datamicroarray');
#save.image('~/Dropbox/RA/singh.RData');

# load dataset
load("~/Dropbox/RA/singh.RData");
# 52 prostate tumor patients and 50 prostate normal patients
group <- singh$y;
# microarrary data: 102 patients, 12600 gene expressions
data_raw <- singh$x;

# 5000 genes with the largest absolute values of t-statistics are used
tstat <- apply(
            X = data_raw,
            MARGIN = 2,
            FUN = function(x) {
                abs(t.test(x[1:52], x[53:102])$statistic)
              }
            );
data <- data_raw[, order(tstat, decreasing = TRUE)[1:5000]];

### test null: Sigma_tumor = Sigma_normal #########################################################
P <- ncol(data);
N <- nrow(data);
# tumor data
data_tumor <- as.matrix(data[1: 52, ]);
# normal data
data_normal <- as.matrix(data[53: 102, ]);

sigma_tumor <- ((52-1)/52)*cov(data_tumor);
sigma_normal <- ((50-1)/50)*cov(data_normal);

# compute theta
###################################################################################################
# FUNCTION: get.theta
# PARAMETERS: 
#           i: row index; 
#           j: colum index; 
#           data: dataset you use: tumor or normal
#
# OUTPUTS:    return the theta_ij which is defined in the paper
###################################################################################################

get.theta <- function(i = i, j = j, data = data) {
  
  n <- nrow(data);
  X_bar_i <- mean(data[, i]);
  X_bar_j <- mean(data[, j]);
  if (52 == n) { sigma <- sigma_tumor[i,j] } else { sigma <- sigma_normal[i,j] };
  
  out <- 0;
  
  for(ii in 1:n) {
    temp <- ((data[ii, i] - X_bar_i)*(data[ii, j] - X_bar_j) - sigma)^2;
    out <- out + temp;
  }
  
  return(out/n);
  
}

# compute M_n
get.Mn <- function(i = i, j = j) {
  #if(i > j) { stop("i should be no larger than j"); }
  sigma_t <- sigma_tumor[i,j];
  sigma_n <- sigma_normal[i,j];
  theta_tumor <- get.theta(i = i, j = j, data = data_tumor);
  theta_normal <- get.theta(i = i, j = j, data = data_normal);
  M_ij <- (sigma_t - sigma_n)^2/(theta_tumor/52 + theta_normal/50);
  
  return(M_ij);

}
get.j <- function(i) { return(seq(i, P, 1)); }

Mn <- max(
        unlist(
          lapply(
            1:P,
            FUN = function(z) { 
                get.Mn(i = z, j = seq(z, P, 1) )
                }
            )
          )
        );
# 76.43568

# type one extreme value distribution
alpha <- 0.05;
q_alpha <- -log(8*pi) -2*log( log (1/(1-alpha) ));
q_alpha + 4*log(P)-log(log(P));
# 34.6429

#library(QRM);
#qGumbel(p = 0.95, mu = -log(8*pi), sigma = 2)

### support recovery ##############################################################################
# the first and the second method
diag <- lapply(
              1:P,
              FUN = function(z) {
                      get.Mn(i = z, j = z) > 2*log(P)
                      }
              );

#diag[diag == TRUE];
# count is 21; matched with results in the paper

# show number of row
names(diag) <- 1:P;
as.numeric(names(diag[diag == TRUE]));

support_rec <- function(i = i) {
  #i!=j
  out <- get.Mn(i = i, j = seq(i, P, 1)) >= 4*log(P);
  return(any(out));
}

off_diag <- lapply(
                1:P,
                FUN = support_rec
                );
names(off_diag) <- 1:P;
length(as.numeric(names(off_diag[off_diag == TRUE])));
# 736

# results
union(as.numeric(names(diag[diag == TRUE])), as.numeric(names(off_diag[off_diag == TRUE])));
# 745

# the third method: row by row
alpha <- 0.1;
q_alpha <- -log(8*pi) -2*log( log (1/(1-alpha) ));

get.Mi <- function(i = i) {
  #i!=j
  out <- max(get.Mn(i = i, j = seq(1:P)[-i])) >= q_alpha + 4*log(P)-log(log(P));  
  return(out);
  
}


row_by_row <- function(i = i) {
  out <- any(get.Mi(i), (get.Mn(i = i, j = i) >= 2*log(P)) );

  return(as.numeric(out));
}

rr <- lapply(
  1:P,
  FUN = row_by_row
  );

mean(unlist(rr))*P;
# 1124


