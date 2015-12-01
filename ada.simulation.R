	# Author: Honglei Xie
	# This is the script for reproducing results in 
	# "Adaptive Thresholding for Sparse Covariance Matrix Estimation", T.Cai, W.Liu (2011)
	# Section 5.1: Simulation study
	###################################################################################################
	  
	library(mvtnorm);
	library(parallel);
  # need some functions
  #dir.path <- '~/Dropbox/RA';
  dir.path <- '~/MyDocuments';
	
	source(file.path(dir.path, 'sim.R'));
  source(file.path(dir.path, 'generate.model.data.R'));
  
  p <- 200;
  
  loss_hard <- mclapply(
                  1:100, # 100 replications
                  mc.cores = detectCores() -1,
                  FUN = function(x) {
                        sim(p = p, n = 100, thresholding = 'hard', model = 'm1')
                        }
                  );
    
  loss_hard <- matrix(unlist(loss_hard), ncol = 2, byrow = TRUE);
  write.table(loss_hard, file.path(dir.path, 'loss_hard.txt'));
  
# 	loss_al <- mclapply(
# 	  1:100, # 100 replications
# 	  mc.cores = detectCores() -1,
# 	  FUN = function(x) {
# 	    sim(p = p, n = 100, thresholding = 'al', model = 'm1')
# 	    }
# 	  );
# 	
# 	loss_al <- matrix(unlist(loss_al), ncol = 2, byrow = TRUE);
# 	write.table(loss_al, file.path(dir.path, 'loss_al.txt'));
  
  
	results <- function(loss) {
	  
	  if('loss_hard' == loss) {
	    loss_table <- read.table(file.path(dir.path, 'loss_hard.txt'), header = TRUE);
	  }else {
	    loss_table <- read.table(file.path(dir.path, 'loss_al.txt'), header = TRUE);
	  }
	  
	  return(
	    list(  
	      mean_operator = mean(loss_table[, 1]),
	      mean_f = mean(loss_table[, 2]),
	      std_operator = 0.1*sqrt(var(loss_table[, 1])),
	      std_f = 0.1*sqrt(var(loss_table[, 2]))
	      )
	    );
	  
	} 
	
	results('loss_hard');
	#results('loss_al');
  