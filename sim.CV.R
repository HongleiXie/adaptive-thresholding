library(parallel);
library(mvtnorm);

# set path to dir
#dir.path <- '~/Dropbox/RA';
dir.path <- '~/MyDocuments';
source(file.path(dir.path, 'sim.find.delta.R'));
source(file.path(dir.path, 'sim.R'));
source(file.path(dir.path, 'generate.model.data.R'));

# set parameter
p <- 200;
measure <- '2';

# N is the fixed integer, refer paper for details
N <- 10;

deltagrid <- expand.grid(fold.id = 1:5, 
                         delta = seq(0, 4, 1/N), 
                         thresholding = c("hard","al"),
                         p = p,
                         n = 100,
                         model = 'm1',
                         measure = measure
                         );
find_delta <- deltagrid;

tos <- Sys.time();
find_delta$loss <- mcmapply(
  sim.find.delta, 
  delta = deltagrid$delta, 
  fold.id = deltagrid$fold.id, 
  thresholding = deltagrid$thresholding,
  n = deltagrid$n,
  p = deltagrid$p,
  model = deltagrid$model,
  measure = deltagrid$measure,
  mc.cores = detectCores() -1
  );
toe <- Sys.time();
print(toe-tos);

# save results
write.table(find_delta, file = file.path(dir.path, 'find_delta.txt'));

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

# print optimal delta
message('the best delta in hard thresholding is ', delta_opt_hard, ' ','under norm', ' ',measure);
message('the best delta in al thresholding is ', delta_opt_al,' ','under norm', ' ', measure);

loss_hard <- mclapply(
  1:100, 
  mc.cores = detectCores()-1,
  FUN = function(x) {
    sim(p = p, n = 100, delta = delta_opt_hard, thresholding = 'hard', model = 'm1')
    }
  );

loss_hard <- matrix(unlist(loss_hard), ncol = 2, byrow = TRUE);


loss_al <- mclapply(
  1:100, 
  mc.cores = detectCores()-1,
  FUN = function(x) {
    sim(p = p, n = 100, delta = delta_opt_al, thresholding = 'al', model = 'm1')
    }
  );

loss_al <- matrix(unlist(loss_al), ncol = 2, byrow = TRUE);

write.table(loss_al, file.path(dir.path, 'CV_delta_loss_al.txt'));
write.table(loss_hard, file.path(dir.path, 'CV_delta_loss_hard.txt'));

results <- function(loss) {
  
  if('loss_hard' == loss) {
  loss_table <- read.table(file.path(dir.path, 'CV_delta_loss_hard.txt'), header = TRUE);
  }else {
    loss_table <- read.table(file.path(dir.path, 'CV_delta_loss_al.txt'), header = TRUE);
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
results('loss_al');
