# _______________________________________________________
# kernel search on Air passenger data with GP regression.
# The search is performed with base the kernels: SE, PER, 
# LIN and RQ.
# _______________________________________________________

# rm(list=ls(all=TRUE))
# install.packages("kernlab")
library(kernlab)
# install.packages("datasets")
library(datasets)
?AirPassengers# data set description
# Air passenger data set
mat <- data.frame(y=as.matrix(AirPassengers),
                  x=time(AirPassengers))
x <- as.numeric(mat$x); y <- mat$y
?AirPassengers# data set description
# variables for plot 
xtest <- seq(1949.1, 1961 + 3, 0.01)
x_graph <- seq(1949, 1970, 0.08)
colors_gpr <- c("tomato2", "steelblue", "steelblue1", "steelblue2", "steelblue3",
                "steelblue4")
plot_data <- function(){
  plot(x, y, xlim = c(1949.1 + 0.3, 1961 + 1.85 + 0.2),
       ylim = c(104,790), ylab = "", xlab="")
}

# 
# base kernels functions
# 

# three squared exponential (SE) base kernel functions:
k_se_0 <- function(x, y){
  (sigma_var)^2*exp(-0.5*(1/(ell_se_0)^2)*sum((x - y)^2))
}
class(k_se_0) <- "kernel"
k_se_1 <- function(x, y){
  (sigma_var)^2*exp(-0.5*(1/(ell_se_1)^2)*sum((x - y)^2))
}
class(k_se_1) <- "kernel"
k_se_2 <- function(x, y){
  (sigma_var)^2*exp(-0.5*(1/(ell_se_2)^2)*sum((x - y)^2))
}
class(k_se_2) <- "kernel"
# three periodic (PER) base kernel functions: 
k_per_0 <- function(x, y){
  (sigma_var)^2*exp(-2*(sin(pi*sum(x-y)/period_0)^2)/ell_per_0^2)
}
class(k_per_0) <- "kernel"
k_per_1 <- function(x, y){
  (sigma_var)^2*exp(-2*(sin(pi*sum(x-y)/period_1)^2)/ell_per_1^2)
}
class(k_per_1) <- "kernel"
k_per_2 <- function(x, y){
  (sigma_var)^2*exp(-2*(sin(pi*sum(x-y)/period_2)^2)/ell_per_2^2)
}
class(k_per_2) <- "kernel"
# three linear (LIN) base kernel functions:
k_lin_0 <- function(x, y){
  (sigma_var)^2*(sum((x-ell_lin_0)*(y-ell_lin_0)))
}
class(k_lin_0) <- "kernel"
k_lin_1 <- function(x, y){
  (sigma_var)^2*(sum((x-ell_lin_1)*(y-ell_lin_1)))
}
class(k_lin_1) <- "kernel"
k_lin_2 <- function(x, y){
  (sigma_var)^2*(sum((x-ell_lin_2)*(y-ell_lin_2)))
}
class(k_lin_2) <- "kernel"
# three rational quadratic (RQ) base kernel functions:
k_rq_0 <- function(x, y){
  (sigma_var)^2*(1+sum(x-y)^2/(2*alpha_0*ell_rq_0^2))^(-alpha_0)
}
class(k_rq_0) <- "kernel"
k_rq_1 <- function(x, y){
  (sigma_var)^2*(1+sum(x-y)^2/(2*alpha_1*ell_rq_1^2))^(-alpha_1)
}
class(k_rq_1) <- "kernel"
k_rq_2 <- function(x, y){
  (sigma_var)^2*(1+sum(x-y)^2/(2*alpha_2*ell_rq_2^2))^(-alpha_2)
}
class(k_rq_2) <- "kernel"
# default kernel, multiplying by one
k_none <- function(x, y){
  1
}
class(k_none) <- "kernel"

# 
# arrangements of list, variables and functions used in the search
# 
kernels_0 <- c(`k_se_0`, `k_per_0`, `k_lin_0`, `k_rq_0`, `k_none`)
kernels_1 <- c(`k_se_1`, `k_per_1`, `k_lin_1`, `k_rq_1`, `k_none`)
kernels_2 <- c(`k_se_2`, `k_per_2`, `k_lin_2`, `k_rq_2`, `k_none`)
# variables used in the kernel functions
kernel_names <- c("k_se", "k_per", "k_lin", "k_rq")# names of the base kernels
k_se <- c("ell_se"); k_per <- c("ell_per", "period")# names of the hyperparameters
k_lin <- c("ell_lin"); k_rq <- c("ell_rq", "alpha")
ell_se <- c(); ell_per <- c(); period <- c()
ell_lin <- c(); ell_rq <- c(); alpha <- c()
for(i in 1:3){
  ell_se[i] <- paste0("ell_se_", i-1)
  assign(paste0("ell_se_", i-1), 1)
}
for(i in 1:3){
  ell_per[i] <- paste0("ell_per_", i-1)
  assign(paste0("ell_per_", i-1), 1)
}
for(i in 1:3){
  period[i] <- paste0("period_", i-1)
  assign(paste0("period_", i-1), 1)
}
for(i in 1:3){
  ell_lin[i] <- paste0("ell_lin_", i-1)
  assign(paste0("ell_lin_", i-1), 1)
}
for(i in 1:3){
  ell_rq[i] <- paste0("ell_rq_", i-1)
  assign(paste0("ell_rq_", i-1), 1)
}
for(i in 1:3){
  alpha[i] <- paste0("alpha_", i-1)
  assign(paste0("alpha_", i-1), 1)
}
# 
for(i in 1:2){# only products as operations
  assign(paste0("o_", i), 2)
}
for(i in 1:3){# base kernels are set to k_none
  assign(paste0("k_", i-1), 5)
}
# search operators
operations_1 <- c(`+`, `*`)
operations_2 <- c(`+`, `*`)
operation_sign <- c("+", "×")
# nested kernel functions, 
# the first base kernel is multiplied or summed with the second base kernel
k_composition_1 <- function(x, y){
  u = operations_1[[o_1]](kernels_0[[k_0]](x, y), kernels_1[[k_1]](x, y))
  return(u)  
}
class(k_composition_1) <- "kernel"
# main kernel function which is the argument for the gausspr function:
k_composition_2 <- function(x, y){
  u = operations_2[[o_2]](k_composition_1(x, y), kernels_2[[k_2]](x, y))
  return(u)  
}
class(k_composition_2) <- "kernel"
# f. to plot different models
plot_model <- function(col_nr, line_t){
  model <- gausspr(x[1:115], y[1:115], kernel = k_composition_2,
                   var = var_noise)
  pred <- predict(model, x_graph)
  lines(x_graph, pred, col = colors_gpr[col_nr],
        lwd = 3, lty = line_t)
}
# model evaluation function used in the k.search
model_eval <- function(){
  table <- c()
  set.seed(1202)
  model <- gausspr(x[1:115], y[1:115], kernel = k_composition_2,
                   var = var_noise, cross = 2)
  res_1 <- cross(model)
  pred_model <- predict(model, x[116:144])
  MSE <- mean((y[116:144] - pred_model)^2)
  res_2 <- MSE
  table <- c(round(res_1, 3), round(res_2, 3))
  return(table)
}

# 
# kernel search stage I
# 
unique_b_k <- combn(1:4, 2)# all combinations of the base kernels, first stage
same_b_k <- matrix(c(1:4), c(1:4), nrow = 2, ncol = 4)# identical base kernel pairs
k_comb <- cbind(unique_b_k, same_b_k)
k_comb <- as.data.frame(t(k_comb))
k_comb[4,] <- c(k_comb[4,2], k_comb[4,1])
names(k_comb) <- c("k_0", "k_1")

# default variables for all GPR models 
var_noise <- 0.01
sigma_var <- 1
# interval for grid search h.param. tuning
interval_gpr <- seq(0.55, 3.25, 0.45)# defining the possible h.param. values
# interval_gpr <- seq(0.55, 3.25, 0.45*4)# interval for test run
two_h_param <- expand.grid(interval_gpr, interval_gpr)# table with values 
# for a kernel which has two hyperparameters
three_h_param <- expand.grid(interval_gpr, interval_gpr, interval_gpr)
four_h_param <- expand.grid(interval_gpr, interval_gpr, interval_gpr, interval_gpr)
five_h_param <- expand.grid(interval_gpr, interval_gpr,
                            interval_gpr, interval_gpr, interval_gpr)
six_h_param <- expand.grid(interval_gpr, interval_gpr, interval_gpr, interval_gpr,
                           interval_gpr, interval_gpr)# six hyperparameters
search_lengths <- c("two_h_param", "three_h_param", "four_h_param",
                    "five_h_param", "six_h_param")
# print table with results
res_table_I <- data.frame(matrix(NA, nrow = 9, ncol = 20))
names(res_table_I) <- c(1:20)
row.names(res_table_I) <- c("CV", "MSE test","k_0", "h_param_1", "h_param_2",  
                            "o_1", "k_1", "h_param_3", "h_param_4")
# k.search for-loop
count <- 1
for(g in 1:2){# iterates through search operations: sums and products
  o_1 <- g
  for(h in 1:dim(k_comb)[1]){# iterates through table with kernel combinations
    k_0 <- k_comb[h, 1] 
    k_1 <- k_comb[h, 2]
    
    two_h_p <- c(F, F)# variable to identify 2-h.param. kernel 
    two_h_p[1] <- ifelse((k_0 == 2)|(k_0 == 4), T, F)
    two_h_p[2] <- ifelse((k_1 == 2)|(k_1 == 4), T, F)
    
    h_param_k_0 <- get(kernel_names[k_0])
    h_param_k_1 <- get(kernel_names[k_1])
    all_h_p <- c(h_param_k_0, h_param_k_1)
    tot_nr_h_p <- length(c(h_param_k_0, h_param_k_1))
    h_param_grid <- get(search_lengths[tot_nr_h_p-1])
    l_h_param_grid <- dim(h_param_grid)[1]
    
    cv_tr <- c()
    mse_te <- c()
    for(i in 1:l_h_param_grid){# grid search for-loop
      h_p_count <- 1
      for(j in 1:2){# iterates through the nr. of base kernels
        variable <- get(all_h_p[h_p_count])[j]
        assign(variable, h_param_grid[i, h_p_count])
        h_p_count <- h_p_count + 1
        if(two_h_p[j]){# TRUE when the jth base kernel has two hyperparameters
          variable <- get(all_h_p[h_p_count])[j]
          assign(variable, h_param_grid[i, h_p_count])
          h_p_count <- h_p_count + 1
        }
      } 
      eval_res <- model_eval()
      cv_tr[i] <- eval_res[1]
      mse_te[i] <- eval_res[2]
    }# grid search end
    nr <- which.min((cv_tr + mse_te)/2)# choosing the best model
    res_table_I[1, count] <- cv_tr[nr]# storing the CV from the best model
    res_table_I[2, count] <- mse_te[nr]# storing the MSE from the best model
    # grid search results inserted into table:
    col_index <- 1
    variable <- get(all_h_p[col_index])[1]
    assign(variable, h_param_grid[nr, col_index])
    res_table_I[4, count] <- get(variable)# storing one hyperparameter value
    col_index <- col_index + 1
    if(two_h_p[1]){
      variable <- get(all_h_p[col_index])[1]
      assign(variable, h_param_grid[nr, col_index])
      res_table_I[5, count] <-  get(variable)
      col_index <- col_index + 1
    }
    
    variable <- get(all_h_p[col_index])[2]
    assign(variable, h_param_grid[nr, col_index])
    res_table_I[8, count] <-  get(variable)
    col_index <- col_index + 1
    if(two_h_p[2]){
      variable <- get(all_h_p[col_index])[2]
      assign(variable, h_param_grid[nr, col_index])
      res_table_I[9, count] <-  get(variable)
    }
    
    name_1 <- kernel_names[k_0]
    name_2 <- kernel_names[k_1]
    res_table_I[c(3, 7), count] <- c(name_1, name_2)
    res_table_I[6, count] <- operation_sign[o_1]
    count <- count + 1# next column of the table
  }   
}# kernel search first stage end
res_table_I[is.na(res_table_I)] <- "-"
res_table_I[1:10]
res_table_I[11:20]

var_1_I <- as.numeric(res_table_I[1,])
var_2_I <- as.numeric(res_table_I[2,])
# choosing the model with the best kernel structure:
var_3_I <- (var_1_I + var_2_I)/2# the mean of the models in table
var_4_I <- sort(var_3_I)
first_I <- match(var_4_I, var_3_I)[1]
res_table_I[,first_I]
# > res_table_I[,first_I]
# [1] "617.386"  "2320.717" "k_se"     "3.25"     "-"        "+"        "k_per"    "0.55"    
# [9] "1.45"   

k_0 <- which(res_table_I[,first_I][3] == kernel_names)
o_1 <- which(res_table_I[,first_I][6] == operation_sign)
k_1 <- which(res_table_I[,first_I][7] == kernel_names)
# Plot res_I model
h_param_k_0 <- get(kernel_names[k_0])
h_param_k_1 <- get(kernel_names[k_1])
assign(get(h_param_k_0[1])[1],
       as.numeric(res_table_I[,first_I][4]), envir = .GlobalEnv)
if(length(h_param_k_0)==2){
  assign(get(h_param_k_0[2])[1],
         as.numeric(res_table_I[,first_I][5]), envir = .GlobalEnv)
}
assign(get(h_param_k_1[1])[2],
       as.numeric(res_table_I[,first_I][8]), envir = .GlobalEnv)
if(length(h_param_k_1)==2){
  assign(get(h_param_k_1[2])[2],
         as.numeric(res_table_I[,first_I][9]), envir = .GlobalEnv)
}
# plot of models in stage I
par(mfrow=c(1,2), mgp = c(1.3, 0.4, 0), mar = c(1.7, 1.7, 0.2, 0.2))
x_graph <- seq(1949, 1970, 0.01)
plot_model2 <- function(col_nr, line_t, bred){
  model <- gausspr(x[1:115], y[1:115], kernel = k_composition_2,
                   var = var_noise)
  pred <- predict(model, x_graph)
  lines(x_graph, pred, col = colors_gpr[col_nr],
        lwd = bred, lty = line_t)
}
de_beste <- match(var_4_I, var_3_I)

# six of the best models from stage I
plot_data()
for(i in 6:1){
  linje_t <- 1
  bredde <- 2
  i_en <- de_beste[i]
  k_0 <- which(res_table_I[,i_en][3] == kernel_names)
  o_1 <- which(res_table_I[,i_en][6] == operation_sign)
  k_1 <- which(res_table_I[,i_en][7] == kernel_names)
  
  h_param_k_0 <- get(kernel_names[k_0])
  h_param_k_1 <- get(kernel_names[k_1])
  
  assign(get(h_param_k_0[1])[1],
         as.numeric(res_table_I[,i_en][4]), envir = .GlobalEnv)
  if(length(h_param_k_0)==2){
    assign(get(h_param_k_0[2])[1],
           as.numeric(res_table_I[,i_en][5]), envir = .GlobalEnv)
  }
  assign(get(h_param_k_1[1])[2],
         as.numeric(res_table_I[,i_en][8]), envir = .GlobalEnv)
  if(length(h_param_k_1)==2){
    assign(get(h_param_k_1[2])[2],
           as.numeric(res_table_I[,i_en][9]), envir = .GlobalEnv)
  }
  if(i == 1){
    linje_t <- 1
    bredde <- 3.5
  }
  plot_model2(i, linje_t, bredde)
}
abline(v = 1960.917, lwd = 2, lty = 2, col = "grey50")

# 
# kernel search stage II
# 
res_table_II <- data.frame(matrix(NA, nrow = 11, ncol = 8))
names(res_table_II) <- c(1:8)
row.names(res_table_II) <- c("CV", "MSE test","k_0{+,×}k_1", "h_param_1", 
                             "h_param_2", "h_param_3", "h_param_4", 
                             "o_2", "k_2", "h_param_5", "h_param_6")
# k.search for-loop
count <- 1
for(g in 1:2){# iterates through search operations: sums and products
  o_1 <- which(res_table_I[,first_I][6] == operation_sign)
  o_2 <- g
  for(h in 1:4){# iterates through base kernels to expand on kernel
    k_0 <- which(res_table_I[,first_I][3] == kernel_names)
    k_1 <- which(res_table_I[,first_I][7] == kernel_names)
    k_2 <- h
    
    two_h_p <- c(F, F, F)# variable to identify 2-h.param. kernel 
    two_h_p[1] <- ifelse((k_0 == 2)|(k_0 == 4), T, F)
    two_h_p[2] <- ifelse((k_1 == 2)|(k_1 == 4), T, F)
    two_h_p[3] <- ifelse((k_2 == 2)|(k_2 == 4), T, F)
    
    h_param_k_0 <- get(kernel_names[k_0])
    h_param_k_1 <- get(kernel_names[k_1])
    h_param_k_2 <- get(kernel_names[k_2])
    all_h_p <- c(h_param_k_0, h_param_k_1, h_param_k_2)
    tot_nr_h_p <- length(c(h_param_k_0, h_param_k_1,
                           h_param_k_2))
    h_param_grid <- get(search_lengths[tot_nr_h_p-1])
    l_h_param_grid <- dim(h_param_grid)[1]
    
    cv_tr <- c()
    mse_te <- c()
    for(i in 1:l_h_param_grid){# grid search for-loop
      h_p_count <- 1
      for(j in 1:3){# iterates through the nr. of base kernels
        variable <- get(all_h_p[h_p_count])[j]
        assign(variable, h_param_grid[i, h_p_count])
        h_p_count <- h_p_count + 1
        if(two_h_p[j] == TRUE){# TRUE when the jth base kernel has two hyperparameters
          variable <- get(all_h_p[h_p_count])[j]
          assign(variable, h_param_grid[i, h_p_count])
          h_p_count <- h_p_count + 1
        }
      }
      eval_res <- model_eval()
      cv_tr[i] <- eval_res[1]
      mse_te[i] <- eval_res[2]
    }# grid search end
    nr <- which.min((cv_tr + mse_te)/2)# choosing the best model
    res_table_II[1, count] <- cv_tr[nr]# storing the CV from the best model
    res_table_II[2, count] <- mse_te[nr]# storing the MSE from the best model
    # grid search results inserted into table:
    col_index <- 1
    variable <- get(all_h_p[col_index])[1]
    assign(variable, h_param_grid[nr, col_index])
    res_table_II[4, count] <- get(variable)
    col_index <- col_index + 1
    if(two_h_p[1]){
      variable <- get(all_h_p[col_index])[1]
      assign(variable, h_param_grid[nr, col_index])
      res_table_II[5, count] <-  get(variable)
      col_index <- col_index + 1
    }
    
    variable <- get(all_h_p[col_index])[2]
    assign(variable, h_param_grid[nr, col_index])
    res_table_II[6, count] <-  get(variable)
    col_index <- col_index + 1
    if(two_h_p[2]){
      variable <- get(all_h_p[col_index])[2]
      assign(variable, h_param_grid[nr, col_index])
      res_table_II[7, count] <-  get(variable)
      col_index <- col_index + 1
    }
    
    variable <- get(all_h_p[col_index])[3]
    assign(variable, h_param_grid[nr, col_index])
    res_table_II[10, count] <-  get(variable)
    col_index <- col_index + 1
    if(two_h_p[3]){
      variable <- get(all_h_p[col_index])[3]
      assign(variable, h_param_grid[nr, col_index])
      res_table_II[11, count] <-  get(variable)
    }
    
    name_1 <- kernel_names[k_0]
    name_2 <- kernel_names[k_1]
    k_comp_1_name <- paste0(name_1, operation_sign[o_1], name_2)
    res_table_II[3, count] <- paste0("(", k_comp_1_name, ")")
    res_table_II[8, count] <- operation_sign[o_2]
    name_3 <- kernel_names[k_2]
    res_table_II[9, count] <- name_3
    count <- count + 1# next column of the table
  }
}# kernel search second stage end
res_table_II[is.na(res_table_II)] <- "-"
res_table_II

a5 <- as.numeric(res_table_II[1,])
a6 <- as.numeric(res_table_II[2,])
a7 <- (a5 + a6)/2
a8 <- sort(a7)
first_II <- match(a8, a7)[1]
res_table_II[,first_II]
# > res_table_II[,first_II]
# [1] "377.338"      "950.497"      "(k_se+k_per)" "3.25"         "-"           
# [6] "0.55"         "0.55"         "×"            "k_per"        "1.45"        
# [11] "1"  

k_0 <- which(res_table_I[,first_I][3] == kernel_names)
o_1 <- which(res_table_I[,first_I][6] == operation_sign)
k_1 <- which(res_table_I[,first_I][7] == kernel_names)

o_2 <- which(res_table_II[,first_II][8] == operation_sign)
k_2 <- which(res_table_II[,first_II][9] == kernel_names)

# 
# plot res_II model
# 
h_param_k_0 <- get(kernel_names[k_0])
h_param_k_1 <- get(kernel_names[k_1])
h_param_k_2 <- get(kernel_names[k_2])
current_h_par <- c()
col_ind <- 1
assign(get(h_param_k_0[1])[1],
       as.numeric(res_table_II[,first_II][4]), envir = .GlobalEnv)
current_h_par[col_ind] <- get(h_param_k_0[1])[1]
col_ind <- col_ind + 1
if(length(h_param_k_0) == 2){
  assign(get(h_param_k_0[2])[1],
         as.numeric(res_table_II[,first_II][5]), envir = .GlobalEnv)
  current_h_par[col_ind] <- get(h_param_k_0[2])[1]
  col_ind <- col_ind + 1
}
assign(get(h_param_k_1[1])[2],
       as.numeric(res_table_II[,first_II][6]), envir = .GlobalEnv)
current_h_par[col_ind] <- get(h_param_k_1[1])[2]
col_ind <- col_ind + 1
if(length(h_param_k_1) == 2){
  assign(get(h_param_k_1[2])[2],
         as.numeric(res_table_II[,first_II][7]), envir = .GlobalEnv)
  current_h_par[col_ind] <- get(h_param_k_1[2])[2]
  col_ind <- col_ind + 1
}
assign(get(h_param_k_2[1])[3],
       as.numeric(res_table_II[,first_II][10]), envir = .GlobalEnv)
current_h_par[col_ind] <- get(h_param_k_2[1])[3]
col_ind <- col_ind + 1
if(length(h_param_k_2) == 2){
  assign(get(h_param_k_2[2])[3],
         as.numeric(res_table_II[,first_II][11]), envir = .GlobalEnv)
  current_h_par[col_ind] <- get(h_param_k_2[2])[3]
}
current_h_par <- c(current_h_par, "sigma_var")
# > current_h_par
# [1] "ell_se_0"  "ell_per_1" "period_1"  "ell_per_2" "period_2"  "sigma_var"
plot_data()
plot_model(4, 1)
model_eval()
# > model_eval()
# [1] 377.338 950.497

# 
# Hyperparameter tuning function
# 
search_operator <- c(`+`, `-`)
adj_hyperparameter <- function(s, adj, itr, hyperparameter){
  improved <- 0
  eval_res <- model_eval()
  CV_0 <- (eval_res[1] + eval_res[2])/2
  print(paste("The initial score is", CV_0))
  for(g in 1:itr){
    assign(hyperparameter,
           search_operator[[s]](get(hyperparameter), adj),
           envir = .GlobalEnv)
    if(get(hyperparameter) <= 0){
      assign(hyperparameter,
             search_operator[[s%%2+1]](get(hyperparameter),
                                       adj), envir = .GlobalEnv)
      print(paste("Error", hyperparameter, "is negative"))
      break
    }
    eval_res <- model_eval()
    CV_1 <- (eval_res[1] + eval_res[2])/2
    if(CV_1 >= CV_0){
      assign(hyperparameter,
             search_operator[[s%%2+1]](get(hyperparameter),
                                       adj), envir = .GlobalEnv)
      print(paste("the score", CV_1, "is not an improvement"))
      break
    }
    print("improved")
    improved <- 1
    CV_0 <- CV_1
  }
  print(paste("Resulting score is", CV_0, "with hyperparameter",
              hyperparameter, "=", get(hyperparameter)))
  table <- list(CV_0, hyperparameter, improved)
  return(table)
}
# 
# While-loop
# 
nr_c_h <- length(current_h_par)
iter <- 1
stop_adj <- rep(F, nr_c_h)
no_success <- 0
while(iter<9999){
  add_sub <- 1
  attempt <- 0
  first_attempt <- 0
  for(i in 1:5){
    index <- (iter+nr_c_h-1)%%nr_c_h+1
    continue <- adj_hyperparameter(add_sub, 0.001, 5,
                                   current_h_par[index])[3]
    if(continue == 1){
      first_attempt <- 1
      stop_adj[index] <- F
    }
    if(continue == 0 && attempt == 1){
      if(first_attempt == 0){
        stop_adj[index] <- T
      }
      break
    }
    if(continue == 0){
      add_sub <- add_sub%%2 + 1
      attempt <- attempt + 1
    }
  }
  if(all(stop_adj)){
    break
  }
  iter <- iter + 1
}

# 
# final results
# 
plot_model(1,1)
abline(v = 1960.917, lwd = 2, lty = 2, col = "grey50")
res_table_II[c(1,2),first_II]
# > res_table_II[c(1,2),first_II]
# [1] "377.338" "950.497"
model_eval()
# > model_eval()
# [1] 384.539 263.163
new_h_p_values <- t(sapply(current_h_par, function(x) get(x)))
res_table_II[c(4:7,10:11),first_II]
# > res_table_II[c(4:7,10:11),first_II]
# [1] "3.25" "-"    "0.55" "0.55" "1.45" "1"  
new_h_p_values
# > new_h_p_values
#      ell_se_0 ell_per_1 period_1 ell_per_2 period_2 sigma_var
# [1,]   14.824     0.715    0.541     1.376    1.061     1.799



