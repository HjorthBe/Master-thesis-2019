# __________________________________________________________
# kernel search on CO2 data for support vector regression.
# The search is performed with the base kernels: SE and PER.  
# __________________________________________________________

# rm(list=ls(all=TRUE))
# install.packages("kernlab")
library(kernlab)
# install.packages("datasets")
library(datasets)
?CO2# data set description
# require(graphics)
# CO2 data set
mat <- data.frame(y_0 = as.matrix(co2), x_0 = time(co2))
x <- mat$x_0; y <- mat$y_0
# variables for plot 
x_graph <- seq(1959, 2012.5, 0.01) 
colors_svr <- c("tomato", "steelblue", "chartreuse","tomato2", "steelblue1",
                "steelblue2","chartreuse2","tomato3", "steelblue3", "tomato4",
                "chartreuse4", "steelblue4")
plot_data <- function(){
  plot(x, y, xlim = c(1960.5, 2010), ylim = c(314, 385),
       ylab = expression("Atmospheric concentration of CO"[2]),
       xlab = "Time")
}

# 
# Base kernels functions
# 

# SE kernels:
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
# PER kernels:
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
# default kernel, multiplying by one
k_none <- function(x, y){
  1
}
class(k_none) <- "kernel"

# 
# Arrangements of lists, variables and functions used in the search
# 
kernels_0 <- c(`k_se_0`, `k_per_0`, `k_none`)
kernels_1 <- c(`k_se_1`, `k_per_1`, `k_none`)
kernels_2 <- c(`k_se_2`, `k_per_2`, `k_none`)
# variables used in the kernel functions
kernel_names <- c("k_se", "k_per")
k_se <- c("ell_se"); k_per <- c("ell_per", "period")
ell_se <- c(); ell_per <- c(); period <- c();
ell_lin <- c(); ell_rq <- c(); alpha <- c();
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
#
for(i in 1:2){
  assign(paste0("o_", i), 2)
}
for(i in 1:3){
  assign(paste0("k_", i), 3)
}
# search operators
operations_1 <- c(`+`, `*`)
operations_2 <- c(`+`, `*`)
operation_sign <- c("+", "×")
# compositional kernel functions
k_composition_1 <- function(x, y){
  u = operations_1[[o_1]](kernels_0[[k_0]](x, y), kernels_1[[k_1]](x, y))
  return(u)  
}
class(k_composition_1) <- "kernel"
k_composition_2 <- function(x, y){
  u = operations_2[[o_2]](k_composition_1(x, y), kernels_2[[k_2]](x, y))
  return(u)  
}
class(k_composition_2) <- "kernel"
# function to plot different models
plot_model <- function(col_nr, line_w){
  model <- ksvm(x[1:374], y[1:374],
                kernel = k_composition_2, C = cost, 
                epsilon = epsilon_value)
  pred <- predict(model, x_graph)
  lines(x_graph, pred, col = colors_svr[col_nr],
        lwd = line_w, lty = 1)
}
# model evaluation function used in the k.search
model_eval <- function(){
  table <- c()
  set.seed(1202)
  model <- ksvm(x[1:374], y[1:374], epsilon = epsilon_value, C = cost,
                kernel = k_composition_2)
  pred_model <- predict(model, x[1:374])
  mse <- mean((pred_model - y[1:374])^2)
  res_1 <- sqrt(mse)
  pred_model <- predict(model, x[375:468])
  res_2 <- mean((pred_model - y[375:468])^2)
  res_2 <- sqrt(res_2)
  table <- c(round(res_1, 3), round(res_2, 3))
  return(table)
} 

# 
# Kernel search stage I
# 
k_comb <- combn(2, 2)
k_comb <- expand.grid(1:2, 1:2)
k_comb <- k_comb[-2,]
names(k_comb) <- c("k_0", "k_1")
# variables for all SVR models 
epsilon_value <- 0.01
cost <- 1.5 
sigma_var <- 1
# interval for grid search h.param. tuning
interval_svr <- seq(0.1, 8.1, 1.8)
# interval_svr <- seq(0.1, 8.1, 1.8*3)# interval for test run
two_h_param <- expand.grid(interval_svr, interval_svr)
three_h_param <- expand.grid(interval_svr, interval_svr, interval_svr)
four_h_param <- expand.grid(interval_svr, interval_svr, interval_svr, interval_svr)
five_h_param <- expand.grid(interval_svr, interval_svr,
                            interval_svr, interval_svr, interval_svr)
six_h_param <- expand.grid(interval_svr, interval_svr, interval_svr, interval_svr,
                           interval_svr, interval_svr)
search_lengths <- c("two_h_param", "three_h_param", "four_h_param", "five_h_param",
                    "six_h_param")
# table with results
res_table_I <- data.frame(matrix(NA, nrow = 10, ncol = 6))
names(res_table_I) <- c(1:6)
row.names(res_table_I) <- c("RMSE train", "RMSE test","k_0", "h_param_1", "h_param_2",  
                            "o_1", "k_1", "h_param_3", "h_param_4", "C")
res_table_I[10,] <- cost
# k.search for-loop
count <- 1
for(g in 1:2){# iterates through search operations: sums and products
  o_1 <- g
  for(h in 1:dim(k_comb)[1]){# iterates through table with kernel combinations
    k_0 <- k_comb[h, 1] 
    k_1 <- k_comb[h, 2]
    
    two_h_p <- c(F, F)# variable to identify 2-h.param. kernel 
    two_h_p[1] <- ifelse(k_0 == 2, T, F)
    two_h_p[2] <- ifelse(k_1 == 2, T, F)
    
    h_param_k_0 <- get(kernel_names[k_0])
    h_param_k_1 <- get(kernel_names[k_1])
    all_h_p <- c(h_param_k_0, h_param_k_1)
    tot_nr_h_p <- length(c(h_param_k_0, h_param_k_1))
    h_param_grid <- get(search_lengths[tot_nr_h_p-1])
    l_h_param_grid <- dim(h_param_grid)[1]
    
    rmse_tr <- c()
    rmse_te <- c()
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
      rmse_tr[i] <- eval_res[1]
      rmse_te[i] <- eval_res[2]
    }# grid search end
    nr <- which.min((rmse_tr + rmse_te)/2)# choosing the best model
    res_table_I[1, count] <- rmse_tr[nr]
    res_table_I[2, count] <- rmse_te[nr]
    # grid search results inserted into table:
    col_index <- 1
    variable <- get(all_h_p[col_index])[1]
    assign(variable, h_param_grid[nr, col_index])
    res_table_I[4, count] <- get(variable)
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
    count <- count + 1
  }
}# kernel search first stage end
res_table_I[is.na(res_table_I)] <- "-"
res_table_I

variable_1 <- as.numeric(res_table_I[1,])
variable_2 <- as.numeric(res_table_I[2,])
variable_3 <- (variable_1 + variable_2)/2
variable_4 <- sort(variable_3)
first_I <- match(variable_4, variable_3)[1]
res_table_I[,first_I]
# > res_table_I[,first_I]
# [1] "0.931" "2.176" "k_se"  "3.7"   "-"     "+"     "k_per" "0.1"   "1.9"   "1.5"  

k_0 <- which(res_table_I[,first_I][3] == kernel_names)
o_1 <- which(res_table_I[,first_I][6] == operation_sign)
k_1 <- which(res_table_I[,first_I][7] == kernel_names)

# 
# kernel search stage II
# 
res_table_II <- data.frame(matrix(NA, nrow = 12, ncol = 4))
names(res_table_II) <- c(1:4)
row.names(res_table_II) <- c("RMSE train", "RMSE test","k_0{+,×}k_1", "h_param_1", "h_param_2", 
                             "h_param_3", "h_param_4", "o_2", "k_2", "h_param_5",
                             "h_param_6", "C")
res_table_II[12,] <- cost
# k.search for-loop
count <- 1
for(g in 1:2){# iterates through search operations: sums and products
   o_1 <- which(res_table_I[,first_I][6] == operation_sign)
   o_2 <- g
  for(h in 1:2){# iterates through base kernels to expand on kernel
    k_0 <- which(res_table_I[,first_I][3] == kernel_names)
    k_1 <- which(res_table_I[,first_I][7] == kernel_names)
    k_2 <- h
    
    two_h_p <- c(F, F, F)# variable to identify 2-h.param. kernel 
    two_h_p[1] <- ifelse(k_0 == 2, T, F)
    two_h_p[2] <- ifelse(k_1 == 2, T, F)
    two_h_p[3] <- ifelse(k_2 == 2, T, F)
    
    h_param_k_0 <- get(kernel_names[k_0])
    h_param_k_1 <- get(kernel_names[k_1])
    h_param_k_2 <- get(kernel_names[k_2])
    all_h_p <- c(h_param_k_0, h_param_k_1,
                 h_param_k_2)
    tot_nr_h_p <- length(c(h_param_k_0, h_param_k_1,
                           h_param_k_2))
    h_param_grid <- get(search_lengths[tot_nr_h_p-1])
    l_h_param_grid <- dim(h_param_grid)[1]
    
    rmse_tr <- c()
    rmse_te <- c()
    for(i in 1:l_h_param_grid){# grid search for-loop
      h_p_count <- 1
      for(j in 1:3){# iterates through the nr. of base kernels
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
      rmse_tr[i] <- eval_res[1]
      rmse_te[i] <- eval_res[2]
    }# grid search end
    nr <- which.min((rmse_tr + rmse_te)/2)# choosing the best model
    res_table_II[1, count] <- rmse_tr[nr]
    res_table_II[2, count] <- rmse_te[nr]
    # grid search results inserted into table:
    col_index <- 1
    variable <- get(all_h_p[col_index])[1]
    assign(variable, h_param_grid[nr, col_index])
    res_table_II[4, count] <- get(variable)
    col_index <- col_index + 1
    if(two_h_p[1]){
      variable <- get(all_h_p[col_index])[1]
      assign(variable, h_param_grid[nr, col_index])
      res_table_II[5, count] <- get(variable)
      col_index <- col_index + 1
    }
    
    variable <- get(all_h_p[col_index])[2]
    assign(variable, h_param_grid[nr, col_index])
    res_table_II[6, count] <- get(variable)
    col_index <- col_index + 1
    if(two_h_p[2]){
      variable <- get(all_h_p[col_index])[2]
      assign(variable, h_param_grid[nr, col_index])
      res_table_II[7, count] <- get(variable)
      col_index <- col_index + 1
    }
    
    variable <- get(all_h_p[col_index])[3]
    assign(variable, h_param_grid[nr, col_index])
    res_table_II[10, count] <- get(variable)
    col_index <- col_index + 1
    if(two_h_p[3]){
      variable <- get(all_h_p[col_index])[3]
      assign(variable, h_param_grid[nr, col_index])
      res_table_II[11, count] <- get(variable)
    }
    
    name_1 <- kernel_names[k_0]
    name_2 <- kernel_names[k_1]
    k_comp_1_name <- c(paste0(name_1, operation_sign[o_1], 
                              name_2))
    res_table_II[3, count] <- paste0("(", k_comp_1_name, ")")
    res_table_II[8, count] <- operation_sign[o_2]
    res_table_II[9, count] <- kernel_names[k_2]
    count <- count + 1
  }
}# kernel search second stage end
res_table_II[is.na(res_table_II)] <- "-"
res_table_II

var_II_1 <- as.numeric(res_table_II[1,])
var_II_2 <- as.numeric(res_table_II[2,])
var_II_3 <- (var_II_1 + var_II_2)/2
var_II_4 <- sort(var_II_3)
first_II <- match(var_II_4, var_II_3)[1]
res_table_II[,first_II]
# > res_table_II[,first_II]
# [1] "0.93"         "2.153"        "(k_se+k_per)" "3.7"          "-"           
# [6] "0.1"          "1.9"          "+"            "k_se"         "3.7"         
# [11] "-"            "1.5"         

k_0 <- which(res_table_I[,first_I][3] == kernel_names)
o_1 <- which(res_table_I[,first_I][6] == operation_sign)
k_1 <- which(res_table_I[,first_I][7] == kernel_names)

o_2 <- which(res_table_II[,first_II][8] == operation_sign)
k_2 <- which(res_table_II[,first_II][9] == kernel_names)
h_param_k_0 <- get(kernel_names[k_0])
h_param_k_1 <- get(kernel_names[k_1])
h_param_k_2 <- get(kernel_names[k_2])

# rerun of grid search on the best kernel to restore variables:
two_h_p <- c(F, F, F) 
two_h_p[1] <- ifelse(k_0 == 2, T, F)
two_h_p[2] <- ifelse(k_1 == 2, T, F)
two_h_p[3] <- ifelse(k_2 == 2, T, F)

tot_nr_h_P <- length(c(h_param_k_0, h_param_k_1,
                       h_param_k_2))
h_param_grid <- get(search_lengths[tot_nr_h_P-1])
l_h_param_grid <- dim(h_param_grid)[1]
all_h_p <- c(h_param_k_0, h_param_k_1,
             h_param_k_2)

rmse_tr <- c()
rmse_te <- c()
for(h in 1:l_h_param_grid){# grid search for-loop
  h_p_count <- 1
  for(i in 1:3){# iterates through the nr. of base kernels
    variable <- get(all_h_p[h_p_count])[i]
    assign(variable, h_param_grid[h, h_p_count])
    h_p_count <- h_p_count + 1
    if(two_h_p[i] == TRUE){# TRUE when the jth base kernel has two hyperparameters
      variable <- get(all_h_p[h_p_count])[i]
      assign(variable, h_param_grid[h, h_p_count])
      h_p_count <- h_p_count + 1
    }
  }
  eval_res <- model_eval()
  rmse_tr[h] <- eval_res[1]
  rmse_te[h] <- eval_res[2]
}# grid search end
nr <- which.min((rmse_tr+rmse_te)/2)

# 
# Function to assign the first-twelfth best models from grid search
# 
assign_values <- function(ordinal_nr, all_h_p){
  two_h_p <- c(F, F, F) 
  two_h_p[1] <- ifelse(k_0 == 2, T, F)
  two_h_p[2] <- ifelse(k_1 == 2, T, F)
  two_h_p[3] <- ifelse(k_2 == 2, T, F)
  h_p_count <- 1
  for(i in 1:3){# iterates through the nr. of base kernels
    variable <- get(all_h_p[h_p_count])[i]
    assign(variable, h_param_grid[ordinal_nr, h_p_count], envir = .GlobalEnv)
    h_p_count <- h_p_count + 1
    if(two_h_p[i] == TRUE){# TRUE when the jth base kernel has two hyperparameters
      variable <- get(all_h_p[h_p_count])[i]
      assign(variable, h_param_grid[ordinal_nr, h_p_count], envir = .GlobalEnv)
      h_p_count <- h_p_count + 1
    }
  }
}

# 
# Collecting the twelve best models from the grid search:
# 
ordinal_numbers <- c("first", "second", "third", "fourth", "fifth",
                     "sixth", "seventh", "eighth", "ninth", "tenth",
                     "eleventh", "twelfth") 
for(i in 1:12){
  score_val_s <- sort((rmse_tr+rmse_te)/2)
  assign(ordinal_numbers[i], match(unique(score_val_s)[i],
                                   (rmse_tr+rmse_te)/2))
} 
ord_n_val <- t(sapply(ordinal_numbers, function(x) get(x)))
ord_n_val
# > ord_n_val
#      first second third fourth fifth sixth seventh eighth ninth tenth eleventh twelfth
# [1,]   278    280   279    530   405   404     153    152   154   155      274     264

# 
# collecting the twelve models with lowest RMSE from the grid search:
# 
ordinal_numbers2 <- c("first2", "second2", "third2", "fourth2", "fifth2",
                     "sixth2", "seventh2", "eighth2", "ninth2", "tenth2",
                     "eleventh2", "twelfth2") 
for(i in 1:12){
  sort_rmse_tr <- sort(rmse_tr)
  assign(ordinal_numbers2[i], match(unique(sort_rmse_tr), rmse_tr)[i])
} 
ord_n_val2 <- t(sapply(ordinal_numbers2, function(x) get(x)))
ord_n_val2
# > ord_n_val2
#      first2 second2 third2 fourth2 fifth2 sixth2 seventh2 eighth2 ninth2 tenth2 eleventh2
# [1,]     26      27    279     278    152    153      154     404    155    405       530
#      twelfth2
# [1,]      177

# 
# Plots used in thesis
# 
par(mfrow=c(1,1), mgp = c(1.3, 0.4, 0), mar = c(2.8, 2.8, 0.2, 0.2))
plot_data()
abline(v = 1997.917, lwd = 2, lty = 2, col = "grey50")
grid_rmse_train <- c()
grid_rmse_test <- c()
for(i in 12:1){
  ordinal_nr <- get(ordinal_numbers[i])
  assign_values(ordinal_nr, all_h_p)
  plot_model(i%%12+1, 2)
  nr <- ord_n_val[i]
# storing the corresponding twelve RMSE train and RMSE test values
  grid_rmse_train[i] <- rmse_tr[nr]
  grid_rmse_test[i] <- rmse_te[nr]
}
# storing values of lowest RMSE train models:
grid_rmse_train2 <- c()
grid_rmse_test2 <- c()
for(i in 1:12){
  ordinal_nr2 <- get(ordinal_numbers2[i])
  # assign_values(ordinal_nr2, all_h_p)
  # plot_model(i%%9+1, 2)
  nr <- ord_n_val2[i]
  grid_rmse_train2[i] <- rmse_tr[nr]
  grid_rmse_test2[i] <- rmse_te[nr]
}

# plot for comparison on model evaluation
par(mfrow=c(1,2), mgp = c(1.3, 0.4, 0), mar = c(2.8, 1.8, 2.6, 0.2))
plot(1:12, grid_rmse_train, col = "black", 
     pch = 1, cex = 1.5, ylim = c(0.9,27), xlab = "Models", ylab = "", 
     main = "Models sorted by the lowest values from taking the mean \n of the RMSE train and RMSE test", font.lab= 2)#2.8
points(1:12, grid_rmse_test, col = "black", 
       pch = 20, cex = 1.6)
plot(1:12, grid_rmse_train2, col = "black", 
     pch = 1, cex = 1.5, ylim = c(0.9,27), xlab = "Models", ylab = "",
     main = "Models sorted by \n the lowest RMSE train values", font.lab= 2)
points(1:12, grid_rmse_test2, col = "black", 
       pch = 20, cex = 1.6)

# 
# Plot of the best model(mean evaluation):
# 
par(mfrow=c(1,1), mgp = c(1.3, 0.4, 0), mar = c(2.8, 2.8, 0.2, 0.2))
assign_values(first, all_h_p)#sets h.param. values of the best model
plot_data()
abline(v = 1997.917, lwd = 2, lty = 2, col = "grey50")
plot_model(6,2)
model_eval()
# > model_eval()
# [1] 0.930 2.153
mean(model_eval())
# > mean(model_eval())
# [1] 1.5415

# 
# Hyperparameter tuning function
# 
search_operator <- c(`+`, `-`)
adj_hyperparameter <- function(s, adj, itr, hyperparameter){
  improved <- 0
  eval_res <- model_eval()
  rmse_tr <- eval_res[1]
  rmse_te <- eval_res[2]
  rmse_0 <- (rmse_tr+rmse_te)/2
  print(paste("The initial score is", rmse_0))
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
    rmse_tr <- eval_res[1]
    rmse_te <- eval_res[2]
    rmse_1 <- (rmse_tr+rmse_te)/2
    if(rmse_1 > rmse_0){
      assign(hyperparameter,
             search_operator[[s%%2+1]](get(hyperparameter),
                                       adj), envir = .GlobalEnv)
      print(paste("the score", rmse_1, "is not an improvement"))
      break
    }
    print("improved")
    improved <- 1
    rmse_0 <- rmse_1
  }
  print(paste("Resulting score is", rmse_0, "with hyperparameter",
              hyperparameter, "=", get(hyperparameter)))
  table <- list(rmse_0, hyperparameter, improved)
  return(table)
}
# h.parameters for tuning.
two_h_p <- c(F, F, F) 
two_h_p[1] <- ifelse(k_0 == 2, T, F)
two_h_p[2] <- ifelse(k_1 == 2, T, F)
two_h_p[3] <- ifelse(k_2 == 2, T, F)
current_h_par <- c()
h_p_count <- 1
for(j in 1:3){
  current_h_par[h_p_count] <- get(all_h_p[h_p_count])[j]
  h_p_count <- h_p_count + 1
  if(two_h_p[j]){
    current_h_par[h_p_count] <- get(all_h_p[h_p_count])[j]
    h_p_count <- h_p_count + 1
  }
}
current_h_par <- c(current_h_par, "sigma_var")
# > current_h_par
# [1] "ell_se_0"  "ell_per_1" "period_1"  "ell_se_2"  "sigma_var"

# 
# While-loop
# 
nr_c_h <- length(current_h_par)
iter <- 1
stop_adj <- rep(F,nr_c_h)
no_success <- 0
while(iter<999){
  add_sub <- 2
  attempt <- 0
  first_attempt <- 0
  for(i in 1:5){
    index <- (iter+nr_c_h-1)%%nr_c_h+1
    continue <- adj_hyperparameter(add_sub, 0.008, 5,
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
# Final results
# 

# final plot:
plot_model(4, 2)

res_table_II[c(1,2),first_II]
# > res_table_II[c(1,2),first_II]
# [1] "0.93"  "2.153"
model_eval()
# > model_eval()
# [1] 0.382 0.618
mean(model_eval())
# > mean(model_eval())
# [1] 0.5

res_table_II[c(4:7,10:11),first_II]
# > res_table_II[c(4:7,10:11),first_II]
# [1] "3.7" "-"   "0.1" "1.9" "3.7" "-" 
new_h_p_values <- t(sapply(current_h_par, function(x) get(x)))
new_h_p_values
# > new_h_p_values
#      ell_se_0 ell_per_1 period_1 ell_se_2 sigma_var
# [1,]    3.692     0.076    1.884    3.708      1.24

