# ______________________________________________________________
# kernel search on Pima Indians diabetes data with GP classifier.
# The search is performed with the base kernels: SE.  
# ______________________________________________________________

# rm(list=ls(all=TRUE))
# install.packages("kernlab")
library(kernlab)
# install.packages("DMwR")
library(DMwR)
# install.packages("caret")
library(caret)
# install.packages("mlbench")
library(mlbench)
data(PimaIndiansDiabetes2)
?PimaIndiansDiabetes2# data set description
data <- PimaIndiansDiabetes2
colSums(is.na(data))# missing values
# KNN imputation with k=5 to fill in NA values:
data[,c(-8,-9)] <- knnImputation(data[,c(-8,-9)], k = 5)
colSums(is.na(data))
# 
nrows <- NROW(data)
set.seed(8017)                          
index <- sample(nrows, 0.8*nrows)   
train_data <- data[index,]
test_data <- data[-index,]

# squared exponential (SE) base kernel functions:
k_se_1 <- function(x, y){
  (sigma_var)^2*exp(-0.5*(1/(ell_se_1)^2)*sum((x - y)^2))
}
class(k_se_1) <- "kernel"
k_se_2 <- function(x, y){
  (sigma_var)^2*exp(-0.5*(1/(ell_se_2)^2)*sum((x - y)^2))
}
class(k_se_2) <- "kernel"
k_se_3 <- function(x, y){
  (sigma_var)^2*exp(-0.5*(1/(ell_se_3)^2)*sum((x - y)^2))
}
class(k_se_3) <- "kernel"
k_se_4 <- function(x, y){
  (sigma_var)^2*exp(-0.5*(1/(ell_se_4)^2)*sum((x - y)^2))
}
class(k_se_4) <- "kernel"
# LIN, not used
k_lin <- function(x, y){
  (sigma_var)^2*(sum((x)*(y)))
}
class(k_lin) <- "kernel"
# default kernel, multiplying by one
k_none <- function(x, y){
  1
}
class(k_none) <- "kernel"

# 
# arrangements of list, variables and functions used in the search
# 
kern_1 <- c(`k_se_1`, `k_lin`, `k_none`)
kern_2 <- c(`k_se_2`, `k_lin`, `k_none`)
kern_3 <- c(`k_se_3`, `k_lin`, `k_none`)
kern_4 <- c(`k_se_4`, `k_lin`, `k_none`)
# variables used in the kernel functions
o <- c()
for(i in 1:4){
  o[i] <- paste0("o_", i)
  assign(paste0("o_", i), 2)
}
k <- c()
for(i in 1:4){
  k[i] <- paste0("k_", i)
  assign(paste0("k_", i), 3)
}
ell_se <- c()
for(i in 1:4){
  ell_se[i] <- paste0("ell_se_", i)
  assign(paste0("ell_se_", i), 1)
}
# search operators
operations <- c(`+`, `*`)
operation_sign <- c("+", "×")
kernel_names <- c("SE", "LIN")
SE <- 1
LIN <- 0

# 
k_comp_1 <- function(x, y){
  u = operations[[o_1]](kern_1[[k_1]](x, y), kern_2[[k_2]](x, y))
  return(u)  
}
class(k_comp_1) <- "kernel"
# main kernel function which is the argument for the gausspr function:
k_comp_2 <- function(x, y){
  u = operations[[o_2]](k_comp_1(x, y), kern_3[[k_3]](x, y))
  return(u)  
}
class(k_comp_2) <- "kernel"
k_comp_3 <- function(x, y){
  u = operations[[o_3]](k_comp_2(x, y), kern_4[[k_4]](x, y))
  return(u)  
}
class(k_comp_3) <- "kernel"
# 

# model evaluation function used in the k.search 
metric_eval <- function(){
  table <- c()
  model <- gausspr(diabetes ~ ., data = train_data, kernel = k_comp_3)
  res_1 <- confusionMatrix(predict(model, train_data[,-9]),
                           train_data$diabetes, positive = "pos")$byClass[1]# Sens. train
  res_2 <- confusionMatrix(predict(model, test_data[,-9]),
                           test_data$diabetes, positive = "pos")$byClass[1]# Sens. test
  res_3 <- confusionMatrix(predict(model, test_data[,-9]),
                           test_data$diabetes, positive = "pos")$overall[1]# Accuracy
  table <- c(round(res_1, 3), round(res_2, 3), round(res_3, 3))
  return(table)
}
# 
# Kernel search stage I
# 
res_table_I <- data.frame(matrix(nrow = 8, ncol = 2))
row.names(res_table_I) <- c("train sens.", "test sens.", "test acc.", "k_1", 
                            "ell_1", "o_1", "k_2", "ell_2")
interval <- c(0.1, 0.5, 1, 1.5, 2, 2.5)
interval <- c(0.1, 0.5, 1, 1.5)
# interval <- c(0.1, 1)# interval for test run
names(res_table_I) <- 1:2
# default variable
sigma_var <- 1

one_h_param <- expand.grid(interval)
two_h_param <- expand.grid(interval, interval)
three_h_param <- expand.grid(interval, interval, interval)
four_h_param <- expand.grid(interval, interval,
                            interval, interval)
search_lengths <- c("one_h_param", "two_h_param", 
                    "three_h_param", "four_h_param")
score_value <- c()
count <- 1
for(h in 1:2){# iterates through search operations: sums and products
  o_1 <- h
  k_1 <- 1
  k_2 <- 1
  res_table_I[6, count] <- operation_sign[o_1]
  res_table_I[4, count] <- kernel_names[k_1]
  res_table_I[7, count] <- kernel_names[k_2]
  nr_h_p <- 2
  r1 <- c()
  r2 <- c()
  r3 <- c()
  if(nr_h_p > 0){
    h_param_grid <- get(search_lengths[nr_h_p])
    h_p_length <- dim(h_param_grid)[1]
    for(j in 1:h_p_length){# grid search for-loop
      assign("ell_se_1", h_param_grid[j, 1], envir = .GlobalEnv)
      assign("ell_se_2", h_param_grid[j, 2], envir = .GlobalEnv)
      eval_res <- metric_eval()
      r1[j] <- eval_res[1]
      r2[j] <- eval_res[2]
      r3[j] <- eval_res[3]
    }# grid search end
  }
  nr <- which.max(r1+r2+r3)# choosing the best model
  score_value[count] <- r1[nr]+r2[nr]+r3[nr]
  assign("ell_se_1", h_param_grid[nr, 1], envir = .GlobalEnv)
  res_table_I[5, count] <- ell_se_1
  assign("ell_se_2", h_param_grid[nr, 2], envir = .GlobalEnv)
  res_table_I[8, count] <- ell_se_2
  res_table_I[1, count] <- r1[nr]
  res_table_I[2, count] <- r2[nr]
  res_table_I[3, count] <- r3[nr]
  count <- count + 1
}
res_table_I[is.na(res_table_I)] <- "-"
res_table_I

choice_I <- which.max(score_value)
res_table_I[,choice_I]
# > res_table_I[,choice_I]
# [1] "0.886" "0.776" "0.838" "SE"    "1"     "+"     "SE"    "0.1"  

o_1 <- which(res_table_I[,choice_I][6] == operation_sign)
# 
# Kernel search stage II
# 
res_table_II <- data.frame(matrix(nrow = 9, ncol = 2))
row.names(res_table_II) <- c("train sens.", "test sens.", "test acc.", "kernel", 
                             "ell_1", "ell_2", "o_2", "k_3", "ell_3")
names(res_table_II) <- 1:2
k_comp_1_name <- c(paste0("(",kernel_names[k_1], operation_sign[o_1], 
                          kernel_names[k_2], ")"))
res_table_II[4,] <- k_comp_1_name

score_value_II <- c()
count <- 1
for(h in 1:2){# iterates through search operations: sums and products
  o_2 <- h
  res_table_II[7, count] <- operation_sign[o_2]
  k_3 <- 1
  res_table_II[8, count] <- kernel_names[k_3]
  nr_h_p <- 3
  r1 <- c()
  r2 <- c()
  r3 <- c()
  if(nr_h_p > 0){
    h_param_grid <- get(search_lengths[nr_h_p])
    h_p_length <- dim(h_param_grid)[1]
    for(j in 1:h_p_length){# grid search for-loop
      assign("ell_se_1", h_param_grid[j, 1], envir = .GlobalEnv)
      assign("ell_se_2", h_param_grid[j, 2], envir = .GlobalEnv)
      assign("ell_se_3", h_param_grid[j, 3], envir = .GlobalEnv)
      eval_res <- metric_eval()
      r1[j] <- eval_res[1]
      r2[j] <- eval_res[2]
      r3[j] <- eval_res[3]
    }
  }
  nr <- which.max(r1+r2+r3)
  score_value_II[count] <- r1[nr]+r2[nr]+r3[nr]
  col_index <- 1
  assign("ell_se_1", h_param_grid[nr, col_index], envir = .GlobalEnv)
  res_table_II[5, count] <- ell_se_1
  col_index <- col_index + 1
  assign("ell_se_2", h_param_grid[nr, col_index], envir = .GlobalEnv)
  res_table_II[6, count] <- ell_se_2
  col_index <- col_index + 1
  assign("ell_se_3", h_param_grid[nr, col_index], envir = .GlobalEnv)
  res_table_II[9, count] <- ell_se_3
  res_table_II[1, count] <- r1[nr]
  res_table_II[2, count] <- r2[nr]
  res_table_II[3, count] <- r3[nr]
  count <- count + 1
}
res_table_II[is.na(res_table_II)] <- "-"
res_table_II

choice_II <- which.max(score_value_II)
res_table_II[,choice_II]
# > res_table_II[,choice_II]
# [1] "0.963"   "0.776"   "0.838"   "(SE+SE)" "1"       "0.5"     "+"       "SE"     
# [9] "0.1"  

o_2 <- which(res_table_II[,choice_II][7] == operation_sign)
# 
# Kernel search stage III
# 
res_table_III <- data.frame(matrix(nrow = 10, ncol = 2))
row.names(res_table_III) <- c("train sens.", "test sens.", "test acc.", "kernel", 
                              "ell_1", "ell_2", "ell_3", "o_3", "k_4", "ell_4")
names(res_table_III) <- 1:2
k_comp_2_name <- c(paste0("(", k_comp_1_name, operation_sign[o_2],
                          kernel_names[k_3],")"))
res_table_III[4,] <- k_comp_2_name

score_value_III <- c()
count <- 1
for(h in 1:2){
  o_3 <- h
  res_table_III[8, count] <- operation_sign[o_3]
  k_4 <- 1
  res_table_III[9, count] <- kernel_names[k_4]
  nr_h_p <- 4
  r1 <- c()
  r2 <- c()
  r3 <- c()
  h_param_grid <- get(search_lengths[nr_h_p])
  h_p_length <- dim(h_param_grid)[1]
  for(j in 1:h_p_length){
    assign("ell_se_1", h_param_grid[j, 1], envir = .GlobalEnv)
    assign("ell_se_2", h_param_grid[j, 2], envir = .GlobalEnv)
    assign("ell_se_3", h_param_grid[j, 3], envir = .GlobalEnv)
    assign("ell_se_4", h_param_grid[j, 4], envir = .GlobalEnv)
    eval_res <- metric_eval()
    r1[j] <- eval_res[1]
    r2[j] <- eval_res[2]
    r3[j] <- eval_res[3]
  }
  
  nr <- which.max(r1+r2+r3)
  score_value_III[count] <- r1[nr]+r2[nr]+r3[nr]
  col_index <- 1
  assign("ell_se_1", h_param_grid[nr, col_index], envir = .GlobalEnv)
  res_table_III[5, count] <- ell_se_1
  col_index <- col_index + 1
  assign("ell_se_2", h_param_grid[nr, col_index], envir = .GlobalEnv)
  res_table_III[6, count] <- ell_se_2
  col_index <- col_index + 1
  assign("ell_se_3", h_param_grid[nr, col_index], envir = .GlobalEnv)
  res_table_III[7, count] <- ell_se_3
  col_index <- col_index + 1
  assign("ell_se_4", h_param_grid[nr, col_index], envir = .GlobalEnv)
  res_table_III[10, count] <- ell_se_3
  col_index <- col_index + 1
  res_table_III[1, count] <- r1[nr]
  res_table_III[2, count] <- r2[nr]
  res_table_III[3, count] <- r3[nr]
  count <- count + 1
}
res_table_III[is.na(res_table_III)] <- "-"
res_table_III

choice_III <- which.max(score_value_III)
res_table_III[,choice_III]
# > res_table_III[,choice_III]
# [1] "0.995"        "0.776"        "0.844"        "((SE+SE)+SE)" "1"           
# [6] "0.5"          "0.1"          "+"            "SE"           "0.1"

# best model from the search at stage III:
o_3 <- which(res_table_III[,choice_III][8] == operation_sign)
assign("ell_se_1", as.numeric(res_table_III[c(5),choice_III]))
assign("ell_se_2", as.numeric(res_table_III[c(6),choice_III]))
assign("ell_se_3", as.numeric(res_table_III[c(7),choice_III]))
assign("ell_se_4", as.numeric(res_table_III[c(10),choice_III]))


metric_eval()
# > metric_eval()
# Sensitivity Sensitivity    Accuracy 
# 0.995       0.776          0.844 
sum(metric_eval())
# > sum(metric_eval())
# [1] 2.615

# 
# Confusion matrix
# 
model <- gausspr(diabetes ~ ., data = train_data, kernel = k_comp_3)
conf_tr <- confusionMatrix(predict(model, train_data[,-9]),
                           train_data$diabetes, positive = "pos")
conf_te <- confusionMatrix(predict(model, test_data[,-9]),
                           test_data$diabetes, positive = "pos")
conf_tr$table
# > conf_tr$table
#           Reference
# Prediction neg pos
# neg        395   1
# pos        0   218
conf_te$table
# > conf_te$table
#           Reference
# Prediction neg pos
# neg         92  11
# pos         13  38

# 
# Comparing with RBF kernel
# 
model_RBF <- gausspr(diabetes ~ ., data = data)
conf_tr_RBF <- confusionMatrix(predict(model_RBF, train_data[,-9]),
                              train_data$diabetes, positive = "pos")
conf_te_RBF <- confusionMatrix(predict(model_RBF, test_data[,-9]),
                              test_data$diabetes, positive = "pos")
conf_tr_RBF$table
# > conf_tr_RBF$table
# Reference
# Prediction neg pos
# neg 348  83
# pos  47 136

conf_te_RBF$table
# > conf_te_RBF$table
# Reference
# Prediction neg pos
# neg  93  11
# pos  12  38

