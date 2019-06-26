# __________________________________________________________
# Kernel search on Heart disease data with SVM.
# The search is performed with the base kernels: SE and LIN.  
# __________________________________________________________

# rm(list=ls(all=TRUE))
# install.packages("kernlab")
library(kernlab)
# install.packages("caret")
library(caret)
# install.packages("ROCR")
library(ROCR)

# ***
# Following data preparation is from Brigitte Mueller: 
# https://rpubs.com/mbbrigitte/heartdisease
# ***
# install.packages("data.table")
library(data.table)
heart.data <- read.csv(paste0("https://archive.ics.uci.edu/ml/machine-learning-databases/",
                              "heart-disease/processed.cleveland.data"), header=FALSE, 
                       sep=",",na.strings = '?')
names(heart.data) <- c("age", "sex", "cp", "trestbps", "chol", "fbs", "restecg", 
                       "thalach","exang", "oldpeak","slope", "ca", "thal", "num")
heart.data$num[heart.data$num > 0] <- 1
chclass <- c("numeric","factor","factor","numeric","numeric",
             "factor","factor","numeric","factor","numeric",
             "factor","factor","factor","factor")
convert.magic <- function(obj,types){
  for (i in 1:length(obj)){
    FUN <- switch(types[i],character = as.character, 
                  numeric = as.numeric, 
                  factor = as.factor)
    obj[,i] <- FUN(obj[,i])
  }
  obj
}
heart.data <- convert.magic(heart.data,chclass)
heart = heart.data #add labels only for plot
# missing data:
s = sum(is.na(heart.data))
dim(heart.data)
# [1] 303  14
heart.data <- na.omit(heart.data)
dim(heart.data)
# [1] 297  14

# training and test set:
set.seed(10)
inTrainRows <- createDataPartition(heart.data$num, p=0.7, list=FALSE)
train_data <- heart.data[inTrainRows,]
test_data <-  heart.data[-inTrainRows,]

# 
# base kernels functions
# 

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
# LIN kernel without h.parameter
k_lin_0 <- function(x, y){
  (sigma_var)^2*(sum((x)*(y)))
}
class(k_lin_0) <- "kernel"
# default kernel, multiplying by one
k_none <- function(x, y){
  1
}
class(k_none) <- "kernel"
# 

# 
# arrangements of list, variables and functions used in the search
# 
kern_1 <- c(`k_se_1`, `k_lin_0`, `k_none`)
kern_2 <- c(`k_se_2`, `k_lin_0`, `k_none`)
kern_3 <- c(`k_se_3`, `k_lin_0`, `k_none`)
kern_4 <- c(`k_se_4`, `k_lin_0`, `k_none`)
# 
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
all_h_p <- get("ell_se")
kernel_names <- c("k_se", "k_lin")
operations <- c(`+`, `*`)
operation_sign <- c("+", "×")
k_se <- c("ell_se")
k_lin <- 0

# infix functions. Sums of products grammar 
`%op_I%` <- function(a, b){
  operations[[o_1]](a, b)
}
`%op_II%` <- function(a, b){
  operations[[o_2]](a, b)
}
`%op_III%` <- function(a, b){
  operations[[o_3]](a, b)
}

# main kernel function which is the argument for the ksvm function:
k_comp <- function(x, y){
  res <- kern_1[[k_1]](x, y)%op_I%kern_2[[k_2]](x, y)%op_II%
    kern_3[[k_3]](x, y)%op_III%kern_4[[k_4]](x, y)
  return(res)
}
class(k_comp) <- "kernel"
# 
sigma_var <- 1
cost <- 1
# model evaluation function used in the k.search
metric_eval <- function(){
  table <- c()
  model <- ksvm(num ~ ., data = train_data, C = cost, kernel = k_comp)
  res_1 <- confusionMatrix(predict(model, train_data[,-14]),
                           train_data$num, positive = "1")$byClass[1]# Sens. train
  res_2 <- confusionMatrix(predict(model, test_data[,-14]),
                           test_data$num, positive = "1")$byClass[1]# Sens. test
  res_3 <- confusionMatrix(predict(model, test_data[,-14]),
                           test_data$num, positive = "1")$overall[1]# Accuracy
  table <- c(round(res_1, 3), round(res_2, 3), round(res_3, 3))
  return(table)
}

# 
# kernel search stage I
# 
res_table_I <- data.frame(matrix(nrow = 9, ncol = 6))
row.names(res_table_I) <- c("train sens.", "test sens.", "test acc.", "k_1", 
                            "ell_1", "o_1", "k_2", "ell_2", "C")
res_table_I[9,] <- cost
interval <- c(0.1, 0.5, 1, 1.5, 2, 2.5, 3)
# interval <- c(0.1, 3)# interval for test run
k_comb <- expand.grid(c(1,2),c(1,2))
k_comb <- k_comb[-2,]
names(k_comb) <- c("k_1", "k_2")
names(res_table_I) <- 1:6

one_h_param <- expand.grid(interval)
two_h_param <- expand.grid(interval, interval)
three_h_param <- expand.grid(interval, interval, interval)
four_h_param <- expand.grid(interval, interval,
                            interval, interval)
search_lengths <- c("one_h_param", "two_h_param", 
                    "three_h_param", "four_h_param")
score_value <- c()
# k.search for-loop
count <- 1
for(h in 1:2){# iterates through search operations: sums and products
  o_1 <- h
  for(i in 1:dim(k_comb)[1]){# iterates through table with kernel combinations
    k_1 <- k_comb[i, 1]
    k_2 <- k_comb[i, 2]
    
    one_h_p <- c(F, F)# variable to identify 2-h.param. kernel 
    one_h_p[1] <- ifelse(k_1 == 1, T, F)
    one_h_p[2] <- ifelse(k_2 == 1, T, F)
    
    tot_nr_h_p <- sum(one_h_p)#number of TRUE in vector
    if(tot_nr_h_p > 0){
      h_param_grid <- get(search_lengths[tot_nr_h_p])
      l_h_param_grid <- dim(h_param_grid)[1]
    }else{
      l_h_param_grid <- 1 #no grid search when no h.param.
    }
    r1 <- c()
    r2 <- c()
    r3 <- c()
    for(i in 1:l_h_param_grid){# grid search for-loop
      h_p_count <- 1
      for(j in 1:2){
        if(one_h_p[j] == TRUE){
          assign(all_h_p[j], h_param_grid[i, h_p_count])
          h_p_count <- h_p_count + 1
        }
      }
      eval_res <- metric_eval()
      r1[i] <- eval_res[1]
      r2[i] <- eval_res[2]
      r3[i] <- eval_res[3]
      
    }# grid search end
    nr <- which.max(r1+r2+r3)# choosing the best model
    score_value[count] <- r1[nr]+r2[nr]+r3[nr]
    # grid search results inserted into table:
    res_table_I[1, count] <- r1[nr]
    res_table_I[2, count] <- r2[nr]
    res_table_I[3, count] <- r3[nr]
    col_index <- 1
    if(one_h_p[1]){
      assign(all_h_p[1], h_param_grid[nr, col_index])
      res_table_I[5, count] <-  get(all_h_p[1])
      col_index <- col_index + 1
    }
    if(one_h_p[2]){
      assign(all_h_p[2], h_param_grid[nr, col_index])
      res_table_I[8, count] <-  get(all_h_p[2])
    }
    
    name_1 <- kernel_names[k_1]
    name_2 <- kernel_names[k_2]
    res_table_I[c(4, 7), count] <- c(name_1, name_2)
    res_table_I[6, count] <- operation_sign[o_1]
    count <- count + 1# next column of the table
  }
}# kernel search first stage end
res_table_I[is.na(res_table_I)] <- "-"
res_table_I

choice_I <- which.max(score_value)
res_table_I[,choice_I]

# 
# kernel search stage II
# 
res_table_II <- data.frame(matrix(nrow = 10, ncol = 4))
row.names(res_table_II) <- c("train sens.", "test sens.", "test acc.", "kernel", 
                             "ell_1", "ell_2", "o_2", "k_3", "ell_3", "C")
res_table_II[10,] <- cost
names(res_table_II) <- 1:4

score_value_II <- c()
count <- 1
for(h in 1:2){
  o_1 <- which(res_table_I[,choice_I][6] == operation_sign)
  o_2 <- h
  for(i in 1:2){
    k_1 <- which(res_table_I[,choice_I][4] == kernel_names)
    k_2 <- which(res_table_I[,choice_I][7] == kernel_names)
    k_3 <- i
    
    one_h_p <- c(F, F, F)
    one_h_p[1] <- ifelse(k_1 == 1, T, F)
    one_h_p[2] <- ifelse(k_2 == 1, T, F)
    one_h_p[3] <- ifelse(k_3 == 1, T, F)
    
    tot_nr_h_p <- sum(one_h_p)#number of TRUE in vector
    if(tot_nr_h_p > 0){
      h_param_grid <- get(search_lengths[tot_nr_h_p])
      l_h_param_grid <- dim(h_param_grid)[1]
    }else{
      l_h_param_grid <- 1
    }
    
    r1 <- c()
    r2 <- c()
    r3 <- c()
    for(i in 1:l_h_param_grid){
      h_p_count <- 1
      for(j in 1:3){
        if(one_h_p[j] == TRUE){
          assign(all_h_p[j], h_param_grid[i, h_p_count])
          h_p_count <- h_p_count + 1
        }
      }
      eval_res <- metric_eval()
      r1[i] <- eval_res[1]
      r2[i] <- eval_res[2]
      r3[i] <- eval_res[3]
    }
    nr <- which.max(r1+r2+r3)
    score_value_II[count] <- r1[nr]+r2[nr]+r3[nr]
    
    res_table_II[1, count] <- r1[nr]
    res_table_II[2, count] <- r2[nr]
    res_table_II[3, count] <- r3[nr]
    
    col_index <- 1
    if(one_h_p[1]){
      assign(all_h_p[1], h_param_grid[nr, col_index])
      res_table_II[5, count] <-  get(all_h_p[1])
      col_index <- col_index + 1
    }
    if(one_h_p[2]){
      assign(all_h_p[2], h_param_grid[nr, col_index])
      res_table_II[6, count] <-  get(all_h_p[2])
      col_index <- col_index + 1
    }
    if(one_h_p[3]){
      assign(all_h_p[3], h_param_grid[nr, col_index])
      res_table_II[9, count] <-  get(all_h_p[3])
    }
    
    name_1 <- kernel_names[k_1]
    name_2 <- kernel_names[k_2]
    k_comp_1_name <- c(paste0(kernel_names[k_1], operation_sign[o_1], 
                              kernel_names[k_2]))
    res_table_II[4,] <- k_comp_1_name
    res_table_II[7, count] <- operation_sign[o_2]
    res_table_II[8, count] <- kernel_names[k_3]
    count <- count + 1
  }
}
res_table_II[is.na(res_table_II)] <- "-"
res_table_II

choice_II <- which.max(score_value_II)
res_table_II[,choice_II]

# 
# kernel search stage III
# 
res_table_III <- data.frame(matrix(nrow = 10, ncol = 4))
row.names(res_table_III) <- c("train sens.", "test sens.", "test acc.", "kernel", 
                              "ell_1", "ell_2", "ell_3", "o_3", "k_4", "ell_4")
names(res_table_III) <- 1:4

score_value_III <- c()
count <- 1
for(h in 1:2){
  o_2 <- which(res_table_II[,choice_II][7] == operation_sign)
  o_3 <- h
  for(i in 1:2){
    k_3 <- which(res_table_II[,choice_II][8] == kernel_names)
    k_4 <- i
    
    one_h_p <- c(F, F, F, F)
    one_h_p[1] <- ifelse(k_1 == 1, T, F)
    one_h_p[2] <- ifelse(k_2 == 1, T, F)
    one_h_p[3] <- ifelse(k_3 == 1, T, F)
    one_h_p[4] <- ifelse(k_4 == 1, T, F)
    
    tot_nr_h_p <- sum(one_h_p)#number of TRUE in vector
    if(tot_nr_h_p > 0){
      h_param_grid <- get(search_lengths[tot_nr_h_p])
      l_h_param_grid <- dim(h_param_grid)[1]
    }else{
      l_h_param_grid <- 1
    }
    r1 <- c()
    r2 <- c()
    r3 <- c()
    for(i in 1:l_h_param_grid){
      h_p_count <- 1
      for(j in 1:4){
        if(one_h_p[j] == TRUE){
          assign(all_h_p[j], h_param_grid[i, h_p_count])
          h_p_count <- h_p_count + 1
        }
      }
      eval_res <- metric_eval()
      r1[i] <- eval_res[1]
      r2[i] <- eval_res[2]
      r3[i] <- eval_res[3]
    }
    nr <- which.max(r1+r2+r3)
    score_value_III[count] <- r1[nr]+r2[nr]+r3[nr]
    res_table_III[1, count] <- r1[nr]
    res_table_III[2, count] <- r2[nr]
    res_table_III[3, count] <- r3[nr]
    
    col_index <- 1
    if(one_h_p[1]){
      assign(all_h_p[1], h_param_grid[nr, col_index])
      res_table_III[5, count] <-  get(all_h_p[1])
      col_index <- col_index + 1
    }
    if(one_h_p[2]){
      assign(all_h_p[2], h_param_grid[nr, col_index])
      res_table_III[6, count] <-  get(all_h_p[2])
      col_index <- col_index + 1
    }
    if(one_h_p[3]){
      assign(all_h_p[3], h_param_grid[nr, col_index])
      res_table_III[7, count] <-  get(all_h_p[3])
      col_index <- col_index + 1
    }
    if(one_h_p[4]){
      assign(all_h_p[4], h_param_grid[nr, col_index])
      res_table_III[10, count] <-  get(all_h_p[4])
    }
    
    k_comp_2_name <- c(paste0(kernel_names[k_1], operation_sign[o_1], 
                              kernel_names[k_2], operation_sign[o_2],
                              kernel_names[k_3]))
    res_table_III[4,] <- k_comp_2_name
    res_table_III[8, count] <- operation_sign[o_3]
    res_table_III[9,count] <- kernel_names[k_4]
    count <- count + 1
  }
}
res_table_III[is.na(res_table_III)] <- "-"
res_table_III

choice_III <- which.max(score_value_III)
res_table_III[,choice_III]

o_3 <- which(res_table_III[,choice_III][8] == operation_sign)
k_4 <- which(res_table_III[,choice_III][9] == kernel_names)

# 
# Kernel search result from table III:
# 
k_1 <- 1; k_2 <- 2; k_3 <- 1; k_4 <- 1
ell_se_1 <- 0.1; ell_se_3 <- 2; ell_se_4 <- 0.5
o_1 <- 1; o_2 <- 2; o_3 <- 1
metric_eval()
# > metric_eval()
# Sensitivity Sensitivity    Accuracy 
# 0.958       0.902       0.820 


# Plot(thesis):
par(mfrow=c(1,2), mgp = c(1.3, 0.4, 0), mar = c(2.8, 2.8, 3.2, 0.2))
# 
# ROC TRAIN DATA:
#

# RBF kernel:
set.seed(392)
model0 <- ksvm(num ~ ., data = heart.data, C = cost)
ypredscore = predict(model0, train_data, type= "decision")
pred0 <- prediction(ypredscore, train_data[,14])
perf0 <- performance(pred0, measure = "tpr", x.measure = "fpr")
plot(perf0, col='steelblue', main = "Train data set", lwd = 2)
auc <- performance(pred0,measure="auc")@y.values
auc
# > auc
# [[1]]
# [1] 0.9617746

# constructed kernel:
model <- ksvm(num ~ ., data = train_data, C = cost, kernel = k_comp)
ypredscore = predict(model, train_data, type= "decision")
pred <- prediction(ypredscore, train_data[,14])
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf, col = 'red', add = 1, lwd = 2)
abline(a = 0, b = 1, lty = 2)
auc <- performance(pred,measure="auc")@y.values
auc
# > auc
# [[1]]
# [1] 0.9861421
legend("bottomright", c("kernel search", "kernlab RBF"), lty=1, lwd = 2,
       col = c("red", "steelblue"), bty="n", cex = 0.8, inset=c(0.05,0.05))

# 
# ROC TEST DATA:
#

# RBF kernel:
set.seed(392)
model <- ksvm(num ~ ., C = cost, 
              data = heart.data)
ypredscore = predict(model, test_data, type= "decision")
pred <- prediction(ypredscore, test_data[,14])
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf, col='steelblue', main = "Test data set", lwd = 2)
legend("bottomright", c("kernel search", "kernlab RBF"), lty=1, 
       col = c("red", "steelblue"), bty="n", cex = 0.8,
       inset=c(0.05,0.05), lwd = 2)
auc <- performance(pred,measure="auc")@y.values
auc
# > auc
# [[1]]
# [1] 0.9563008

# constructed kernel:
model <- ksvm(num ~ ., data = train_data,
              C = cost, kernel = k_comp)
ypredscore = predict(model, test_data, type= "decision")
pred <- prediction(ypredscore, test_data[,14])
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf, col='red', add = 1, lwd = 2)
abline(a=0,b=1, lty = 2)
auc <- performance(pred,measure="auc")@y.values
auc
# > auc
# [[1]]
# [1] 0.8922764

