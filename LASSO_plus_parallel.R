library(ggpubr)
library(corrplot)
library(glmnet)
library(caret)
library(CBCgrps)
library(tidyverse)
library(rms)
library(pROC)
library(readxl)
library(tableone)
library(this.path)
library(parallel)
library(foreach)
library(doParallel)
library(writexl)
library(doSNOW)
library(readr)
library(dplyr)
library(ggplot2)
currdir = dirname(this.path())
setwd(currdir)

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 确保提供了参数
if (length(args) != 0) {
  train_num <- as.numeric(args[1])
}else{
  train_num <- 1
}
print('##################################')
print(train_num)
print('##################################')


############### function for data ###############
data_extract=function(dataX0,dataY0,randseed,sbjnum_ratio) {
  set.seed(randseed)
  
  sbjtype1=which(dataY0==1)
  sbjtype0=which(dataY0==0)
  
  sbjmin=min(c(length(sbjtype1),length(sbjtype0)))            # 取类型数量相对较少的一组的sample数
  
  sbjtype1sel=sample(sbjtype1,round(sbjmin*sbjnum_ratio))     # 取部分sample作为训练集
  sbjtype0sel=sample(sbjtype0,round(sbjmin*sbjnum_ratio))     # 取部分sample作为训练集
  # 剩余样本，为training样本。剩余的5人，为test样本（在特征筛选阶段，test样本不起作用）
  
  sbjsel=c(sbjtype0sel,sbjtype1sel)
  
  dataXtrain=dataX0[sbjsel,]
  dataYtrain=dataY0[sbjsel]
  
  dataYtest_pool=dataY0[-sbjsel]
  dataXtest_pool=dataX0[-sbjsel,]
  
  sbjtype1test=which(dataYtest_pool==1)
  sbjtype0test=which(dataYtest_pool==0)
  sbjmin_test=min(c(length(sbjtype1test),length(sbjtype0test)))
  
  sbjtype1sel_test=sample(sbjtype1test,sbjmin_test)     # 取部分sample作为训练集
  sbjtype0sel_test=sample(sbjtype0test,sbjmin_test)     # 取部分sample作为训练集
  
  dataXtest=dataXtest_pool[c(sbjtype1sel_test,sbjtype0sel_test),]
  dataYtest=dataYtest_pool[c(sbjtype1sel_test,sbjtype0sel_test)]
  
  return(list(dataXtrain,dataYtrain,dataXtest,dataYtest))
}


# 提取特征表
get_feature_table <- function(dataX, sorted_indices) {
  feature_table <- dataX[, sorted_indices]
  # feature_table[] <- lapply(feature_table, as.numeric)
  return(feature_table)
}

############# sub function for lasso #############
trylasso=function(dataX0, dataY0, randseed, hyper) {

  set.seed(randseed)
  sbjtype1=which(dataY0==1)
  sbjtype0=which(dataY0==0)
  sbjtype1sel=sample(sbjtype1,length(sbjtype1)-5)     # 治愈组随机减去5人
  sbjtype0sel=sample(sbjtype0,length(sbjtype0)-5)     # 非治愈组随机减去5人
  # 剩余样本，为training样本。剩余的5人，为test样本（在特征筛选阶段，test样本不起作用）
  sbjsel=c(sbjtype0sel,sbjtype1sel);
  dataX=dataX0[sbjsel,]
  dataY=dataY0[sbjsel]
  dataXtest=dataX0[-sbjsel,]
  dataYtest=dataY0[-sbjsel]

  fit<- glmnet(dataX,dataY,family = "binomial",lambda = hyper,alpha = 1)    # training model
  roilog=which(coef(fit)!=0)-1         # index of selected feature (-1 means index of intercept is removed)

###临时，占位用的
  accmax = 0
  accmax_test = 0

  res=list(roilog,accmax,accmax_test)
  return(res)

}


#逻辑
logistic_regression = function( dataX0, dataY0, randseed, feature_table) {
  test_accuracy = 0
  # 使用 data_extract 函数进行数据处理
  data_splits = data_extract(feature_table, dataY0, randseed, 0.85)
  dataXtrain = data_splits[[1]]
  dataYtrain = data_splits[[2]]
  dataXtest = data_splits[[3]]
  dataYtest = data_splits[[4]]
  # 将提取的特征数据转换为数据框
  feature_table_scaled_train <- as.data.frame(dataXtrain)
  feature_table_scaled_test <- as.data.frame(dataXtest)

  # 逻辑回归模型训练
  logistic_model <- glm(dataYtrain ~ ., data = feature_table_scaled_train, family = binomial)
  # 模型预测
  test_predictions <- predict(logistic_model, newdata = feature_table_scaled_test, type = "response")
  test_pred_classes <- ifelse(test_predictions > 0.5, 1, 0)
  # 计算测试集上的准确率
  test_accuracy <- sum(test_pred_classes == dataYtest) / length(dataYtest)
  return(list(model = logistic_model, accuracy = test_accuracy))
}




data <- read.csv("musk.csv")

roimaxpred = 50   # 最终预测时的最大特征数
randseedbase = 100
sample_train = 100          # size for training set
sbjnum_ratio = 0.85
sample_test = 200           # size for testing set

x <- as.matrix(data[, -ncol(data)])  # 排除目标变量列
y <- as.factor(data[, ncol(data)])  # 目标变量

set.seed(1234)

# extract data training pool and testing pool
set.seed(randseedbase + train_num * 10)
while (TRUE) {
  ind_train = sample(1:round(nrow(x) / 3), sample_train)
  if (sum(y[ind_train] == 1) >= 0.3 * sample_train & sum(y[ind_train] == 1) <= 0.7 * sample_train) {
    break
  }
}
print(y)
x_train = x[ind_train, ]
y_train = y[ind_train]
# 测试集为原始数据减去训练集的部分
ind_test = setdiff(1:nrow(x), ind_train)
x_test = x[ind_test, ]
y_test = y[ind_test]


print('##################################')
print("Training data:")
print(x_train)
print(y_train)
print("Testing data:")
print(x_test)
print(y_test)


# 将训练数据合并成一个数据框
train_data <- data.frame(x_train)
train_data$Label <- y_train

# 动态生成文件名
train_data_name <- paste0("train_data_",train_num, ".csv")

# 保存训练数据到CSV文件
write.csv(train_data, train_data_name, row.names = FALSE)



dataX0=x_train
dataY0=y_train


grid <-  10^seq(2, -4, length = 100)
randseed <- 1235673
numCores <- detectCores()
print(numCores)
cl <- makeCluster(numCores)
registerDoSNOW(cl)
pb <- txtProgressBar(min = 0, max = length(grid), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

foreach_result <- foreach(fixed_lambda = grid, .combine = rbind, .packages = c('glmnet'),.inorder=TRUE,.options.snow=opts) %dopar% {
  best_feature_table = c()
  best_accuracy = 0.1
  local_results <- data.frame()
  m_accuracy = 0
  # 初始化存储特征选择次数的向量
  num_features <- ncol(dataX0)
  feature_counts <- integer(num_features)
  # 循环1000次
  for (i in 1:1000) {
    res <- trylasso(dataX0, dataY0, randseed + i, fixed_lambda)
    roilog <- res[[1]]
    # 更新特征选择次数
    if (length(roilog) > 0) {
      feature_counts[roilog] <- feature_counts[roilog] + 1
    }
  }
  # 对特征选择次数进行降序排序
  sorted_indices <- order(feature_counts, decreasing = TRUE)
  # 去除出现次数为0的特征索引
  sorted_indices <- sorted_indices[feature_counts[sorted_indices] > 0]
  #排除特征数量小于2的，不晓得为什么小于2会出问题，2似乎没问题
  if (length(sorted_indices) >= 2){
  #进行切片
    ll = length(sorted_indices)
    if (length(sorted_indices) > 50){
      ll = 50
    }
    for (i in 2:ll) {
      #进行切片
      spilted_indices = head(sorted_indices, i)
      # 获取特征表
      feature_table <- get_feature_table(dataX0, spilted_indices)
      if (length(feature_table)!=0){
      # #逻辑回归
        a_accuracy = 0
        nn = 0
        for (j in 1:1000) {
          result = logistic_regression(dataX0, dataY0,randseed+j,feature_table)
          if ( length(result$accuracy) > 0 && !is.na(result$accuracy)){
            nn = nn + 1
            a_accuracy = result$accuracy + a_accuracy
            m_accuracy = a_accuracy/nn
          }
        }
        min_selected_count <- min(feature_counts[spilted_indices])
      }
      local_results <- rbind(local_results, data.frame(FixedLambda = fixed_lambda, Indices = paste(spilted_indices, collapse = ","), Accuracy = m_accuracy,freq = min_selected_count))
    }
  }
  return(local_results)
}
print("结束啦")
results <- foreach_result
# 保存所有特征索引及其准确率

# 动态生成文件名
feature_selection_results_name <- paste0("feature_selection_results_",train_num, ".csv")

write.csv(results, feature_selection_results_name, row.names = FALSE)
stopCluster(cl)


# 读取结果数据
results <- read.csv(feature_selection_results_name)

# 过滤出频率大于500的结果
filtered_results <- results %>% filter(freq > 500)

# 找到最高的准确率
max_accuracy <- max(filtered_results$Accuracy)

# 筛选出最高准确率的结果
best_results <- filtered_results %>% filter(Accuracy == max_accuracy)

# 如果出现多个准确率相同的结果，则选择特征数量最少的组合
min_feature_count <- min(nchar(best_results$Indices))
best_results_min_features <- best_results %>% filter(nchar(Indices) == min_feature_count)

# 如果特征数量和准确率相同，则选择频率最高的组合
final_result <- best_results_min_features %>% filter(freq == max(freq))

# 将最终选择的特征组合的Indices赋值给feature_list
indices_string <- final_result$Indices[1]
indices_char <- unlist(strsplit(indices_string, ","))
feature_list <- as.numeric(indices_char)

# 显示或保存最终选择的特征组合
print(final_result)
print(paste("Selected Feature List: ", feature_list))

# 保存最终结果
# 动态生成文件名
best_feature_selection_name <- paste0("best_feature_selection_",train_num, ".csv")
write.csv(final_result, best_feature_selection_name, row.names = FALSE)

####训练俩个模型####
###传统lasso###
cat("使用Lasso进行特征选择\n")

# 使用Lasso进行特征选择并进行交叉验证
lasso_model <- glmnet(dataX0, dataY0, family = "binomial", nlambda = 1000, alpha = 1)
cv_fit <- cv.glmnet(dataX0, dataY0, family = "binomial", alpha = 1)

# 获取最佳lambda值
best_lambda <- cv_fit$lambda.min

# 使用最佳lambda值提取特征
best_model <- glmnet(dataX0, dataY0, family = "binomial", alpha = 1, lambda = best_lambda)
selected_features <- which(coef(best_model) != 0)[-1]-1
dataX0_selected_glm <- dataX0[, selected_features]
logistic_model <- glm(dataY0 ~ ., data = data.frame(dataX0_selected_glm), family = binomial)

dataX0_selected_lasso <- dataX0[, feature_list]
my_model <- glm(dataY0 ~ ., data = data.frame(dataX0_selected_lasso), family = binomial)
##进行1000次test

cl <- makeCluster(numCores)
registerDoSNOW(cl)
pb <- txtProgressBar(min = 0, max = 1000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

foreach_result_2 <- foreach(id = 1:1000, .combine = rbind, .packages = c('glmnet'),.inorder=TRUE,.options.snow=opts) %dopar% {
    random_seed <- 1234 + id
    set.seed(random_seed)

    sbjtype1=which(y_test==1)
    sbjtype0=which(y_test==0)
    sbjtype1sel=sample(sbjtype1,100)     # 1取100个
    sbjtype0sel=sample(sbjtype0,100)     # 0取100个

    sbjsel=c(sbjtype0sel,sbjtype1sel);
    test_dataX = x_test[sbjsel,]
    test_dataY = y_test[sbjsel]

    test_dataX_selected_glm = test_dataX[, selected_features]
    test_dataX_selected_lasso = test_dataX[, feature_list]

    predictions_glm <- predict(logistic_model, newdata = data.frame(test_dataX_selected_glm), type = "response")
    predicted_classes_glm <- ifelse(predictions_glm > 0.5, 1, 0)  # 使用0.5作为阈值进行分类
    glm_accuracy <- sum(predicted_classes_glm == test_dataY) / length(test_dataY)

    predictions_lasso <- predict(my_model, newdata = data.frame(test_dataX_selected_lasso), type = "response")
    predicted_classes_lasso <- ifelse(predictions_lasso > 0.5, 1, 0)  # 使用0.5作为阈值进行分类
    lasso_accuracy <- sum(predicted_classes_lasso == test_dataY) / length(test_dataY)

  return(c(lasso_accuracy, glm_accuracy,random_seed))
}

stopCluster(cl)
results_df <- as.data.frame(foreach_result_2)
colnames(results_df) <- c("LASSO_Plus_Accuracy", "GLM_Accuracy", "random_seed")

results_name <- paste0("results_",train_num, ".csv")

write.csv(results_df, file = results_name, row.names = FALSE)
cat("Results saved successfully.\n")

# 读取 CSV 文件
data <- read.csv(results_name)

# 绘制 y1 和 y2 的密度图
density_plot <- ggplot(data) + 
  geom_density(aes(x = LASSO_Plus_Accuracy, color = "LASSO_Plus_Accuracy"), size = 1) +
  geom_density(aes(x = GLM_Accuracy, color = "GLM_Accuracy"), size = 1) +
  labs(x = "Value", y = "Density", color = "Legend") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))

# # 显示图像
# print(density_plot)

# 保存图像
imgname <- paste0("plot_",train_num, ".png")
ggsave(imgname, plot = density_plot, width = 8, height = 6)
