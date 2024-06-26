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

currdir = dirname(this.path())
setwd(currdir)

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
  set.seed(randseed)
  sbjtype1=which(dataY0==1)
  sbjtype0=which(dataY0==0)
  sbjtype1sel=sample(sbjtype1,length(sbjtype1)-5)     # 治愈组随机减去5人
  sbjtype0sel=sample(sbjtype0,length(sbjtype0)-5)     # 非治愈组随机减去5人
  # 剩余样本，为training样本。剩余的5人，为test样本（在特征筛选阶段，test样本不起作用）

  sbjsel=c(sbjtype0sel,sbjtype1sel);

  dataX_glm=feature_table[sbjsel,]
  dataY_glm=dataY0[sbjsel]

  dataXtest_glm=feature_table[-sbjsel,]
  dataYtest_glm=dataY0[-sbjsel]

  feature_table_scaled <- as.data.frame(dataX_glm)
  feature_table_scaled_test <- as.data.frame(dataXtest_glm)

  # 逻辑回归
  logistic_model <- glm(dataY_glm ~ ., data = feature_table_scaled, family = binomial)
  test_predictions <- predict(logistic_model, newdata = feature_table_scaled_test, type = "response")
  test_pred_classes <- ifelse(test_predictions > 0.5, 1, 0)
  test_accuracy <- sum(test_pred_classes == dataYtest_glm) / length(dataYtest_glm)
  return(list(model = logistic_model, accuracy = test_accuracy))
}



data <- read.csv("TUANDROMD.csv")

x <- as.matrix(data[, -ncol(data)])  # 排除目标变量列
y <- as.factor(data[, ncol(data)])  # 目标变量

set.seed(123)
train_indices <- sample(1:nrow(x), 150)
x_train <- x[train_indices, ]
y_train <- y[train_indices]

x_test <- x[-train_indices, ]
y_test <- y[-train_indices]


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
write.csv(results, "feature_selection_results.csv", row.names = FALSE)
stopCluster(cl)


# 读取结果数据
results <- read.csv("feature_selection_results.csv")

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
indices_string <- final_result$Indices
indices_char <- unlist(strsplit(indices_string, ","))
feature_list <- as.numeric(indices_char)

# 显示或保存最终选择的特征组合
print(final_result)
print(paste("Selected Feature List: ", feature_list))

# 保存最终结果
write.csv(final_result, "best_feature_selection.csv", row.names = FALSE)

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
selected_features <- which(coef(best_model) != 0)[-1]
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
    set.seed(Sys.time() + id)

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

  return(c(lasso_accuracy, glm_accuracy))
}

stopCluster(cl)
results_df <- as.data.frame(foreach_result_2)
colnames(results_df) <- c("LASSO_Plus_Accuracy", "GLM_Accuracy")
write.csv(results_df, file = "results.csv", row.names = FALSE)
cat("Results saved successfully.\n")
