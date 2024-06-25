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

# currdir = dirname(this.path())
# setwd(currdir)

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

  #dataY_pred=predict(fit,dataX,lambda=hyper,type='class')
  #accmax=sum(dataY_pred==dataY)/length(dataY)

  #dataYtest_pred=predict(fit,dataXtest,lambda=hyper,type='class')
  #accmax_test=sum(dataYtest_pred==dataYtest)/length(dataYtest)
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
  # feature_table_scaled[] <- lapply(feature_table_scaled, as.numeric)   # 转换为数值型数据
  # feature_table_scaled_test[] <- lapply(feature_table_scaled_test, as.numeric)

  # 逻辑回归
  logistic_model <- glm(dataY_glm ~ ., data = feature_table_scaled, family = binomial)
  # #训练集
  # train_predictions <- predict(logistic_model, newdata = feature_table_scaled, type = "response")
  # train_pred_classes <- ifelse(train_predictions > 0.5, 1, 0)
  # train_accuracy <- sum(train_pred_classes == dataY_glm) / length(dataY_glm)
  # #cat("训练集准确率:", train_accuracy, "\n")
  # #测试集
  test_predictions <- predict(logistic_model, newdata = feature_table_scaled_test, type = "response")
  test_pred_classes <- ifelse(test_predictions > 0.5, 1, 0)
  test_accuracy <- sum(test_pred_classes == dataYtest_glm) / length(dataYtest_glm)
  #cat("测试集准确率:", test_accuracy, "\n")
  return(list(model = logistic_model, accuracy = test_accuracy))
}

# =============== load data
loaddata=function(OCDfile,sheetname,sbjsel_flag) {
  # ============== load data ================
  data0=read_xlsx(OCDfile,'缓解率')
  data1=read_xlsx(OCDfile,sheetname)
  tmpdata_rTMS=matrix(1,31,1)
  tmpdata_cTBS=matrix(2,48,1)
  tmpdata_TMS=as.data.frame(rbind(tmpdata_rTMS,tmpdata_cTBS))
  tmpdata0=cbind(tmpdata_TMS,data0[,c(4:6,8,10:12)])
  extralist=c('TMStype','sex','age',
              'HAMApre','HAMDpre',
              'ThoughtPre','BehavePre','Y_BOCSpre')
  names(tmpdata0)[1:length(extralist)]=extralist
  
  if (sbjsel_flag==1) {
    sbjsel=tmpdata0$ThoughtPre>=10 | tmpdata0$BehavePre>=10 | tmpdata0$Y_BOCSpre>=15
    dataX=cbind(tmpdata0[sbjsel,],data1[sbjsel,2:(roimax+1)])
    dataY=as.factor(data1$response[sbjsel])
  } else {
    dataX=cbind(tmpdata0,data1[,2:(roimax+1)])
    dataY=as.factor(data1$response)
  }
  dataX[]=lapply(dataX,as.numeric)
  dataX=scale(dataX)     # SVM中需要标准化
  return(list(dataX, dataY))
}



# # 下载数据集文件并读取
# url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/arrhythmia/arrhythmia.data"
# data <- read_csv(url, col_names = FALSE)

data <- read.csv("TUNADROMD.csv")

# # 设置列名（特征列）
# feature_names <- paste0("V", 1:279)
# colnames(data) <- c(feature_names, "Class")

# # 将 Class 列转换为二分类：1 为正常心律，其他为不正常心律
# data <- data %>%
#   mutate(Class = ifelse(Class == 1, 1, 0))

# # 选取了1到279列，特征量
# selected_features <- feature_names[1:279]
# filtered_data <- data %>% select(all_of(selected_features), Class)

# # 随机打乱数据集
# set.seed(123) # 为了可重复性设置随机种子
# shuffled_data <- filtered_data %>% sample_frac()

# # 随机选择 150 个样本作为训练集
# # train_data <- shuffled_data %>% slice(1:150)

# # 将剩余的样本作为测试集
# test_data <- shuffled_data %>% slice(151:n())

# # 提取训练集特征和标签
# x_train <- as.matrix(train_data %>% select(-Class))
# y_train <- train_data$Class

X <- as.matrix(data[, -ncol(data)])  # 排除目标变量列
y <- as.factor(data[, ncol(data)])  # 目标变量

set.seed(123)
train_indices <- sample(1:nrow(X), 150)
X_train <- X[train_indices, ]
y_train <- y[train_indices]
X_test <- X[-train_indices, ]
y_test <- y[-train_indices]


dataX0=x_train
dataY0=y_train

# 固定的 lambda 值
#fixed_lambda <- 0.1
grid <-  10^seq(2, -4, length = 100)
randseed <- 1235673
numCores <- detectCores()
print(numCores)
cl <- makeCluster(numCores)
#registerDoParallel(cl)
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
  #   next
  # }
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
          #cat(i,"-accuracy",accuracy,"m_accuracy",a_accuracy,"\n")
        }
        #cat("测试集平均准确率:", m_accuracy, "\n")

        #我也不晓得为什么会出现长度为0
        # if (length(m_accuracy)!=0 && length(best_accuracy)!=0 ){
        #   next
        # }
        min_selected_count <- min(feature_counts[spilted_indices])
        
        
      }
      local_results <- rbind(local_results, data.frame(FixedLambda = fixed_lambda, Indices = paste(spilted_indices, collapse = ","), Accuracy = m_accuracy,freq = min_selected_count))
    }
  }
  return(local_results)
}
print("结束啦")
#print(foreach_result)
results <- foreach_result

# 保存所有特征索引及其准确率
write.csv(results, "feature_selection_results.csv", row.names = FALSE)

#write_xlsx(results, "feature_selection_results.xlsx")

stopCluster(cl)
