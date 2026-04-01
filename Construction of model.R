wd<- "D:/16s-metagenome抗性基因预测/预测/genus"
setwd(wd)

#Loading R package
library(caret)
library(randomForest)
library(e1071)
library(ggplot2)
library(patchwork)
library(class)
library(ggpubr)
library(ggExtra)
library(ggbreak)

#Loading data
data <- read.csv("data_genusnolog.csv",row.names = 1)
head(data)

#Extract column name
col_names <- colnames(data)
selected_col_names <- col_names[1:156]

#The risk orabundance of Risk is determined by species abundance
formula <- as.formula(paste("Risk ~", paste(selected_col_names, collapse = "+")))
impVars <- data[, selected_col_names]

# Create Resample (80 % of input data and 20 Resamples)
set.seed(321)
train <- createDataPartition(data$Risk, p=.80,list = FALSE)
data_train <- data[train,]
data_test <- data[-train,]

#10-fold cross-validation(There were ten ten-fold cross-checks to be more accurate,
#If the data set is lRiske, it is recommended to reduce the number of iterations)
tr <- trainControl(method = "cv",number = 10)

###############################################
########################     Random Forest      ################################
###############################################
#树的数量（ntree）：增加树的数量可以提高模型的稳定性和准确性，
#但也会增加计算时间。一般来说，增加树的数量直到模型性能趋于稳定为止。
#特征数（mtry）：mtry参数控制每个决策树在分裂节点时随机选择的特征数。
#较小的mtry值会增加树之间的差异性，但可能会降低模型的准确性。
#较大的mtry值会增加模型的稳定性，但可能会导致模型过度拟合。
#一般推荐使用默认值sqrt(p)，其中p是特征的总数。
#nodesize通常用于避免生成具有非常小叶子节点的树，
#因为这可能会导致模型过拟合，对训练数据过于敏感，而泛化性能不佳。
#一般来说，你可以根据你的数据集的大小和性质来选择合适的nodesize值。
#较大的nodesize值会导致树的深度减小，模型更加简化，
#对训练数据的拟合会减弱，但泛化性能可能更好。
#较小的nodesize值会导致树的深度增加，模型更复杂，
#对训练数据的拟合更好，但泛化性能可能受到影响。
custom_RF <- list(
  type = "Regression",
  library = "randomForest",
  loop = NULL,
  parameters = data.frame(parameter = c("mtry", "ntree", "nodesize"), 
                          class = rep("numeric", 3), 
                          label = c("mtry", "ntree", "nodesize")),
  
  grid = function(x, y, len = NULL, search = "grid") {
    if(search == "grid") {
      out <- expand.grid(mtry = caret::var_seq(p = ncol(x),
                                               classification = FALSE,
                                               len = len),
                         ntree = c(500, 800, 1000, 1200, 1500, 2000),
                         nodesize = c(3, 4, 5, 6, 7, 8, 9, 10))  # add nodesize parameter
    } else {
      out <- data.frame(mtry = unique(sample(1:ncol(x), size = len, replace = TRUE)),
                        ntree = unique(sample(c(500, 700, 900, 1000, 1500), 
                                              size = len, replace = TRUE)),
                        nodesize = unique(sample(c(3, 4, 5, 6, 7, 8, 9, 10), 
                                                 size = len, replace = TRUE)))  # add nodesize parameter
    }
  },
  fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    randomForest(x, y, 
                 mtry = param$mtry, 
                 ntree = param$ntree, 
                 nodesize = param$nodesize, ...)
  },
  predict = function(modelFit, 
                     newdata, 
                     preProc = NULL, 
                     submodels = NULL) {
    predict(modelFit, newdata)
  },
  prob = NULL,
  sort = NULL
)

# train model
if(file.exists('rf_fit.rda')) {
  rf_fit <- readRDS("rf_fit.rda")
} else {
  tunegrid <- expand.grid(mtry = c(10, 20, 30, 40, 50, 100), 
                          ntree = c(500, 800, 1000, 1200, 1500, 2000),
                          nodesize = c(3, 4, 5, 6, 7, 8, 9, 10))
  
  rf_fit <- train(Risk ~ .,
                  data = data_train,
                  method = custom_RF,  # This is the regression model
                  metric = "Rsquared",
                  tuneGrid = tunegrid,
                  trControl = tr)
  saveRDS(rf_fit, "rf_fit.rda")
}
#output result of Randomforest
result_rf <- as.data.frame(rf_fit$result)
write.csv(result_rf,"Randomforest_fit.csv")
result_rf$ntree <- factor(result_rf$ntree)
result_rf$nodesize <- factor(result_rf$nodesize)

# Draw a line chart
P1 <-ggplot(data=result_rf, aes(x=mtry,y= Rsquared, group=ntree, color=ntree)) +
  geom_line(size = 0.5)+
  geom_point(size = 2)+
  scale_color_manual(values = c("#4169B2","#B1A4C0","#479E9B",
                                         "#BB2BA0","#DDA0DD","#BC8F8F"))+
                                           facet_wrap(~nodesize)+ 
  theme_bw(base_size = 15)
P1
##########################################
###################     Support vector machine      ###########################
##########################################

svm_grid <- expand.grid(
  sigma = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1),
  C = c(1, 2, 3, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 5))

svm_fit <- train(
  Risk ~ .,
  data = data_train,
  method = "svmRadial",
  tuneGrid = svm_grid,
  trControl = tr)

print(svm_fit)
saveRDS(svm_fit, "svm_fit.rda")

#output result of SVM
result_svm <- as.data.frame(svm_fit$results)
write.csv(result_svm,"SVM_fit.csv")
result_svm$sigma <- factor(result_svm$sigma)

# Draw a line chart
P2 <-ggplot(data=result_svm, aes(x=C,y=Rsquared, group=sigma, color=sigma)) +
  geom_line(size = 0.5)+
  geom_point(size = 2)+
  scale_color_manual(values = c("#4169B2","#B1A4C0","#479E9B",
                                         "#BB2BA0","#DDA0DD","#BC8F8F",
                                         "#FFDAB9","#B4EEB4","#99CCFF",
                                         "#AC9179","#CDD7CB","#594B69",
                                         "#DBCD9D","#73844F","#D499A4")) +
                                           theme_bw(base_size = 15)
P2

##########################################
################   XGBoost（eXtreme Gradient Boosting   #####################
##########################################

xgb_grid <- expand.grid(
  nrounds = c(300),
  # Iteration rounds
  max_depth = c(2, 4, 6, 8, 10, 12), 
  # maximum depth of a tree
  eta = c(0.1), 
  # control the learning rate: scale the contribution of each tree by
  # a factor of 0 < eta < 1 when it is added to the current approximation. 
  # Used to prevent overfitting by making the boosting process more conservative.
  # Lower value for eta implies lRisker value for nrounds: 
  # low eta value means model more robust to overfitting but slower to compute.
  gamma = c(0.1, 0.2, 0.3, 0.4, 0.5,0.6), 
  # minimum loss reduction required to make a further partition on a leaf node of the tree. 
  # the lRisker, the more conservative the algorithm will be.
  colsample_bytree = c(0.8),
  # subsample ratio of columns when constructing each tree.
  min_child_weight = c(1, 2, 3, 4, 5, 6), 
  # minimum sum of instance weight (hessian) needed in a child. 
  # If the tree partition step results in a leaf node with the sum of instance 
  # weight less than min_child_weight, then the building process will give up 
  # further partitioning. In linear regression mode, this simply corresponds to 
  # minimum number of instances needed to be in each node. 
  # The lRisker, the more conservative the algorithm will be.（min_child_weight）
  subsample = c(0.8))
# subsample ratio of the training instance. Setting it to 0.5 means that
# xgboost randomly collected half of the data instances to grow trees and this
# will prevent overfitting. It makes computation shorter (because less data to analyse). 
# It is advised to use this parameter with eta and increase nrounds. 

xgb_fit <- train(
  Risk ~ .,
  data = data_train,
  method = "xgbTree",
  tuneGrid = xgb_grid,
  trControl = tr)
xgb_fit
saveRDS(xgb_fit, "xgb_fit.rda")

#output result of XGB
result_xgb <- as.data.frame(xgb_fit$results)
write.csv(result_xgb,"XGB_fit.csv")
result_xgb$gamma <- factor(result_xgb$gamma)
result_xgb$min_child_weight <- factor(result_xgb$min_child_weight)
# Draw a line chart
P3 <-ggplot(data=result_xgb, aes(x=max_depth,y=Rsquared, group=gamma, color=gamma)) +
  geom_line(size = 0.5)+
  geom_point(size = 2)+
  scale_color_manual(values = c("#4169B2","#B1A4C0","#479E9B",
                                         "#BB2BA0","#DDA0DD","#BC8F8F"))+
                                           facet_wrap(~min_child_weight)+ 
  theme_bw(base_size = 15)
P3

##########################################
################        KNN  k-nearest neighbor        #####################
##########################################
# k值是一个重要的参数，代表着在进行分类或回归任务时，用于确定最近邻居数量的超参数。
# KNN算法的基本思想是根据输入样本的特征，在训练数据集中找到距离该样本最近的k个数据点，
# 然后根据这k个数据点的类别（对于分类任务）或平均值（对于回归任务）来预测该样本的类别或值

knn_grid <- expand.grid(k = c(2, 4, 8, 16, 32, 64, 128, 256, 500))

knn_fit_triangular <- train(
  Risk ~ .,
  data = data_train,
  method = "knn",
  tuneGrid = knn_grid,
  trControl = tr,
  kernel = "triangular")
# Kernel to use. Possible choices are "rectangular" (which is standard unweighted knn), 
# "triangular", "epanechnikov" (or beta(2,2)), "biweight" (or beta(3,3)), 
# "triweight" (or beta(4,4)), "cos", "inv", "gaussian" and "optimal".

print(knn_fit_epanechnikov)
saveRDS(knn_fit_rectangular, "knn_fit_rectangular.rda")
saveRDS(knn_fit_triangular, "knn_fit_triangular.rda")
saveRDS(knn_fit_epanechnikov, "knn_fit_epanechnikov.rda")
saveRDS(knn_fit_gaussian, "knn_fit_gaussian.rda")
saveRDS(knn_fit_optimal, "knn_fit_optimal.rda")
saveRDS(knn_fit_triweight, "knn_fit_triweight.rda")

#output result of knn
knn_fit_rectangular <- as.data.frame(knn_fit_rectangular$results)
write.csv(knn_fit_rectangular,"knn_fit_rectangular.csv")

knn_fit_triangular <- as.data.frame(knn_fit_triangular$results)
write.csv(knn_fit_triangular,"knn_fit_triangular.csv")

knn_fit_epanechnikov <- as.data.frame(knn_fit_epanechnikov$results)
write.csv(knn_fit_epanechnikov,"knn_fit_epanechnikov.csv")

knn_fit_gaussian <- as.data.frame(knn_fit_gaussian$results)
write.csv(knn_fit_gaussian,"knn_fit_gaussian.csv")

knn_fit_optimal <- as.data.frame(knn_fit_optimal$results)
write.csv(knn_fit_optimal,"knn_fit_optimal.csv")

result_knn <- as.data.frame(knn_fit_triweight$results)
write.csv(result_knn,"knn_fit_triweight.csv")

knn_fit <- read.csv("knn_fit.csv")
knn_fit$kernel <- factor(knn_fit$kernel)
# Draw a line chart
P4 <-ggplot(data=knn_fit, aes(x=k,y=Rsquared, group=kernel, color=kernel)) +
  geom_line(size = 0.5)+
  geom_point(size = 2)+
  scale_color_manual(values = c("#4169B2","#B1A4C0","#479E9B",
                                         "#BB2BA0","#DDA0DD","#BC8F8F",
                                         "#FFDAB9","#B4EEB4","#99CCFF",
                                         "#AC9179","#CDD7CB","#594B69",
                                         "#DBCD9D","#73844F","#D499A4")) +
                                           theme_bw(base_size = 15)
P4

#################################################
############################         Lasso         #############################
#################################################

Lasso_grid <- expand.grid(
  alpha = c(1,0),
  lambda = c(1, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10,15,20))
# glmnet()中alpha=0表示岭回归，alpha=1表示lasso回归

Lasso_fit_poisson <- train(
  Risk ~ .,
  data = data_train,
  method = "glmnet",
  tuneGrid = Lasso_grid,
  trControl = tr,
  family = "poisson")
#"gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"
print(Lasso_fit_gaussian)
print(Lasso_fit_poisson)
saveRDS(Lasso_fit_gaussian, "Lasso_fit_gaussian.rda")
saveRDS(Lasso_fit_poisson, "Lasso_fit_poisson.rda")

#output result of Lasso
result_Lasso <- as.data.frame(Lasso_fit_poisson$results)
write.csv(result_Lasso,"Lasso_fit_poisson.csv")

result_Lasso <- as.data.frame(Lasso_fit_gaussian$results)
write.csv(result_Lasso,"Lasso_fit_gaussian.csv")

Lasso_fit <- read.csv("Lasso_fit.csv")
Lasso_fit$family <- factor(Lasso_fit$family)
Lasso_fit$alpha <- factor(Lasso_fit$alpha)

# Draw a line chart
P5 <-ggplot(data=Lasso_fit, aes(x=lambda,y=Rsquared, group=family, color=family)) +
  geom_line(size = 0.5)+
  geom_point(size = 2)+
  scale_color_manual(values = c("#4169B2","#B1A4C0","#479E9B",
                                         "#BB2BA0","#DDA0DD","#BC8F8F",
                                         "#FFDAB9","#B4EEB4","#99CCFF",
                                         "#AC9179","#CDD7CB","#594B69",
                                         "#DBCD9D","#73844F","#D499A4"))+
                                           theme_bw(base_size = 15)
P5

##########################################
#################            Neural Network            #####################
##########################################
#size表示隐藏层的节点个数，通常为自变量个数的1.2-1.5倍，
#size为0时表示无隐藏层；rang表示初始权重范围[-rang,rang]，
#一般rang与x的绝对值中的最大值的乘积大约等于1；decay表示模型权重值的衰减精度；

NN_grid <- expand.grid(
  size = c( 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,2),
  decay = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1))

NN_fit <- train(
  Risk ~ .,
  data = data_train,
  method = "nnet",
  tuneGrid = NN_grid,
  trControl = tr)

print(NN_fit)
saveRDS(NN_fit, "NN_fit.rda")
#output result of NN
result_NN <- as.data.frame(NN_fit$results)
# 删除包含 NaN 值的行
#result_NN <- subset(result_NN, !is.nan(Rsquared))
write.csv(result_NN,"NN_fit.csv")
result_NN <- read.csv("NN_fit.csv")
result_NN$decay <- factor(result_NN$decay)

# Draw a line chart
P6 <-ggplot(data=result_NN, aes(x=size,y=Rsquared, group=decay, color=decay)) +
  geom_line(size = 0.5)+
  geom_point(size = 2)+
  scale_color_manual(values = c("#4169B2","#B1A4C0","#479E9B",
                                         "#BB2BA0","#DDA0DD","#BC8F8F",
                                         "#FFDAB9","#B4EEB4","#99CCFF",
                                         "#AC9179","#CDD7CB","#594B69",
                                         "#DBCD9D","#73844F","#D499A4")) +
                                           theme_bw(base_size = 15)
P6

#Merge all model tuning result plots
p_fit <- (P1+P2+P3)/(P4+P5+P6)
p_fit

#Model comparison under gaussian parameters

################Random forest#########################

rf_grid <- expand.grid(mtry = 50)

rf_best <- train(Risk ~ .,
                 data = data,
                 method = "rf",
                 metric = "Rsquared",  
                 trControl = tr,
                 tuneGrid = rf_grid,
                 ntree = 800,
                 nodesize = 7)

rf_best$results
print(rf_best)

results_rf <- data.frame(Variables = integer(0), R_squared = numeric(0))

# 从一个自变量开始逐步增加自变量数量
for (num_vars in 1:(ncol(data) - 1)) {
  # 选择自变量
  predictors <- names(data)[1:num_vars]
  formula <- as.formula(paste("Risk ~", paste(predictors, collapse = " + ")))
  rf_grid <- expand.grid(mtry = 50)
  # 构建线性回归模型
  model_rf <- train(formula, data = data, method = "rf", trControl = tr,
                    ntree = 800,
                    tuneGrid = rf_grid,
                    nodesize = 7)
  
  # 计算模型性能
  r_squared_rf <- model_rf$results$Rsquared
  
  # 存储结果
  results_rf <- rbind(results_rf, data.frame(Variables = num_vars, R_squared = r_squared_rf))
}

# 打印或可视化结果
print(results_rf)
write.csv(results_rf,"results_rf.csv")

# 创建一个新的数据框，添加起始点 (0, 0)
start_point <- data.frame(Variables = 0, R_squared = 0)

P_RF <-ggplot(data = rbind(start_point, results_rf), aes(x=Variables,y=R_squared)) +
  geom_line(size=1,color = "#4169B2")+
  # geom_smooth(method = "glm", color = "#4169B2")+
  geom_point(size = 4,color = "#dc3023",alpha = 0.6) +
  labs(x = 'Number of variables', y = '10-fold cross validation R2')+
  theme_bw(base_size = 15)
P_RF

################ Support vector machine#########################
svm_grid <- expand.grid(
  sigma = c(0.01),
  C = c(3.5))

svm_best <- train(
  Risk ~ .,
  data = data,
  method = "svmRadial",
  trControl = tr,
  tuneGrid = svm_grid)

svm_best$results
print(svm_best)

results_svm <- data.frame(Variables = integer(0), R_squared = numeric(0))

# 从一个自变量开始逐步增加自变量数量
for (num_vars in 1:(ncol(data) - 1)) {
  # 选择自变量
  predictors <- names(data)[1:num_vars]
  formula <- as.formula(paste("Risk ~", paste(predictors, collapse = " + ")))
  # 构建线性回归模型
  model_svm <- train(formula, data = data, method = "svmRadial", trControl = tr,
                     tuneGrid = svm_grid)
  
  # 计算模型性能
  r_squared_svm <- model_svm$results$Rsquared
  
  # 存储结果
  results_svm <- rbind(results_svm, data.frame(Variables = num_vars, R_squared = r_squared_svm))
}

# 打印或可视化结果
print(results_svm)
write.csv(results_svm,"results_svm.csv")

# 创建一个新的数据框，添加起始点 (0, 0)
start_point <- data.frame(Variables = 0, R_squared = 0)

P_svm <-ggplot(data = rbind(start_point, results_svm), aes(x=Variables,y=R_squared)) +
  geom_line(size=1,color = "#4169B2")+
  # geom_smooth(method = "glm", color = "#4169B2")+
  geom_point(size = 4,color = "#dc3023",alpha = 0.6) +
  labs(x = 'Number of variables', y = '10-fold cross validation R2',title = "SVM R2:0.3698279")+
  theme_bw(base_size = 15)
P_svm

################ XGBoost（eXtreme Gradient Boosting #########################

xgb_grid <- expand.grid(
  nrounds = c(300), # 迭代轮数（nrounds）
  max_depth = c(6), # 最大树深度（max_depth）
  eta = c(0.1), # 学习率（eta）
  gamma = c(0.3), # 树分裂所需的最小损失减少值
  colsample_bytree = c(0.8), # 特征子采样比例（colsample_bytree）
  min_child_weight = c(3), # 叶子节点的最小权重和（min_child_weight）
  subsample = c(0.8)) # 和样本子采样比例（subsample）

xgb_best <- train(
  Risk ~ .,
  data = data,
  method = "xgbTree",
  tuneGrid = xgb_grid,
  trControl = tr)
xgb_best

xgb_best$results
print(xgb_best)

results_xgb <- data.frame(Variables = integer(0), R_squared = numeric(0))

# 从一个自变量开始逐步增加自变量数量
for (num_vars in 1:(ncol(data) - 1)) {
  # 选择自变量
  predictors <- names(data)[1:num_vars]
  formula <- as.formula(paste("Risk ~", paste(predictors, collapse = " + ")))
  rf_grid <- expand.grid(mtry = 10)
  # 构建线性回归模型
  model_xgb <- train(formula, data = data, method = "xgbTree", trControl = tr,
                     tuneGrid = xgb_grid)
  
  # 计算模型性能
  r_squared_xgb <- model_xgb$results$Rsquared
  
  # 存储结果
  results_xgb <- rbind(results_xgb, data.frame(Variables = num_vars, R_squared = r_squared_xgb))
}

# 打印或可视化结果
print(results_xgb)
write.csv(results_xgb,"results_xgb.csv")

# 创建一个新的数据框，添加起始点 (0, 0)
start_point <- data.frame(Variables = 0, R_squared = 0)

P_xgb <-ggplot(data = rbind(start_point, results_xgb), aes(x=Variables,y=R_squared)) +
  geom_line(size=1,color = "#4169B2")+
  # geom_smooth(method = "glm", color = "#4169B2")+
  geom_point(size = 2,color = "#dc3023",alpha = 0.6) +
  labs(x = 'Number of variables', y = '10-fold cross validation R2')+
  theme_bw(base_size = 15)
P_xgb

################ KNN  k-nearest neighbor #########################

knn_grid <- expand.grid(k = c(16))

knn_best <- train(
  Risk ~ .,
  data = data,
  method = "knn",
  tuneGrid = knn_grid,
  trControl = tr,
  kernel = "triangular")

print(knn_best)

results_knn <- data.frame(Variables = integer(0), R_squared = numeric(0))

# 从一个自变量开始逐步增加自变量数量
for (num_vars in 1:(ncol(data_train) - 1)) {
  # 选择自变量
  predictors <- names(data_train)[1:num_vars]
  formula <- as.formula(paste("Risk ~", paste(predictors, collapse = " + ")))
  # 构建线性回归模型
  model_knn <- train(formula, data = data_train, method = "knn", trControl = tr,
                     tuneGrid = knn_grid)
  
  # 计算模型性能
  r_squared_knn <- model_knn$results$Rsquared
  
  # 存储结果
  results_knn <- rbind(results_knn, data.frame(Variables = num_vars, R_squared = r_squared_knn))
}

# 打印或可视化结果
print(results_knn)

# 创建一个新的数据框，添加起始点 (0, 0)
start_point <- data.frame(Variables = 0, R_squared = 0)

P_knn <-ggplot(data = rbind(start_point, results_knn), aes(x=Variables,y=R_squared)) +
  geom_line(size=1,color = "#4169B2")+
  # geom_smooth(method = "glm", color = "#4169B2")+
  geom_point(size = 4,color = "#dc3023",alpha = 0.6) +
  labs(x = 'Number of variables', y = '10-fold cross validation R2',title = "knn R2:0.3577008")+
  theme_bw(base_size = 15)
P_knn

############################### Lasso ###################################

Lasso_grid <- expand.grid(
  alpha = c(0),
  lambda = c(20))

Lasso_best <- train(
  Risk ~ .,
  data = data,
  method = "glmnet",
  family = "poisson",
  tuneGrid = Lasso_grid,
  trControl = tr)

print(Lasso_best)

results_Lasso <- data.frame(Variables = integer(0), R_squared = numeric(0))
# 设置初始的 num_vars 为2
num_vars <- 2
# 开始逐步增加自变量数量的循环
for (num_vars in num_vars:(ncol(data_train) - 1)) {
  # 选择自变量
  predictors <- names(data_train)[1:num_vars]
  formula <- as.formula(paste("Risk ~", paste(predictors, collapse = " + ")))
  # 构建线性回归模型
  model_Lasso <- train(formula, data = data_train, method = "glmnet", trControl = tr,
                       tuneGrid = Lasso_grid)
  
  # 计算模型性能
  r_squared_Lasso <- model_Lasso$results$Rsquared
  
  # 存储结果
  results_Lasso <- rbind(results_Lasso, data.frame(Variables = num_vars, R_squared = r_squared_Lasso))
}

# 打印或可视化结果
print(results_Lasso)

# 创建一个新的数据框，添加起始点 (0, 0)
start_point <- data.frame(Variables = 0, R_squared = 0)

P_Lasso <-ggplot(data = rbind(start_point, results_Lasso), aes(x=Variables,y=R_squared)) +
  geom_line(size=1,color = "#4169B2")+
  # geom_smooth(method = "glm", color = "#4169B2")+
  geom_point(size = 4,color = "#dc3023",alpha = 0.6) +
  labs(x = 'Number of variables', y = '10-fold cross validation R2',title = "Lasso R2:0.2018684")+
  theme_bw(base_size = 15)
P_Lasso

########################### Neural Network  ###################################

NN_grid <- expand.grid(
  size = c(2),
  decay = c(0.01))

NN_best <- train(
  Risk ~ .,
  data = data,
  method = "nnet",
  tuneGrid = NN_grid,
  trControl = tr)

print(NN_best)

results_NN <- data.frame(Variables = integer(0), R_squared = numeric(0))

# 开始逐步增加自变量数量的循环
for (num_vars in num_vars:(ncol(data_train) - 1)) {
  # 选择自变量
  predictors <- names(data_train)[1:num_vars]
  formula <- as.formula(paste("Risk ~", paste(predictors, collapse = " + ")))
  # 构建线性回归模型
  model_NN <- train(Risk ~ ., data = data_train, method = "nnet", tuneGrid = NN_grid,
                    trControl = tr)
  
  # 计算模型性能
  r_squared_NN <- model_NN$results$Rsquared
  
  # 存储结果
  results_NN <- rbind(results_NN, data.frame(Variables = num_vars, R_squared = r_squared_NN))
}


# 打印或可视化结果
print(results_NN)

# 创建一个新的数据框，添加起始点 (0, 0)
start_point <- data.frame(Variables = 0, R_squared = 0)

P_NN <-ggplot(data = rbind(start_point, results_NN), aes(x=Variables,y=R_squared)) +
  geom_line(size=1,color = "#4169B2")+
  # geom_smooth(method = "glm", color = "#4169B2")+
  geom_point(size = 4,color = "#dc3023",alpha = 0.6) +
  labs(x = 'Number of variables', y = '10-fold cross validation R2',title = "NN R2 :0.01390247")+
  theme_bw(base_size = 15)
P_NN

#合并上述表格
new_row <- c("1","0")
results_Lasso_new <- rbind(results_Lasso,new_row)

dat_cbind<-cbind(results_rf, results_svm, results_knn, results_Lasso_new, results_xgb, results_NN)
write.csv(dat_cbind, "Genus_10-fold_R2-variables_best-model-parameters.csv")

#------------------------------------------------------------------------------#
P_cv_best_fit <- (P_RF + P_svm + P_xgb)/(P_knn + P_Lasso + P_NN)
P_cv_best_fit
######################################预测######################################
result_10_fold <- read.csv("results_10_fold.csv",row.names = 1)

P_10_fold <-ggplot(data = result_10_fold, aes(x=Variables,y=R_squared, group = group)) +
  geom_line(size=1,color = "#4169B2")+
  # geom_smooth(method = "glm", color = "#4169B2")+
  geom_point(size = 4,aes(color = group),alpha = 0.7) +
  scale_color_manual(values=c("#007EB1","#00A074","#821193"))+
  labs(x = 'Number of variables', y = '10-fold cross validation R2')+
  annotate('text', label = 'RF: R2 = 0.5134, RMSE =0.4064', x =100, y = 0.05, size =4,color="black")+
  annotate('text', label = 'XGB: R2 = 0.5008, RMSE = 0.4077', x =100, y = 0.1, size =4,color="black")+
  annotate('text', label = 'SVM: R2 = 0.4462, RMSE =0.4405', x =100, y = 0.15, size =4,color="black")+
  theme(panel.background =element_rect(color ='black',fill ='White'),
        panel.grid =element_blank(),
        legend.position ='none')
P_10_fold

#在测试集检验
pred_test_rf <- predict(rf_best, data_test)
pred_test_rf <- data.frame(pred = pred_test_rf, Risk = data_test$Risk)

#测试集合线性回归
lm_test_rf <- lm(Risk ~ pred, data = pred_test_rf)
summary(lm_test_rf)

#绘制结果
p_test_rf <- ggplot(pred_test_rf, aes(Risk, pred))+
  geom_point(size=5, col="#2171b5",alpha=0.6)+ 
  geom_smooth(method="lm",se = T,size=1.5) + 
  geom_abline(intercept = 0, slope = 1, size=1.5, col="red")+
  scale_y_continuous(expand = c(0, 0), limit = c(5.5,9.5))+
  scale_x_continuous(expand = c(0, 0), limit = c(5.5,9.5))+
  annotate('text', label = 'RF:R2 = 0.9566, P < 2.2e-16, RSE = 0.1242', x =6.5, y = 9, size =4,color="black")+
  labs(x = "Observed ARGs risk", y = "Predicted ARGs risk")+
  theme_bw()
p_test_rf

#在测试集检验
pred_test_svm <- predict(svm_best, data_test)
pred_test_svm <- data.frame(pred = pred_test_svm, Risk = data_test$Risk)

#测试集合线性回归
lm_test_svm <- lm(Risk ~ pred, data = pred_test_svm)
summary(lm_test_svm)

#绘制结果
p_test_svm <- ggplot(pred_test_svm, aes(Risk, pred))+
  geom_point(size=5, col="#2171b5",alpha=0.6)+ 
  geom_smooth(method="lm",se = T,size=1.5) + 
  geom_abline(intercept = 0, slope = 1, size=1.5, col="red")+
  scale_y_continuous(expand = c(0, 0), limit = c(5.5,9.5))+
  scale_x_continuous(expand = c(0, 0), limit = c(5.5,9.5))+
  annotate('text', label = 'SVM:R2 = 0.883, P < 2.2e-16, RSE = 0.2038', x =6.5, y = 9, size =4,color="black")+
  labs(x = "Observed ARGs risk", y = "Predicted ARGs risk")+
  theme_bw()
p_test_svm

#在测试集检验
pred_test_xgb <- predict(xgb_best, data_test)
pred_test_xgb <- data.frame(pred = pred_test_xgb, Risk = data_test$Risk)

#测试集合线性回归
lm_test_xgb <- lm(Risk ~ pred, data = pred_test_xgb)
summary(lm_test_xgb)

#绘制结果
p_test_xgb <- ggplot(pred_test_xgb, aes(Risk, pred))+
  geom_point(size=5, col="#2171b5",alpha=0.6)+ 
  geom_smooth(method="lm",se = T,size=1.5) + 
  geom_abline(intercept = 0, slope = 1, size=1.5, col="red")+
  scale_y_continuous(expand = c(0, 0), limit = c(5.7,9.2))+
  scale_x_continuous(expand = c(0, 0), limit = c(5.7,9.2))+
  annotate('text', label = 'XGB:R2 = 0.9547, P < 2.2e-16, RSE = 0.1268', x =6.5, y = 9, size =4,color="black")+
  labs(x = "Observed ARGs risk", y = "Predicted ARGs risk")+
  theme_bw()
p_test_xgb
