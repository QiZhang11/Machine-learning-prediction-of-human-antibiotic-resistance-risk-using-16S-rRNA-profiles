# 设置工作目录
wd <- "C:/Users/zhangqi/Desktop/"
setwd(wd)

library(caret)
library(randomForest)

# 读取模型
rf_fit <- readRDS("rf_fit.rda")

# 读取新数据
new_data <- read.csv("prediction.csv", row.names = 1, check.names = FALSE)

# 查看训练模型需要哪些变量
required_features <- rf_fit$finalModel$xNames
print(required_features)

# 如果新数据缺少某些特征，补 0
missing_features <- setdiff(required_features, colnames(new_data))
if(length(missing_features) > 0){
  for(f in missing_features){
    new_data[[f]] <- 0
  }
}

# 如果新数据有多余特征，去掉
new_data <- new_data[, required_features, drop = FALSE]

# 开始预测
pred_risk <- predict(rf_fit, newdata = new_data)

# 输出结果
pred_result <- data.frame(
  SampleID = rownames(new_data),
  Predicted_Risk = pred_risk
)

print(head(pred_result))

# 保存结果
write.csv(pred_result, "predicted_risk_from_rf_fit.csv", row.names = FALSE)