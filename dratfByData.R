# 1. data struction
# y2, y1, x1, .....x10
# 1   1
# 3   1
# 2   2
# where 3 is the swing group, y2 is the real classification, y1 is the present classification
##  label
y1 = c(rep(1,200),rep(2,100))
y2 = c(rep(1,100), rep(3,100), rep(2,100))
## x1-x6 1=3
x1 = c(rnorm(200,0,0.1),  rnorm(100, 1,0.05))
x2 = c(rnorm(200,0.5,0.2),  rnorm(100,0,0.1))
x3 = c(rnorm(200,0.7,0.2),  rnorm(100,0.5,0.1))
x4 = c(rnorm(200,0.5,0.1),rnorm(100,0.8,0.1))
x5 = c(rnorm(200,0.2,0.1),  rnorm(100,0.22,0.2))
x6 = c(rnorm(200,0.15,0.2),rnorm(100,0.17,0.2))

## x7-x8 1=2
x7 = c(rnorm(100, 0.7,0.1), rnorm(100,0.5,0.04), rnorm(100, 0.7,0.1))
x8 = c(rnorm(100,0.3,0.2), rnorm(100,0.2,0.1),rnorm(100, 0.3,0.2))

## x9-x10 2=3
x9 = c(rnorm(100, 0.1,0.1), rnorm(200,0.6,0.1))
x10 = c(rnorm(100, 0.13,0.2), rnorm(200, 0.1,0.2))

data = data.frame(cbind(y1,y2,x1, x2, x3, x4, x5, x6, x7, x8, x9, x10))
data[,1] = as.factor(data[,1])
data[,2] = as.factor(data[,2])


# 2. tree-classify with present label
library(randomForest)
tree = randomForest(y1~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = data, mtry = 4, importance=TRUE, proximity = TRUE)
print(tree$importance)# 结果显示并不能很好的筛选出决定性特征

# 3. CASGI(ours)
## select the feature
x = data[,3:12]
mNum = 4
TimesTry = 1000
id_selected = c(sample(c(seq(length(x))),size=mNum))
x_selected = c()
for(i in id_selected){
  x_selected = cbind(x_selected, x[,i])
}

data_selected = data.matrix(cbind(y1, x_selected))

treeOld<-rpart(y1~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data, method = "class", parms = list(split="information"))
treeNew<-rpart(y2~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data, method = "class", parms = list(split="information"))
rpart.plot(treeNew)
rpart.plot(treeOld)
## get the center of different groups
temp = c(as.numeric(names(table(data$y1)))) # 标签类别提取,返回结果为出现的标签
center = data.frame(matrix(0, length(temp), mNum+1))
names(center) = c("y", id_selected)
for(i in 1:length(temp)){
  center[i,1] = temp[i]
  for (j in 1:mNum){
    center[i, j+1] = mean(data_selected[,j+1])
  }
}


## 计算所有点到每个类中心的距离，返回一个矩阵，列名称为标签，行名称为点,最后一列为序号
distToCenter = matrix(0, dim(data)[1], length(temp)+1)
distToCenter[,-1] = seq(1:dim(data[1]))
for (i in 1:length(temp)){
  for(j in 1:dim(data)[1]){
    distToCenter[j, i] = sqrt(sum(center[i, 2:(mNum+1)]-data_selected[j, 2:(mNum+1)])^2)
  }
}

QuantileRank = ceiling(0.75*dim(data)[1])
NewCenter







