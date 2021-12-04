WeightingDistance<-function(x, center, weighting){
  # input:
  # x为一个列向量，代表每个点
  # center为一个列向量，表示当前的类中心
  # w为一个行向量，表示权重
  # output：
  # 加权后的距离
  return(sqrt(sum(weighting%*% t(x-center)^2 )))
}


calIC<-function(X, labels){
  m = dim(X)[2]
  n = dim(X)[1]
  labels_n <- length(unique(labels))
  samples_n <- nrow(X)
  X_mean <- apply(X,2,mean)
  in_disp <- c()
  for (i in c(1:labels_n)) {
    cluster_k <- X[which(labels==i),]
    mean_k <- apply(cluster_k,2,mean)
    a1 <- nrow(cluster_k)*sum((mean_k-X_mean)^2)
    a2 <- sum((t(t(cluster_k)-mean_k))^2)
    in_disp <- c(in_disp,a2)
  }
  D<- sum(in_disp)
  k = length(unique(labels))
  AIC = D+2*m*k
  BIC = D+log(n)*m*k
  return(c(AIC, BIC))
}



calCH <- function(X,labels){
  # 根据X和标签计算CH值
  # input:
  #   X
  #   labels: 标签
  # output:
  #   (numeric):CH值
  labels_n <- length(unique(labels))
  samples_n <- nrow(X)
  X_mean <- apply(X,2,mean)
  ex_disp <- c()
  in_disp <- c()
  for (i in c(1:labels_n)) {
    cluster_k <- X[which(labels==i),]
    mean_k <- apply(cluster_k,2,mean)
    a1 <- nrow(cluster_k)*sum((mean_k-X_mean)^2)
    ex_disp <- c(ex_disp,a1)
    a2 <- sum((t(t(cluster_k)-mean_k))^2)
    in_disp <- c(in_disp,a2)
  }
  k1<- sum(ex_disp)
  k2<- sum(in_disp)
  if(k2==0)
  {
    return(1)
  }
  else
  {
    return((k1*(samples_n-labels_n))/(k2*(labels_n-1)))
  }
}

getWeightingMatrix<-function(p){
  # input:
  # X
  # Y
  # output:
  # WeightingMatrix
  # 返回一个01矩阵（2^p-1列，p行），0表示当前变量没有被选中，1表示选中当前变量
  set <- 0:(2^p-1)
  rst <- matrix(0,ncol = p,nrow = 2^p)
  for (i in 1:p){
    rst[, i] = ifelse((set-rowSums(rst*rep(c(2^((p-1):0)), each=2^p)))/(2^(p-i))>=1, 1, 0)
  }
  rst = rst[2:2^p,]
  return(rst)
}




CWW_kmeans<-function(X, Y, k, center, weighting, max_iter=1000, max_eps){
  # input:
  #   X
  #   Y
  #   k:表示当前新增的k
  #   max_iter:最多迭代次数
  #   max_eps:对类中心变化容忍度
  # return:
  #   cl(vector):输入数据的新标签
  max_iter = 1
  p = dim(X)[2]
  XY = data.frame(Y,X)
  class = c(as.numeric(names(table(Y)))) # 标签元素
  class_len = length(class) # 标签元素个数
  center = matrix(0, (k+class_len), p)
  #先把已知的类中心放在最后
  for (i in 1:class_len){
        group = XY[which(XY$Y == i), ][,2:dim(XY)[2]]
        center[k+i,] = apply(group, 2, mean)
  }
  
  # 对于每次新生成的类，采取类似于kmeans++的方法
  for (i in 1:k){
    centernow = center[(rowSums(center)) != 0, ]
    DisMatrix = matrix(0, dim(centernow)[1], dim(X)[1])
    for(j in 1:dim(centernow)[1]){DisMatrix[j, ]=apply(X, 1, WeightingDistance, weighting = weighting, center  = centernow[j, ])}
    idx = order(apply(DisMatrix, 2, sum), decreasing = FALSE)[ceiling(0.75*length(DisMatrix[i,]))]
    center[i, ] = as.numeric(X[idx, ])
  }
  
  # 下面进行聚类
  iter = 1
  
  repeat{
    eps = 0 # 误差清零
    DisMatrix = matrix(0, dim(center)[1], dim(X)[1])
    for (i in 1:dim(center)[1]){
      DisMatrix[i, ] = apply(X, 1, WeightingDistance, weighting = weighting, center = center[i, ])
    }
    cl = apply(DisMatrix, 2, which.min)
    clX = data.frame(cl,X)
    newcenter = center
    for(j in 1:dim(center)[1]){
      newcenter[j, ] = as.numeric(apply(clX[clX$cl == j,], 2, mean)[2:(dim(center)[2]+1)])
    }
    for(j in 1:dim(center)[1]){eps = eps + sum(abs(newcenter[j, ]-center[j, ]))}
    eps = eps/dim(center)[1]
    center = newcenter # 更新新重心
    iter = iter+1
    if(iter<max_iter | eps < max_eps) break
  }
  return(clX$cl)
}

CASGI<-function(X,Y){
  # input:
  #   X和Y
  # output:
  #   cl(vector):根据当前数据生成的新标签
  CASGIres = matrix(0, 2, dim(X)[1])   
  XY = data.frame(Y,X)
  p = dim(X)[2]
  WeightingMatrix = getWeightingMatrix(p)
  res_id = 1
  res = matrix(0, (p*2^(p-1)), dim(X)[1])
  for (i in 1:dim(WeightingMatrix)[1]){
    weighting = WeightingMatrix[i, ]
    kmax = sum(weighting)
    for(k in 1:kmax){
      res[res_id, ] = CWW_kmeans(X, Y, k, weighting, max_eps = 0.01)
      res_id = res_id+1
    }
  }
  CH = c(rep(0, (p*2^(p-1))))
  AIC = c(rep(0, (p*2^(p-1))))
  BIC = c(rep(0, (p*2^(p-1))))
  for(i in 1:(p*2^(p-1))){
    CH[i] = calCH(X, res[i,])
    AIC[i] = calIC(X, res[i, ])[1]
    BIC[i] = calIC(X, res[i, ])[2]
    }
  idxBest  = which.min(temp) # 表示模型最好的那个
  CASGIcl = res[idxBest, ]
  
  # 找到最好的模型所对应的权重,
  idx_weightingInf = 1
  idx_weightingSup = 0
   
  for(i in 1:dim(WeightingMatrix)[1]){
    idx=i
    idx_weightingSup = idx_weightingInf +factorial(sum(WeightingMatrix[i,]))-1
    if(idxBest>idx_weightingInf & idxBest<idx_weightingSup) break
    idx_weightingInf = idx_weightingSup + 1
  }
  CASGIrweighting = WeightingMatrix[idx,]
  return(list(CASGIcl, CASGIrweighting))
}




randomCASGI<-function(X,Y){
  # 采取抽样的方式减少计算量
  
  pall = dim(X)[2]
  pickNum = ceiling(pall/3)
  sam = getWeightingMatrix(pall)
  sam[rowSums(sam)==pickNum,]->sam
  
  
  cl = matrix(0,  dim(sam)[1],dim(X)[1])
  WeightingList = matrix(0,  dim(sam)[1], dim(X)[2])
  # 对每个小组只获得最优的那个
  for (i in 1:dim(sam)[1]){
    temp = sam[i,]
    Xselect = X[,which(temp != 0)]
    XselectId = c(which(temp!=0))
    CASGIres = CASGI(Xselect, Y)
    cl[i, ]  = CASGIres[[1]]
    itemid = 1
    for(item in XselectId){
      WeightingList[i, item] = CASGIres[[2]][itemid]
      itemid = itemid+1
    }
    
  }
  return(list(cl, WeightingList))
}


X = as.matrix(data[,3:12])
test<-CASGI(X, Y)
CHall = c(rep(0, dim(test[[1]])[1]))
AICall = c(rep(0, dim(test[[1]])[1]))
BIC = c(rep(0, dim(test[[1]])[1]))


for (i in 1:length(CHall)){
  CHall[i] = calCH(X, c(test[[1]][i,]))
}
which.min(CHall)
clnew=test[[1]][which.min(CHall), ] # 新的标签
cleSelected=test[[2]][which.min(CHall), ] # 被选择的变量
num = choose(   (dim(X)[2]-1) , (ceiling (dim(X)[2]/4) -1) ) # 每个变量应该被抽到的次数

