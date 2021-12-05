require(nnet)
require(stringr)

CH<-function(cl, X){
  p = dim(X)[2]
  
  clNum = length(names(table(cl))) # number of clusterings
  container = matrix(0, clNum, 2) 
  # for each item in list, corresponding MSA and MSE to cl is included, eg list[[1]] includes MSA and MSE of cl 1  
  ## calulate center of each clustering
  center = matrix(0, clNum, p)
  class = c(as.numeric(names(table(cl))))
  i=1
  for(item in class){
    idx = which(cl == item)
    temp = as.matrix(X[idx ,])
    center[i, ] = as.vector(apply(temp, 2, getmode))
    i=i+1
  }
  
  ## calculate the SSA & SSE of each clustering and then save them in the matrix where the first col is SSA and the second col is SSE
  i=1
  for(item in class){
    idx = which(cl == item)
    temp = as.matrix(X[idx ,]) 
    
    SSA =  sum(apply(center, 1, itemDis, center[i, ]))/(clNum-1)# sum of squares between groups
    SSE = sum(apply(temp, 1, itemDis, center[i, ]))/ dim(temp)[1]# sum of squares within groups
    container[i,]  = c(SSA, SSE)
    i=i+1
  }
  CH = ( sum(container[,1])/ sum(container[,2]) ) - log(clNum)
  return(CH)
}

IC<-function(cl, X){
  p = dim(X)[2]
  
  clNum = length(names(table(cl))) # number of clusterings
  container = matrix(0, clNum, 2) 
  # for each item in list, corresponding MSA and MSE to cl is included, eg list[[1]] includes MSA and MSE of cl 1  
  ## calulate center of each clustering
  center = matrix(0, clNum, p)
  class = c(as.numeric(names(table(cl))))
  
  i=1
  for(item in class){
    idx = which(cl == item)
    temp = as.matrix(X[idx ,])
    center[i, ] = as.vector(apply(temp, 2, getmode))
    i=i+1
  }
  
  ## calculate the SSA & SSE of each clustering and then save them in the matrix where the first col is SSA and the second col is SSE
  i=1
  for(item in class){
    idx = which(cl == item)
    temp = as.matrix(X[idx ,]) 
    
    SSA =  sum(apply(center, 1, itemDis, center[i, ]))/(clNum-1)# sum of squares between groups
    SSE = sum(apply(temp, 1, itemDis, center[i, ]))/ dim(temp)[1]# sum of squares within groups
    container[i,]  = c(SSA, SSE)
    i=i+1
  }
  ic = 2*p - ( sum(container[,1])/ sum(container[,2]) )
  return(ic)
}

itemDis<-function(a, b){
  return(sum(abs(a-b)))
}

onehot<-function(x){
  p = dim(x)[2]
  res=class.ind(x[,1])[,1:dim(class.ind(x[,1]))[2]-1]
  for(i in 2:p){
    res = cbind(res, class.ind(x[,i])[,1:dim(class.ind(x[,i]))[2]-1])
  }
  return(res)
}


getmode<-function(x){
  # input 
  #   x: vector
  # return
  #   mode:numeric
  uniqx = unique(x)
  return(uniqx[which.max(tabulate(match(x, uniqx)))])
}

getCenter<-function(x){
  # input
  #   x: matrix n*p
  # return
  #   center:vector length=p
  return(apply(x, 2, getmode))
}

kmodespp<-function(X, Y, kNew, max_iter=1000){
  # input: 
  #   X(matrix): observations
  #   Y(vector): origin categories
  #   max_iter(numeric): max iteration number
  #   kNew(numeric): number of merging new categories
  # return:
  #   cl(vecter): final clustering results
  # max_iter=1
  p = dim(X)[2]
  XY = data.frame(Y,X)
  class = c(as.numeric(names(table(Y)))) # label elements
  class_len = length(class) # number of label elements
  center = matrix(0, (kNew+class_len), p)
  #put the known categories at the very final lines of matrix
  for (i in 1:class_len){
    group = XY[which(XY$Y == i), ][,2:dim(XY)[2]]
    center[kNew+i,] = getCenter(group)
  }
  
  # for each new merging label, adopt the methods like kmeans++
  for (i in 1:kNew){
    centernow = center[(rowSums(center)) != 0, ] # centernow:  center of known clusterings
    DisMatrix = matrix(0, dim(centernow)[1], dim(X)[1])
    for(j in 1:dim(centernow)[1]){DisMatrix[j, ]=apply(X, 1, itemDis, b  = centernow[j, ])}
    idx = order(apply(DisMatrix, 2, sum), decreasing = FALSE)[ceiling(0.75*length(DisMatrix[i,]))]
    center[i, ] = as.numeric(X[idx, ])
  }
  
  # process of kmodes++
  iter = 1
  
  repeat{
    eps = 0 #  clear eps
    DisMatrix = matrix(0, dim(center)[1], dim(X)[1])
    for (i in 1:dim(center)[1]){
      DisMatrix[i, ] = apply(X, 1, itemDis, b  = center[i, ])
    }
    cl = apply(DisMatrix, 2, which.min)
    clX = data.frame(cl,X)
    newcenter = center # renew newcenter to fit the dim of center
    for(j in 1:dim(center)[1]){
      newcenter[j, ] = as.numeric(apply(clX[clX$cl == j,], 2, mean)[2:(dim(center)[2]+1)])
    }
    for(j in 1:dim(center)[1]){eps = eps + sum(abs(newcenter[j, ]-center[j, ]))} 
    eps = eps/dim(center)[1]
    center = newcenter # renew center
    iter = iter+1
    if(iter<max_iter | eps >0) break
  }
  return(clX$cl)
  
}


chooseK<-function(X, Y, maxiter = 1000){
  p = dim(X)[2] # num of variabels
  NumY = length(unique(Y))
  n = dim(X)[1]
  Res = matrix(0, (p-NumY), n)
  for(kNew in 1:(p-NumY)){
    Res[kNew, ]=kmodespp(X,Y, kNew = kNew)
    for(i in 1:n){
      clnow = Res[kNew, j]
      for(j in 1:kNew+NumY){}
    }
  }
}








