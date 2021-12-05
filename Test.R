require(nnet)
require(stringr)
source("./kmodes++.R")

# soybean
data = read.table("./soybean-small.txt", sep=",")
data = as.matrix(data)
data = data[,1:35]
Y = data[,36]
t = c()
for(i in 1:35){if(length(table(data[,i]))==1){t = append(t, i)}}
data=subset(data, select=-t)
Y <- str_replace(Y, "D","")
YReal = Y
Y = gsub("4", "1", Y )
X = onehot(data)



