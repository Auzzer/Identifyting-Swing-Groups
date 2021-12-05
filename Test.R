data = read.table("./soybean-small.txt", sep=",")
X = data[,1:35]
Y = data[,36]
Y <- str_replace(Y, "D","")
YReal = Y
Y = gsub("4", "1", Y )
X = onehot(X)