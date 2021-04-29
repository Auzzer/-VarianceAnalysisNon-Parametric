data = c( 23.1, 57.6, 10.5, 23.6, 11.9, 54.6, 21.0, 20.3,
          22.7, 53.2,  9.7, 19.6, 13.8, 47.1, 13.6, 23.6,
          22.5, 53.7, 10.8, 21.1, 13.7, 39.2, 13.7, 16.3,
          22.6, 53.1,  8.3, 21.6, 13.3, 37.0, 14.8, 14.8)
treat = c(rep(1,8),rep(2,8),rep(3,8),rep(4,8))
block = c(rep(c(1:8),4))
x = cbind.data.frame(data,treat,block)
b = length(unique(block)) # 区组数目
k = length(unique(treat)) # 处理数目
plot = 0 # initial plot
method = "mean"
R = matrix(0,k,b)
for (i in 1:b){
  R[,i]=rank(subset(x, block == i)$data)
}
Rj = c(0, k)
for (i in 1:k){
  Rj[i] = sum(R[i,])
}
Aligen = c(0,b)
if (method == "mean"){
  
  for (i in 1:b){
    Aligen[i] = mean(subset(x, block == i)$data)
  }
}if(method == "median"){
  for (i in 1:b){
    Aligen[i] = median(subset(x, block == i)$data)
  }
}else{print("Please choose the right way to aligen observation")}
AligenedObservation <- c()
for (i in 1:b){
  tmp <- subset(x, block == i)$data-Aligen[i]
  AligenedObservation <- cbind(AligenedObservation, tmp)
}
rm(tmp)

R_adj <- matrix(rank(AligenedObservation),k,b,byrow = FALSE)

Rj_adj <- c(0,k)
Ri_adj <- c(0,b)
for (i in 1:k){
  Ri_adj[i] = sum(R_adj[i,])
}


for (i in 1:b){
  Rj_adj[i] = sum(R_adj[,i])
}

for (i in 1:k){
  if (max(R_adj[i, ]) != 1){
    plot = 1
  }
}


if (method == "mean"){
  for (i in 1:b){
    Aligen[i] = mean(subset(x, block == i)$data)
  }
}else if(method == "median"){
  for (i in 1:b){
    Aligen[i] = median(subset(x, block == i)$data)
  }
}else{
  print("Please choose the right way to aligen observation")
}
  
