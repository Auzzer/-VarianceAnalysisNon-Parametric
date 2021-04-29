
###结的计算
a = c(1,1,1,13,3,3,3,3,2,2,2,2,2,5,2,2)
rank(a)
jie = c(0,length(a))
for (i in 1:length(a)){
  for (j in 1:(length(a)-i))
    if(a[i] == a[i+j]){jie[i] = jie[i]+1}
}

#########################Brown-Mood median test########################################################

BrownMoodMedian<-function(data,treat,alpha, ...){
  xy = cbind.data.frame(data,treat)
  x = xy[which(xy$treat == 1), 1]
  y = xy[which(xy$treat == 2), 1]
  xy_mid = median(xy[,1])
  A = sum(x>xy_mid)
  B = sum(y>xy_mid)
  C = sum(x<xy_mid)
  D = sum(y<xy_mid)
  t = sum(A,B)
  m = sum(A,C)
  n = sum(B,D)
  if(alt == "Small"){
    p <- phyper(A,m,n,t)
  }else if(alt == "Large"){
    z = (A-m*t/(m+n)) / sqrt(m*n*t*(m+n-t)/((m+n)^3))
    p = pnorm(z)
  }else{print("No defined parameter, Please check the note")}
  Contingency_table <- matrix(c(A,B,t,C,D,(m+n)-(A+B),m,n,m+n),
                              3,3)
  col_name <- c("X","Y","X+Y")
  row_name <- c(">M_XY","<M_XY","Total")
  dimnames(Contingency_table) = list(col_name,row_name)
  result<-list(ContingencyTable=Contingency_table, p.value=p)
  return(result)
}

####################mood方差检验############################################################
#：调用mood.test()函数，统一数据接口：
md.test <- function(data,treat,...){
  xy = cbind.data.frame(data,treat)
  x = xy[which(xy$treat == 1), 1]
  y = xy[which(xy$treat == 2), 1]
  result <- mood.test(x,y)
  return(result)
} 





####################Moses 方差检验#############################################################

data = c(8.2, 10.7, 7.5, 14.6, 6.3, 9.2, 11.9, 5.6, 12.8, 5.2, 4.9, 13.5,
         4.7, 6.3, 5.2, 6.8, 5.6, 4.2, 6.0, 7.4, 8.1, 6.5)
treat = c(rep(1,12), rep(2, 10))
m1=4
m2=3
k=3
MosesVar<-function(data, treat, m1, m2, k,...){
  xy = cbind.data.frame(data,treat)
  x = xy[xy$treat == 1, 1]
  y = xy[xy$treat == 2, 1]
  n1=m1*k
  n2=m2*k
  Sample1 = matrix(sample(x,n1), nrow = m1)
  Sample2 = matrix(sample(y,n2), nrow = m2)
  SSA=c(rep(0,m1))
  SSB=c(rep(0,m2))
  for (i in 1:m1){SSA[i] = sum((Sample1[i,]-mean(Sample1[i,]))^2)}
  for (i in 1:m2){SSB[i] = sum((Sample2[i,]-mean(Sample2[i,]))^2)}
  SSR = c(SSA,SSB)
  R = rank(SSR)
  R1 = R[1:m1]
  R2 = R[(m1+1):(m1+m2)]
  S = min(sum(R1),sum(R2))
  # 有一个比较的步骤，看S是来源于哪里
  if(sum(R1)<sum(R2)){m=m1}else{m=m2}
  T_M = S-m*(m+1)/2
}







################################################################################
################################################################################
##################Part 2: 多组数据位置推断##############################################################

# 在使用Kruskal-wallis检验后，可以使用dunn检验进行两两之间的比较
dunn.test<-function(data,treat,alpha,...){
  #首先是准备工作，分别需要定义处理内平均秩；结点出现时对MST、SE修正
  x = data
  n = length(x) ##数据总个数
  k = length(unique(treat)) ##类别数
  # 若abs(d_ij) >= Z_{1-alphastar}，则表示第i和第j处理见有显著差异
  alphastar = alpha/(choose(k,2)*2)
  Z = qnorm(1-alphastar)
  
  #接下来计算秩均值
  # 每一次计算都要初始化
  rank<- NULL
  m<- NULL
  rank_mean <- NULL
  for (i in 1:k){
    rank[i]=sum(subset(rank(x),treat==i))##处理内秩和
    m[i]=sum(treat==i)##类内个数
    rank_mean[i] <- rank[i]/m[i]
  }
  # 当数据中存在结点时
  t = as.vector(table(rank(x))) ##计算结点长度
  #当结点存在时，需要对MST,SE进行修正,首先定义一个修正项:correct
  correct = sum(t^3-t)/(12*(n-1))
  mst = n*(n+1)/12-correct
  
  
  #进行比较
  rowname <- NULL
  diff_ij <- NULL
  se_ij <- NULL
  d_ij <-NULL
  p <- NULL
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      diff <- abs(rank_mean[i]-rank_mean[j])
      se <- sqrt(mst*(1/m[i]+1/m[j])) # se = sqrt(MST(1/n_i+1/n_j))
      d = diff/se
      #将结果保存
      diff_ij <- cbind(diff_ij,diff)
      se_ij <- cbind(se_ij,se)
      d_ij <- cbind(d_ij,d)
      p <- cbind(p,1-pnorm(d))
      rowname <- cbind(rowname,paste(i,"VS", j))
      
    }
  }
  
  result_table <- matrix(c(diff_ij,se_ij,d_ij,p),choose(k,2),4,byrow = FALSE)
  rownames(result_table) <- rowname
  colnames(result_table) <- c("|\bar{R_.i}-\bar{R.j}|","SE","d_ij","p")
  result <- list(Table = result_table, alphastar = alphastar,Z=Z)
  return(result)
  
}

########################### Jonckheere–Terpstra检验######################################################

JT.test<-function(data, treat, ...){
  x = data
  N = length(x) ##查看样本容量
  k = length(unique(treat)) ## 计算组数
  #计算每个组的样本容量
  n<-NULL
  for (i in 1:k){
    n[i] = length(subset(x,treat==i))
  }
  jie = 0 ##先初始化结，当jie = 0时，表示没有结存在，在这里先定义为没有
  # 计算秩
  w_ij = 0
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      xi = subset(x, treat==i)
      xj = subset(x, treat==j)
      for (a in 1:n[i]){
        for (b in 1:n[j]){
          if (xi[a]<xj[b]){
            w_ij = w_ij+1
          }else if(xi[a]==xj[b]){ #处理结存在时的情况
            w_ij = w_ij+1/2
            jie = 1 # 更新结的状态
          }
        }
      }
    }
  }
  J <- w_ij
  # 计算J的期望方差
  ## 期望
  EJ = (N^2 - sum(n^2))/4
  ## 方差
  if (jie==0){
    VarJ = (N^2(2*N+3)-sum(2*(n^3)+3*(n^2)))/72
  }else{
    t = as.vector(table(rank(x)))
    VarJ = 
      ( (N*(N-1)*(2*N+5))  - sum(n*(n-1)*(2*n+5)) - sum(t*(t-1)*(2*t+5)) )/72#第一项 
    +
      1 / (36*N*(N-1)*(N-2)) * (sum(n*(n-1)*(n-2))) * (sum(t*(t-1)*(t-2)))#第二项
    +
      1 / (8*N*(N-1)) * (sum(n*(n-1))) * (sum(t*(t-1)))
  }
  ## 依据J的近似分布计算p值
  
  Z = (J-EJ)/sqrt(VarJ)
  p <- 1-pnorm(Z)
  result <- list(J.Value=J, EJ=EJ, SD=sqrt(VarJ), Z.value=Z, p.value=p)
  return(result)
}


###########################Friedman.test()######################################################

PekingFish <-c(85,82,82,79,
               87,75,86,82,
               90,81,80,76,
               80,75,81,75)
treat.PekingFish <- c(rep(1:4, time=4))
block.PekingFish <- c(rep(1:4, each=4))
friedman.test(PekingFish, treat.PekingFish, block.PekingFish)
# 在Friedman检验之后，如果结果存在差异时，可以使用Hollander-wolfe 两两比较



##################################Hollander-Wolfe 两处理间比较###############################################

HollanderWolfe.test<- function(data, treat, block,...){
  x = cbind.data.frame(data,treat,block)
  # 计算block数量
  b = length(unique(block)) # 对应书里的b
  # 计算treat数量
  k = length(unique(treat))
  # 初始化秩矩阵
  R = matrix(0,b,k)
  for (i in 1:k){
    R[i,]=rank(subset(x, treat ==i)$data)
  }
  Rj <- NULL
  for (i in 1:b){
    Rj<- cbind(Rj, sum(R[,i]))
  }
    
  # 判定是否存在结
  tau = c(0,n)
  library(plyr)
  for(i in 1:n){tau[i]=max(count(R[i,])$freq)}#tau 表示结 
  
  if (sum(tau)==length(tau)){# sum(tau)=length(tau) 代表没有结
    print("no plot")
    SE = sqrt(k*(k+1)/6)
  }else{
    print("plots exit")
    SE = sqrt((k*(k+1)/6)-n*(sum(tau^3-tau))/(6*(k-1)))
  }
  
  diff_ij = c()
  rowname <- NULL
  for (i in 1:(l-1)){
    for ( j in 1:(l-i)){
      diff = abs(Rj[i]-Rj[(i+j)])
      # print("asd")
      diff_ij <- cbind(diff_ij,diff)
      rowname <- cbind(rowname,paste(i,"VS", (i+j)))
    }
  }
  
  D_ij <- NULL
  for (i in 1:length(diff_ij)) {
    D = diff_ij[i]/SE
    D_ij <- cbind(D_ij,D)
  }
  dimNum = choose(length(Rj),2)
  se = rep(SE, dimNum)
  alpha.star = 0.05/(choose(length(Rj),2))
  z = qnorm(1-alpha.star)
  Z = rep(z, choose(length(Rj),2))
  result_table <- matrix(c(diff_ij,se,D_ij,Z),choose(length(Rj),2),4,byrow = FALSE)
  rownames(result_table) <- rowname
  colnames(result_table) <- c("|\bar{R_.i}-\bar{R.j}|","SE","D_ij","Z")
  result <- list(Table = result_table, alphastar = alpha.star,Z=Z)
  print("if D_ij is larger than Z, it shows that these two have connections.")
return(result)
  
  
  
}



#################调整秩和检验，亦称为Hodges-Lehmmann 检验（HL检验）################################################################

data = c(23.1, 57.6, 10.5, 23.6, 11.9, 54.6, 21.0, 20.3,
         22.7, 53.2,  9.7, 19.6, 13.8, 47.1, 13.6, 23.6,
         22.5, 53.7, 10.8, 21.1, 13.7, 39.2, 13.7, 16.3,
         22.6, 53.1,  8.3, 21.6, 13.3, 37.0, 14.8, 14.8)
treat = c(rep(1,8),rep(2,8),rep(3,8),rep(4,8))
block = c(rep(c(1:8),4))
HL.test <- function(data,treat,block,method,...){
  #H0:theat_1=theta_2=...=theta_k 
  #H1: exits i,j in 1,2,..k, s.t. theta_i neq theta_j
  x = cbind.data.frame(data,treat,block)
  b = length(unique(block)) # 区组数目
  k = length(unique(treat)) # 处理数目
  jie = 0 # initial plot
  
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
    
  
  if (jie == 0){
    Q <- (k-1)(sum(Rj_adj)-k*b*b*(k*b+1)*(k*b+1)/4) /
      (k*k*(k*b+1)*(2*k*b+1)/6 - sum(Ri_adj)/k)
  }else{
    Q <- (k-1)*(sum(Rj_adj)-k*b*b*(k*b+1)*(k*b+1)/4) / 
      (sum(R_adj)-sum(Ri_adj)/k)
  }
  
  p = 1-pchisq(Q, df=(k-1))
  result <- list(Q.Stats = Q, P.value=p)
return(result)
  
}





############################# Cochran检验####################################################


cochran.test<-function(data,treat,block,...){
  b = length(unique(block)) #区数
  k = length(unique(treat)) #组数
  y = matrix(NA,k,b)
  for (i in 1:k){
    for (j in 1:b){
      y[i,j]=subset(data,(treat==i)&(block==j))
    }
  }
  nj = apply(y, 1, sum)
  ni = apply(y,2,sum)
  # 计算Q统计量
  Q = ((k-1)*(sum(nj^2)-(sum(nj)*sum(nj))/k))  /  (sum(ni)-sum(ni^2)/k)
  p=1-pchisq(Q,k-1)
  result <- list(Q.value=Q,P.value=p)
  return(result)
}




################Durbin 检验##########################################

data = c(73, 74, NA, 71, 
         NA, 75, 67, 72, 
         74, 75, 68, NA,
         75, NA, 72, 75)
block = c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
treat = c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
Durbin.test <- function(data,treat, block,...){
  x = cbind.data.frame(data,treat, block)
  k = length(unique(treat))
  b = length(unique(block))
  
  R = matrix(0, k, b)
  for (i in 1:b){
    R[, i] <- rank(subset(x, block ==i)$data)
  }
  
  Rj = c(0, k)
  for (j in 1:k){
    Rj[j] <- sum(R[j,])
  }
  Rj<- Rj-b
  ##in the complete group design, t=k,r=b
  t = length(na.omit(subset(x, treat == 1)$data))
  r = length(na.omit(subset(x, block == 1)$data))
  #计算D统计量
  D = 12*(k-1)*sum(Rj^2)/(r*k*(t^2-1)) - 3*r*(k-1)*(t+1)/(t-1)
  
  p = 1-pchisq(D, df=(k-1))
  result <- list(D.Stats = D, P.value=p)
return(result)
}





























