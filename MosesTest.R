data = c(8.2, 10.7, 7.5, 14.6, 6.3, 9.2, 11.9, 5.6, 12.8, 5.2, 4.9, 13.5,
         4.7, 6.3, 5.2, 6.8, 5.6, 4.2, 6.0, 7.4, 8.1, 6.5)
treat = c(rep(1,12), rep(2, 10))
m1=4
m2=3
k=3
MWSta <- function()

MosesTreat<-function(data, treat, m1, m2, k, 
                     alpha=c(0.001,0.005,0.01,0.025,0.05,0.1),...){
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





