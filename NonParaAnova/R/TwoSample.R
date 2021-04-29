# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

# This part is aimed for inferring the location and scale of the two samples
#it includes: Brown-Mood test  Wilcoxon-Mann-Whitney test which infers the mean of two samples
#and Mood test & Moses test which infers the variance of two samples.
# For more details of the tests below, maybe you can check the referrence in the end of doc.


# Brown-mood test
BM.test <- function(data,treat,alt,...){
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
  result<-structure(list(ContingencyTable=Contingency_table, p.value=p))
  return(result)
}

# Wilcoxon-Mann-Whitney test
# this part makes use of the wilcox.test() which is already inside R,
#while to unified the input of data,we make some adaptions.

WMW.test <- function(data,treat){
  xy = cbind.data.frame(data,treat)
  x = xy[which(xy$treat == 1), 1]
  y = xy[which(xy$treat == 2), 1]
  result <- wilcox.test(x,y)
  return(result)
}

# Mood test
# this part makes use of the mood.test() which is already inside R,
# while to unified the input of data,we make some adaptions.

md.test <- function(data,treat,...){
  xy = cbind.data.frame(data,treat)
  x = xy[which(xy$treat == 1), 1]
  y = xy[which(xy$treat == 2), 1]
  result <- mood.test(x,y)
  return(result)
}


# This function differs from the MoseTest in R package "DescTools" which tests the extreme, while this tests the
# variance of two samples. For more infor, please check the the reference[4]
MoseRank.test <- function(data, treat, m1, m2, k,...){
  xy = cbind.data.frame(data, treat)
  x = xy[which(xy$treat == 1), 1]
  y = xy[which(xy$treat == 2), 1]
  n1 = k*m1
  n2 = k*m2

  Sample1 = sample(x, n1)
  Sample2 = sample(y, n2)
  SSA = matrix(0, m1, k+2)
  SSB = matrix(0, m2, k+2)
  SSA = matrix(sample(x, n1))
  SSB = matrix(sample(x, n1))
  # calculate SSR of each part
  for (i in 1:m1){
    SSA[i, k+1] = sum(SSA[i, 1:k]-sum(SSA[i, 1:k]))
  }
  for (i in 1:m2){
    SSA[i, k+1] = sum(SSB[i, 1:k]-sum(SSB[i, 1:k]))
  }
  SSR = cbind(SSA[, k+1], SSB[, k+1])
  R = rank(SSR)
  SSA[,k+2] = R[1:m1]
  SSB[,K+2] = R[m1+1:(m1+m2)]



}


MoseRank.test <- function(data,treat,...){
  xy = cbind.data.frame(data,treat)
  x = xy[which(xy$treat == 1), 1]
  y = xy[which(xy$treat == 2), 1]



}


# Reference:
# [1]Brown-Mood Test:
# [2]
# [3]
# [4]Mose Test: Lincoln E. Moses. "Rank Tests of Dispersion." Ann. Math. Statist. 34 (3) 973 - 983, September, 1963. https://doi.org/10.1214/aoms/1177704020


