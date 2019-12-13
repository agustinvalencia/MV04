library(expm)
library(knitr)
library(CCP)

data <- read.table("./Data/P10-16.DAT")

#number of observations (patients)
n <- 46
#number of primary variables
p <- 3
#number of secondary variables
q <- 2

#separating the variance-covariance matrix
sigma11 <- as.matrix(data[1:3,1:3])
sigma22 <- as.matrix(data[4:5, 4:5])
sigma21 <- as.matrix(data[4:5, 1:3])
sigma12 <- as.matrix(data[1:3, 4:5])


A <- sqrtm(solve(sigma11)) %*% sigma12 %*% solve(sigma22) %*% sigma21 %*% sqrtm(solve(sigma11))
B <- sqrtm(solve(sigma22)) %*% sigma21 %*% solve(sigma11) %*% sigma12 %*% sqrtm(solve(sigma22))

e <- eigen(A)$vectors[,1:2]
f <- eigen(B)$vectors

a <- t(e) %*% sqrtm(solve(sigma11))
a <- t(a)
b <- t(f) %*% sqrtm(solve(sigma22))
b <- t(b)
a
b


ro1 <- sqrt(eigen(A)$values) 
ro1[3] <- 0

ro2 <- sqrt(eigen(B)$values)


#Hypothesis testing
alpha <- 0.05
#critical value
crit <- qchisq(p = (1-alpha), df = p*q)
#test statistic
test1 <- -(n-1-(0.5*(p+q+1))) * log(prod(1-ro1^2))  #13.74948
test2 <- -(n-1-(0.5*(p+q+1))) * log(prod(1-ro2^2))  #13.74948 they are the same

#test1 > critical value >>>>> We reject H0

