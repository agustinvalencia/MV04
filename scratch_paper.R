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





#################################################################################



data <- read.table("./Data/P10-16.DAT")
data <- cov2cor(as.matrix(data))

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

e <- eigen(A)$vectors
f <- eigen(B)$vectors

a <- t(e) %*% sqrtm(solve(sigma11))
a <- t(a)
b <- t(f) %*% sqrtm(solve(sigma22))
b <- t(b)
a
b


ro1 <- sqrt(eigen(A)$values) 


ro2 <- sqrt(eigen(B)$values)

# # Inverse matrices
# a_ <- solve(a)
# b_ <- solve(b)
# 
# zeros <- rep(0, ncol(a_)*ncol(a_))
# az <- matrix(zeros, ncol = ncol(a_),nrow = ncol(a_) )
# for(i in 1:ncol(a_)) {
#     az <- az + a_[,i] %*% t(a_[,i])
# }
# az
# 
# corr_az <- sigma11 - az
# res_a <- 1/p * sum(diag(corr_az))
# quality_a <- 1 - res_a
# 
# 
# zeros <- rep(0, ncol(b_)*ncol(b_))
# bz <- matrix(zeros, ncol = ncol(b_),nrow = ncol(b_) )
# for(i in 1:ncol(b_)) {
#     bz <- bz + b_[,i] %*% t(b_[,i])
# }
# bz
# 
# corr_bz <- sigma22 - bz
# res_b <- 1/q * sum(diag(corr_bz))
# quality_b <- 1 - res_b
# quality_b













