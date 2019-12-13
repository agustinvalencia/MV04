library(expm)
library(knitr)

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

#get correlation
R11 <- cor(sigma11)
R22 <- cor(sigma22)
R12 <- cor(sigma12)
R21 <- cor(sigma21)


e <- sqrtm(sigma11) %*% sigma12 %*% solve(sigma22) %*% sigma21 %*% sqrtm(sigma11)

f <- sqrtm(sigma22) %*% sigma21 %*% solve(sigma11) %*% sigma12 %*% sqrtm(sigma22)


ro1 <- sqrt(eigen(e)$values)
ro2 <- sqrt(eigen(f)$values)


#Hypothesis testing
alpha <- 0.05
#critical value
crit <- qchisq(p = (1-alpha), df = p*q)
#test statistic
test1 <- -(n-1-(0.5*(p+q+1))) * log(prod(1-(ro1)^2))
