#############################################################################
# Box's M 

# k is the number of covar matrices testing
# p is the number of groups (populations)
# n is the sample size of one group (must be equal sample sizes)
# N is the sample size of both groups (2*n)

# Box's M Function
Box.M <- function(grp1, grp2){
	k <- 2; p <- length(grp1[1,]); n <- length(grp1[,1]); N <- n*k;
	
	cov.1 <- cov(grp1); cov.2 <- cov(grp2);
	cov.pool <- (1/(2*n-2))*((n-1)*cov.1+(n-1)*cov.2)
	
	df <-  0.5*(k-1)*p*(p+1);
	C <- ((2*(p^2)+3*p-1)*(k+1))/(6*(p+1)*k*(n-1)); 
	M <- 0.5*((n-1)*log(det(cov.1))+(n-1)*log(det(cov.2)))-0.5*(2*n-2)*log(det(cov.pool))
	return(pchisq(-2*(1-C)*M, df, lower.tail=F));
}  

# Example from using Beall 1945 Dataset
# y1 - 
# y2 -
# y3 - 
# y4 - 
# y5 - Sex of individual (Grouping Variable)

Beall <- read.table('~/Desktop/Beall 1945.txt', header =T)
Beall
Males <- Beall[Beall$y5 == 'M',1:4]
Females <- Beall[Beall$y5 == 'F',1:4]
Box.M(Males,Females)

# > Beall <- read.table('~/Desktop/Beall 1945.txt', header =T)
# > Beall
#    y1 y2 y3 y4 y5
# 1  15 17 24 14  M
# 2  17 15 32 26  M
# 3  15 14 29 23  M
# 4  13 12 10 16  M
# 5  20 17 26 28  M
# 6  15 21 26 21  M
# 7  15 13 26 22  M
# 8  13  5 22 22  M
# 9  14  7 30 17  M
# 10 17 15 30 27  M
# 11 17 17 26 20  M
# 12 17 20 28 24  M
# 13 15 15 29 24  M
# 14 18 19 32 28  M
# 15 18 18 31 27  M
# 16 15 14 26 21  M
# 17 18 17 33 26  M
# 18 10 14 19 17  M
# 19 18 21 30 29  M
# 20 18 21 34 26  M
# 21 13 17 30 24  M
# 22 16 16 16 16  M
# 23 11 15 25 23  M
# 24 16 13 26 16  M
# 25 16 13 23 21  M
# 26 18 18 34 24  M
# 27 16 15 28 27  M
# 28 15 16 29 24  M
# 29 18 19 32 23  M
# 30 18 16 33 23  M
# 31 17 20 21 21  M
# 32 19 19 30 28  M
# 33 13 14 12 21  F
# 34 14 12 14 26  F
# 35 12 19 21 21  F
# 36 12 13 10 16  F
# 37 11 20 16 16  F
# 38 12  9 14 18  F
# 39 10 13 18 24  F
# 40 10  8 13 23  F
# 41 12 20 19 23  F
# 42 11 10 11 27  F
# 43 12 18 25 25  F
# 44 14 18 13 26  F
# 45 14 10 25 28  F
# 46 13 16  8 14  F
# 47 14  8 13 25  F
# 48 13 16 23 28  F
# 49 16 21 26 26  F
# 50 14 17 14 14  F
# 51 16 16 15 23  F
# 52 13 16 23 24  F
# 53  2  6 16 21  F
# 54 14 16 22 26  F
# 55 17 17 22 28  F
# 56 16 13 16 14  F
# 57 15 14 20 26  F
# 58 12 10 12  9  F
# 59 14 17 24 23  F
# 60 13 15 18 20  F
# 61 11 16 18 28  F
# 62  7  7 19 18  F
# 63 12 15  7 28  F
# 64  6  5  6 13  F
# > Males <- Beall[Beall$y5 == 'M',1:4]
# > Females <- Beall[Beall$y5 == 'F',1:4]
# > Box.M(Males,Females)
# [1] "Returns p-value"
# [1] 0.1944866








