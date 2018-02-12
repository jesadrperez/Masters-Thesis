## EFFECTS OF DIFFERENT VARCOVARIANCE STRUCTURES ON GAMMA ##
# This program will calculate the gamma of 5 varcovariance structures to see if the structure has any effect on it.

# mnormt contains the function for generating random multivariate values
library('mnormt')

# VarCovariance Structures to be Tested
# 3x3 VarCovariance Structures will be used
covA <- matrix(c(c(1,0,0), c(0,1,0), c(0,0,1)), ncol=3, nrow=3)
covB <- matrix(c(c(1,0.5,0.25), c(0.5,1,0.5), c(0.25,0.5,1)), ncol=3, nrow=3)
covC <- matrix(c(c(2,0,0), c(0,4,0), c(0,0,8)), ncol=3, nrow=3)
covD <- matrix(c(c(2,0.5,0.25), c(0.5,4,0.5), c(0.25,0.5,8)), ncol=3, nrow=3)
covE <- matrix(c(c(10,0,0), c(0,1,0), c(0,0,1)), ncol=3, nrow=3)

# Creates a list of test varcov matrices
VarCovs <- list(covA, covB, covC, covD, covE)
label <- c('covA', 'covB', 'covC', 'covD', 'covE')

# Means, constant for all matrices
mean = rep(0,3)

# Hawkins' Test 
# Returns Z1 and Z2
# Gamma = abs(Z1.1 - Z1.2)
Hawkins <- function(grp1, grp2) {
	p <- length(grp1[1,]);
	g <- 2
	n <- length(grp1[,1]);
	m <- length(grp2[,1]);
	N <- n + m
	v <- N - g - 1;
	
	# Calculates sample covariance matrices for each grp and pooled sample 
	# covariance matix. 
	cov.1 <- cov(grp1); cov.2 <- cov(grp2);
	cov.pool <- (1/(N-g))*((n-1)*cov.1+(m-1)*cov.2)
	
	# Calculates mahalanobis distance. This is Vij in the paper.	
	maha.1 <- mahalanobis(grp1, colMeans(grp1), cov.pool)
	maha.2 <- mahalanobis(grp2, colMeans(grp2), cov.pool)
	
	# This is eq 2  
	hotel.1 <- ((N-g-p)*n*maha.1)/(p*((n-1)*(N-g)-n*maha.1))
	hotel.2 <- ((N-g-p)*m*maha.2)/(p*((m-1)*(N-g)-m*maha.2))
	
	# Gives Pr[F>hotel] with p and v-p+1 df	
	a.1 <- pf(hotel.1, p, v-p+1, lower.tail=F)
	a.2 <- pf(hotel.2, p, v-p+1, lower.tail=F)
	#a.T <- c(a.1, a.2)
	
	# return(ad.test(a.T, punif)$p.value)
	
	# This is that eq w/ the Legendre Polynomials
	# (I got the formula from another book)	
	Z1.1 <- -sqrt(3/n)*sum(2*a.1-1)
	Z1.2 <- -sqrt(3/m)*sum(2*a.2-1)
	#Z1.T <- -sqrt(3/(n+m))*sum(2*a.T-1)

	# diff <- abs(Z1.1 - Z1.2)
	# p.value <-  2*pnorm(diff, mean = 0, sd = 1.444983, lower.tail = FALSE)

	# return(diff)
	return(c(Z1.1, Z1.2))
	# print(Z1.1)
	# print(Z1.2)
	# print(Z1.T)	
}	

# Number of simulations to run for each matrix
sim <- 10000
Var.Data <- NULL
system.time(
for (i in 1:length(VarCovs)) {
	VarCov <- VarCovs[[i]];
	Var.Temp <- data.frame(NULL)
	print(cat("Working on", label[i],"\t",sep=" "));
	
	for (j in 1:sim) {
	 	grp1 <- rmnorm(30, mean, VarCov);
		grp2 <- rmnorm(30, mean, VarCov);	
		zs <- Hawkins(grp1,grp2);
		z.1 <- zs[1];
		z.2 <- zs[2];
		Var.Temp[j,1] <- label[i];
		Var.Temp[j,2] <- z.1;
		Var.Temp[j,3] <- z.2;
		if (j%%1000 == 0) {
			print(cat("Finished:", format(j), "out of", format(sim), "\t" ,sep=" "));
		}
	}
	Var.Data <- rbind(Var.Data, Var.Temp);
}		
)


