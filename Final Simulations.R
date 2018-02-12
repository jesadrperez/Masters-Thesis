####### Final Simulations #########

# Loads required packages
require('mnormt')
require('adk')
require('Hmisc')

# Loads identity matrix and mean vector
cov <- matrix(c(c(1,0,0), c(0,1,0), c(0,0,1)), ncol=3, nrow=3)
mean = rep(0,3)

# Loads Statistical Test 
# non-parametric test
require(adk)
Adk <- function(grp1,grp2){	
	p <- length(grp1[1,]);
	g <- 2
	n <- length(grp1[,1]);
	m <- length(grp2[,1]);
	N <- n + m
	v <- N - g - 1;
	
	cov.1 <- cov(grp1); cov.2 <- cov(grp2);
	cov.pool <- (1/(N-g))*((n-1)*cov.1+(m-1)*cov.2)
		
	maha.1 <- mahalanobis(grp1, colMeans(grp1), cov.pool)
	maha.2 <- mahalanobis(grp2, colMeans(grp2), cov.pool)
	
	hotel.1 <- ((N-g-p)*n*maha.1)/(p*((n-1)*(N-g)-n*maha.1))
	hotel.2 <- ((N-g-p)*m*maha.2)/(p*((m-1)*(N-g)-m*maha.2))
	
	return(adk.test(list(hotel.1,hotel.2))$adk[1,2])
}	
# Box's M
Box.M <- function(grp1, grp2){
	k <- 2; p <- length(grp1[1,]); n <- length(grp1[,1]); N <- n*k;
	
	cov.1 <- cov(grp1); cov.2 <- cov(grp2);
	cov.pool <- (1/(2*n-2))*((n-1)*cov.1+(n-1)*cov.2)
	
	df <-  0.5*(k-1)*p*(p+1);
	C <- ((2*(p^2)+3*p-1)*(k+1))/(6*(p+1)*k*(n-1)); 
	M <- 0.5*((n-1)*log(det(cov.1))+(n-1)*log(det(cov.2)))-0.5*(2*n-2)*log(det(cov.pool))
	return(pchisq(-2*(1-C)*M, df, lower.tail=F));
}  
# Hawkins' Test
Hawkins <- function(grp1, grp2) {
	p <- length(grp1[1,]);
	g <- 2
	n <- length(grp1[,1]);
	m <- length(grp2[,1]);
	N <- n + m
	v <- N - g - 1;
	
	cov.1 <- cov(grp1); cov.2 <- cov(grp2);
	cov.pool <- (1/(N-g))*((n-1)*cov.1+(m-1)*cov.2)
		
	maha.1 <- mahalanobis(grp1, colMeans(grp1), cov.pool)
	maha.2 <- mahalanobis(grp2, colMeans(grp2), cov.pool)
	 
	hotel.1 <- ((N-g-p)*n*maha.1)/(p*((n-1)*(N-g)-n*maha.1))
	hotel.2 <- ((N-g-p)*m*maha.2)/(p*((m-1)*(N-g)-m*maha.2))
	
	a.1 <- pf(hotel.1, p, v-p+1, lower.tail=F)
	a.2 <- pf(hotel.2, p, v-p+1, lower.tail=F)

	Z1.1 <- -sqrt(3/n)*sum(2*a.1-1)
	Z1.2 <- -sqrt(3/m)*sum(2*a.2-1)

	diff <- abs(Z1.1 - Z1.2)
	p.value <-  2*pnorm(diff, mean = 0, sd = 1.444983, lower.tail = FALSE)

	return(p.value)
}	

# Creates function for joining datasets
combine <- function(x,y) {merge(x,y, all = TRUE, sort = FALSE)}

# Simulation Preferences
sim <- 20000 # Num of Sumlation Runs
N <- c(10, 20, 40, 100) # Different Sample Sizes

# Starts Simulation
# nulls p.value vectors
p.values <- results05 <- results10 <- NULL
counter <- 1
for (n in N) {
	for (c in seq(0,1,0.25)) {
		mat <- matrix(c(c(2^(1*c),0,0), c(0,2^(2*c),0), c(0,0,2^(3*c))), ncol=3, nrow=3)
		Box.p <- Adk.p <- Hawkins.p <- NULL
		for (i in 1:sim) {
			# generates random MN data 
			grp1 <- rmnorm(n, mean, cov);
			grp2 <- rmnorm(n, mean, mat);
			# perfroms test on random data
			Box.p <- append(Box.p, Box.M(grp1, grp2));
			Adk.p <- append(Adk.p, Adk(grp1, grp2));
			Hawkins.p <- append(Hawkins.p, Hawkins(grp1, grp2));
			}
		# Helps deal with NULL database problems	
		# creates dataframe of p.values 
		cov.mat <- rep('covB', sim)
		constant <- rep(c, sim)	
		sam.size <- rep(n, sim)	
		if (counter == 1) {
			p.values <- data.frame(cov.mat, constant, n, Box.p, Adk.p, Hawkins.p)
		}	
		else {
			results <- data.frame(cov.mat, constant, n, Box.p, Adk.p, Hawkins.p)
			p.values <- combine(list(p.values), list(results))
		}
		
		# creates dataframe of the number of rejections at 0.05
		alpha <- 0.05
		sam.size <- n
		test.matrix <- 'covB'
		Box.rejects <- length(subset(Box.p, Box.p <= alpha))
		Adk.rejects <- length(subset(Adk.p, Adk.p <= alpha))
		Hawkins.rejects <- length(subset(Hawkins.p, Hawkins.p <= alpha))
		if (counter == 1) {
			results05 <- data.frame(test.matrix, c, sam.size, alpha, Box.rejects, Adk.rejects, Hawkins.rejects)
		}
		else {	
			rejects05 <- data.frame(test.matrix, c, sam.size, alpha, Box.rejects, Adk.rejects, Hawkins.rejects)
			results05 <- combine(list(results05), list(rejects05))
		}
		
		# creates dataframe of the number of rejections at 0.10
		alpha <- 0.10
		Box.rejects <- length(subset(Box.p, Box.p <= alpha))
		Adk.rejects <- length(subset(Adk.p, Adk.p <= alpha))
		Hawkins.rejects <- length(subset(Hawkins.p, Hawkins.p <= alpha))
		if (counter == 1) {
			results10 <- data.frame(test.matrix, c, sam.size, alpha, Box.rejects, Adk.rejects, Hawkins.rejects)
		}
		else {	
			rejects10 <- data.frame(test.matrix, c, sam.size, alpha, Box.rejects, Adk.rejects, Hawkins.rejects)
			results10 <- combine(list(results10), list(rejects10))
		}
		counter <- counter + 1	
	}
}	


# Exporting files
write.table(results05, '~/Desktop/Results/Results05 - covB.txt', sep = "\t", row.names = FALSE)
write.table(results10, '~/Desktop/Results/Results10 - covB.txt', sep = "\t", row.names = FALSE)
write.table(p.values, '~/Desktop/Results/p.values - covB.txt', sep = "\t", row.names = FALSE)

write.table(Type1, '~/Desktop/Results/Type1.txt', sep = "\t", row.names = FALSE)

results05 <- read.table('~/Desktop/Results/Results05 - covB.txt', sep = "\t", header = TRUE)
results10 <- read.table('~/Desktop/Results/Results10 - covB.txt', sep = "\t", header = TRUE)
p.values <- read.table('~/Desktop/Results/p.values - covB.txt', sep = "\t", header = TRUE)

require(Hmisc)

#Confidence Intervals for alpha = 0.05
Box.CI05 <- binconf(results05$Box.rejects, 20000)
row.names(Box.CI05) <- 1:20; 
Box.CI05 <- as.data.frame(Box.CI05)
names(Box.CI05) <- c("Box.PointEst", "Box.Lower", "Box.Upper")

Adk.CI05 <- binconf(results05$Adk.rejects, 20000)
row.names(Adk.CI05) <- 1:20; 
Adk.CI05 <- as.data.frame(Adk.CI05)
names(Adk.CI05) <- c("Adk.PointEst", "Adk.Lower", "Adk.Upper")

Hawkins.CI05 <- binconf(results05$Hawkins.rejects, 20000)
row.names(Hawkins.CI05) <- 1:20; 
Hawkins.CI05 <- as.data.frame(Hawkins.CI05)
names(Hawkins.CI05) <- c("Hawkins.PointEst", "Hawkins.Lower", "Hawkins.Upper")

percents05 <- data.frame(results05[,1:4], Box.CI05, Adk.CI05, Hawkins.CI05)

Box.CI05 <- data.frame(results05[,1:4], Box.CI05)
Adk.CI05 <- data.frame(results05[,1:4], Adk.CI05)
Hawkins.CI05 <- data.frame(results05[,1:4], Hawkins.CI05)

#Confidence Intervals for alpha = 0.10
Box.CI10 <- binconf(results10$Box.rejects, 20000)
row.names(Box.CI10) <- 1:20; 
Box.CI10 <- as.data.frame(Box.CI10)
names(Box.CI10) <- c("Box.PointEst", "Box.Lower", "Box.Upper")

Adk.CI10 <- binconf(results10$Adk.rejects, 20000)
row.names(Adk.CI10) <- 1:20; 
Adk.CI10 <- as.data.frame(Adk.CI10)
names(Adk.CI10) <- c("Adk.PointEst", "Adk.Lower", "Adk.Upper")

Hawkins.CI10 <- binconf(results10$Hawkins.rejects, 20000)
row.names(Hawkins.CI10) <- 1:20; 
Hawkins.CI10 <- as.data.frame(Hawkins.CI10)
names(Hawkins.CI10) <- c("Hawkins.PointEst", "Hawkins.Lower", "Hawkins.Upper")

percents10 <- data.frame(results10[,1:4], Box.CI10, Adk.CI10, Hawkins.CI10)

Box.CI10 <- data.frame(results10[,1:4], Box.CI10)
Adk.CI10 <- data.frame(results10[,1:4], Adk.CI10)
Hawkins.CI10 <- data.frame(results10[,1:4], Hawkins.CI10)

write.table(percents05, '~/Desktop/Results/percents05 - covA.txt', sep = "\t", row.names = FALSE)
write.table(percents10, '~/Desktop/Results/percents10 - covA.txt', sep = "\t", row.names = FALSE)

write.table(Box.CI05, '~/Desktop/Results/Box.CI05 - covA.txt', sep = "\t", row.names = FALSE)
write.table(Adk.CI05, '~/Desktop/Results/Adk.CI05 - covA.txt', sep = "\t", row.names = FALSE)
write.table(Hawkins.CI05, '~/Desktop/Results/Hawkins.CI05 - covA.txt', sep = "\t", row.names = FALSE)

write.table(Box.CI10, '~/Desktop/Results/Box.CI10 - covA.txt', sep = "\t", row.names = FALSE)
write.table(Adk.CI10, '~/Desktop/Results/Adk.CI10 - covA.txt', sep = "\t", row.names = FALSE)
write.table(Hawkins.CI10, '~/Desktop/Results/Hawkins.CI10 - covA.txt', sep = "\t", row.names = FALSE)

write.table(results05.H0, '~/Desktop/Results/results05.H0.txt', sep = "\t", row.names = FALSE)
write.table(results10.H0, '~/Desktop/Results/results10.H0.txt', sep = "\t", row.names = FALSE)

extract <- function(x) {
	dat <- rbind(
		c(x[1,2], x[1,3], x[1,4], x[1,5] + x[5,5], x[1,6] + x[5,6], x[1,7] + x[5,7]),
		c(x[2,2], x[2,3], x[2,4], x[2,5] + x[6,5], x[2,6] + x[6,6], x[2,7] + x[6,7]),
		c(x[3,2], x[3,3], x[3,4], x[3,5] + x[7,5], x[3,6] + x[7,6], x[3,7] + x[7,7]),
		c(x[4,2], x[4,3], x[4,4], x[4,5] + x[8,5], x[4,6] + x[8,6], x[4,7] + x[8,7]))
	dat <- as.data.frame(dat)	
	names(dat) <- names(x)[2:7]
	dat
}	

Hawkins.05 <- binconf(results05.H0$Hawkins.rejects, 40000)
Hawkins.10 <- binconf(results10.H0$Hawkins.rejects, 40000)
row.names(Hawkins.05) <- row.names(Hawkins.10) <- 1:4 
Hawkins.05 <- as.data.frame(Hawkins.05)
Hawkins.10 <- as.data.frame(Hawkins.10)
names(Hawkins.05) <- names(Hawkins.10) <- c("Hawkins.PointEst", "Hawkins.Lower", "Hawkins.Upper")

Adk.05 <- binconf(results05.H0$Adk.rejects, 40000)
Adk.10 <- binconf(results10.H0$Adk.rejects, 40000)
row.names(Adk.05) <- row.names(Adk.10) <- 1:4 
Adk.05 <- as.data.frame(Adk.05)
Adk.10 <- as.data.frame(Adk.10)
names(Adk.05) <- names(Adk.10) <- c("Adk.PointEst", "Adk.Lower", "Adk.Upper")

Box.05 <- binconf(results05.H0$Box.rejects, 40000)
Box.10 <- binconf(results10.H0$Box.rejects, 40000)
row.names(Box.05) <- row.names(Box.10) <- 1:4 
Box.05 <- as.data.frame(Box.05)
Box.10 <- as.data.frame(Box.10)
names(Box.05) <- names(Box.10) <- c("Box.PointEst", "Box.Lower", "Box.Upper")

CI05.H0 <- data.frame(results05.H0[1:3], Box.05, Hawkins.05, Adk.05)
CI10.H0 <- data.frame(results10.H0[1:3], Box.10, Hawkins.10, Adk.10)

write.table(CI05.H0, '~/Desktop/Results/CI05.H0.txt', sep='\t', row.names = FALSE)
write.table(CI10.H0, '~/Desktop/Results/CI10.H0.txt', sep='\t', row.names = FALSE)
