# Loads required packages
require('MASS');
require('MBESS');

# Loads Hawkins' Test
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

# Outer Loop - Saves files, clears data frame
for (i in 1:40){	
	# Number of in simulations (25000*40 = 10^6)
	sim <- 25000;
	# resets inner loop
	j <- 1;
	# nullifes data frames
	Results <- data.frame(NULL);
	Coefficents <- data.frame(NULL);
	# prints simulation outer loop progress
	print(paste('Started:', i, 'out of', '14', sep =' '));
	# Inner Loop - Performs simulations till it reachs sim amount
	while (j <= sim) {
		# Creates Correlation Matrix
		# generates diagionals, all equal to 1
		a11 <- a22 <- a33 <- 1;
		# generates off-diagionals from U(0,1)
		a12 <- a21 <- runif(1, 0, 1);
		a13 <- a31 <- runif(1, 0, 1);
		a23 <- a32 <- runif(1, 0, 1);
		# creates cor matrix from values
		cor <- matrix(c(a11, a21, a31, a12, a22, a23, a31, a32, a33), nrow = 3, ncol = 3);
		# generates sd from chisq
		sds <- rchisq(3, 5);
		# converts cor to cov
		cov <- cor2cov(cor, sds);
		# Starts grp data simulation and Hawkins' Test
		# checks for positive definative matrices 
		if (det(cov) > 0) {	
			# saves coefficents
			Coefficents <- rbind(Coefficents, data.frame(t(as.vector(cov))));
			# generates grp data
			grp1 <- rmnorm(30, c(0,0,0), cov);
			grp2 <- rmnorm(30, c(0,0,0), cov);
			# performs Hawkins' Test
			zs <- Hawkins(grp1, grp2);
			# save results from Hawkins' Test
			Results <- rbind(Results, zs);
			# controls iterations 
			j <- j+1;
			# prints inner loop progress
			if (j%%5000 == 0) {
				print(paste("Finished:", format(j), "out of", format(sim),sep=" "));
			}# end print if
		}# end grp data simulation if
	}#end while (inner) loop
	# labels Cofficent's rows with cov mat indexies 	
	names(Coefficents) <- c('a11', 'a12', 'a13', 'a21', 'a22', 'a23', 'a31', 'a32', 'a33');
	# labels Results's row with Z values
	names(Results) <-  c('z1', 'z2');
	# Combines Results and Coefficents into one data frame
	Results <- cbind(Coefficents, Results);
	# Saves Results data frame
	write.table(Results, paste('~/Desktop/Results/Results_', i, '.txt', sep=''));	
}

# Loads simulation files and creates one file
# Results <- NULL;
# for (i in 1:40) {
	# Temp <- NULL
	# Temp <- read.table(paste('~/Desktop/Results/Results_', i, '.txt', sep=''));
	# Results <- rbind(Results, Temp);
	# write.table(Results, '~/Desktop/Results.txt');
# }

Results <- read.table(paste('~/Dropbox/School Work/Thesis/Results/Reg Results/Results.txt', sep=''));

# Creates Gamma and Performs Regression
gamma <- abs(Results$z1 - Results$z2);
Results <- cbind(Results, gamma);
summary(lm(gamma~a11 + a22 + a33 + a12 + a13 + a23, data=Results))

