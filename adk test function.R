library(adk)
Adk <- function(grp1,grp2){	
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
	
	#print(hotel.1)
	#print(hotel.2)
	return(adk.test(list(hotel.1,hotel.2))$adk[1,2])
}	

grp1 <- rmnorm(20, mean, covA);
grp2 <- rmnorm(20, mean, covA);	
adktest(grp1,grp2)

stats <- NULL
for (i in 1:10000) {
	grp1 <- rmnorm(20, mean, covA);
	grp2 <- rmnorm(20, mean, covC);	
	stats <- append(stats, adktest(grp1,grp2))
}
length(subset(stats, stats <=0.05))/length(stats)