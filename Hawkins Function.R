# 0 is H0 is False
# 1 is H0 is True

hawkins <- function(grp1, grp2, threshold){	
	p <- length(grp1[1,]);
	n <- length(grp2[,1]);
	g <- 2
	v <- n - g - p;
	
	cov.1 <- cov(grp1); cov.2 <- cov(grp2);
	cov.pool <- (1/(2*n-2))*((n-1)*cov.1+(n-1)*cov.2)
	
	maha.1 <- mahalanobis(grp1, colMeans(grp1), cov.pool)
	maha.2 <- mahalanobis(grp2, colMeans(grp2), cov.pool)
	
	six.1 <- ((v*n/p)*maha.1)/((v+p)*(n-1)-n*maha.1)
	six.2 <- ((v*n/p)*maha.2)/((v+p)*(n-1)-n*maha.2)
	
	a.1 <- pf(six.1, p, v, lower.tail = F)
	a.2 <- pf(six.2, p, v, lower.tail = F)
	
	a.T <- c(a.1, a.2)
	
	Z1.1 <- -sqrt(3/n)*sum(2*a.1-1)
	Z1.2 <- -sqrt(3/n)*sum(2*a.2-1)
	Z1.T <- -sqrt(3/(2*n))*sum(2*a.T-1)
	
	#Z2.1 <- -sqrt(5/n)*sum(0.5*(3*(2*a.1-1)^2)-1)
	#Z2.1 <- -sqrt(5/n)*sum(0.5*(3*(2*a.2-1)^2)-1)
	#Z2.T <- -sqrt(5/(2*n))*sum(0.5*(3*(2*a.T-1)^2)-1)
	
	if (sign(Z1.1)*sign(Z1.2) > 0){
		return(1)
	}
	if (abs(abs(Z1.1)-abs(Z1.2)) < threshold) {
		return(1)
	}
	else 
		return(0)	
}	

Hawks <- NULL
for (i in 1:10000){
	grp1 <- rmnorm(30, mean, covA);
	grp2 <- rmnorm(30, mean, covA);
	Hawks <- append(Hawks, hawkins(grp1,grp2, 1))
}
1-mean(Hawks)


stats1 <- stats2 <- value <- NULL;
for (i in 1:10000) {
	grp1 <- rmnorm(30, mean, covA);
	grp2 <- rmnorm(30, mean, covA);
	stats1 <- append(stats1, hawkins(grp1, grp2)[[1]]);
	stats2 <- append(stats2, hawkins(grp1, grp2)[[2]]);
}	
value <- append(value, sign(stats1)*sign(stats2));
length(subset(value, value<0))/length(value)

neg <- value == -1
diff <- abs(stats1[neg])-abs(stats2[neg])
(length(diff[abs(diff)<1])+length(subset(value, value == 1)))/length(value)