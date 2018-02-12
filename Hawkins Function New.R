############################## HAWKINS TEST ##################################
# p is the number of variables measured
# g is the number of groups (fixed to 2)
# n is the number of obs in group 1
# m is the number of obs in group 2
# N is the total number of obs (n+m)
# v is a value needed for dfs calculations

# 0 is H0 is True
# 1 is H0 is False


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

	diff <- abs(Z1.1 - Z1.2)
	p.value <-  2*pnorm(diff, mean = 0, sd = 1.444983, lower.tail = FALSE)

	return(p.value)
	# return(c(Z1.1, Z1.2))
	# print(Z1.1)
	# print(Z1.2)
	# print(Z1.T)	
}	

p.value <- NULL
for (i in 1:100000) {
	grp1 <- rmnorm(30, mean, covA);
	grp2 <- rmnorm(30, mean, covA);	
	p.value <- append(p.value, Hawkins(grp1, grp2))
	}

length(subset(p.value, p.value <= 0.025))/length(p.value)
length(subset(p.value, p.value <= 0.05))/length(p.value)


loc <- '~/Desktop/Aitchison.txt'
Aitchison <- read.table(loc, header=T)
Aitchison.A <- split(Aitchison, Aitchison$Group)$A
Aitchison.B <- split(Aitchison, Aitchison$Group)$B

hawkins(Aitchison.A[,3:8],Aitchison.B[,3:8])

stats <- sign <- diff <- add <- abs.diff <- NULL
for (i in 1:100000) {
	grp1 <- rmnorm(30, mean, covB);
	grp2 <- rmnorm(30, mean, covC);	
	stats <- rbind(stats, hawkins(grp1, grp2))
	#sign <- append(sign, sign(stats[i,1])*sign(stats[i,2]))
	#diff <- append(diff, stats[i,1]-stats[i,2])
	#add <- append(diff, stats[i,1]+stats[i,2])
	#abs.diff <- append(diff, abs(abs(stats[i,1])-abs(stats[i,2])))
}	
stats <- cbind(stats) #, sign, diff, add, abs diff)
stats <- as.data.frame(stats)
names(stats) <- c('Z1', 'Z2')
#reject <- subset(stats, stats$sign < 0 & stats$diff > 1)
#length(reject[,1])/length(stats[,1])

#locate <- function(x) {
#	reject <- subset(stats, stats$sign < 0 & stats$diff > x)
#	return(length(reject[,1])/length(stats[,1]))
#}


probs <- quantile(probs=seq(0,1,0.005),x=stats$diff)
names <- seq(0, 1, .005)
order <- order(names, decreasing = TRUE)
plot(probs, names[order], type = 'l', xlab = "Threshold Value", ylab = "Alpha")
abline(h=0.05, col = 'red')
abline(h=0.10, col = 'blue')


####################
# Aplha covA 	covB	covC	covD	covE
# 0.10	0.7333	0.7406	0.7358	0.7340	0.7365
# 0.05	0.9075	0.9082	0.9093	0.9040	0.9104

quantile(probs=c(0.95),x=diff)

hist(stats$diff, freq = FALSE, border = '2')
curve(dnorm(x = stats$diff, mean = 0, sd = 0.675062346), col = 1, lty = 2, lwd = 2, add = TRUE)
	
total$sign <-  sign(total$Z1)*sign(total$Z2)
length(subset(total, total$sign > 0)[,1])/length(total[,1])
total$diff <- total$Z1 - total$Z2
total$sum <- total$Z1 + total$Z2
total$abs <- abs(total$Z1) - abs(total$Z2)

quantile(probs = c(0.003, 0.05, 0.32, 0.5, 0.68, 0.95, 0.997), x = total$diff)

png(file = "~/Desktop/Figure 2.%d.png", height = 5.53, width = 8.32, units = "in", res = 90)

hist(total$diff, freq = FALSE, main = expression(paste(bold("Figure 2.1"), " - Histogram of γ from Simulations")), xlab = "Value of γ ")

hist(total$diff, freq = FALSE, main = expression(paste(bold("Figure 2.2"), " - Histogram of γ from Simulations with N(0,2.087976)")), xlab = "Value of γ ", breaks = 100)
x <- total$diff
curve(dnorm(x, mean = 0,  sd = 1.444983), col = 2, lty = 2, lwd = 2, add = TRUE)
legend(x = "topright", "N(0,2.087976)", col = "red", lty = 2, lwd = 2, box.lwd = 0)

qqnorm(sort(total$diff), type = "l", main = expression(paste(bold("Figure 2.3"), " - Normal Q-Q Plot")))
qqline(sort(total$diff), col=2, lty = 2)
legend(x = "bottomright", "Q-Q Line", col="red", lty = 2, box.lwd = 0)

dev.off()

stats$diff <- stats$Z1 - stats$Z2
stats$abs <-  abs(stats$Z1) - abs(stats$Z2)
stats$p.value <- pnorm(stats$diff, mean = 0, sd = 1.444983, lower.tail = TRUE)

