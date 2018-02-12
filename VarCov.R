##################################################################################
# Possible Covariance Variance Matrices
# covA is the identity matrix
# covB is the identity matrix with decreasing covariance by half 
# covC is a matrix with increasing variances by power
# mean is the mean vector 

library('mnormt')
library('adk')

covA <- matrix(c(c(1,0,0), c(0,1,0), c(0,0,1)), ncol=3, nrow=3)
covB <- matrix(c(c(1,0.5,0.25), c(0.5,1,0.5), c(0.25,0.5,1)), ncol=3, nrow=3)
covC <- matrix(c(c(2,0,0), c(0,4,0), c(0,0,8)), ncol=3, nrow=3)
covD <- matrix(c(c(2,0.5,0.25), c(0.5,4,0.5), c(0.25,0.5,8)), ncol=3, nrow=3)
covEd<- matrix(c(c(10,0,0), c(0,1,0), c(0,0,1)), ncol=3, nrow=3)
mean = rep(0,3)

N <- c(10, 20, 40, 100) 

for (c in seq(0,1,0.25)) {
	mat <- matrix(c(c(1,0.5*c,0.25*c), c(0.5*c,1,0.5*c), c(0.25*c,0.5*c,1)), ncol=3, nrow=3)
	print(mat)
}

for (c in seq(0,1,0.25)) {
	mat <- matrix(c(c(2^(1*c),0,0), c(0,2^(2*c),0), c(0,0,2^(3*c))), ncol=3, nrow=3)
	print(mat)
}
 
Box <- Adk <- Hawk <- Con <- NULL
Box.data <- Adk.data <- Hawk.data <- NULL
for (c in seq(0,1,0.25)) {
	mat <- matrix(c(c(1,0.5*c,0.25*c), c(0.5*c,1,0.50*c), c(0.25*c,0.5*c,1)), ncol=3, nrow=3)
	Box.p <- Adk.p <- Hawk.p <- NULL
	for (i in 1:10000) {
		grp1 <- rmnorm(30, mean, covA);
		grp2 <- rmnorm(30, mean, mat);
		
		Box.p <- append(Box.p, Box.M(grp1, grp2));
		Adk.p <- append(Adk.p, adktest(grp1, grp2)); 
		Hawk.p <- append(Hawk.p, hawkins(grp1, grp2, 1));		
	}
	Box.data <- cbind(Box.data, Box.p)
	Adk.data <- cbind(Adk.data, Adk.p)
	Hawk.data <- cbind(Hawk.data, Hawk.p)
			
	Con <- append(Con, c)
	Box <- append(Box, length(subset(Box.p, Box.p<=0.05))/length(Box.p))
	Adk <- append(Adk, length(subset(Adk.p, Adk.p<=0.05))/length(Adk.p))
	Hawk <- append(Hawk, 1-mean(Hawk.p))
}
data.frame(Con, Box, Adk, Hawk)


Boxs <- Adks <- Hawks <- Cons <- Sizes <- NULL
Box.data <- Adk.data <- Hawk.data <- NULL

for(n in N) {
	Box <- Adk <- Hawk <- Con <- Size <- NULL
	for (c in seq(0,1,0.25)) {
		mat <- matrix(c(c(1,0.5*c,0.25*c), c(0.5*c,1,0.50*c), c(0.25*c,0.5*c,1)), ncol=3, nrow=3)
		Box.p <- Adk.p <- Hawk.p <- NULL
		for (i in 1:10000) {
			grp1 <- rmnorm(n, mean, covA);
			grp2 <- rmnorm(n, mean, mat);
			
			Box.p <- append(Box.p, Box.M(grp1, grp2));
			Adk.p <- append(Adk.p, adktest(grp1, grp2)); 
			Hawk.p <- append(Hawk.p, hawkins(grp1, grp2, 0.8));		
		}
		Box.data <- cbind(Box.data, Box.p)
		Adk.data <- cbind(Adk.data, Adk.p)
		Hawk.data <- cbind(Hawk.data, Hawk.p)
				
		Con <- append(Con, c)
		Size <- append(Size, n)
		Box <- append(Box, length(subset(Box.p, Box.p<=0.1))/length(Box.p))
		Adk <- append(Adk, length(subset(Adk.p, Adk.p<=0.1))/length(Adk.p))
		Hawk <- append(Hawk, 1-mean(Hawk.p))
	}
	data.frame(Con, Size, Box, Adk, Hawk)
	
	Cons <- append(Cons, Con)
	Sizes <- append(Sizes, Size)
	Boxs <- append(Boxs, Box)
	Adks <- append(Adks, Adk)
	Hawks <- append(Hawks, Hawk)
	
}	
data.frame(Cons, Sizes, Boxs, Adks, Hawks)