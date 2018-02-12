# Example using iris data

virginca <- split(iris, iris$Species)$virginica[,1:4]
versicolor <- split(iris, iris$Species)$versicolor[,1:4]
setosa <- split(iris, iris$Species)$setosa[,1:4]

################   Box Plots  ###############################

constant <- rep(seq(0,1,0.25),4)
sample.size <- rep(c(10,20,40,100),c(5,5,5,5))

#### Box.M

Box.data <- as.data.frame(Box.data)
# By constant, fixed sample size
names(Box.data) <- constants
# n = 10
boxplot(Box.data[,1:5], ylim = c(0,1))
# n = 20
boxplot(Box.data[,6:10], ylim = c(0,1))
# n = 40
boxplot(Box.data[,11:15], ylim = c(0,1))
# n = 100
boxplot(Box.data[,16:20], range = 0, ylim = c(0,1))

# By sample size, fixed constant
names(Box.data) <- sample.size
# c = 0
boxplot(Box.data[,c(1,6,11,16)], range = 0, ylim = c(0,1))
# c = 0.25
boxplot(Box.data[,c(2,7,12,17)], range = 0, ylim = c(0,1))
# c = 0.5
boxplot(Box.data[,c(3,8,13,18)], range = 0, ylim = c(0,1))
# c = 0.75
boxplot(Box.data[,c(4,9,14,19)], range = 0, ylim = c(0,1))
# c = 1.0
boxplot(Box.data[,c(5,10,15,20)], range = 0, ylim = c(0,1))
#### Adk

Adk.data <- as.data.frame(Adk.data)
# By constant, fixed sample size
names(Adk.data) <- constants
# n = 10
boxplot(Adk.data[,1:5],  range = 0, ylim = c(0,1))
# n = 20
boxplot(Adk.data[,6:10],  range = 0, ylim = c(0,1))
# n = 40
boxplot(Adk.data[,11:15],  range = 0, ylim = c(0,1))
# n = 100
boxplot(Adk.data[,16:20], range = 0, ylim = c(0,1))

# By sample size, fixed constant
names(Adk.data) <- sample.size
# c = 0
boxplot(Adk.data[,c(1,6,11,16)], range = 0, ylim = c(0,1))
# c = 0.25
boxplot(Adk.data[,c(2,7,12,17)], range = 0, ylim = c(0,1))
# c = 0.5
boxplot(Adk.data[,c(3,8,13,18)], range = 0, ylim = c(0,1))
# c = 0.75
boxplot(Adk.data[,c(4,9,14,19)], range = 0, ylim = c(0,1))
# c = 1.0
boxplot(Adk.data[,c(5,10,15,20)], range = 0, ylim = c(0,1))