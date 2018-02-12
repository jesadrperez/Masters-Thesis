######## NULL HYPOTHESIS PLOTS ########## 
p.values <- read.table('~/Desktop/Results/Type1.txt', sep = "\t", header = TRUE)
p.values <- p.values[,c(1:3,4,6,5)]

require(vioplot)

#### By Test
vio.test <- function(loc, subtitle) {
	vioplot(p.values$Box.p[loc], p.values$Hawkins.p[loc], p.values$Adk.p[loc], names = x.values, col = "white", drawRect = F, wex = 0.75, range = 0, ylim = c(0,1))
boxplot(p.values$Box.p[loc], p.values$Hawkins.p[loc], p.values$Adk.p[loc], names = x.values, border = "black", pars = list(boxwex = .25), add = T, outline = F)
	title(main = title, sub = subtitle,	ylab = y.title,	xlab = x.title)
	abline(h=0.05, col = "red", lty = 5)
	abline(h=0.1, col = "red", lty = 5)
	text(x=0.65, y = 0.07, "α = 0.05")
	text(x=0.65, y = 0.12, "α = 0.10")
}

title <- expression(paste(bold("Figure 4.1.?")," - Violin Plots for P-Values When H", scriptstyle(0), " is True"))

x.values <- c("Box's M", "Hawkins", "Non-Parametric")
x.title <- "Test"
y.title <- "P-Value"
d <- 1
png(file="~/Desktop/Plots/Null/Test/vioplot 4.1.%d.png", height = 5.53, width = 8.32, units = "in", res = 90)

vio.test(c(1:20000,80001:100000), "(n=10)")
vio.test(c(20001:40000,100001:120000), "(n=20)")
vio.test(c(40001:60000,120001:140000), "(n=40)")
vio.test(c(60001:80000,140001:160000), "(n=100)")

dev.off()


### Violin Plots by Constant for Cov A #############

require(vioplot)
p.values <- read.table('~/Desktop/Results/p.values - covA.txt', sep = "\t", header = TRUE)
p.values <- p.values[,c(1:3,4,6,5)]

vio.constant <- function(one, two, three, four, subtitle) {
	vioplot(one, two, three, four, names = x.values, col = "white", drawRect = F, wex = 0.75, range = 0, ylim = c(0,1))
boxplot(one, two, three, four, names = x.values, border = "black", pars = list(boxwex = .25), add = T, outline = F)
	title(main = title, sub = subtitle,	ylab = y.title,	xlab = x.title)
	abline(h=0.05, col = "red", lty = 5)
	abline(h=0.1, col = "red", lty = 5)
	text(x=0.65, y = 0.07, "α = 0.05")
	text(x=0.65, y = 0.12, "α = 0.10")
}

title <- expression(paste("Violin Plots for P-Values When H", scriptstyle(0), " is False for ", Sigma, scriptstyle(A)))
x.values <- c('0.25', '0.50', '0.75', '1.00')
x.title <- "Constant Value"
y.title <- "P-Value"
sub <- NULL

png(file="~/Desktop/Plots/Cov A/CV/vioplot%d.png", height = 5.53, width = 8.32, units = "in", res = 90)

vio.constant(p.values$Box.p[20001:40000], p.values$Box.p[40001:60000], p.values$Box.p[60001:80000], p.values$Box.p[80001:100000], "Box's M (n=10)")

vio.constant(p.values$Hawkins.p[20001:40000], p.values$Hawkins.p[40001:60000], p.values$Hawkins.p[60001:180000], p.values$Hawkins.p[80001:100000], "Hawkins' Test (n=10)")

vio.constant(p.values$Adk.p[20001:40000], p.values$Adk.p[40001:60000], p.values$Adk.p[60001:80000], p.values$Adk.p[80001:100000], "Non-Parametric Test (n=10)")

vio.constant(p.values$Box.p[120001:140000], p.values$Box.p[140001:160000], p.values$Box.p[160001:180000], p.values$Box.p[180001:200000], "Box's M (n=20)")

vio.constant(p.values$Hawkins.p[120001:140000], p.values$Hawkins.p[140001:160000], p.values$Hawkins.p[160001:180000], p.values$Hawkins.p[180001:200000], "Hawkins' Test (n=20)")

vio.constant(p.values$Adk.p[120001:140000], p.values$Adk.p[140001:160000], p.values$Adk.p[160001:180000], p.values$Adk.p[180001:200000], "Non-Parametric Test (n=20)")

vio.constant(p.values$Box.p[220001:240000], p.values$Box.p[240001:260000], p.values$Box.p[260001:280000], p.values$Box.p[280001:300000], "Box's M (n=40)")

vio.constant(p.values$Hawkins.p[220001:240000], p.values$Hawkins.p[240001:260000], p.values$Hawkins.p[260001:280000], p.values$Hawkins.p[280001:300000], "Hawkins' Test (n=40)")

vio.constant(p.values$Adk.p[220001:240000], p.values$Adk.p[240001:260000], p.values$Adk.p[260001:280000], p.values$Adk.p[280001:300000], "Non-Parametric Test (n=40)")

vio.constant(p.values$Box.p[320001:340000], p.values$Box.p[340001:360000], p.values$Box.p[360001:380000], p.values$Box.p[380001:400000], "Box's M (n=100)")

vio.constant(p.values$Hawkins.p[320001:340000], p.values$Hawkins.p[340001:360000], p.values$Hawkins.p[360001:380000], p.values$Hawkins.p[380001:400000], "Hawkins' Test (n=100)")

vio.constant(p.values$Adk.p[320001:340000], p.values$Adk.p[340001:360000], p.values$Adk.p[360001:380000], p.values$Adk.p[380001:400000], "Non-Parametric Test (n=100)")

dev.off()

######## Violin Plots by Test for Cov A ########## 

vio.test <- function(loc, subtitle) {
	vioplot(p.values$Box.p[loc], p.values$Hawkins.p[loc], p.values$Adk.p[loc], names = x.values, col = "white", drawRect = F, wex = 0.75, range = 0, ylim = c(0,1))
boxplot(p.values$Box.p[loc], p.values$Hawkins.p[loc], p.values$Adk.p[loc], names = x.values, border = "black", pars = list(boxwex = .25), add = T, outline = F)
	title(main = title, sub = subtitle,	ylab = y.title,	xlab = x.title)
	abline(h=0.05, col = "red", lty = 5)
	abline(h=0.1, col = "red", lty = 5)
	text(x=0.65, y = 0.07, "α = 0.05")
	text(x=0.65, y = 0.12, "α = 0.10")
}

title <- expression(paste(bold("Figure 4.?? "), "- Violin Plots for P-Values When H", scriptstyle(0), " is False for ", Sigma, scriptstyle(A), "."))
x.values <- c("Box's M", "Hawkins", "Non-Parametric")
x.title <- "Test"
y.title <- "P-Value"

png(file="~/Desktop/Plots/Cov A/Test/vioplot%d.png", height = 5.53, width = 8.32, units = "in", res = 90)

vio.test(20001:40000, "Constant Value 0.25 (n=10)")
vio.test(40001:60000, "Constant Value 0.50 (n=10)")
vio.test(60001:80000, "Constant Value 0.75 (n=10)")
vio.test(80001:100000, "Constant Value 1.00 (n=10)")

vio.test(120001:140000, "Constant Value 0.25 (n=20)")
vio.test(140001:160000, "Constant Value 0.50 (n=20)")
vio.test(160001:180000, "Constant Value 0.75 (n=20)")
vio.test(180001:200000, "Constant Value 1.00 (n=20)")

vio.test(220001:140000, "Constant Value 0.25 (n=40)")
vio.test(240001:260000, "Constant Value 0.50 (n=40)")
vio.test(260001:280000, "Constant Value 0.75 (n=40)")
vio.test(280001:300000, "Constant Value 1.00 (n=40)")

vio.test(320001:340000, "Constant Value 0.25 (n=100)")
vio.test(340001:360000, "Constant Value 0.50 (n=100)")
vio.test(360001:380000, "Constant Value 0.75 (n=100)")
vio.test(380001:400000, "Constant Value 1.00 (n=100)")

dev.off()

################ COV B Vio Plots #############################


p.values <- read.table('~/Desktop/Results/p.values - covB.txt', sep = "\t", header = TRUE)
p.values <- p.values[,c(1:3,4,6,5)]

#### By Constant Value

vio.constant <- function(one, two, three, four, subtitle) {
	vioplot(one, two, three, four, names = x.values, col = "white", drawRect = F, wex = 0.75, range = 0, ylim = c(0,1))
boxplot(one, two, three, four, names = x.values, border = "black", pars = list(boxwex = .25), add = T, outline = F)
	title(main = title, sub = subtitle,	ylab = y.title,	xlab = x.title)
	abline(h=0.05, col = "red", lty = 5)
	abline(h=0.1, col = "red", lty = 5)
	text(x=0.65, y = 0.07, "α = 0.05")
	text(x=0.65, y = 0.12, "α = 0.10")
}

title <- expression(paste("Violin Plots for P-Values When H", scriptstyle(0), " is False for ", Sigma, scriptstyle(B)))
x.values <- c('0.25', '0.50', '0.75', '1.00')
x.title <- "Constant Value"
y.title <- "P-Value"
sub <- NULL

png(file="~/Desktop/Plots/Cov B/CV/vioplot%d.png", height = 5.53, width = 8.32, units = "in", res = 90)

vio.constant(p.values$Box.p[20001:40000], p.values$Box.p[40001:60000], p.values$Box.p[60001:80000], p.values$Box.p[80001:100000], "Box's M (n=10)")

vio.constant(p.values$Hawkins.p[20001:40000], p.values$Hawkins.p[40001:60000], p.values$Hawkins.p[60001:180000], p.values$Hawkins.p[80001:100000], "Hawkins' Test (n=10)")

vio.constant(p.values$Adk.p[20001:40000], p.values$Adk.p[40001:60000], p.values$Adk.p[60001:80000], p.values$Adk.p[80001:100000], "Non-Parametric Test (n=10)")

vio.constant(p.values$Box.p[120001:140000], p.values$Box.p[140001:160000], p.values$Box.p[160001:180000], p.values$Box.p[180001:200000], "Box's M (n=20)")

vio.constant(p.values$Hawkins.p[120001:140000], p.values$Hawkins.p[140001:160000], p.values$Hawkins.p[160001:180000], p.values$Hawkins.p[180001:200000], "Hawkins' Test (n=20)")

vio.constant(p.values$Adk.p[120001:140000], p.values$Adk.p[140001:160000], p.values$Adk.p[160001:180000], p.values$Adk.p[180001:200000], "Non-Parametric Test (n=20)")

vio.constant(p.values$Box.p[220001:240000], p.values$Box.p[240001:260000], p.values$Box.p[260001:280000], p.values$Box.p[280001:300000], "Box's M (n=40)")

vio.constant(p.values$Hawkins.p[220001:240000], p.values$Hawkins.p[240001:260000], p.values$Hawkins.p[260001:280000], p.values$Hawkins.p[280001:300000], "Hawkins' Test (n=40)")

vio.constant(p.values$Adk.p[220001:240000], p.values$Adk.p[240001:260000], p.values$Adk.p[260001:280000], p.values$Adk.p[280001:300000], "Non-Parametric Test (n=40)")

vio.constant(p.values$Box.p[320001:340000], p.values$Box.p[340001:360000], p.values$Box.p[360001:380000], p.values$Box.p[380001:400000], "Box's M (n=100)")

vio.constant(p.values$Hawkins.p[320001:340000], p.values$Hawkins.p[340001:360000], p.values$Hawkins.p[360001:380000], p.values$Hawkins.p[380001:400000], "Hawkins' Test (n=100)")

vio.constant(p.values$Adk.p[320001:340000], p.values$Adk.p[340001:360000], p.values$Adk.p[360001:380000], p.values$Adk.p[360001:380000], "Non-Parametric' Test (n=100)")

dev.off()

###### By Statistical Test

vio.test <- function(loc, subtitle) {
	vioplot(p.values$Box.p[loc], p.values$Hawkins.p[loc], p.values$Adk.p[loc], names = x.values, col = "white", drawRect = F, wex = 0.75, range = 0, ylim = c(0,1))
boxplot(p.values$Box.p[loc], p.values$Hawkins.p[loc], p.values$Adk.p[loc], names = x.values, border = "black", pars = list(boxwex = .25), add = T, outline = F)
	title(main = title, sub = subtitle,	ylab = y.title,	xlab = x.title)
	abline(h=0.05, col = "red", lty = 5)
	abline(h=0.1, col = "red", lty = 5)
	text(x=0.65, y = 0.07, "α = 0.05")
	text(x=0.65, y = 0.12, "α = 0.10")
}

title <- expression(paste("Violin Plots for P-Values When H", scriptstyle(0), " is False for ", Sigma, scriptstyle(B)))
x.values <- c("Box's M", "Hawkins", "Non-Parametric")
x.title <- "Test"
y.title <- "P-Value"

#png(file="~/Desktop/Plots/Cov A/Test/vioplot%d.png", height = 4.67, width = 7.07, units = "in", res = 72)

png(file="~/Desktop/Plots/Cov B/Test/vioplot%d.png", height = 5.53, width = 8.32, units = "in", res = 90)

vio.test(20001:40000, "Constant Value 0.25 (n=10)")
vio.test(40001:60000, "Constant Value 0.50 (n=10)")
vio.test(60001:80000, "Constant Value 0.75 (n=10)")
vio.test(80001:100000, "Constant Value 1.00 (n=10)")

vio.test(120001:140000, "Constant Value 0.25 (n=20)")
vio.test(140001:160000, "Constant Value 0.50 (n=20)")
vio.test(160001:180000, "Constant Value 0.75 (n=20)")
vio.test(180001:200000, "Constant Value 1.00 (n=20)")

vio.test(220001:140000, "Constant Value 0.25 (n=40)")
vio.test(240001:260000, "Constant Value 0.50 (n=40)")
vio.test(260001:280000, "Constant Value 0.75 (n=40)")
vio.test(280001:300000, "Constant Value 1.00 (n=40)")

vio.test(320001:340000, "Constant Value 0.25 (n=100)")
vio.test(340001:360000, "Constant Value 0.50 (n=100)")
vio.test(360001:380000, "Constant Value 0.75 (n=100)")
vio.test(360001:380000, "Constant Value 1.00 (n=100)")

dev.off()

########## OTHER PLOTS ######################

cov.mat <- rep(p.values.H0$cov.mat,3)
constant <- rep(p.values.H0$constant,3)
n <- rep(p.values.H0$n,3)
test <- rep(c("Box's M", " Hawkins' Test", "Non-Parametric"), rep(length(p.values.H0$Box.p),3))
results <- c(p.values.H0$Box.p, p.values.H0$Hawkins.p, p.values.H0$Adk.p)
rp.values <- data.frame(cov.mat, constant, n, test, results)

sam.size.10 <- subset(rp.values, rp.values$n == 10)

bwplot(test ~ results, sam.size.10,
       panel = function(..., box.ratio) {
           panel.violin(..., col = "transparent",
                        varwidth = T, box.ratio = box.ratio)
           panel.bwplot(..., fill = NULL, box.ratio = .1)
       } )

