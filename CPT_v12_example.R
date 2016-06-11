rm(list = ls())

# Set working directory
setwd("/Users/wfinsing/Documents/R_packages/CHAR_Changepoint/CPT_v12")

# Load packages
require(changepoint)

# Load cpt and pretreat.full source files
source('cpt_v12.r')
source('pretreatment_full.r')

# Load a random charcoal dataset already formatted ####
# for the pretreat.full function
pip <- read.csv("random_CharConc.csv", header=T)


# Interpolate records with pretreat.full function ####
dat <- pip
dat.i <- pretreat.full(params=dat[ ,c(1:5)], serie=dat$CharConc, yrInterp=40, plotit=T,Int=T)

dat.i <- as.data.frame(dat.i)
dat.i <- dat.i[which(complete.cases(dat.i)), ] # omits rows with NAs at bottom of data.frame

# Check interpolated CHAR record with a plot
x.lim <- c(max(dat.i$AgeI), min(dat.i$AgeI))
par(mfrow=c(3,1))
par(mar = c(0.5,5,0.5,1))
par(oma = c(5,1,1,1), cex=0.7)
plot(dat.i$AgeI, dat.i$ConcI, xlim=x.lim, type="s", xaxt="n", ylab="random char.\nconc. interp.")
plot(dat.i$AgeI, dat.i$ArI, xlim=x.lim, type="s", xaxt="n", ylab="random CHAR\ninterpolated")
plot(dat.i$AgeI, dat.i$sedArI, xlim=x.lim, type="s", ylab="sediment AR\n(cm/yr)")

# Analyse the record with the change-point analysis function "proxy.cpt()" ####
proxy.cpt(serie=dat.i, Name="Pip_NoBoot")
proxy.cpt(serie=dat.i, Name="Pip_NoBoot_Penalty3log", bootstrap = F, penalty = "3*log(n)")

toRm = c("dat", "dat.i", "pip", "x.lim")
rm(list = toRm)




# Creates a random charcoal concentration record ("rand.Conc") for an age-depth model record ####
# Here start with a known age model, and known min and max charcoal concentrations
SeriesAge = read.table("Ech_interpolated_ages.txt", header=T)
SeriesAge <- SeriesAge[ ,-c(2:3)]
EchSar = 1/SeriesAge$accrate

Data.length = length(SeriesAge[,1])
Conc.min = 0
Conc.max = 0.6 * 60    #mode of max values * yr/cm as from Rius et al. (VHA)
rand.Conc = data.frame(seq(1:(100*Data.length)), runif((100*Data.length), Conc.min, Conc.max))
colnames(rand.Conc) [2] = "RandConc"

# Samples randomly with replacement from rand.Conc
Conc.rand <- sample(rand.Conc$RandConc, size=Data.length, replace=T)

# Merge random char dataset with age model
pip1 <- data.frame(SeriesAge, Conc.rand)

# Formats data.frame for pretreat.full function
pip1 <- pip1[ ,-3]
colnames(pip1) <- c("CmTop", "AgeTop",  "CharConc") # rename columns
n.rows <- length(pip1$CmTop)
pip1$CmBot <- pip1$CmTop + 1
age.bot <- pip1[2:n.rows, 2]
pip1$AgeBot <- c(age.bot, NA)
pip1$Vol <- 1
pip1 <- pip1[ ,c(1,4,2,5,6,3)] # order columns for paleofire package
pip1 <- pip1[which(complete.cases(pip1)), ] # removes last row where AgeBot = NA


# Interpolate records with pretreatment function ####
dat <- pip1

dat.i <- pretreat.full(params=dat[ ,c(1:5)], serie=dat$CharConc, yrInterp=40, plotit=T,Int=T)

dat.i <- as.data.frame(dat.i)
dat.i <- dat.i[which(complete.cases(dat.i)), ] # omits rows with NAs

# Check interpolated CHAR record with a plot
x.lim <- c(max(dat.i$AgeI), min(dat.i$AgeI))
par(mfrow=c(3,1))
par(mar = c(0.5,5,0.5,1))
par(oma = c(5,1,1,1), cex=0.7)
plot(dat.i$AgeI, dat.i$ConcI, xlim=x.lim, type="s", xaxt="n", ylab="random\nchar conc.")
plot(dat.i$AgeI, dat.i$ArI, xlim=x.lim, type="s", xaxt="n", ylab="random CHAR\ninterpolated")
plot(dat.i$AgeI, dat.i$sedArI, xlim=x.lim, type="s", ylab="sediment AR\n(cm/yr)")


# Analyse the record with the change-point analysis function "proxy.cpt()" ####
proxy.cpt(serie=dat.i, Name="Pip1_NoBoot")
proxy.cpt(serie=dat.i, Name="Pip1_NoBoot_Penalty3log", bootstrap = F, penalty = "3*log(n)")

