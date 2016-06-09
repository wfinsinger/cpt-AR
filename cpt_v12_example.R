rm(list = ls())

# Set working directory
setwd("/Users/wfinsing/Documents/R_packages/CHAR_Changepoint/CPT_v12")

# Load cpt source file
source('cpt_v12.r')
source('pretreatment_full.r')

# Load a random charcoal dataset ####
pip <- read.csv("random_charcconc.csv", header=T)

# Creates random charcoal conc record for an age-depth model record
# Here start with a known age model, and known min and max charcoal concentrations

EchAge = read.table("Ech_interpolated_ages.txt", header=T)
EchAge <- EchAge[ ,-c(2:3)]
EchSar = 1/EchAge$accrate

Data.length = length(EchAge[,1])
Conc.min = 0
Conc.max = 0.6 * 60    #mode of max values * yr/cm as from Rius et al. (VHA)
rand.Conc = data.frame(seq(1:(100*Data.length)), runif((100*Data.length), Conc.min, Conc.max))
colnames(rand.Conc) [2] = "RandConc"

# Samples random data
Conc.rand <- sample(rand.Conc$RandConc, size=Data.length, replace=T)

CHAR.rand <- Conc.rand * EchSar

par(mfrow=c(3,1))
par(mar = c(0.5,5,0.5,1))
par(oma = c(5,1,1,1), cex=0.7)
plot(EchAge$best, Conc.rand, type="l", xaxt="n", ylab="random\nchar conc.")
plot(EchAge$best, CHAR.rand, type="l", xaxt="n", ylab="random\nCHAR")
plot(EchAge$best, EchSar, type="l", ylab="Sed. AR")

# Clean environment
toRm = c("rand.Conc", "CHAR.rand", "Conc.max", "Conc.min", "EchSar")
rm(list = toRm)

#write()

# Merge random char dataset with age model
pip1 <- data.frame(EchAge, Conc.rand)

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

# Modify pip in order to be able to work with the pretreatment function from Phil
pip <- pip[ ,-c(3,5)] # delete useless columns
colnames(pip) <- c("CmTop", "AgeTop", "CharConc") # rename columns

n.rows <- length(pip$CmTop)
pip$CmBot <- pip$CmTop + 1
age.bot <- pip[2:n.rows, 2]
pip$AgeBot <- c(age.bot, NA)
pip$Vol <- 1
pip <- pip[ ,c(1,4,2,5,6,3)] # order columns for paleofire package
pip <- pip[-n.rows, ] # removes last row where AgeBot = NA

# Clean environment
toRm <- c("EchAge", "age.bot", "Conc.rand", "Data.length", "n.rows")
rm(list = toRm)



# Interpolate records with pretreatment function ####

dat <- pip

dat.i <- pretreat.full(params=dat[ ,c(1:5)], serie=dat$CharConc, yrInterp=40, plotit=T,Int=T)

dat.i <- as.data.frame(dat.i)
dat.i <- dat.i[which(complete.cases(dat.i)), ] # omits rows with NAs

# Check interpolated CHAR record with a plot
x.lim <- c(max(dat.i$AgeI), min(dat.i$AgeI))
par(mfrow=c(3,1))
par(mar = c(0.5,5,0.5,1))
par(oma = c(5,1,1,1), cex=0.7)
plot(dat.i$AgeI, dat.i$ConcI, xlim=x.lim, type="s", xaxt="n", ylab="random char.\nconc. interp.")
plot(dat.i$AgeI, dat.i$ArI, xlim=x.lim, type="s", xaxt="n", ylab="random CHAR\ninterpolated")
plot(dat.i$AgeI, dat.i$sedArI, xlim=x.lim, type="s", ylab="sediment AR\n(cm/yr)")

proxy.cpt(serie=dat.i, Name="Pip_NoBoot")
proxy.cpt(serie=dat.i, Name="Pip_NoBoot_Penalty3log", bootstrap = F, penalty = "3*log(n)")

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

proxy.cpt(serie=dat.i, Name="Pip1_NoBoot")
proxy.cpt(serie=dat.i, Name="Pip1_NoBoot_Penalty3log", bootstrap = F, penalty = "3*log(n)")


