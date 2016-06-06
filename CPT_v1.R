#----------------------------------------------------------------------------------------#
#  Determines zone boundaries for single influx records using the change-point analysis  #
#  as described in Finsinger et al. (2016).                                              #
#----------------------------------------------------------------------------------------#
# The function requires one input file with five columns:
#    Column 1: depth          =   sample depth
#    Column 2: age            =   sample age (as cal yrs BP)
#    Column 3: sed.AR.Data    =   sediment accumulation rate as from age-depth model (as cm year-1)
#    Column 4: Conc           =   concentration (as pieces cm-3)
#    Column 5: AR             =   influx (or accumulation rate) value (as pieces cm-2 yr-1)
#
# NB: the data has to be interpolated to a constant temporal resolution!!
#
# The function also requires following additional parameters
#  (which can be left at default values)
#    Name      =   Site name
#    bootstrap =   if FALSE the random dataset is generated with the runif() function
#    q         =   number of random datasets generated to determine 
#    n.Q       =   maximum number of change points
#    n.screen  =   a change point in the random datasets is validated if it occurs in more than
#                   n.screen datasets. By default n.screen = q * 0.025 (thus with q=1000 this equals
#                   2.5% chance of occurrence) 
#
#
# Suggested citation: Finsinger W., Magyari E.K., Fevre J., Orban I., Pal I., Vincze I., Hubay K,
#                     Birks H.H., Braun M., Toth M.  (2016) – Holocene fire regimes near the treeline
#                     in the Retezat Mts. (Southern Carpathians). Quaternary International.
#                     doi: 10.1016/j.quaint.2016.04.029. In press
#

# ------ Defines FUNCTION ------------
proxy.cpt = function(Data, Name, bootstrap=F, q=1000, n.Q=10, n.screen=q*0.025,
                     output.dir=file.path(".","cpt_output"))
{
  

# ------------ SETUP ------------ #
  
  # -------- Load required libraries ----------- #
  require(changepoint)

  # -------- Create output directory
  dir.create(output.dir)
  
  
# ---------- Defines default parameters
#   bootstrap = F
#   resol = 40
#   q = 1000
#   n.Q = 10
#   n.screen = q*0.025

# ---------- Settings giving relatively robust results with Lia and Brazi CHARc and 'runif' function
  test=cpt.meanvar
  meth.cpt="BinSeg"
  pen="Manual"
  pen.val="4*log(n)"
  t.stat="Normal"


# --------- First sets some values useful in the following

Data.Conc = Data[,4]
sed.AR.Data = Data[,3]
Data.length = length(Data.Conc)
Data.min = min(Data.Conc)
Data.max = max(Data.Conc)

if(bootstrap) {
  
  ## CPT analysis with bootstrapped Concentration dataset ####
  
  cpts.rand.constant = matrix() # stores cpts from rand series
  cpts.rand.model = matrix() # stores cpts from rand series with age-depth model from Data
  
  ## Loop cpt analyses through rand datasets
  for (i in 1:q) {
    Conc.rand = sample(Data.Conc, size=Data.length, replace=T)
    
    ### Merge rand dataset with Data
    rand = cbind(Data, Conc.rand)
    # plot(rand[,4], rand[,6], type="l", main=paste(Name, "- randised Conc"), xlab="Data Conc", ylab="rand Conc")
    
    ## Calculate randAR with linear age-depth model and with records' age-depth model
    
    sed.AR.median = median(rand$sed.AR.Data)
    
    rand$AR.rand.constant = rand$Conc.rand * sed.AR.median
    rand$AR.rand.model = rand$Conc.rand * sed.AR.Data
    
    #plot(rand$age, rand$AR.rand.constant, type="l")
    #lines(rand$age, rand$AR.rand.model, col="red")
    #lines(rand$age, rand$AR, col="blue")
    
    ### Run change-point analysis with rand AR datasets
    #cpt.AR.rand = cpt.mean(rand$AR.rand, method=meth.cpt)
    cpt.AR.rand.constant = test(rand$AR.rand.constant, method=meth.cpt, Q=n.Q, test.stat=t.stat, penalty=pen, pen.value=pen.val)
    cpts.rand.constant = c(cpts.rand.constant, cpts(cpt.AR.rand.constant))
    
    #cpt.AR.Data.rand = cpt.mean(rand$AR.Data.rand, method=meth.cpt)
    cpt.AR.rand.model = test(rand$AR.rand.model, method=meth.cpt, Q=n.Q, test.stat=t.stat, penalty=pen, pen.value=pen.val)
    cpts.rand.model = c(cpts.rand.model, cpts(cpt.AR.rand.model))
  }
  
  ### Run cpt analysis with true CHAR record
  cpt.AR.Data = test(rand$AR, method=meth.cpt, Q=n.Q, test.stat=t.stat, penalty=pen, pen.value=pen.val)
  cpts(cpt.AR.Data)
  
  # screens cpts from rand series
  cpts.rand.screen = data.frame(tabulate(cpts.rand.constant), 1.2*max(rand$AR))
  colnames(cpts.rand.screen) = c("cpt.count", "X1")
  cpts.rand.screen$id = seq(1:(length(cpts.rand.screen[,1])))
  rand.screen = cpts.rand.screen[ which(cpts.rand.screen[,1]>n.screen), ]
  dummy = matrix(data=NA, nrow=1, ncol=3, dimnames=list("1", c("cpt.count", "X1", "id"))) # adds a 0 line to plot in case of no cpt occurrence
  rand.screen = rbind(rand.screen, dummy)
  
  cpts.Data.rand.screen = data.frame(tabulate(cpts.rand.model), 1.2*max(rand$AR))
  colnames(cpts.Data.rand.screen) = c("cpt.count", "X1")
  cpts.Data.rand.screen$id = seq(1:(length(cpts.Data.rand.screen[,1])))
  Data.rand.screen = cpts.Data.rand.screen[ which(cpts.Data.rand.screen[,1]>n.screen), ]
  dummy = matrix(data=NA, nrow=1, ncol=3, dimnames=list("1", c("cpt.count", "X1", "id"))) # adds a 0 line to plot in case of no cpt occurrence
  Data.rand.screen = rbind(Data.rand.screen, dummy)
  
  
  ## Diagnostic plot with true CHAR record and cpts
  pdf(file.path(output.dir, paste(Name,".pdf", sep="")))
  par(mfrow=c(2,1))
  par(mar = c(0.5,5,0.5,1))
  par(oma = c(3,1,2,1), cex=0.8)
  x.lim.cpt = rev(range(c(1, Data.length)))
  x.lim.deptime = rev(range(c(min(rand$age), max(rand$age))))
  y.lim=c(0, 1.4*max(rand$AR))
  plot(cpt.AR.Data, type="s", ylab="CHAR", xaxt="n", xlim=x.lim.cpt, main=Name)
  abline(v=c(cpts(cpt.AR.Data)), lwd=2, col="blue")
  par(new=T)
  plot(rand.screen$id, rand.screen$X1, ylab="", yaxt="n", xaxt="n", pch=1, xlim=x.lim.cpt, ylim=y.lim)
  par(new=T)
  plot(Data.rand.screen$id, Data.rand.screen$X1, ylab="", yaxt="n", xaxt="n", pch=1, xlim=x.lim.cpt, ylim=y.lim)
  plot(rand$age, rand$sed.AR.Data, type="l", ylab="Sediment-accumulation\nrate (cm/yr)", xlim=x.lim.deptime)
  mtext(("Age (cal yrs BP)"), side = 1, line = 2.5, cex=0.9)
  dev.off()
  
  } else {  
  
  
  ## CPT Analysis with rand CharConcentration dataset with 'runif' function ####
  rand.Data = data.frame(seq(1:(100*Data.length)), runif((100*Data.length), Data.min, Data.max))
  colnames(rand.Data) [2] = "Conc.rand"
  
  ## Loop cpt analyses through rand datasets
  cpts.rand.constant = matrix() # stores cpts from rand series
  cpts.rand.model = matrix() # stores cpts from rand series with age-depth model from Data
  
  for (i in 1:q) {
    Conc.rand = sample(rand.Data$Conc.rand, size=Data.length, replace=T)
    # plot(rand.Data[,1], rand.Data[,2], type="l", main="randised_charData", xlab="Age", ylab="rand CHAR#")
    
    ### Merge rand dataset with Data
    rand = cbind(Data, Conc.rand)
    # plot(rand[,4], rand[,6], type="l", main=paste(Name, "- randised Conc"), xlab="Data Conc", ylab="rand Conc")
    
    ## Calculate randAR with linear age-depth model and with records' age-depth model
    
    sed.AR.median = median(sed.AR.Data)
    
    rand$AR.rand.constant = rand$Conc.rand * sed.AR.median
    rand$AR.rand.model = rand$Conc.rand * sed.AR.Data
    
    #plot(rand$age, rand$AR.rand.constant, type="l")
    #lines(rand$age, rand$AR.rand.model, col="red")
    #lines(rand$age, rand$AR, col="blue")
    
    
    ### Run change-point analysis with rand CHAR datasets
    cpt.AR.rand.constant = test(rand$AR.rand.constant, method=meth.cpt, Q=n.Q, test.stat=t.stat, penalty=pen, pen.value=pen.val)
    cpts.rand.constant = c(cpts.rand.constant, cpts(cpt.AR.rand.constant))
    
    #cpt.AR.Data.rand = cpt.mean(rand$AR.Data.rand, method=meth.cpt)
    cpt.AR.rand.model = test(rand$AR.rand.model, method=meth.cpt, Q=n.Q, test.stat=t.stat, penalty=pen, pen.value=pen.val)
    cpts.rand.model = c(cpts.rand.model, cpts(cpt.AR.rand.model))
  }
  
  ### Run cpt analysis with true CHAR record
  cpt.AR.Data = test(rand$AR, method=meth.cpt, Q=n.Q, test.stat=t.stat, penalty=pen, pen.value=pen.val)
  cpts(cpt.AR.Data)
  
  # screens cpts from rand series
  cpts.rand.screen = data.frame(tabulate(cpts.rand.constant), 1.2*max(rand$AR))
  colnames(cpts.rand.screen) = c("cpt.count", "X1")
  cpts.rand.screen$id = seq(1:(length(cpts.rand.screen[,1])))
  rand.screen = cpts.rand.screen[ which(cpts.rand.screen[,1]>n.screen), ]
  dummy = matrix(data=NA, nrow=1, ncol=3, dimnames=list("1", c("cpt.count", "X1", "id"))) # adds a 0 line to plot in case of no cpt occurrence
  rand.screen = rbind(rand.screen, dummy)
  
  cpts.Data.rand.screen = data.frame(tabulate(cpts.rand.model), 1.2*max(rand$AR))
  colnames(cpts.Data.rand.screen) = c("cpt.count", "X1")
  cpts.Data.rand.screen$id = seq(1:(length(cpts.Data.rand.screen[,1])))
  Data.rand.screen = cpts.Data.rand.screen[ which(cpts.Data.rand.screen[,1]>n.screen), ]
  dummy = matrix(data=NA, nrow=1, ncol=3, dimnames=list("1", c("cpt.count", "X1", "id"))) # adds a 0 line to plot in case of no cpt occurrence
  Data.rand.screen = rbind(Data.rand.screen, dummy)
  
  
  ## Diagnostic plot  with true CHAR record and cpts
  pdf(file.path(output.dir, paste(Name,".pdf", sep="")))
  par(mfrow=c(2,1))
  par(mar = c(0.5,5,0.5,1))
  par(oma = c(5,1,1,1), cex=0.7)
  x.lim.cpt = rev(range(c(1, Data.length)))
  x.lim.age = rev(range(c(min(rand$age), max(rand$age))))
  y.lim=c(0, 1.4*max(rand$AR))
  plot(cpt.AR.Data, type="s", ylab="AR", xaxt="n", xlim=x.lim.cpt, main=Name)
  abline(v=c(cpts(cpt.AR.Data)), lwd=2, col="blue")
  par(new=T)
  plot(rand.screen$id, rand.screen$X1, ylab="", yaxt="n", xaxt="n", pch=1, xlim=x.lim.cpt, ylim=y.lim)
  par(new=T)
  plot(Data.rand.screen$id, Data.rand.screen$X1, ylab="", yaxt="n", xaxt="n", pch=1, xlim=x.lim.cpt, ylim=y.lim)
  plot(rand$age, rand$sed.AR.Data, xlim=x.lim.age, type="l", ylab="Sediment-accumulation\nrate (cm/yr)", xaxt="n")
  at = seq(from=0, to=max(x.lim.age), by=1000)
  axis(side=1, at=at, las=2, hadj=0.9)
  mtext(("Age (cal yrs BP)"), side = 1, line = 3.5, cex=0.8)
  dev.off()
  }

}


# Welcome
cat("Hi there, welcome to Change-point analysis with influx records (v.1.0)\n")
cat(" \n")
cat("The function requires one input file (interpolated to a constant temporal resolution)\n with five columns:\n")
cat("  Column 1: depth          =   sample depth\n")
cat("  Column 2: age            =   sample age (as cal yrs BP)\n")
cat("  Column 3: sed.AR.Data    =   sediment accumulation rate as from age-depth model (as cm year-1)\n")
cat("  Column 4: Conc           =   concentration (as pieces cm-3)\n")
cat("  Column 5: AR             =   influx (or accumulation rate) value (as pieces cm-2 yr-1)\n")