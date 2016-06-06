# Load cpt source file
source('CPT_v1.r')

# Load record
xample = read.csv("CPT_example.csv", header=T)


# Run Change-point analysis function
proxy.cpt(Data=xample, Name="xample")
