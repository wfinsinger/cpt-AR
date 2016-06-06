# cpt-AR (Change-point analysis with single influx records)

-------------------------------------------------------------------------------------
# Determines zone boundaries for single influx records using the change-point analysis
as described in Finsinger et al. (2016).

--
The function requires one input file with five columns:

Column 1: depth          =   sample depth

Column 2: age            =   sample age (as cal yrs BP)

Column 3: sed.AR.Data    =   sediment accumulation rate as from age-depth model (as cm year-1)

Column 4: Conc           =   concentration (as pieces cm-3)

Column 5: AR             =   influx (or accumulation rate) value (as pieces cm-2 yr-1)

NB: the data has to be interpolated to a constant temporal resolution!!

--
The function also requires following additional parameters
(which can be left at default values)
 Name      =   Site name
 
 bootstrap =   if FALSE the random dataset is generated with the runif() function
 
 q         =   number of random datasets generated to determine 
 
 n.Q       =   maximum number of change points
 
 n.screen  =   a change point in the random datasets is validated if it occurs in more than
                n.screen datasets. By default n.screen = q * 0.025 (thus with q=1000 this equals
                2.5% chance of occurrence) 

--
See the example.csv file as a template for data input.

-------------------------------------------------------------------------------------
Suggested citation: Finsinger W., Magyari E.K., Fevre J., Orban I., Pal I., Vincze I., Hubay K,
                     Birks H.H., Braun M., Toth M.  (2016) – Holocene fire regimes near the treeline
                     in the Retezat Mts. (Southern Carpathians). Quaternary International.
                     doi: 10.1016/j.quaint.2016.04.029. In press
