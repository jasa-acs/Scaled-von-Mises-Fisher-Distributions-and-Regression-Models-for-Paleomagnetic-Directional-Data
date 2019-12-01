# Scaled von Mises–Fisher Distributions and Regression Models for Paleomagnetic Directional Data


## Data

The datafile EIFdata3 is a subset of the archeomagnetic data found in the GEOMAGIA50.v3 online database. For further details see Brown et al. (2015) which is a cited reference in the paper. We extracted the data from the database and the EIFdata3 dataset contains observations on and near the Eifel maars lakes in Germany. 
File format is a csv file. The data items are age in years, dec which is declination angle, inc which is inclination angle, lat which is latitude and lon which is longitude.

## Code

The code calculates parameter estimates in the IID case (see section 4.1 in the paper) and also estimated regression parameters (see section 4.2 in the paper). We also include code for the Kent model.

The code is in R and we include all relevant functions. 


The R file functions contains R functions. This file needs to be run first because all the other files below use these functions.

The R file IID_models calculates parameter estimates in the IID model case. It can be used to create Figure 2 and the case 1 data parameter estimates (see section 5 and 6 of the paper). At the end of the file we also include a section titled Simulation which generates one simulated sample from the fitted a1=1 SvMF model and calculates estimates for this sample. By repeating this simulation code 1000 times the simulation results in Section 6 can be obtained for the a1=1 SvMF model.

The R file regression_models calculates the regression model parameter estimates for the case 2 data. The code also produces Figure 3 and calculates the observed information standard errors for the regression coefficients (see Table 1).  With minor adjustments to the code the case 3 data estimates can be produced. The bootstrap code is not given here, because it is similar to the above code, but with minor adjustments to do repeated sampling. 

The R file figure 1 produces the Figure 1 plot.
