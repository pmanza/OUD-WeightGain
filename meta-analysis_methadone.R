library("robumeta")
library("metafor")
library("dplyr")

# this code was adapted from examples from Daniel S. Quintana, from the following paper:
# Quintana, D. S. (2015). From pre-registration to publication: a non-technical primer for conducting a meta-analysis to synthesize correlational data. Frontiers in psychology, 1549.


# NOTE: Make sure to change your working directory to wherever this script and the .csv file are located! 

# import data
dat <- read.csv("FinalDataMethadone.csv", header=TRUE)  

# calculate effect sizes for the studies that include Means and SDs for change in BMI
dat2 <- escalc(measure="SMD", m1i=Mean_BMI_Change_Female, sd1i=SD_BMI_Change_Female, n1i=N_Female,
               m2i=Mean_BMI_Change_Male, sd2i=SD_BMI_Change_Male, n2i=N_Male, data=dat)

# and calculate it for the one study that didn't have SDs but did have p-value
dat2$tVal <- replmiss(dat2$tVal, with(dat2, qt(pVal/2, df=N_Male+N_Female-2, lower.tail=FALSE)))

# calculate cohens d
dat2$dVal <- replmiss(dat2$dVal, with(dat2, tVal * sqrt(1/N_Male + 1/N_Female)))

# and now create standardized effect size based on Cohen's D (with bias correction)
dat2$yi <- replmiss(dat2$yi, with(dat2, (1 - 3/(4*(N_Male+N_Female-2) - 1)) * dVal))

#  Any missing values for the vi variable (the sampling variances) can now be replaced with: 
dat2$vi <- replmiss(dat2$vi, with(dat2, 1/N_Male+ 1/N_Female + yi^2/(2*(N_Male+N_Female))))


# perform the meta-analysis using a random-effects model. The following commands will print out the data and also calculates and print the confidence interval for the amount of heterogeneity (I^2).

res <- rma(yi, vi, data=dat2) 
res 
predict(res, digits=3, transf=transf.ztor)
confint(res)  



# While the Q-statistic and I^2 can provide evidence for heterogeneity, they do not provide information on which studies may be influencing to overall heterogeneity. If there is evidence of overall heterogeneity, construction of a Bajaut plot can illustrate studies that are contribute to overall heterogeneity and the overall result. Study IDs are used to identify studies

b_res <- rma(yi, vi, data=dat2, slab=StudyID)  # New meta-analysis with study ID identifier  

# The next command will plot a Baujat plot.

baujat(b_res)

# Studies that fall to the top right quadrant of the Baujat plot contribute most to both these factors. 
# A set of diagnostics are also available to identify potential outliers and influential cases.

inf <- influence(res)
print(inf)
plot(inf)

# Now we visualize the meta-analysis with a forest plot. 

forest(res, slab=paste(dat2$Author, dat2$Year, sep=", "), xlim = c(-16, 8),
       ilab=cbind(dat2$Age, dat2$Dose..mg, dat2$Treatment.Duration..Days.),
       ilab.xpos=c(-8,-5,-2),
       digits=c(2,1), cex=.75, header = TRUE)

### add text with Q-value, dfs, p-value, and I^2 statistic
text(-13, -.98, pos=4, cex=0.75, bquote(paste("(Q(", .(res$k - res$p), ") = ", .(formatC(res$QE, digits=2, format="f")),
                                            ", p = ", .(formatC(res$QEp, digits=3, format="f")), "; ", I^2, " = ",
                                            .(formatC(res$I2, digits=1, format="f")), "%)")))


