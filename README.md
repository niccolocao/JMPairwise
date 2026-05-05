# JMPairwise
The R-scripts are related to the under-review article "Joint modeling of high-dimensional mixed-type longitudinal and survival data: An application to an orthopedic dataset with dynamic predictions"

In particular, the code implements the pairwise-fitting estimation (Fieuws and Verbeke, 2006) for joint modelling ordinal, continuous, and time-to-event outcomes with dynamic predictions: 
    (i) the ordinal outcomes are modelled through a cumulative probit mixed model, 
    (ii) the continuous outcome by a Gaussian linear mixed model, 
    (iii) the time-to-event by a Weibull proportional hazards model with a univariate log-normal frailty. 
    Note that the longitudinal outcomes have correlated errors across time points.

The code uses GLMMadaptive package for estimation which can be downloaded for the used version by `remotes::install_github("drizopoulos/GLMMadaptive@36b1b26")`

# Quick description
### R-scripts for fitting bivariate models:
 - ord+ord.R: functions for fitting (in parallel) the pairs of ordinal outcomes
 - ord+cont.R: functions for fitting (in parallel) the pairs of one ordinal outcome and one continuous outcome
 - ord+time.R: functions for fitting (in parallel) the pairs of one ordinal outcome and one time-to-event outcome
 - cont+time.R: functions for fitting (in parallel) the pairs of one continuous outcome and one time-to-event outcome
 - univariateFuns.R: functions for unvariate GLMMs models used for starting values in pairwise-fitting estimation
   
### R-scripts for the analysis:
 - Application.R: code used in the application
 - pairwiseUtilities.R: functions for extracting the bivariate model estiamtes, averaging the overlapping estiamtes, obtaining variance-covariance matrix of the estiamtes, and applying delta method
 - genericUtilities.R: functions for Wald tests, tables, and plots
 - dynpred.R: functions for dynamic predictions
 - dynplot.R: functions for plotting dynamic predictions
