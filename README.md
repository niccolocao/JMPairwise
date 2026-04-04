# JMPairwise
The R-scripts are related to the under-review article ...
In particular, the code implements the pairwise-fitting estimation (Fieuws and Verbeke, 2006) for joint modelling ordinal, continuous, and time-to-event outcomes: (i) the ordinal outcomes are modelled through a cumulative probit mixed model, (ii) the continuous outcome by a Gaussian linear mixed model, and (iii) the time-to-event by a Weibull proportional hazards model with a univariate log-normal frailty. 

# Quick description
 - ord+ord.R: functions for fitting (in parallel) the pairs of ordinal outcomes
 - ord+cont.R: functions for fitting (in parallel) the pairs of one ordinal outcome and one continuous outcome
 - ord+time.R: functions for fitting (in parallel) the pairs of one ordinal outcome and one time-to-event outcome
 - cont+time.R: functions for fitting (in parallel) the pairs of one continuous outcome and one time-to-event outcome
 - univariateFuns.R: functions for unvariate GLMMs models used for starting values in pairwise-fitting estimation
