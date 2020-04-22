# Chloroquine-concentration-fatality

This code does two things. First it fits a Bayesian logistic regression model to pooled patient data from self-poisoning studies (see pooled_data.csv). The final model uses all the prospectively colected data (n=258).
Second, it uses output Cmax distributions from a standard pharmacokinetic model (written in NONMEM) to estimate the risk of fatal overdose under a set of chloroquine treatment regimens (5 regimens for COVID-19 and one malaria regimen).

