# Chloroquine-concentration-fatality

This code does two things, bundled up as an RMarkdown file. 

First it fits a Bayesian logistic regression model to pooled patient data from self-poisoning studies (see pooled_data.csv). The final model uses all the prospectively studied patients (n=258) as there is significant bias in the retrospectively studies patients. It also uses whole blood chloroquine concentrations measured on admission (no change when using observed peak concentrations in individuals with multiple measurements).

Second, it uses output Cmax distributions from a standard pharmacokinetic model (written in NONMEM) to estimate the risk of fatal overdose under a set of chloroquine treatment regimens (5 regimens for COVID-19 and one malaria regimen).

## Sensitivity analyses

We've done a couple of sensitivity analyses to both the PK and PD models.
