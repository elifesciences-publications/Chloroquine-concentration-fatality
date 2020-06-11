---
title: "Chloroquine concentration-fatality curve"
author: "James Watson"
date: "4/20/2020"
output:
  html_document:
    df_print: paged
    keep_md: yes
  pdf_document: default
---



## Data from self-poisoning cohorts and healthy volunteers

Data from Riou were extracted from published graph (see Figure 3 in NEJM Riou et al, 1988) using WebPlotDigitiser


```r
pooled_data = read.csv('Pooled_QRS_data.csv')
```


```r
ind_clem = pooled_data$study==1
ys = rep(0, nrow(pooled_data)); ys[ind_clem] = rnorm(sum(ind_clem),mean = 0,sd = 1)
par(las=1)
plot(log10(pooled_data$CQ_uMol), ys+pooled_data$QRS, xaxt='n',
     col = pooled_data$study, ylim=c(60,250), ylab='QRS', 
     xlab='Whole blood chloroquine concentration (umol)', pch = pooled_data$died+1)
axis(1, at = seq(-1, 2, by = 0.5), labels = round(10^seq(-1, 2, by = 0.5),1))
legend('topleft', legend = c('Healthy volunteers (malaria dose)','Clemessy','Megarbane'),
       col = c(3,1,2), lwd=2, inset = 0.03)
abline(h=150, v=1)
```

![](Analysis_QRS_files/figure-html/Sfig_raw_QRS_data-1.png)<!-- -->


## Stan model

Emax type model


```r
Conc_QRS_Emax = "
functions {
  real sigmoid(real log10_conc, real ed50, real log_slope, real max_effect, real min_effect){
    return min_effect + (max_effect-min_effect)/(1 + exp(exp(log_slope)*(log10_conc-ed50)));
  }
}

data {
  int<lower=0> N;
  real log10_conc[N];
  real QRS[N];
  int<lower=0,upper=1> study[N];
  real ed50_mu;
  real max_effect_prior_mu;
  real max_effect_prior_sd;
  real min_effect_prior_mu;
  real min_effect_prior_sd;
  real log_slope_prior_mu;
  real log_slope_prior_sd;
}

parameters {
  real log_slope;
  real min_effect;  
  real max_effect;
  real bias_term;
  real ed50;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
}

model {
  bias_term ~ normal(-20,10);
  log_slope ~ normal(log_slope_prior_mu, log_slope_prior_sd);
  ed50 ~ normal(ed50_mu, 1);
  max_effect ~ normal(max_effect_prior_mu,max_effect_prior_sd);
  min_effect ~ normal(min_effect_prior_mu,min_effect_prior_sd);
  
  sigma1 ~ normal(10,2);
  sigma2 ~ normal(20,5);
  
  for (j in 1:N){
    real QRS_pred;
    real sigma;
    if(study[j] == 3) { sigma = sigma1; } else { sigma = sigma2; }
    QRS_pred = sigmoid(log10_conc[j], ed50, log_slope, max_effect, min_effect);
    QRS[j] ~ normal(QRS_pred + bias_term*study[j], sigma);
  }
}
"
if(RUN_MODELS) conc_QRS_mod = stan_model(model_code = Conc_QRS_Emax)
```

## Fit main model to prospective data





```r
sigmoid = function(log10_conc, ed50, log_slope, max_effect, min_effect){
  return (min_effect + (max_effect-min_effect)/(1 + exp(exp(log_slope)*(log10_conc-ed50))))
}
load('mod_QRS_full.stanout')
thetas = extract(mod_QRS_full)


xs = seq(-1,2,length.out = 100)
ys = array(dim = c(100,length(thetas$log_slope)))
for(i in 1:length(thetas$log_slope)){
  ys[,i]=sigmoid(log10_conc = xs, ed50 = (thetas$ed50[i]),log_slope = (thetas$log_slope[i]),
        max_effect = (thetas$max_effect[i]), min_effect = (thetas$min_effect[i]))
}

col_study = RColorBrewer::brewer.pal(n = 3, name = 'Set1')
par(las=1, bty='n', family='serif')
jitter_QRS = rnorm(nrow(pooled_data), mean = , sd = 1*as.numeric(pooled_data$study==1))
plot(x = log10(pooled_data$CQ_uMol), panel.first = grid(),
    y = pooled_data$QRS - mean(thetas$bias_term)*as.numeric(pooled_data$study==1) + jitter_QRS, 
     xaxt='n', col = col_study[pooled_data$study], ylab='QRS (msec)', 
     xlab=expression(paste('Whole blood chloroquine concentration (',mu,'mol/L)')), 
    pch = pooled_data$died+1)
axis(1, at = seq(-1, 2, by = 0.5), labels = round(10^seq(-1, 2, by = 0.5),1))
legend('topleft',  pch = c(1,1,1,1,2),
       legend = c('Healthy volunteers (malaria dose)','Clemessy',
                  'Megarbane','Survivors','Fatal cases'),
       col = c(col_study[c(3,1,2)],'black','black'), inset = 0.03)
qs=quantile(pooled_data$QRS[pooled_data$CQ_uMol==0], probs = c(0.025,0.975,0.5))
polygon(x = c(-3,3,3,-3),y = c(qs[1],qs[1],qs[2],qs[2]), 
        col = adjustcolor('grey',alpha.f = .4),border = NA)

abline(h=150, v=1, lty=2)
lines(xs,rowMeans(ys),lwd=3)
lines(xs,apply(ys,1,quantile,probs=0.025),lwd=1,lty=2)
lines(xs,apply(ys,1,quantile,probs=0.975),lwd=1,lty=2)


abline(h=qs[3], lty=1, col='grey',lwd=3)

points(log10(pooled_data$CQ_uMol), 
     pooled_data$QRS - mean(thetas$bias_term)*as.numeric(pooled_data$study==1)+jitter_QRS, 
      col = col_study[pooled_data$study],  pch = pooled_data$died+1)
```

![](Analysis_QRS_files/figure-html/QRS_fit-1.png)<!-- -->



```r
par(mfrow=c(3,2),las=1, cex.lab=1.5)
hist(thetas$max_effect, freq = F,breaks = 50,main = '', 
     xlab = 'Mean lower QRS value', col = 'grey',ylab='', yaxt='n')
xs = seq(min(thetas$max_effect), max(thetas$max_effect), length.out=500)
lines(xs, dnorm(xs, mean = prior_params$max_effect_prior_mu, 
                sd = prior_params$max_effect_prior_sd), lwd=3, col='red')

hist(thetas$min_effect, freq = F,breaks = 50,main = '', 
     xlab = 'Mean upper QRS value', col = 'grey',ylab='', yaxt='n')
xs = seq(min(thetas$min_effect), max(thetas$min_effect), length.out=500)
lines(xs, dnorm(xs, mean = prior_params$min_effect_prior_mu, 
                sd = prior_params$min_effect_prior_sd), lwd=3, col='red')

hist(thetas$bias_term,freq=F,breaks = 50,main = '', xlab = 'Bias term (Clemessy)', 
     col = 'grey',ylab='', yaxt='n')
xs = seq(min(thetas$bias_term), max(thetas$bias_term), length.out=500)
lines(xs, dnorm(xs, mean = -20, sd = 10), lwd=3, col='red')

hist(thetas$ed50, freq=F,breaks = 50,main = '', xlab = 'ED_50', col = 'grey',
     ylab='', yaxt='n')
xs=seq(0,2,length.out = 2000); 
lines(xs, dnorm(xs,mean = prior_params$ed50_mu,sd = 1),col='red',lwd=3)

hist(thetas$sigma1, freq=F,breaks = 50,main = '', xlab = 'Sigma 1', col = 'grey',
     ylab='', yaxt='n')
xs=seq(0,1,length.out = 500); lines(xs, dnorm(xs,mean=0.7, sd=0.065),col='red',lwd=3)

hist(thetas$sigma2, freq=F,breaks = 50,main = '', xlab = 'Sigma 2', col = 'grey',
     ylab='', yaxt='n')
xs=seq(0,1,length.out = 500); lines(xs, dnorm(xs,mean=0.7, sd=0.065),col='red',lwd=3)
```

![](Analysis_QRS_files/figure-html/priors_vs_posteriors-1.png)<!-- -->

