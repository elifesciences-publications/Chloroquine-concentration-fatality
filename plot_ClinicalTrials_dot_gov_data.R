dat = readxl::read_excel('../../Data/trial_doses.xlsx')

pdf('~/Dropbox/Apps/Overleaf/Chloroquine dose-response/ClinicalTrialgov.pdf')
par(las=1, family = 'serif', bty='n')
plot(dat$`Total Duration`,dat$`Total dose (base)`/1000, col = as.numeric(dat$Drug=='HCQ')+1,
     xlab='Duration of treatment (days)', ylab='Total dose in base equivalent (g)', pch=16)
legend('bottomright', col=1:2, pch=16,legend = c('Chloroquine','Hydroxychloroquine'),inset=0.03,bty='n')

text(x = 10, y=12.1,labels = 'Borba et al')
text(x = 14, y=12.6,labels = 'PATCH',col=2)

text(x = 10, y=8.9,labels = 'NCT04351620',col=2)
text(x = 10+1.6, y=7.13,labels = 'RECOVERY',col=2)

abline(h = 6.82, lty=2)
dev.off()
