dat = readxl::read_excel('trial_doses.xlsx')

pdf('ClinicalTrialgov.pdf')

par(las=1, family = 'serif', bty='n')
plot(dat$`Total Duration`,dat$`Total dose (base)`/1000, 
     col = as.numeric(dat$Drug=='HCQ')+1, yaxt='n',
     xlab='Duration of treatment (days)', xaxt='n',
     ylab='Total dose in base equivalent (g)', pch=16)
legend('bottomright', col=c(1,2,1), pch=c(16,16,0),
       legend = c('Chloroquine','Hydroxychloroquine','Simulated flat regimens'),
       inset=0.03)
axis(2, at = seq(1,13,by=2), tick = T)
axis(1, at = seq(3, 21, by=3), tick = T)

xs = c(10, 10 , 7, 3)
ys = c(620*2*10, 620*2 + 310*2*9, 620*2 + 310*2*6, 620*2+310)/1000
points(xs, ys, pch = 0, cex=1.5)

text(x = 10, y=12.1,labels = 'Borba et al')
text(x = 14, y=12.6,labels = 'PATCH',col=2)

text(x = 10, y=8.9,labels = 'NCT04351620',col=2)
text(x = 10+2, y=7.44,labels = 'RECOVERY',col=2)

text(x = 10+2, y=6.82,labels = 'SOLIDARITY',col=2)

dev.off()
