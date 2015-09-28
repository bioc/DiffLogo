
library(MotifDb)
library(seqLogo)
library(DiffLogo)

hitIndeces <- grep ('CTCF', values (MotifDb)$geneSymbol, ignore.case=TRUE)
list    <- as.list(MotifDb[hitIndeces])

pwm1 <- list$"Hsapiens-JASPAR_CORE-CTCF-MA0139.1"[4:1, 18:2] ## trim and reverse complement
pwm2 <- list$"Hsapiens-jolma2013-CTCF"

par(mfrow=c(2,1), pin=c(3, 1), mar = c(2, 4, 1, 1))
seqLogo::seqLogo(pwm = pwm1)
seqLogo::seqLogo(pwm = pwm2)
par(mfrow=c(1,1), pin=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1))

par(mfrow=c(2,1), pin=c(3, 1), mar = c(2, 4, 1, 1))
seqLogo(pwm = pwm1)
seqLogo(pwm = pwm2, stackHeight = sumProbabilities)
par(mfrow=c(1,1), pin=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1))

#par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
diffLogoFromPwm(pwm1, pwm2)
