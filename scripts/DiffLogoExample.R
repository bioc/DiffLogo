
###########################################################################################################
## built
## 
## cba
## knitr
## 


###########################################################################################################
## tests
seqLogo(pwm = CTCF[["1"]], sparse = FALSE)

###########################################################################################################
## toy example

pwm11 <- unlist(list(1, 0, 0, 0))
pwm12 <- unlist(list(0, 0, 0, 1))

pwm21 <- unlist(list(0.5, 0, 0, 0.5))
pwm22 <- unlist(list(0, 0.5, 0.5, 0))

pwm31 <- unlist(list(0.5, 0, 0, 0.5))
pwm32 <- unlist(list(0.5, 0.5, 0, 0))

e <- 0.2145
pwm41 <- unlist(list(1 - e, 0, 0, e))
pwm42 <- unlist(list(e, 0, 0, 1 - e))

pwm51 <- unlist(list(0.7, 0.2, 0.1, 0.0))
pwm52 <- unlist(list(0.7, 0.2, 0.1, 0.0))

#pwm1 <- unlist(list(0.9, 0, 0, 0.1))
#pwm2 <- unlist(list(0, 0.1, 0.9, 0))

## parameters
widthToHeightRatio = 16/10;
size = 1.1
resolution <- 300
width <- size * widthToHeightRatio
height <- size

## export DiffLogo table as png image
fileName <- "SeqLogo_11.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
par(mar=c(0.5,1.5,0.5,0.5))
seqLogo(pwm = pwm11, sparse = TRUE, alphabet = DNA)
dev.off()
fileName <- "SeqLogo_12.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
par(mar=c(0.5,1.5,0.5,0.5))
seqLogo(pwm = pwm12, sparse = TRUE, alphabet = DNA)
dev.off()
fileName <- "SeqLogo_21.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
par(mar=c(0.5,1.5,0.5,0.5))
seqLogo(pwm = pwm21, sparse = TRUE, alphabet = DNA)
dev.off()
fileName <- "SeqLogo_22.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
par(mar=c(0.5,1.5,0.5,0.5))
seqLogo(pwm = pwm22, sparse = TRUE, alphabet = DNA)
dev.off()
fileName <- "SeqLogo_31.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
par(mar=c(0.5,1.5,0.5,0.5))
seqLogo(pwm = pwm31, sparse = TRUE, alphabet = DNA)
dev.off()
fileName <- "SeqLogo_32.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
par(mar=c(0.5,1.5,0.5,0.5))
seqLogo(pwm = pwm32, sparse = TRUE, alphabet = DNA)
dev.off()
fileName <- "SeqLogo_41.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
par(mar=c(0.5,1.5,0.5,0.5))
seqLogo(pwm = pwm41, sparse = TRUE, alphabet = DNA)
dev.off()
fileName <- "SeqLogo_42.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
par(mar=c(0.5,1.5,0.5,0.5))
seqLogo(pwm = pwm42, sparse = TRUE, alphabet = DNA)
dev.off()
fileName <- "SeqLogo_51.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
par(mar=c(0.5,1.5,0.5,0.5))
seqLogo(pwm = pwm51, sparse = TRUE, alphabet = DNA)
dev.off()
fileName <- "SeqLogo_52.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
par(mar=c(0.5,1.5,0.5,0.5))
seqLogo(pwm = pwm52, sparse = TRUE, alphabet = DNA)
dev.off()

diffLogoObj = createDiffLogoObject(pwm1 = pwm11, pwm2 = pwm12, alphabet = DNA)
diffLogoObj$heights
diffLogoObj = createDiffLogoObject(pwm1 = pwm21, pwm2 = pwm22, alphabet = DNA)
diffLogoObj$heights
diffLogoObj = createDiffLogoObject(pwm1 = pwm31, pwm2 = pwm32, alphabet = DNA)
diffLogoObj$heights
diffLogoObj = createDiffLogoObject(pwm1 = pwm41, pwm2 = pwm42, alphabet = DNA)
diffLogoObj$heights
diffLogoObj = createDiffLogoObject(pwm1 = pwm51, pwm2 = pwm52, alphabet = DNA)
diffLogoObj$heights

diffLogo(diffLogoObj)

e <- 0.2145
p <- unlist(list(1 - e, 0, 0, e))
q <- unlist(list(e, 0, 0, 1 - e))
m = (p + q) / 2
d <- 0.5*sum( p * (log2(p) - log2(m)),na.rm=T) + 0.5*sum( q * (log2(q) - log2(m)),na.rm=T)
d

###########################################################################################################
## vignette
library(DiffLogo)

## import pwms
motif_folder = "pwm";
motif_names = c("HepG2","MCF7","HUVEC","ProgFib");
motifs = list();
for (name in motif_names) {
  file = system.file(paste(motif_folder,"/",name,".txt",sep=""), package = "DiffLogo")
  motifs[[name]] = as.matrix(read.delim(file,header=F));
}

## plot classic sequence logo
pwm1 = motifs[[motif_names[[1]]]]
pwm2 = motifs[[motif_names[[6]]]]

par(mfrow=c(1,2), pin=c(4, 2))
seqLogo(pwm = pwm1)
seqLogo(pwm = pwm2)

## plot DiffLogo
#diffLogoFromPwm(pwm2 = motifs[["eucarya"]], pwm1 = motifs[["archaea"]], alphabet = ASN)
diffLogoFromPwm(pwm1 = pwm1, pwm2 = pwm2)

## diffLogoFromPwm is a convenience function for
diffLogoObj = createDiffLogoObject(pwm1 = pwm1, pwm2 = pwm2);
diffLogo(diffLogoObj)

## plot DiffLogo table
diffLogoTable(motifs);

###########################################################################################################
## tool comparison
library(DiffLogo)

## import pwms
motif_folder = "pwm";
motif_names = c("GM12878", "H1-hESC", "HeLa-S3", "HepG2", "HUVEC", "K562", "MCF7", "NHEK","ProgFib");
motifs = list();
for (name in motif_names) {
  print(name)
  file = system.file(paste(motif_folder,"/",name,".txt",sep=""), package = "DiffLogo")
  print(file)
  motifs[[name]] = as.matrix(read.delim(file,header=F));
}

dim = length(motifs);
similarities = matrix(0,dim,dim);
for ( i in 1:dim) {
  for ( k in 1:dim) {
    similarities[i,k] = 0
    if( i != k ) {
      diffLogoObj = createDiffLogoObject(motifs[[ i ]],motifs[[ k ]],stackHeight=shannonDivergence, baseDistribution=normalizedDifferenceOfProbabilities, alphabet=DNA);
      similarities[i,k] = diffLogoObj$distance;
      
      print(paste(i, k, motif_names[[i]], motif_names[[k]], diffLogoObj$distance))
    }
  }
}


## which(similarities == min(similarities), arr.ind=TRUE)
## min 1 6 GM12878 K562 0.00510014885465492
## 
## which(similarities == max(similarities), arr.ind=TRUE)
## max 2 5 H1-hESC HUVEC 0.360831234802894"
## 
## + seqLogo
## - iceLogo --> down
## - MotifStack --> not installable
## + STAMP --> no tree because all distances = zero
## + Two Sample Logo --> default but DNA / RNA alphabet; 300 dpi resolution
## 
pwmSimilar1 = motifs[[motif_names[[1]]]]
pwmSimilar2 = motifs[[motif_names[[6]]]]
pwmDissimilar1 = motifs[[motif_names[[2]]]]
pwmDissimilar2 = motifs[[motif_names[[5]]]]
pwmSimilar1 = getPwmFromAlignment(readLines("/home/htreutle/Downloads/DiffLogo/Similar_250_GM12878.txt"), DNA, 1)
pwmSimilar2 = getPwmFromAlignment(readLines("/home/htreutle/Downloads/DiffLogo/Similar_250_K562.txt"), DNA, 1)
pwmDissimilar1 = getPwmFromAlignment(readLines("/home/htreutle/Downloads/DiffLogo/Dissimilar_250_H1-hESC.txt"), DNA, 1)
pwmDissimilar2 = getPwmFromAlignment(readLines("/home/htreutle/Downloads/DiffLogo/Dissimilar_250_HUVEC.txt"), DNA, 1)

## parameters
widthToHeightRatio = 3/1;
size = 2
resolution <- 300
width <- size * widthToHeightRatio
height <- size

## seqLogo
fileName <- "/home/htreutle/Downloads/DiffLogo/Similar_250_GM12878_seqLogo.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
par(mar=c(0.5,1.5,0.5,0.5))
seqLogo(pwm = pwmSimilar1, sparse = TRUE, alphabet = DNA)
dev.off()
fileName <- "/home/htreutle/Downloads/DiffLogo/Similar_250_K562_seqLogo.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
par(mar=c(0.5,1.5,0.5,0.5))
seqLogo(pwm = pwmSimilar2, sparse = TRUE, alphabet = DNA)
dev.off()
fileName <- "/home/htreutle/Downloads/DiffLogo/Similar_250_H1-hESC_seqLogo.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
par(mar=c(0.5,1.5,0.5,0.5))
seqLogo(pwm = pwmDissimilar1, sparse = TRUE, alphabet = DNA)
dev.off()
fileName <- "/home/htreutle/Downloads/DiffLogo/Similar_250_HUVEC_seqLogo.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
par(mar=c(0.5,1.5,0.5,0.5))
seqLogo(pwm = pwmDissimilar2, sparse = TRUE, alphabet = DNA)
dev.off()

## generate sequences
motifNamesHere <- c(motif_names[[1]], motif_names[[6]], motif_names[[2]], motif_names[[5]])
motifsHere <- list(
  "GM12878" = pwmSimilar1,
  "K562" = pwmSimilar2,
  "H1-hESC" = pwmDissimilar1,
  "HUVEC" = pwmDissimilar2
)

numberOfSequences <- 250
for(idx in 1:length(motifsHere)){
  motifNameHere <- motifNamesHere[[idx]]
  motifHere <- motifsHere[[idx]]
  
  seqs <- list()
  for(idx2 in 1:numberOfSequences){
    seq <- list()
    for(idx3 in 1:ncol(motifHere))
      seq[[length(seq) + 1]] <- sample(x = c("A", "C", "G", "T"), prob = motifHere[, idx3], replace = TRUE, size = 1)
    seq <- paste(unlist(seq), sep = "", collapse = "")
    seqs[[length(seqs) + 1]] <- seq
  }
  file <- paste("/home/htreutle/Downloads/", motifNameHere, ".txt", sep = "")
  write.table(x = unlist(seqs), file = file, append = FALSE, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

motif1 = getPwmFromAlignment(readLines("/home/htreutle/Downloads/DiffLogo/Similar_250_GM12878.txt"), DNA, 1)
motif2 = getPwmFromAlignment(readLines("/home/htreutle/Downloads/DiffLogo/Similar_1000_GM12878.txt"), DNA, 1)
motif1 = getPwmFromAlignment(readLines("/home/htreutle/Downloads/DiffLogo/Similar_250_K562.txt"), DNA, 1)
motif2 = getPwmFromAlignment(readLines("/home/htreutle/Downloads/DiffLogo/Similar_1000_K562.txt"), DNA, 1)

seqLogo(pwm = motif1)
seqLogo(pwm = motif2)

## motifStack
#source("http://bioconductor.org/biocLite.R")
#biocLite("motifStack")

library(motifStack)
#####Input#####
pcms<-readPCM(file.path(find.package("motifStack"), "extdata"),"pcm$")
motifs<-lapply(pcms,pcm2pfm)
## plot stacks with hierarchical tree
motifStack(motifs, layout="tree")

## DiffLogo table
widthToHeightRatio = 16/10;
size = 10
resolution <- 300
width <- size * widthToHeightRatio
height <- size

fileName <- "/home/htreutle/Downloads/DiffLogo/DiffLogoTable_250_GM12878_K562_H1-hESC_HUVEC.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
diffLogoTable(PWMs = motifsHere, alphabet = DNA)
dev.off()

## Difference logos
widthToHeightRatio = 16/10;
size = 4
resolution <- 300
width <- size * widthToHeightRatio
height <- size

fileName <- "/home/htreutle/Downloads/DiffLogo/DiffLogo_Similar_250_GM12878_vs_K562.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
diffLogoFromPwm(pwm1 = pwmSimilar1, pwm2 = pwmSimilar2)
dev.off()
fileName <- "/home/htreutle/Downloads/DiffLogo/DiffLogo_Similar_250_K562_vs_GM12878.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
diffLogoFromPwm(pwm1 = pwmSimilar2, pwm2 = pwmSimilar1)
dev.off()
fileName <- "/home/htreutle/Downloads/DiffLogo/DiffLogo_Dissimilar_250_H1-hESC_vs_HUVEC.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
diffLogoFromPwm(pwm1 = pwmDissimilar1, pwm2 = pwmDissimilar2)
dev.off()
fileName <- "/home/htreutle/Downloads/DiffLogo/DiffLogo_Dissimilar_250_HUVEC_vs_H1-hESC.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
diffLogoFromPwm(pwm1 = pwmDissimilar2, pwm2 = pwmDissimilar1)
dev.off()


###########################################################################################################
## ASN - Inteins

table <- read.csv(file = "/home/htreutle/Data/Motifs/InBase_BlockA.csv", header = TRUE, sep = "\t")

eucarya <- as.character(table[, 1][table[, 1] != ""])
eubacteria <- as.character(table[, 2][table[, 2] != ""])
archaea <- as.character(table[, 3][table[, 3] != ""])

pseudoCount <- 0
pwmEucarya <- getPwmFromAlignment(eucarya, ASN, pseudoCount)
pwmEubacteria <- getPwmFromAlignment(eubacteria, ASN, pseudoCount)
pwmArchaea <- getPwmFromAlignment(archaea, ASN, pseudoCount)

motifs <- list();
motifs[["eucarya"]] <- pwmEucarya
motifs[["eubacteria"]] <- pwmEubacteria
motifs[["archaea"]] <- pwmArchaea
diffLogoTable(PWMs = motifs, alphabet = ASN)

## export DiffLogo table as png image
fileName <- "/home/htreutle/Downloads/Intein_BlockA_comparison.pdf"
pdf(file = fileName, width = width, height = height)
diffLogoTable(PWMs = motifs, alphabet = ASN)
dev.off()

fileName <- "Intein_BlockA_comparison.png"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
diffLogoTable(PWMs = motifs, alphabet = ASN)
dev.off()

###########################################################################################################
## cluster tree error
aaa <- as.character(c("A", c("A", "A")))
aac <- as.character(c("A", c("A", "C")))
agc <- as.character(c("A", c("G", "C")))

pwmAaa <- getPwmFromAlignment(aaa, ASN, 1)
pwmAac <- getPwmFromAlignment(aac, ASN, 1)
pwmAgc <- getPwmFromAlignment(agc, ASN, 1)

motifs <- list();
motifs[["aaa"]] <- pwmAaa
motifs[["aac"]] <- pwmAac
motifs[["agc"]] <- pwmAgc
diffLogoTable(PWMs = motifs, alphabet = ASN)

###########################################################################################################
## results
motifs[["eucarya"]]
motifs[["eubacteria"]]
motifs[["archaea"]]

pwm <- motifs[["archaea"]]
pwm <- motifs[["eubacteria"]]
pwm <- motifs[["eucarya"]]
for (j in 1:dim(pwm)[2]) {
  column = pwm[, j]
  sh = informationContent(column);
  print(paste(j, sh[1]))
}


###########################################################################################################
## export

## parameters
widthToHeightRatio = 16/10;
size = length(motifs) * 2
resolution <- 300
width <- size * widthToHeightRatio
height <- size

## export single DiffLogo as pdf document
fileName <- "Comparison_of_two_motifs.pdf"
pdf(file = fileName, width = width, height = height)
diffLogoFromPwm(pwm1 = pwm1, pwm2 = pwm2)
dev.off()

## export DiffLogo table as png image
fileName <- "Archaea"
png(filename = fileName, width = width * resolution, height = height * resolution, res = resolution)
diffLogoTable(PWMs = motifs)
dev.off()
