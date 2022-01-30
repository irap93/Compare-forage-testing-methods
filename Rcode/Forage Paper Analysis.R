rm(list = ls())
library(tidyverse)
library(lubridate)
library(formattable)
library(car)
library(betareg)
forage <- read.csv("MSU_Grazingecology/Rawdata/Forage/Forage.csv", 
                   header = TRUE)
str(forage)
forage$NIR.DateCollected <- date(forage$NIR.DateCollected)
length(unique(forage$NIR.DateCollected))

wet.dat <- filter(forage, WetID != "NA")
length(colnames(wet.dat)) # How many columns?
length(unique(wet.dat$SampleID)) # see how many sampling stations are in the dataset

wet.dat <- dplyr::select(wet.dat,-c("NIR.Ca.PDB","NIR.P.PDB","NIR.DM.PAR","NIR.K.PDB","NIR.Mg.PDB",
                                    "NIR.IVTDMD","NIR.Fat.PDB", "NIR.Ash.PDB","NIR.dNDF48.PDB",
                                    "NIR.Starch.PDB","NIR.Fructan.PDB","NIR.WSC.PDB","NIR.CP.PAR"))
colnames(wet.dat)

wet.dat$NIR.DM.prc <- 100-((wet.dat$NIR.Wetwt.g - wet.dat$NIR.Drywt.g)/wet.dat$NIR.Wetwt.g)*100
wet.dat$Wet.DM.prc <- 100 -((wet.dat$Wet.Wetwt.g - wet.dat$Wet.Drywt.g)/wet.dat$Wet.Wetwt.g)*100

wet.dat$complete <- complete.cases(wet.dat) # Check for complete observations
wet.dat <- filter(wet.dat, complete == TRUE)
wet.dat <- wet.dat[,1:24] # Remove complete Cases Column

### Look at the data
par(mfrow=c(4,5))
for (i in 6:24) {
  hist(wet.dat[,i], main = colnames(wet.dat)[i])
}
sm <- summary(wet.dat)
formattable(sm)
### Correlations
cor.dat <- data.frame(cor(wet.dat[,6:22], use = 'complete.obs',method = 'pearson'))
formattable(print(cor.dat[9:17,1:7]))

par.cor <- data.frame(cor(wet.dat[,6:22],use = "pairwise.complete.obs",
                          method = "pearson"))
formattable(par.cor[c(3:7,17),c(12,14:17)])

### Look at differences using pairwise T-test
t.test(wet.dat$NIR.Wetwt.g,wet.dat$Wet.Wetwt.g, paired = T)
t.test(wet.dat$NIR.Drywt.g,wet.dat$Wet.Drywt.g, paired = T)
t.test(wet.dat$NIR.DM.prc,wet.dat$Wet.DM.prc, paired = T)
t.test(wet.dat$NIR.CP.PDB,wet.dat$Wet.CP.PDP, paired = T)
t.test(wet.dat$NIR.ADF.PDB,wet.dat$Wet.ADF.PDB, paired = T)
t.test(wet.dat$NIR.NDF.PDB, wet.dat$Wet.NDF.PDB, paired = T)
t.test(wet.dat$NIR.RFV,wet.dat$Wet.RFV, paired = T)
t.test(wet.dat$NIR.TDN.PDB, wet.dat$Wet.TDN.PDB, paired = T)

#### Graph the differences with a boxplot
par(mfrow=c(2,4))
boxplot(wet.dat$NIR.Wetwt.g,wet.dat$Wet.Wetwt.g, names = c("NIR","Wet"), main = "Wet Weight")
boxplot(wet.dat$NIR.Drywt.g, wet.dat$Wet.Drywt.g,names = c("NIR","Wet"), main = "Dry Weight")
boxplot(wet.dat$NIR.DM.prc,wet.dat$Wet.DM.prc, names = c("NIR","Wet"),main = "Percent Dry Matter")
boxplot(wet.dat$NIR.CP.PDB, wet.dat$Wet.CP.PDP, names = c("NIR", "Wet"), main = "Crude Protein")
boxplot(wet.dat$NIR.ADF.PDB,wet.dat$Wet.ADF.PDB, names = c("NIR","Wet"), main = "Acid Detergent Fiber")
boxplot(wet.dat$NIR.NDF.PDB, wet.dat$Wet.NDF.PDB,names = c("NIR","Wet"), main = "Neutral Detergent Fiber")
boxplot(wet.dat$NIR.RFV,wet.dat$Wet.RFV,names = c("NIR","Wet"), main = "Relative Feed Value")
boxplot(wet.dat$NIR.TDN.PDB,wet.dat$Wet.TDN.PDB,names = c("NIR","Wet"), main = "Total Digestible Nutrients")

#### Plot as continuos variable
plot(wet.dat$NIR.Wetwt.g, wet.dat$Wet.Wetwt.g, xlab = "NIR", ylab = "Wet Chemistry",
     main = "Sample Wet Weight")
plot(wet.dat$NIR.Drywt.g,wet.dat$Wet.Drywt.g, xlab = "NIR", ylab = "Wet Chemistry",
     main = "Sample Dry Weight")
plot(wet.dat$NIR.DM.prc,wet.dat$Wet.DM.prc, 
     xlab = "NIR", ylab = "Wet Chemistry", main = "Percent Dry Matter")
plot(wet.dat$NIR.CP.PDB, wet.dat$Wet.CP.PDP, main = "Crude Protein",
     xlab = "NIR", ylab = "Wet Chemistry")
plot(wet.dat$Wet.ADF.PDB,wet.dat$NIR.ADF.PDB,
     xlab = "NIR", ylab = "Wet Chemistry", main = "Acid Detergent Fiber")
plot(wet.dat$NIR.NDF.PDB, wet.dat$Wet.NDF.PDB,
     xlab = "NIR",ylab = "Wet Chemistry", main = "Neurtral Detergent Fiber")
plot(wet.dat$NIR.RFV,wet.dat$Wet.RFV,
     xlab = "NIR",ylab = "Wet", main = "Relative Feed Value")
plot(wet.dat$NIR.TDN.PDB,wet.dat$Wet.TDN.PDB,
     xlab = "NIR", ylab = "Wet", main = "Total Digestible Nutrients")
par(mfrow=c(1,1))

#### Compare variation and distribution #### 
# Comparing the distributions recieved from various sampling methods
# F test in R assuming normal distribution
s1 <- sd(wet.dat$NIR.Wetwt.g)
s2 <- sd(wet.dat$Wet.Wetwt.g)
F <- s1^2 / s2^2
FL <- qf(0.05,323,323)
FU <- qf(1-0.05,323,323)

# Brown Forsyth test (Brown, M.B., Forsythe A.B. 1974a)
library(onewaytests)
dat <- pivot_longer(wet.dat, cols = c(6:24)) %>% 
  separate(col = name, into = c("Test","Cons","Unit"))
dat$Test <- as.factor(dat$Test)

bf.test(value ~ Test, data = dat[dat$Cons == "Wetwt",], alpha = 0.05, na.rm = TRUE)
bf.test(value ~ Test, data = dat[dat$Cons == "Drywt",], alpha = 0.05, na.rm = TRUE)
bf.test(value ~ Test, data = dat[dat$Cons == "CP",], alpha = 0.05, na.rm = TRUE)
bf.test(value ~ Test, data = dat[dat$Cons == "ADF",], alpha = 0.05, na.rm = TRUE)
bf.test(value ~ Test, data = dat[dat$Cons == "NDF",], alpha = 0.05, na.rm = TRUE)
bf.test(value ~ Test, data = dat[dat$Cons == "RFV",], alpha = 0.05, na.rm = TRUE)
bf.test(value ~ Test, data = dat[dat$Cons == "TDN",], alpha = 0.05, na.rm = TRUE)
bf.test(value ~ Test, data = dat[dat$Cons == "Lignin",], alpha = 0.05, na.rm = TRUE)
bf.test(value ~ Test, data = dat[dat$Cons == "DM",], alpha = 0.05, na.rm = TRUE)

##### GLM model ####
d <- as.data.frame(unique(wet.dat$NIR.DateCollected))
colnames(d) <- "Date"
d$Date <- ymd(d$Date)
d$rank <- rank(d$Date)
dt <- left_join(wet.dat,d, by = c("NIR.DateCollected" = "Date"))

m1 <- betareg(I(NIR.CP.PDB/100) ~ Wet.CP.PDP + NIR.DateCollected, data = dt )
summary(m1)
