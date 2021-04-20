library(brranching)
library(phytools)
library(phylocomr)
library(Hmisc)
library(ape)
library(scales)
library(corrplot)

setwd("your file path PCE.Figures")
alldat=read.csv("PCE.data.csv")
unique(alldat$species)
spnames=alldat$fam.gen.spe
mytree=brranching::phylomatic(spnames, outformat="newick", storedtree="R20120829", get="POST") 

#The ages of the nodes (used to calculate branch lengths) according to: http://www.scielo.br/scielo.php?script=sci_arttext&pid=S1519-69842016000300619
#Note: data from the reference above needs to be obtained before analysis can proced
ages_df=read.csv("ages.from.gastatuer.et.al.csv", header=T) 
mynw=write.tree(mytree)# Turn phylogeny into character format

res <- phylocomr::ph_bladj(ages_df, mynw)#this designates branch lengths
my.new.tree <- read.tree(text = res)

#Remove unerscores from tip labels can capitalizse first letter in string using Hmisc
my.new.tree$tip.label=Hmisc::capitalize(gsub("_", " ", my.new.tree$tip.label))
#insect
class(my.new.tree);str(my.new.tree) 

#Sort the heat toerances so they match the order of tree's tip.labels
alldat2=alldat[match(my.new.tree$tip.label, capitalize(as.character(alldat$species))),]

#Color Palette
ggPal=colorRampPalette(c("dark green","green"))
as.matrix(alldat2$mean.Popt)
traits=as.matrix(alldat2$mean.Popt)
row.names(traits)=alldat2$species
str(traits)

kolors=ggPal(21)[as.numeric(cut((traits[,1]),breaks = 21))]
tree=untangle(my.new.tree, "read.tree")

dat=alldat2$mean.Popt
names(dat)<-alldat2$species
dat<-(dat[my.new.tree$tip.label])
ii<-c(1:floor(0.25*Ntip(my.new.tree)),
      ceiling(0.75*Ntip(my.new.tree)):Ntip(my.new.tree))
jj<-setdiff(1:Ntip(my.new.tree),ii)
hack<-my.new.tree
hack$tip.label[ii]<-paste("  ",my.new.tree$tip.label[ii],sep="")
hack$tip.label[jj]<-paste(my.new.tree$tip.label[jj],"  ",sep="")

#plot the optimum photosynthetic rate for each species. 
par(mfcol=c(1,1))
plotTree.wBars(my.new.tree,(traits[,1]),fsize=.8, fsize=.8, tip.labels=T,  type="fan",  col=kolors, method="plotTree", lwd=2 )


#-----------------------------------------------------------------------------------------------------------------------
# Check for phlyogenetic signal among traits
#-----------------------------------------------------------------------------------------------------------------------
#Resolves dichomomies required for pic() in ape
my.new.tree2=multi2di(my.new.tree)
is.binary(my.new.tree2)
alldat2$species==my.new.tree2$tip.label#check that tips and labels match

#Example of how to use function to calculate phylogenetic signal
phylosig(my.new.tree2, alldat2$tcrit,  method = "K", test=T, nsim=10000)
phylosig(my.new.tree2, alldat2$t50,  method = "K", test=T, nsim=10000)
phylosig(my.new.tree2, alldat2$t100,  method = "K", test=T, nsim=10000)
#There are no phylogenetic signal observed in any of the traits

#-----------------------------------------------------------------------------------------------------------------------
#Correlation bt photosynthetic traits while controlling for phylogeny
#-----------------------------------------------------------------------------------------------------------------------
str(alldat2)
datsp=alldat2[,c(1,3:11)]

datmat=matrix(as.numeric(unlist(datsp[,-1])),nrow=nrow(datsp[,-1]))
colnames(datmat) <- colnames(datsp[,-1])
rownames(datmat) <- datsp$species

datsp$species==my.new.tree2$tip.label#just order one more time

#phylogenetic variance co-variance matrix
obj=phyl.vcv(datmat, vcv(my.new.tree2),1)
rmat=cov2cor(obj$R)
# t-statistic & P-value
tmat<-rmat*sqrt((Ntip(my.new.tree2)-2)/(1-rmat^2))
Pmat<-2*pt(abs(tmat),df=21-2,lower.tail=F)

colnames(rmat) <- c("=T[opt]", "=P[opt]", ":Omega", "=T[opt2]", "=P[opt2]",  "=T[max]", "=T[crit]", "=T[50]", "=T[95]")#
rownames(rmat) <-  c("=T[opt]", "=P[opt]", ":Omega", "=T[opt2]", "=P[opt2]",  "=T[max]", "=T[crit]", "=T[50]", "=T[95]")#
rownames(Pmat) <-  c("=T[opt]", "=P[opt]", ":Omega", "=T[opt2]", "=P[opt2]",  "=T[max]", "=T[crit]", "=T[50]", "=T[95]")#
rownames(Pmat) <-  c("=T[opt]", "=P[opt]", ":Omega", "=T[opt2]", "=P[opt2]",  "=T[max]", "=T[crit]", "=T[50]", "=T[95]")#

#Make a correlation plot all all traits
corrplot(rmat, type = "upper", p.mat = Pmat, sig.level = 0.05, tl.col="black")

#-----------------------------------------------------------------------------------------------------------------------
# Plotting thermal photosynthetic traits
#-----------------------------------------------------------------------------------------------------------------------

fbp=data.frame(temps=c(datsp$mean.Topt, datsp$mean.Tmax, datsp$tcrit, datsp$t50, datsp$t100), 
               type=c(rep("topt", 21), rep("tmax", 21), rep("tcrit", 21), rep("t50", 21), rep("t95", 21)),
               colrz=c(rep("dark green", 21), rep("dark gray", 21), rep("dark blue", 21), rep("dark orange", 21), rep("dark red", 21)))
fbp$type<- factor(fbp$type , levels=c("topt", "tcrit", "tmax", "t50", "t95"))
meds=sapply(datsp[,c(2,4,6,7:10)], mean, na.rm=TRUE)

#Set up plotting areas
zones=matrix( c(c(rep(1,5)), c(rep(2,10)), c(rep(3, 30))), nrow=9, byrow=TRUE) #Zone theory: Take a break and google Tim & ERic Zone Theory
layout(zones, widths = rep.int(0.5, ncol(zones)),
       heights = rep.int(1, nrow(zones)), respect = FALSE)

#Adjust plotting parameters
par( mgp=c(10, 1, 0), oma=c(7, 5, 1, 1), mar=c(0,12,0,0))
plot(NULL, NULL, xlim=c(28,65) ,ylim=c(0,3), yaxt='n', xaxt='n',  ylab="" , xlab="", bty='n')

#Create symbol legend
text(46.5, 2.75, "Dataset means and symbol legend:", cex=1.5)
lines( c(meds[1], meds[4]), c(1,1), lty=2)
points(meds[1], 1, bg=alpha("black", 0.75), pch=21, col="black", cex=2)
text(meds[1], 1.75, expression("T"[opt]), cex=1.5)
library(plotrix)
gradient.rect(meds[5], 1-0.25, meds[6], 1+0.25, gradient="x", 
              col=alpha( smoothColors("dark blue",21,"dark orange"), 0.75) ,border=T)
text(meds[5], 1.75, expression("T"[crit]), cex=1.5)
gradient.rect(meds[6], 1-0.25, meds[7], 1+0.25, gradient="x", 
              col=alpha( smoothColors("dark orange",24,"dark red"), 0.75) ,border=T)
text(meds[6], 1.75, expression("T"[50]), cex=1.5)
text(meds[7], 1.75, expression("T"[95]), cex=1.5)
points(meds[4], 1-0.09, bg=alpha("gray76", 1), pch=24, col="black" ,cex=2)
text(meds[4], .2, expression("T"[max]), cex=1.5, pos=4)

#Means and quantile boxplots of traits thermal traits
boxplot(fbp$temps~fbp$type, horizontal=T, ylim=c(28,65), yaxt='n', ylab='', xaxt='n', xlab='' ,
        col=alpha(c("dark green", "dark blue", "gray76", "dark orange", "dark red"),  0.75))
axis(2, 1:5, labels = c( expression("T"[opt]),  expression("T"[crit]), 
                         expression("T"[max]), expression("T"[50]), expression("T"[95])), 
     font=3, las=2, cex.axis=1.5)
stripchart(temps ~ type,  data = fbp, 
           method = "jitter", add = TRUE, pch=21, 
           col=alpha(c("dark green", "dark blue", "gray76", "dark orange", "dark red"), 0.75))

#Thermal traits for each species
plot(NULL, NULL, xlim=c(28,65) ,ylim=c(1,21), yaxt='n', xaxt='n',  ylab="" , 
     xlab=expression(paste("Temperature ("~degree~"C)"))  ) 
datsp=datsp[order(datsp$mean.Tmax),]
for(i in 1:length(datsp[,1])){
  df=datsp[i,]
  #Tcrit-T50 box
  gradient.rect(df$tcrit, i-0.25, df$t50, i+0.25, gradient="x", 
                col=alpha( smoothColors("dark blue",21,"dark orange"), 0.75) ,border=T)
  #T50-T100 box
  gradient.rect(df$t50, i-0.25, df$t100, i+0.25, gradient="x", 
                col=alpha( smoothColors("dark orange",24,"dark red"), 0.75) ,border=T)
  #Topt
  lines( c(df$mean.Topt, df$mean.Tmax), c(i,i), lty=2)
  points(df$mean.Topt, i-0.01, bg=alpha("dark green", 0.75), pch=21, col="black", cex=2)
  #Tmax
  points(df$mean.Tmax, i-0.09, bg=alpha("gray76", 1), pch=24, col="black" ,cex=2)
  axis(2, at=i, labels = df$species, font=3, las=2, cex.axis=1.5)
}
axis(3, labels = FALSE)
axis(1, cex.axis=1.5)
mtext(expression(paste("Temperature ("~degree~"C)")), 1, 4 ,outer=TRUE, cex=1.15)

#-----------------------------------------------------------------------------------------------------------------------
# Perform Phylogenetic independent contrasts & Create new data frame
#-----------------------------------------------------------------------------------------------------------------------

datsp=datsp[match(my.new.tree2$tip.label, datsp$species),]
datsp$species==my.new.tree2$tip.label
#Perform phylogenetic independent contrasts for the traits of interest
dfpic=data.frame(mean.Topt=unname(pic(datsp$mean.Topt, my.new.tree2)),
                 mean.Popt=unname(pic(datsp$mean.Popt, my.new.tree2)),
                 mean.Topt2=unname(pic(datsp$mean.Topt2, my.new.tree2)),
                 mean.Popt2=unname(pic(datsp$mean.Popt2, my.new.tree2)),
                 mean.Tomeg=unname(pic(datsp$mean.Tomeg, my.new.tree2)),
                 mean.Tmax=unname(pic(datsp$mean.Tmax, my.new.tree2)),
                 Tcrit=unname(pic(datsp$tcrit, my.new.tree2)),
                 T50=unname(pic(datsp$t50, my.new.tree2)),
                 T95=unname(pic(datsp$t100, my.new.tree2)))
#-----------------------------------------------------------------------------------------------------------------------
#Plot results from phylogentic independent contrast for Photosynthetic & vs. PSii heat tolerance traits
#-----------------------------------------------------------------------------------------------------------------------
#All correlations and p-values determined from phylogenetic correlation plots above.

par(mfrow=c(4,3), mar=c(2.4,3.1,0,0), mgp=c(1.5,.25, 0), oma=c(1, 1, .5, .5), cex.lab=1.25, cex.axis=0.75)
min(dfpic$mean.Tomeg);max(dfpic$mean.Tomeg)

#Tmax
plot(dfpic$Tcrit, dfpic$mean.Tmax, ylim=c(-1, 1), xlim=c(-.5, .5), pch=19, xaxt='n', bty="l", col=alpha("black", 0.5),
     xlab='',
     ylab=expression(paste("T"[max]~degree~'C')))
axis(1, outer = F, labels = F)
text(-.5, 1, "A", font=2)
text(-0.5,-0.8,  expression("r = -0.334"), pos=4)
text(-0.5,-1,  expression("p = 0.138"), pos=4)

plot(dfpic$T50, dfpic$mean.Tmax, ylim=c(-1, 1), xlim=c(-.5, .5), pch=19, xaxt='n',  yaxt='n', bty="l",col=alpha("black", 0.5),
     xlab='',
     ylab='')
axis(1, outer = F, labels = F)
axis(2, outer = F, labels = F)
text( -0.5, 1, "B", font=2)
text(-0.5,-0.8,  expression("r = 0.270"), pos=4)
text(-0.5,-1,  expression("p = 0.237"), pos=4)

plot(dfpic$T95, dfpic$mean.Tmax, ylim=c(-1, 1), xlim=c(-1, 1), pch=19, xaxt='n', yaxt='n', bty="l",col=alpha("black", 0.5),
     xlab='',
     ylab='')
axis(1, outer = F, labels = F)
axis(2, outer = F, labels = F)
text(-1, 1, "C", font=2)
text(-1,-0.8,  expression("r = 0.372"), pos=4)
text(-1,-1,  expression("p = 0.256"), pos=4)

#Omega
plot(dfpic$Tcrit, dfpic$mean.Tomeg, ylim=c(-2, 1), xlim=c(-.5, .5), pch=19, xaxt='n', bty="l",col=alpha("black", 0.5),
     xlab='',
     ylab=expression(paste(Omega)~degree~'C') )
axis(1, outer = F, labels = F)
text(-.5, 1, "D", font=2)
text(-0.5,-1.6,  expression("r = -0.190"), pos=4)
text(-0.5,-2,  expression("p = 0.409"), pos=4)


plot(dfpic$T50, dfpic$mean.Tomeg, ylim=c(-2, 1), xlim=c(-.5, .5), pch=19,  xaxt='n',  yaxt='n', bty="l",col=alpha("black", 0.5),
     xlab='',
     ylab='')
axis(1, outer = F, labels = F)
axis(2, outer = F, labels = F)
text(-.5, 1, "E", font=2)
summary(lm( dfpic$mean.Tomeg~dfpic$T50))
abline(lm( dfpic$mean.Tomeg~dfpic$T50))
text(-0.5,-1.6,  expression("r = 0.581"), pos=4)
text(-0.5,-2,  expression("p = 0.006"), pos=4)


plot(dfpic$T95, dfpic$mean.Tomeg, ylim=c(-2, 1), xlim=c(-1, 1), pch=19,  yaxt='n', xaxt='n', bty="l",col=alpha("black", 0.5),
     xlab='',
     ylab='')
axis(1, outer = F, labels = F)
axis(2, outer = F, labels = F)
text(-1, 1, "F", font=2)
summary(lm( dfpic$mean.Tomeg~0+dfpic$T95))
abline(lm(dfpic$mean.Tomeg~dfpic$T95))
text(-1,-1.6,  expression("r = 0.590"), pos=4)
text(-1,-2,  expression("p = 0.005"), pos=4)


#Popt
plot(dfpic$Tcrit, dfpic$mean.Popt, ylim=c(-1, 1.5),xlim=c(-.5, .5), xaxt='n',pch=19, bty="l",col=alpha("black", 0.5),
     xlab='',
     ylab=expression("P"[opt]~"("~mu~"mol m"^-2~"s"^-1~")") )
axis(1, outer = F, labels = F)
text(-.5, 1.5, "G", font=2)
text(-0.5,-0.7,  expression("r = 0.211"), pos=4)
text(-0.5,-1,  expression("p = 0.359"), pos=4)


plot(dfpic$T50, dfpic$mean.Popt, ylim=c(-1, 1.5), xlim=c(-.5, .5), pch=19, xaxt='n',yaxt='n',bty="l",col=alpha("black", 0.5),
     xlab='',
     ylab='' )
axis(1, outer = F, labels = F)
axis(2, outer = F, labels = F)
text(-.5, 1.5, "H", font=2)
summary(lm( dfpic$mean.Popt~dfpic$T50))
abline(lm(dfpic$mean.Popt~dfpic$T50), lty=1)
text(-0.5,-0.7,  expression("r = -0.495"), pos=4)
text(-0.5,-1,  expression("p = 0.022"), pos=4)

plot(dfpic$T95, dfpic$mean.Popt, ylim=c(-1, 1.5), xlim=c(-1, 1), pch=19, xaxt='n',yaxt='n',bty="l",col=alpha("black", 0.5),
     xlab='',
     ylab='' )
axis(1, outer = F, labels = F)
axis(2, outer = F, labels = F)
text(-1, 1.5, "I", font=2, cex=1)
summary(lm( dfpic$mean.Popt~dfpic$T95))
abline(lm(dfpic$mean.Popt~dfpic$T95))
text(-1,-0.7,  expression("r = -0.521"), pos=4)
text(-1,-1,  expression("p = 0.015"), pos=4)


#Topt
plot(dfpic$Tcrit, dfpic$mean.Topt, ylim=c(-.5, .5), xlim=c(-.5, .5), pch=19, bty="l", col=alpha("black", 0.5),
     xlab=expression(paste("T"[crit]~degree~'C')),
     ylab=expression(paste("T"[opt]~degree~'C')) )
axis(1, outer = F, labels = F)
axis(2, outer = F, labels = F)
text(-.5, .5, "J", font=2)
text(-0.5,-0.4,  expression("r = 0.193"), pos=4)
text(-0.5,-0.5,  expression("p = 0.401"), pos=4)

plot(dfpic$T50, dfpic$mean.Topt, ylim=c(-.5, .5), xlim=c(-.5, .5), pch=19, yaxt='n',  yaxt='n', bty="l",col=alpha("black", 0.5),
     xlab=expression(paste("T"[50])~degree~'C'),
     ylab='' )
axis(1, outer = F, labels = F)
axis(2, outer = F, labels = F)
text(-.5, .5, "K", font=2)
summary(lm( dfpic$mean.Topt~dfpic$T50))
abline(lm(dfpic$mean.Topt~dfpic$T95), lty=2)
text(-0.5,-0.4,  expression("r = -0.432"), pos=4)
text(-0.5,-0.5,  expression("p = 0.051"), pos=4)

plot(dfpic$T95, dfpic$mean.Topt, ylim=c(-.5, .5), xlim=c(-1, 1), pch=19, yaxt='n',  yaxt='n', bty="l",col=alpha("black", 0.5),
     xlab=expression(paste("T"[95])~degree~'C'),
     ylab='')
axis(1, outer = F, labels = F)
axis(2, outer = F, labels = F)
text(-1, .5, "L", font=2)
summary(lm( dfpic$mean.Topt~dfpic$T95))
abline(lm(dfpic$mean.Topt~dfpic$T95))
text(-1,-0.4,  expression("r = -0.452"), pos=4)
text(-1,-0.5,  expression("p = 0.039"), pos=4)



#-----------------------------------------------------------------------------------------------------------------------
# Phylogentic independent contrast plots for omega & vs. other photosynthetic traits
#-----------------------------------------------------------------------------------------------------------------------
par(mfrow=c(3,1), mar=c(2.4,3.1,0,0), mgp=c(1.5,.25, 0), oma=c(1, 1, .5, .5), cex.lab=1.25, cex.axis=0.75)
min(dfpic$mean.Tomeg);max(dfpic$mean.Tomeg)
min(dfpic$mean.Topt);max(dfpic$mean.Topt)
min(dfpic$mean.Popt);max(dfpic$mean.Popt)
min(dfpic$mean.Tmax);max(dfpic$mean.Tmax)
#Omeg

plot(dfpic$mean.Tomeg, dfpic$mean.Tmax, ylim=c(-1, 1), xlim=c(-2, 1), pch=19,  xaxt='n', bty="l",col=alpha("black", 0.5),
     ylab=expression("T"[max]~degree~'C'),
     xlab="")
axis(1, outer = F, labels = F)
summary(lm( dfpic$mean.Tmax~0+dfpic$mean.Tomeg))
abline(lm(dfpic$mean.Tmax~0+dfpic$mean.Tomeg))
text(-2, 1, "A", font=2)
text(-2,-0.8,  expression("r = 0.567"), pos=4)
text(-2,-1,  expression("p = 0.007"), pos=4)

plot(dfpic$mean.Tomeg, dfpic$mean.Topt, ylim=c(-0.5, 0.5), xlim=c(-2, 1), pch=19,  xaxt="n", bty="l",col=alpha("black", 0.5),
     ylab=expression("T"[opt]~degree~'C'),
     xlab='')
axis(1, outer = F, labels = F)
axis(2, outer = F, labels = F)
summary(lm( dfpic$mean.Topt~dfpic$mean.Tomeg))
abline(lm(dfpic$mean.Topt~dfpic$mean.Tomeg))
text(-2, 0.5, "B", font=2)
text(-2,-0.4,  expression("r = -0.489"), pos=4)
text(-2,-0.5,  expression("p = 0.024"), pos=4)

plot(dfpic$mean.Tomeg, dfpic$mean.Popt, ylim=c(-1, 1.5), xlim=c(-2, 1), pch=19,  bty="l",col=alpha("black", 0.5),
     ylab=expression("P"[opt]~"("~mu~"mol m"^-2~"s"^-1~")"), 
     xlab=expression(paste(Omega)~degree~'C') )
axis(1, outer = F, labels = F)
axis(2, outer = F, labels = F)
text(-2, 1.5, "C", font=2)
text(-2,-0.7,  expression("r = -0.094"), pos=4)
text(-2,-1.0,  expression("p = 0.584"), pos=4)
