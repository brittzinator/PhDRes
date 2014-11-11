rm(list=ls(all=TRUE))

#Model 
b<-100 						#Number of seeds per adult
g<-0.25							#Probability of maturing (or surviving) until adulthood

#Test population growth over a range of habitat qualities and availabilities (these units become the major tiles of the multicolor plots)
R1vec<-c(0.01,0.25,0.5,0.75,0.99)			#Proportion of habitat affected by Q1 			
R2vec<-c(0.01,0.25,0.5,0.75,0.99)			#Proportion of habitat affected by Q2 
Q1vec<-c(0.01,0.25,0.5,0.75,0.99)			#Range of values for Q1 (Penalty in Habitat 1) 
Q2vec<-c(0.01,0.25,0.5,0.75,0.99)			#Range of values for Q2 (Penalty in Habitat 2)

#Generates LOW RESOLUTION PLOTS
D1vec<-c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99)		#Range of values for D1	(Proportion dispersing away from Habitat 1)
D2vec<-c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99)  	#Range of values for D2 (Proportion dispersing away from Habitat 2)

#Generates HIGH RESOLUTION PLOTS (and takes much longer; >24 hours)
#D1vec<-c(0.01,seq(.025,.975,.05),.99)	#Range of values for D1	(Proportion dispersing away from Habitat 1)
#D2vec<-c(0.01,seq(.025,.975,.05),.99)	#Range of values for D2 (Proportion dispersing away from Habitat 2)
z<- matrix(,ncol=8)					#Matrix for results

#Begin six "for loops"
for (R2 in R2vec) {
for (R1 in R1vec) {
for (Q2 in Q2vec) {
for (Q1 in Q1vec) {
for (D2 in D2vec) {
for (D1 in D1vec) {

#Transition matrix describing the metapopulation
M<-matrix(c(
0,	b*Q1*(1-D1)+b*Q1*D1*R1, 	0, 	b*Q2*D2*R1,
g,	0,				0,	0,
0,	b*Q1*D1*R2,			0,	b*Q2*(1-D2)+b*Q2*D2*R2,
0,	0,				g,	0),
nrow=4,byrow=TRUE)

EIGM<-max(Re(eigen(M)$value))			#Dominant eigenvalue for transition matrix
SUMR<-R1+R2					#Sum of R's (not more than 100%)	
z<-rbind(z,c(R1,R2,D1,D2,Q1,Q2,SUMR,EIGM))		#Store results in matrix z
#print(R2)
}}}}}}

z<-na.omit(z)						#Omit NAs
Z<- as.data.frame(z)					#Turn matrix into a dataframe
names(Z)<-c("R1","R2","D1","D2","Q1","Q2","SUMR","EIGM")	#Name the columns

##Only use the code below if you are running a HIGH RES simulation
#setwd("#Name your working drive here")
#write.table(Z,file="largedataset.txt")

##If you want to find the data without running the simulation again:
##Load data from seperate file

#setwd(#Name your working drive here")			#Set active directory

#data<-read.table("largedataset.txt",header=TRUE)	#Read in dataset
#head(data)
#Z<-data

#Figure 2 (Black and White)

par(mfrow=c(1,3))
subs<-subset(Z, Z$R1==0.5 & Z$R2==0.5 & Q2==0.01 & D2==0.5)
subs2<-subset(subs, subs$D1==0.01 | subs$D1==0.5 | subs$D1==0.99)
plot(format(EIGM, digits=5)~Q1, data=subs2[subs2$D1==.01,], col="grey20", pch=3, type="b", ylab="Metapop. growth rate",xlab="",ylim=c(0,1.6),cex.lab=1.5, bty="n", cex=2, las=1)
points(format(EIGM, digits=5)~Q1, data=subs2[subs2$D1==.5,], col="grey50", pch=2, type="b", cex=2)
points(format(EIGM, digits=5)~Q1, data=subs2[subs2$D1==.99,], col="grey80", pch=1, type="b", cex=2)
abline(h=1)
text(.1,1.5, "a", cex=2)

subs<-subset(Z, Z$R1==0.5 & Z$R2==0.5 & Q2==0.5 & D2==0.5)
subs2<-subset(subs, subs$D1==0.01 | subs$D1==0.5 | subs$D1==0.99)

plot(format(EIGM, digits=5)~Q1, data=subs2[subs2$D1==.01,], col="grey20", pch=3, type="b", ylab=" ", ylim=c(0,1.6),xlab="Local habitat quality (Q)",cex.lab=1.5, bty="n", cex=2, las=1)
points(format(EIGM, digits=5)~Q1, data=subs2[subs2$D1==.5,], col="grey50", pch=2, type="b", cex=2)
points(format(EIGM, digits=5)~Q1, data=subs2[subs2$D1==.99,], col="grey80", pch=1, type="b", cex=2)
abline(h=1)
text(.1,1.5, "b", cex=2)

subs<-subset(Z, Z$R1==0.5 & Z$R2==0.5 & Q2==0.99 & D2==0.5)
subs2<-subset(subs, subs$D1==0.01 | subs$D1==0.5 | subs$D1==0.99)

plot(format(EIGM, digits=5)~Q1, data=subs2[subs2$D1==.01,], col="grey20", pch=3, type="b", ylab=" ", xlab="",ylim=c(0,1.6),cex.lab=1.5, bty="n", cex=2, las=1)
points(format(EIGM, digits=5)~Q1, data=subs2[subs2$D1==.5,], col="grey50", pch=2, type="b", cex=2)
points(format(EIGM, digits=5)~Q1, data=subs2[subs2$D1==.99,], col="grey80", pch=1, type="b", cex=2)
abline(h=1)
text(.1,1.5, "c", cex=2)

legend(x=0.55, y=0.9, col=c("grey80", "grey50", "grey20"), pch=c(1,2,3), c("D=0.99","D=0.5","D=0.01"), bty="n", cex=1.2)

#Multicolored plots
require(lattice)     			#Load Lattice package
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))

#Figure 3 
#Areas above the left-diagonal are discarded because they imply >1 total habitat availability
levelplot(format(EIGM, digits=5)~D1*D2|R1*R2, 
          data=Z[Z$Q1==0.5 & Z$Q2==0.5,],
          between=list(x=.5,y=.5),
          #at=seq(0,1.6,.2),contour=T,
          cuts=50, contour=F,
          col.regions=jet.colors(51),
          scales=list(tick.number=5, alternating=0),
          pretty=T, region=T, labels=T, 
          strip=FALSE, colorkey=F)

#Figure 4
#Triplot A
levelplot(format(EIGM, digits=15)~D1*D2|Q1*Q2, 
data=Z[Z$R1==0.5 & Z$R2==0.5,], 
between=list(x=.5,y=.5),
at=c(seq(0,1.6,.2)),contour=T,
#cuts=50, contour=F,
col.regions=jet.colors(51),
scales=list(tick.number=5, alternating=0),
pretty=T, region=T, labels=T, 
strip=FALSE, colorkey=F)

#Triplot B
levelplot(format(EIGM, digits=15)~D1*D2|Q1*Q2, 
data=Z[Z$R1==0.25 & Z$R2==0.25,], 
between=list(x=.5,y=.5),
at=c(seq(0,1.6,.2)),contour=T,
#cuts=50, contour=F,
col.regions=jet.colors(51),
scales=list(tick.number=5, alternating=0),
pretty=T, region=T, labels=T, 
strip=FALSE, colorkey=F)

#Triplot C
levelplot(format(EIGM, digits=5)~D1*D2|Q1*Q2, 
data=Z[Z$R1==0.01 & Z$R2==0.01,], 
between=list(x=.5,y=.5),
at=c(seq(0,1.6,.2)),contour=T,
#cuts=50, contour=F,
col.regions=jet.colors(100),
scales=list(tick.number=5, alternating=0),
pretty=T, region=T, labels=T, 
strip=FALSE, colorkey=F)









