#input data
sphere=as.matrix(read.table(file="EIFdata3.csv",header=T,sep=","))
n=length(sphere[,1])

#axial dipole correction
inc_new=sphere[,3]+(180/pi)*atan(2*tan(50.12*pi/180))-(180/pi)*atan(2*tan(sphere[,4]*pi/180))
#inc_new2=90-inc_new


depth=sphere[,1]


dec=sphere[,2]
dec2=sphere[,2]
for (j in 1:n)
{
	if (dec[j] > 200) {dec2[j]=dec[j]-360}
}
par(mfrow = c(2,1),mar=c(4.5,4.5,0.5,0.5),cex.axis=1.5)
plot(depth,dec2,xlab="Age (in years)",ylab="Declination",cex=0.7,cex.lab=1.5,pch=16)
plot(depth,inc_new,xlab="Age (in years)",ylab="Inclination",cex=0.7,cex.lab=1.5,pch=16)


#dev.off()