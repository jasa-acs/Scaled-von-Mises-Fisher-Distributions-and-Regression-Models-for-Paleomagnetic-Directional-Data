################################################################
################################################################
#############################Input and set up data
###############################################################

sphere=as.matrix(read.table(file="EIFdata3.csv",header=T,sep=","))
n=length(sphere[,1])
sphere_temp=matrix(0,1,5)
for (j in 1:n)
{
	if (sphere[j,1] <= 1500 && sphere[j,1] >= 0 ) {sphere_temp=rbind(sphere_temp,sphere[j,])}

}
n=length(sphere_temp[,1])
sphere=sphere_temp[2:n,]
n=length(sphere[,1])

#axial dipole correction
inc_new=sphere[,3]+(180/pi)*atan(2*tan(50.12*pi/180))-(180/pi)*atan(2*tan(sphere[,4]*pi/180))
inc_new2=90-inc_new

depth=sphere[,1]
depth=(depth-min(depth))/(max(depth)-min(depth))+1
theta=inc_new2*(pi/180)
phi=sphere[,2]*(pi/180)

#response data
y=matrix(0,n,3)
y[,1]=cos(theta)
y[,2]=sin(theta)*cos(phi)
y[,3]=sin(theta)*sin(phi)


#hist(y[,1],nclass=50)
#hist(y[,2],nclass=50)
#hist(y[,3],nclass=50)

##recentre at northpole

p=3

initial=momu(y)
mag=initial$mag
ybar=initial$ybar
S_orig=initial$S
#initial moment estimate of mu
mu=ybar/mag
initialK=momK(y,mu,S_orig,ybar)
K_est=initialK$K


H=diag(1,p)
H[,1]=t(t(mu))
H[1,]=t(mu)
mu_L=t(t(mu[2:p]))
H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))


Gamma=H%*%K_est

#recentred response data is yp
yp=matrix(0,n,p)
for (i in 1:n)
{
	yp[i,]=y[i,]%*%(Gamma)	

}

#plots
#hist(yp[,1])
#hist(yp[,2])
#hist(yp[,3],nclass=100)


#rotate to centre at positive orthant
mu=matrix(sqrt(1/p),p,1)
H=diag(1,p)
H[,1]=t(t(mu))
H[1,]=t(mu)
mu_L=t(t(mu[2:p]))
H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))

Gamma=H

yp2=matrix(0,n,p)
for (i in 1:n)
{
	yp2[i,]=yp[i,]%*%t(Gamma)	

}
#yp2 is the recentred response data (preliminary transformation applied)

#plots
#hist(yp2[,1],nclass=100)
#hist(yp2[,2],nclass=100)
#hist(yp2[,3],nclass=100)

rm(Gamma,H,inc_new,inc_new2,initial,initialK,K_est,mag,mu,mu_L,phi,S_orig,sphere,sphere_temp,theta,ybar)

#Set up the data for regression
					
#q is a p-1 vector containing the number of covariates used to model
#each of the p-1 coordinates associated with the mean direction.
q=matrix(0,p-1,1)
q[1]=5
q[2]=7


#total number of covariates:
Q=sum(q)

#qcum and qcum2 keep track of covariate positions
temp=0
qcum=matrix(0,p-1,1)
qcum2=matrix(0,p-1,1)
qcum2[1]=1
for (j in 1:sum(p,-1))
{
	temp=q[j]+temp
	qcum[j]=temp	
}
for (j in 2:sum(p,-1))
{
	qcum2[j]=qcum[j-1]+1
}


#x is the n by Q matrix of covariates
x=matrix(0,n,Q)

#x[,1]=1
#x[,2]=depth
#x[,3]=depth^2
#x[,4]=depth^3
#x[,5]=depth
#x[,6]=depth^2
#x[,7]=depth^3
#x[,8]=depth^4

x[,1]=1
x[,2]=depth
x[,3]=depth^2
x[,4]=depth^3
x[,5]=depth^4
x[,6]=1
x[,7]=depth
x[,8]=depth^2
x[,9]=depth^3
x[,10]=depth^4
x[,11]=depth^5
x[,12]=depth^6

##creating response matrix y
y=yp2
##covariate needed in covariance structure
vx=depth

#order of bessel functions
nu=sum(p/2,-1)


####################################################################
####################################################################
######compute initial values of parameters
#####################################################################

##initial values of regression coefficients
co=c(rep(0,Q))
##inital estimate of mu for each datapoint
mu_est=mu_calc(co)

###initial values of variance components
sigma1=0.05
sigma2=0.06
c1=0
#delta1=-1
#delta2=-0.2
delta1=0
delta2=0
m=sigma1*sigma2*(1-c1^2)^0.5
delta3=(delta1+delta2)/2
sigma3=sigma1/sigma2
delta4=(delta1-delta2)/2



####################################################################
####################################################################
######Fit SvMF regression model for a1=1 and Kent model
#####################################################################


#Fit the SvMF regression model for a1=1
tol1=0.00001
tol2=0.000001
est1=SvMFreg(co,sigma3,c1,delta4,m,delta3,tol1,tol2,1)

#Fit the Kent regression model
est2=Kentreg(co,sigma3,c1,delta4,m,delta3,tol1,tol2)

#save estimates for SvMF regression model for a1=1
a=est1$a
sigma3=est1$sigma3
c1=est1$c1
delta4=est1$delta4
m=est1$m
delta3=est1$delta3

####################################################################
####################################################################
#####Print regression model a1=1 estimates 
#####################################################################


estimates=rbind(a,sigma3,c1,delta4,m,delta3)
write.table(t(estimates),"estimates_regression_SvMF_a1_1.csv",sep=",",row.names=FALSE,col.names=FALSE)

#estimated regression coefficients
signif(a,3)

#observed information standard errors for regression coefficients
signif(sqrt(diag(solve(est1$info))),3)

#significant parameters are > 2 in absolute value:
a/sqrt(diag(solve(est1$info)))

#variance component estimates (m is sigma4 in paper)
sigma3
c1
delta4
m
delta3

#delta1
delta1=delta3+delta4

#delta2
delta2=delta3-delta4


####################################################################
####################################################################
#####Print Kent model estimates 
#####################################################################


estimates=rbind(est2$a,est2$sigma3,est2$c1,est2$delta4,est2$m,est2$delta3)
write.table(t(estimates),"estimates_regression_Kent.csv",sep=",",row.names=FALSE,col.names=FALSE)

#estimated regression coefficients
signif(est2$a,3)

#observed information standard errors for regression coefficients
signif(sqrt(diag(solve(est2$info))),3)

#significant parameters are > 2 in absolute value:
est2$a/sqrt(diag(solve(est2$info)))

#variance component estimates (m is sigma4 in paper)
est2$sigma3
est2$c1
est2$delta4
est2$m
est2$delta3




####################################################################
####################################################################
######Produce Figure 3 (case 2 a1=1)
#####################################################################

#rotate data back to north pole

mu=matrix(sqrt(1/p),p,1)
H=diag(1,p)
H[,1]=t(t(mu))
H[1,]=t(mu)
mu_L=t(t(mu[2:p]))
H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))

Gamma=H

yp=matrix(0,n,p)
for (i in 1:n)
{
	yp[i,]=y[i,]%*%(Gamma)	

}

#rotate mean back to northpole

mu_est2=matrix(0,n,p)
mu_est=mu_calc(a)
for (i in 1:n)
{
	mu_est2[i,]=mu_est[i,]%*%(Gamma)	

}

#residuals
res=est1$res

par(mfrow = c(3,2),mar=c(4.5,4.5,0.5,0.5),cex.axis=1.5)
plot(depth,yp[,1],cex=0.8,ylab=expression(y[1]),xlab=expression(x),cex.lab=1.5,pch=16)
lines(depth,mu_est2[,1],lwd=2)


plot(depth,yp[,2],cex=0.8,ylab=expression(y[2]),xlab=expression(x),cex.lab=1.5,pch=16)
lines(depth,mu_est2[,2],lwd=2,lty=1)
lines(smooth.spline(depth,yp[,2]),lty=2,lwd=2)

plot(depth,yp[,3],cex=0.8,ylab=expression(y[3]),xlab=expression(x),cex.lab=1.5,pch=16)
lines(depth,mu_est2[,3],lwd=2)
lines(smooth.spline(depth,yp[,3]),lty=2,lwd=2)


plot(depth,res[,2]/vx^(delta1),cex=0.8,ylab=expression(r[2]),xlab=expression(x),cex.lab=1.5,pch=16)
abline(0,0,lwd=2,lty=1)
var1=res[,2]/vx^(delta1)
#var1=res[,2]
lines(smooth.spline(depth,var1),lty=2,lwd=2)

plot(depth,res[,3]/vx^(delta2),cex=0.8,ylab=expression(r[3]),xlab=expression(x),cex.lab=1.5,pch=16)
abline(0,0,lwd=2,lty=1)
var1=res[,3]/vx^(delta2)
#var1=res[,3]
lines(smooth.spline(depth,var1),lty=2,lwd=2)

plot(res[,2]/vx^(delta1),res[,3]/vx^(delta2),cex=0.8,ylab=expression(r[3]),xlab=expression(r[2]),ylim=c(-0.5,0.5),xlim=c(-0.5,0.5),cex.lab=1.5,pch=16)



####################################################################
####################################################################
######Fit SvMF regression model for a1=6 (uses inital values from 
######previous successively)
#####################################################################


#Fit the SvMF regression model for a1=2
est1=SvMFreg(a,sigma3,c1,delta4,m,delta3,tol1,tol2,2)
a=est1$a
sigma3=est1$sigma3
c1=est1$c1
delta4=est1$delta4
m=est1$m
delta3=est1$delta3
#Fit the SvMF regression model for a1=4
est1=SvMFreg(a,sigma3,c1,delta4,m,delta3,tol1,tol2,4)
a=est1$a
sigma3=est1$sigma3
c1=est1$c1
delta4=est1$delta4
m=est1$m
delta3=est1$delta3
#Fit the SvMF regression model for a1=6
est1=SvMFreg(a,sigma3,c1,delta4,m,delta3,tol1,tol2,6)
a=est1$a
sigma3=est1$sigma3
c1=est1$c1
delta4=est1$delta4
m=est1$m
delta3=est1$delta3

estimates=rbind(a,sigma3,c1,delta4,m,delta3)
write.table(t(estimates),"estimates_regression_SvMF_a1_6.csv",sep=",",row.names=FALSE,col.names=FALSE)

#estimated regression coefficients
signif(a,3)

#observed information standard errors for regression coefficients
signif(sqrt(diag(solve(est1$info))),3)

#significant parameters are > 2 in absolute value:
a/sqrt(diag(solve(est1$info)))

#variance component estimates (m is sigma4 in paper)
sigma3
c1
delta4
m
delta3

