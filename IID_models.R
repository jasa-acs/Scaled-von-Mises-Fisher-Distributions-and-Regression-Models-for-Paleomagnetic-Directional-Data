################################################################
################################################################
#############################Input and set up data
###############################################################

sphere=as.matrix(read.table(file="EIFdata3.csv",header=T,sep=","))
n=length(sphere[,1])
sphere_temp=matrix(0,1,5)
for (j in 1:n)
{
	if (sphere[j,1] == 1250) {sphere_temp=rbind(sphere_temp,sphere[j,])}
}
n=length(sphere_temp[,1])
sphere=sphere_temp[2:n,]
n=length(sphere[,1])

#axial dipole correction
inc_new=sphere[,3]+(180/pi)*atan(2*tan(50.12*pi/180))-(180/pi)*atan(2*tan(sphere[,4]*pi/180))
inc_new2=90-inc_new

depth=sphere[,1]
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


##recentre at northpole (preliminary transformation)

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
#hist(yp[,2],nclass=100)
#hist(yp[,3],nclass=100)
#plot(yp[,2],yp[,3],ylim=c(-1,1),xlim=c(-1,1))
#library(scatterplot3d)
#scatterplot3d(yp[,1],yp[,2],yp[,3], main="3D Scatterplot",angle=60)


#y is now the recentred response data
y=yp

#library(scatterplot3d)
#scatterplot3d(y[,1],y[,2],y[,3], main="3D Scatterplot",angle=60)

rm(mu,depth,Gamma,H,i,j,inc_new,inc_new2,initial,initialK,K_est,mag,mu_L,phi,S_orig,sphere,sphere_temp,theta,ybar,yp)

####################################################################
####################################################################
######compute initial values of mu
#####################################################################


initial=momu(y)
mag=initial$mag
ybar=initial$ybar

#save moment estimate of mu
mom_mu=ybar/mag

#save normalised spatial median
library(ICSNP) 
#spatial.median(y)/sqrt(sum(spatial.median(y)^2))
mu_spatial=t(t(spatial.median(y)/sqrt(sum(spatial.median(y)^2))))

rm(initial,mag,ybar)

####################################################################
####################################################################
######Fit IID models
#####################################################################


#Fit the SvMF model for a1=6
est=SvMFmodel(y,6)

#parameter estimates for a1=6 model (true values in Tables 3 and 4):
prML_mu=est$mu
prML_V=est$V
prML_kappa=est$kappav


#Fit the SvMF model for a1=1
est2=SvMFmodel(y,1)

#parameter estimates for a1=1 model (true values in Tables 3 and 4):
prML_mu1=est2$mu
prML_V1=est2$V
prML_kappa1=est2$kappav


#Fit the Kent model (uses saddlepoint aproximation)
Kentest=Kentmodel(y,1)

	
#estimate of mu for the Kent model
KML_mu=Kentest$mu



##################################################################

#save original data
ysphere=y

####################################################################



##################################################################
############density plots (Figure 2 in paper)
##################################################################

n=100000

#simulate sample from SvMF distribution with a1=6
y_project=simProject(prML_kappa,prML_V,prML_mu,6)$y

#simulate sample from Kent distribution:
sims=simKent(Kentest$kappa,Kentest$beta,KML_mu,Kentest$K)
#simulated sample from Kent distribution:
y_kent=sims$y

#simulate sample from SvMF distribution with a1=1
y_project1=simProject(prML_kappa1,prML_V1,prML_mu1,1)$y



par(mfrow = c(2,2),mar=c(4.5,4.5,0.5,0.5),cex.axis=1.5)


##############density plots component 2

mydensity1=density(y_kent[,2],bw=0.03)
mydensity2=density(y_project1[,2],bw=0.03)
mydensity3=density(y_project[,2],bw=0.03)


plot(density(ysphere[,2],bw=0.03),lty=2,lwd=2,main="",cex.lab=1.5)
rug(jitter(ysphere[,2]))
lines(mydensity1,lty=1,lwd=2)
lines(mydensity2,lty=4,lwd=2)
lines(mydensity3,lty=3,lwd=2)

##############density plots component 3


mydensity1=density(y_kent[,3],bw=0.04)
mydensity2=density(y_project1[,3],bw=0.04)
mydensity3=density(y_project[,3],bw=0.04)


plot(density(ysphere[,3],bw=0.04),lty=2,lwd=2,main="",cex.lab=1.5)
rug(jitter(ysphere[,3]))
lines(mydensity1,lty=1,lwd=2)
lines(mydensity2,lty=4,lwd=2)
lines(mydensity3,lty=3,lwd=2)



##############density plots component 1



mydensity1=density(y_kent[,1],bw=0.003)
mydensity2=density(y_project1[,1],bw=0.003)
mydensity3=density(y_project[,1],bw=0.003)


plot(density(ysphere[,1],bw=0.003),lty=2,lwd=2,main="",cex.lab=1.5)
rug(jitter(ysphere[,1]))
lines(mydensity1,lty=1,lwd=2)
lines(mydensity2,lty=4,lwd=2)
lines(mydensity3,lty=3,lwd=2)


######### KS tests


ks.test(y_kent[,2],ysphere[,2])
ks.test(y_project1[,2],ysphere[,2])
ks.test(y_project[,2],ysphere[,2])


ks.test(y_kent[,3],ysphere[,3])
ks.test(y_project1[,3],ysphere[,3])
ks.test(y_project[,3],ysphere[,3])


ks.test(y_kent[,1],ysphere[,1])
ks.test(y_project1[,1],ysphere[,1])
ks.test(y_project[,1],ysphere[,1])

#final plot

plot(ysphere[,2],ysphere[,3],ylim=c(-0.5,0.5),xlim=c(-0.5,0.5),cex=0.7,cex.lab=1.5,pch=16,xlab=expression(y[2]),ylab=expression(y[3]))


##########################################################################
#################Simulation###############################################
##This section generates one simulated sample only (seed=1). 
##To obtain the simulation results in Section 6 of the paper 
##the below code needs 
##to be repeated 1000 times using seeds 1, 2, ... , 1000. 
##Here we simulate from the SvMF model with a1=1 
##(the case a1=6 is similar and is omitted)
#######################################################################

##save true values conditioned on in the sim below:
true_mu=prML_mu1
ptrue_kappa=prML_kappa1
ptrue_V=prML_V1


#####################Repeat estimation below for different seeds:

#seed1=as.matrix(read.table(file="seed.txt",header=F,sep=","))
#set.seed(seed1[1])
set.seed(1)
n=50
#simulate data
y_project1=simProject(ptrue_kappa,ptrue_V,true_mu,1)$y
y=y_project1
#estimation:
initial=momu(y)
mag=initial$mag
ybar=initial$ybar
#save moment estimate of mu
mom_mu=ybar/mag
#save normalised spatial median
library(ICSNP) 
mu_spatial=t(t(spatial.median(y)/sqrt(sum(spatial.median(y)^2))))
rm(initial,mag,ybar)
#Fit the SvMF model for a1=6
est=SvMFmodel(y,6)
#parameter estimates for a1=6 model:
prML_mu=est$mu
prML_V=est$V
prML_kappa=est$kappav
#Fit the SvMF model for a1=1
est2=SvMFmodel(y,1)
#parameter estimates for a1=1 model:
prML_mu1=est2$mu
prML_V1=est2$V
prML_kappa1=est2$kappav
#Fit the Kent model (uses asymptotic normal aproximation)
Kentest=Kentmodel(y,0)
#estimate of mu for the Kent model
KML_mu=Kentest$mu
true=rbind(true_mu,true_mu,true_mu,true_mu,true_mu,ptrue_kappa,t(t(as.vector(ptrue_V))))
estimates=rbind(mom_mu,mu_spatial,prML_mu,prML_mu1,KML_mu,prML_kappa1,t(t(as.vector(prML_V1))))
#output parameter estimates and true values
write.table(t(estimates),"estimates.csv",sep=",",row.names=FALSE,col.names=FALSE)
write.table(t(true),"true.csv",sep=",",row.names=FALSE,col.names=FALSE)


