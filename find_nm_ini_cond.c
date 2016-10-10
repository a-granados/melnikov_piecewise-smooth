# include <stdlib.h>
#include <stdio.h>
#include <math.h>


//This is an extension of find_nm_ini_cond for the conservative case. We implement here a generalization to
//(n,m)-periodic orbits.

/*Our purpose here is to find an initial condition of (n,m)-periodic orbit of the perturbed system near the point
(\bar{y}_0,\bar{t}_0). \bar{y}_0 is given as initial condition of a nT/m-periodic orbit of the unperturbed system, while \bar{t}_0 is given by the Subharmonic Melnikov function. These are given in the file "system.dat" and saved in by0 and bt0, respectively.
*/

double dist(double *x1,double *x2);
void obtain_new_point(double *newpoint,double *oldpoint,int n, int m,double eps,double r,double omega);
void find_impacts(double **impacts,double *inipoint,int m,double eps,double r, double omega);
void get_next_impact(double *finalimpact, double *iniimpact,double eps,double omega);
void calcinvDF(double **invDF,double **impacts,int m,double eps,double r,double omega);
void mmultiply(double **A,double **B,double **C);

int main (int argc, char *argv[])
{
int i,n,m,maxiter;
double by0,bt0,tol1,eps,omega,T,aux,r;
double delta,ratio,te;
double *newpoint,*oldpoint;
double **finalimpacts;
static double pi=3.14159265358979323846;
FILE *fout, *fin;
FILE *fout2;

finalimpacts=(double**)malloc(50*sizeof(double*));
for (i=0;i<50;i++)
	{
	finalimpacts[i]=(double*)malloc(2*sizeof(double));
	}

if (argc==4){
	fin=fopen(argv[1],"r");
	fout=fopen(argv[2],"w");
	fout2=fopen(argv[3],"w");	
}
else if (argc==1){
	fin=fopen("system.dat","r");
	fout=fopen("output.dat","w");
	fout2=fopen("impacts.dat","w");
}
else{
	printf("Wrong number of input files.\n Please provide system data, output data and impacts data files.\n");
	exit(0);
}
tol1=1e-14;
//maxiter=1000000;
maxiter=10000;

newpoint=(double*)malloc(2*sizeof(double));
oldpoint=(double*)malloc(2*sizeof(double));


//Now we have two options depending on what are the values in fin:

//fscanf(fin,"%i %i %lf %lf %lf %lf %lf",&n,&m,&by0,&bt0,&eps,&r,&omega);
///*
fscanf(fin,"%i %i %lf %lf %lf %lf %lf",&n,&m,&by0,&bt0,&delta,&ratio,&omega);
te=1;
eps=te*delta;
//r=pow(1-ratio*te*delta,1.0/(4.0*m));
r=1-ratio*te*delta;
//*/

printf("r=%20.20f eps=%20.20f\n\n",r,eps);
//exit(0);

//Let us compute the initial y0. In principle it is given in the file, but we can recompute it.
T=2*pi/omega;
aux=exp(n*T/(2*m));
//by0=(aux-1)/(1+aux);

//First evaluation
oldpoint[0]=by0;
oldpoint[1]=bt0;
printf("%20.20f %20.20f\n",oldpoint[0],oldpoint[1]);
obtain_new_point(newpoint,oldpoint,n,m,eps,r,omega);
printf("%20.20f %20.20f\n",newpoint[0],newpoint[1]);
i=1;
//we start the newton process to find the zero of F_{m,eps}(y0,t0)=(0,0) near (by0,ty0)
while (dist(oldpoint,newpoint) > tol1 && i<maxiter)
	{
	oldpoint[0]=newpoint[0];
	oldpoint[1]=newpoint[1];
	obtain_new_point(newpoint,oldpoint,n,m,eps,r,omega);
	printf("%20.20f %20.20f\n",newpoint[0],newpoint[1]);
	i++;
	}
if (i>=maxiter)
	{
	printf("Newton method in main diverges\n");
	}
//printf("%20.20f %20.20f %i\n",oldpoint[0]-newpoint[0],oldpoint[1]-newpoint[1],i);


fprintf(fout,"%20.20f %20.20f\n",newpoint[0],newpoint[1]);


find_impacts(finalimpacts,newpoint,m,eps,r,omega);

//printf("Errors at the final point:\n  %g %2g\n",finalimpacts[0][0]-finalimpacts[2*m][0],finalimpacts[0][1]-finalimpacts[2*m][1]+n*2*pi/omega);


for (i=0;i<=2*m;i++){
	fprintf(fout2,"%20.20f %20.20f\n",finalimpacts[i][0],finalimpacts[i][1]);
	}

free(newpoint);
free(oldpoint);
fclose(fout);
fclose(fout2);
fclose(fin);
for (i=0;i<50;i++)
	{
	free(finalimpacts[i]);
	}

free(finalimpacts);
}

double dist(double *x1,double *x2){

return(sqrt(pow(x1[0]-x2[0],2)+pow(x1[1]-x2[1],2)));

}

void obtain_new_point(double *newpoint,double *oldpoint,int n, int m,double eps,double r,double omega){
//This function evaluates the equation F_{m,eps}(y0,t0).

double **impacts,T;
double **invDF;
int i;
static double pi=3.14159265358979323846;

T=2*pi/omega;

impacts=(double**)malloc(100*sizeof(double*));
for (i=0;i<100;i++)
	{
	impacts[i]=(double*)malloc(2*sizeof(double));
	}
invDF=(double**)malloc(2*sizeof(double*));
invDF[0]=(double*)malloc(2*sizeof(double));
invDF[1]=(double*)malloc(2*sizeof(double));

find_impacts(impacts,oldpoint,m,eps,r,omega);

/*
for (i=0;i<=2*m;i++){
	printf("Hi ");
	printf("%lf %lf\n",impacts[i][0],impacts[i][1]);
	}
*/

calcinvDF(invDF,impacts,m,eps,r,omega);

//newpoint[0]=-(invDF[0][0]*impacts[2*m][0]+invDF[0][1]*impacts[2*m][1])+oldpoint[0]-impacts[0][0];
//newpoint[1]=-(invDF[1][0]*impacts[2*m][0]+invDF[1][1]*impacts[2*m][1])+oldpoint[1]-(impacts[0][1]+n*T);

newpoint[0]=-(invDF[0][0]*(impacts[2*m][0]-impacts[0][0])+invDF[0][1]*(impacts[2*m][1]-impacts[0][1]-n*T))+oldpoint[0];
newpoint[1]=-(invDF[1][0]*(impacts[2*m][0]-impacts[0][0])+invDF[1][1]*(impacts[2*m][1]-impacts[0][1]-n*T))+oldpoint[1];


//This is to try to plot something that has to be close to the Melnikov function:
/*
FILE *MelFitx;
double aux;
MelFitx=fopen("MelFitx.dat","a");
aux=-4/(pow(omega,2)+1)*cos(omega*oldpoint[1]);
//aux=-4*cos(omega*oldpoint[1]);
fprintf(MelFitx,"%20.20f %20.20f\n",oldpoint[1],(pow(impacts[2*m][0],2)-pow(impacts[0][0],2))/(2*eps*aux));
//fprintf(MelFitx,"%20.20f %20.20f\n",oldpoint[1],(impacts[2*m][0]-impacts[0][0])/eps);
fclose(MelFitx);
*/

for (i=0;i<100;i++)
	{
	free(impacts[i]);
	}
free(impacts);
free(invDF[0]);
free(invDF[1]);
free(invDF);
}

void find_impacts(double **impacts,double *inipoint,int m,double eps,double r,double omega){

int i;

impacts[0][0]=inipoint[0];
impacts[0][1]=inipoint[1];
get_next_impact(impacts[1],impacts[0],eps,omega);
impacts[1][0]=r*impacts[1][0];

for (i=1;i<2*m;i++)
	{
	get_next_impact(impacts[i+1],impacts[i],eps,omega);
	impacts[i+1][0]=r*impacts[i+1][0];
	}
}


///Use the following function for H_1=x*cos(omega*t)
///*
void get_next_impact(double *finalimpact, double *iniimpact,double eps,double omega){

double tol,err,newtime,oldtime;
double y0,t0,C1p,C1m,C2p,C2m;
double error,df;
int i,maxiter;

tol=1e-15;
maxiter=2000;

//Careful! The perturbation parameter in the explicit solution that we are using here
//is the original one scaled by omega^2+1!!!
eps=eps/(pow(omega,2)+1);

y0=iniimpact[0];
t0=iniimpact[1];
C1p=(y0 - eps*cos(omega*t0) + eps*omega*sin(omega*t0) - 1)/2*exp(-t0);
C1m=(y0 - eps*cos(omega*t0) + eps*omega*sin(omega*t0) + 1)/2*exp(-t0);
C2p=(-y0 -eps*cos(omega*t0) - eps*omega*sin(omega*t0) - 1)/2*exp(t0);
C2m=(-y0 -eps*cos(omega*t0) - eps*omega*sin(omega*t0) + 1)/2*exp(t0);

//Half of the period of the unperturbed orbit going through y0
oldtime=t0+2*log(sqrt(1-pow(y0,2.0))/(1.0-fabs(y0)));


if (y0>0){
	error=C1p*exp(oldtime)+C2p*exp(-oldtime)+eps*cos(omega*oldtime)+1;
	df=C1p*exp(oldtime)-C2p*exp(-oldtime)-eps*omega*sin(omega*oldtime);
	newtime=-error/df+oldtime;
	}
else {
	error=C1m*exp(oldtime)+C2m*exp(-oldtime)+eps*cos(omega*oldtime)-1;
	df=C1m*exp(oldtime)-C2m*exp(-oldtime)-eps*omega*sin(omega*oldtime);
	newtime=-error/df+oldtime;
	}
i=1;
while (fabs(error)>tol && i<maxiter){
	oldtime=newtime;
	if (y0>0){
		error=C1p*exp(oldtime)+C2p*exp(-oldtime)+eps*cos(omega*oldtime)+1;
		df=C1p*exp(oldtime)-C2p*exp(-oldtime)-eps*omega*sin(omega*oldtime);
		newtime=-error/df+oldtime;
		}
	else{
		error=C1m*exp(oldtime)+C2m*exp(-oldtime)+eps*cos(omega*oldtime)-1;
		df=C1m*exp(oldtime)-C2m*exp(-oldtime)-eps*omega*sin(omega*oldtime);
		newtime=-error/df+oldtime;
		}
	i++;
	}

if (i>=maxiter){
	printf("Newton in get_next_impact diverges\n");
	}

//printf("Impact error: %g\n",error);
finalimpact[0]=df;
finalimpact[1]=newtime;

}
//*/


///Use the following function for H_1=x*(cos(omega*t)+cos(k*omega*t))

/*
void get_next_impact(double *finalimpact, double *iniimpact,double eps,double omega){

double tol,err,newtime,oldtime;
double y0,t0,C1p,C1m,C2p,C2m,gamma1,gamma2,k,t;
double error,df;
int i,maxiter;

tol=1e-15;
maxiter=2000;

k=3.0;
gamma1=eps/(1+pow(omega,2));
gamma2=eps/(1+pow(k*omega,2));
//t=oldtime;

y0=iniimpact[0];
t0=iniimpact[1];
C1p=(y0-gamma1*cos(omega*t0)+omega*gamma1*sin(omega*t0)-gamma2*cos(k*omega*t0)+gamma2*k*omega*sin(k*omega*t0)-1)*exp(-t0)/2;
C1m=(y0-gamma1*cos(omega*t0)+omega*gamma1*sin(omega*t0)-gamma2*cos(k*omega*t0)+gamma2*k*omega*sin(k*omega*t0)+1)*exp(-t0)/2;
C2p=(-y0-gamma1*cos(omega*t0)-omega*gamma1*sin(omega*t0)-gamma2*cos(k*omega*t0)-gamma2*k*omega*sin(k*omega*t0)-1)*exp(t0)/2;
C2m=(-y0-gamma1*cos(omega*t0)-omega*gamma1*sin(omega*t0)-gamma2*cos(k*omega*t0)-gamma2*k*omega*sin(k*omega*t0)+1)*exp(t0)/2;

//Half of the period of the unperturbed orbit going through y0
oldtime=t0+2*log(sqrt(1-pow(y0,2.0))/(1.0-fabs(y0)));

if (y0>0){
	error=C1p*exp(oldtime)+C2p*exp(-oldtime)+1+gamma1*cos(omega*oldtime)+gamma2*cos(k*omega*oldtime);
	df=C1p*exp(oldtime)-C2p*exp(-oldtime)-gamma1*omega*sin(omega*oldtime)-gamma2*k*omega*sin(k*omega*oldtime);
	newtime=-error/df+oldtime;
	}
else {
	error=C1m*exp(oldtime)+C2m*exp(-oldtime)-1+gamma1*cos(omega*oldtime)+gamma2*cos(k*omega*oldtime);
	df=C1m*exp(oldtime)-C2m*exp(-oldtime)-gamma1*omega*sin(omega*oldtime)-gamma2*k*omega*sin(k*omega*oldtime);
	newtime=-error/df+oldtime;
	}
i=1;
while (fabs(error)>tol && i<maxiter){
	oldtime=newtime;
	if (y0>0){
		error=C1p*exp(oldtime)+C2p*exp(-oldtime)+1+gamma1*cos(omega*oldtime)+gamma2*cos(k*omega*oldtime);
		df=C1p*exp(oldtime)-C2p*exp(-oldtime)-gamma1*omega*sin(omega*oldtime)-gamma2*k*omega*sin(k*omega*oldtime);
		newtime=-error/df+oldtime;
		}
	else{
		error=C1m*exp(oldtime)+C2m*exp(-oldtime)-1+gamma1*cos(omega*oldtime)+gamma2*cos(k*omega*oldtime);
		df=C1m*exp(oldtime)-C2m*exp(-oldtime)-gamma1*omega*sin(omega*oldtime)-gamma2*k*omega*sin(k*omega*oldtime);
		newtime=-error/df+oldtime;
		}
	i++;
	}

if (i>=maxiter){
	printf("Newton in get_next_impact diverges\n");
	}

//printf("Impact error: %g\n",error);
finalimpact[0]=df;
finalimpact[1]=newtime;

}
*/

////////////EVERYTHING WHAT FOLLOWS HAS BEEN PUT IN differential_matrices.c
//////////IT SHOULD BE REPLACE REPLACED BY SOME GENERAL FUNCTION calcinvDF






//*/

/*
void calcinvDF(double **invDF,double **impacts,int m,double eps){
//We compute here the first order of the differential of F, wher F=0 is the funcion we want to solve.
//In fact, we will use a Whittaker method.
//Then we invert this matrix.

int i;
double s,alphap1,alphap2,sum;

s=0;
for (i=0;i<=m-1;i++)
	{
	alphap1=4/(1-pow(impacts[0][2*i],2));
	alphap2=4/(1-pow(impacts[0][2*i+1],2));
	s=s+(-alphap2+alphap1)/2;
	}

//invDF[0][0]=0;
//invDF[0][1]=1/s;
//invDF[1][0]=1/eps;
//invDF[1][1]=0;

invDF[0][0]=-eps;
invDF[0][1]=0;
invDF[1][0]=0;
invDF[1][1]=-eps;
}
*/

/*
void calcinvDF(double **invDF,double **impacts,int m,double eps,double r,double omega){
//Here we compoute the differential numerically

double delta,*F,det;
double F1d1f,F1d2f,F2d1f,F2d2f;
double F1d1b,F1d2b,F2d1b,F2d2b;
double *dx1f,*dx2f;
double *dx1b,*dx2b;

double **aux,**DF;
int i,j;


aux=(double**)malloc(100*sizeof(double*));
for (i=0;i<100;i++)
	{
	aux[i]=(double*)malloc(2*sizeof(double));
	}

DF=(double**)malloc(2*sizeof(double*));
DF[0]=(double*)malloc(2*sizeof(double));
DF[1]=(double*)malloc(2*sizeof(double));
F=(double*)malloc(2*sizeof(double));
dx1f=(double*)malloc(2*sizeof(double));
dx2f=(double*)malloc(2*sizeof(double));
dx1b=(double*)malloc(2*sizeof(double));
dx2b=(double*)malloc(2*sizeof(double));


delta=1e-10;

j=2*m;//For debugging reasons. Just to choose at which impact we want the differential. Normally this should be equal to 2*m
F[0]=impacts[j][0];
F[1]=impacts[j][1];

dx1f[0]=impacts[0][0]+delta;
dx1f[1]=impacts[0][1];
dx2f[0]=impacts[0][0];
dx2f[1]=impacts[0][1]+delta;
dx1b[0]=impacts[0][0]-delta;
dx1b[1]=impacts[0][1];
dx2b[0]=impacts[0][0];
dx2b[1]=impacts[0][1]-delta;

find_impacts(aux,dx1f,m,eps,r,omega);
F1d1f=aux[j][0];
F2d1f=aux[j][1];
find_impacts(aux,dx2f,m,eps,r,omega);
F1d2f=aux[j][0];
F2d2f=aux[j][1];

find_impacts(aux,dx1b,m,eps,r,omega);
F1d1b=aux[j][0];
F2d1b=aux[j][1];
find_impacts(aux,dx2b,m,eps,r,omega);
F1d2b=aux[j][0];
F2d2b=aux[j][1];


DF[0][0]=(F1d1f-F1d1b)/(2*delta)-1;
DF[0][1]=(F1d2f-F1d2b)/(2*delta);
DF[1][0]=(F2d1f-F2d1b)/(2*delta);
DF[1][1]=(F2d2f-F2d2b)/(2*delta)-1;

//printf("%g %g\n%g %g\n\n",DF[0][0]+1,DF[0][1],DF[1][0],DF[1][1]+1);
//exit(0);

det=DF[0][0]*DF[1][1]-DF[1][0]*DF[0][1];

//printf("%g\n\n",det);

invDF[0][0]=DF[1][1]/det;
invDF[0][1]=-DF[0][1]/det;
invDF[1][0]=-DF[1][0]/det;
invDF[1][1]=DF[0][0]/det;


for (i=0;i<100;i++)
	{
	free(aux[i]);
	}
free(aux);
free(DF[0]);
free(DF[1]);
free(DF);
free(F);
free(dx1f);
free(dx2f);
free(dx1b);
free(dx2b);
}
*/


void calcinvDF(double **invDF,double **impacts,int m,double eps,double r, double omega){
//Here we compute the differential analytcally for the linear case.



///S'ha de multiplicar la diferencial per r on toqui!!!
double det,y0,t0,t;
double **DF,**aux,**DR;
int i;

//Careful! The perturbation parameter in the explicit solution that we are using here
//is the original one scaled by omega^2+1!!!
eps=eps/(pow(omega,2)+1);


DF=(double**)malloc(2*sizeof(double*));
DF[0]=(double*)malloc(2*sizeof(double));
DF[1]=(double*)malloc(2*sizeof(double));
aux=(double**)malloc(2*sizeof(double*));
aux[0]=(double*)malloc(2*sizeof(double));
aux[1]=(double*)malloc(2*sizeof(double));
DR=(double**)malloc(2*sizeof(double*));
DR[0]=(double*)malloc(2*sizeof(double));
DR[1]=(double*)malloc(2*sizeof(double));

DF[0][0]=1;
DF[0][1]=0;
DF[1][0]=0;
DF[1][1]=1;

DR[0][0]=r;
DR[0][1]=0;
DR[1][0]=0;
DR[1][1]=1;

for (i=1;i<=2*m;i++){
	y0=impacts[i-1][0];
	t0=impacts[i-1][1];
	t=impacts[i][1];
	if (y0>0){
//		aux[0][0]= exp(-t0 + t) / 0.2e1 + exp(t0 - t) / 0.2e1;
//		aux[0][1]=exp(-t0 + t) * eps * omega * omega * cos(omega * t0) / 0.2e1 - exp(-t0 + t) * y0 / 0.2e1 + exp(-t + t) * eps * cos(omega * t0) / 0.2e1 + exp(-t0 + t) / 0.2e1 + exp(t0 - t) * eps * omega * omega * cos(omega * t0) / 0.2e1 + exp(t0 - t) * y0 / 0.2e1 + exp(t0 - t) * eps * cos(omega * t0) / 0.2e1 + exp(t0 - t) / 0.2e1;
		aux[0][0]=exp(-t0 + t) / 0.2e1 + exp(t0 - t) / 0.2e1 - (exp(-t0 + t) * y0 / 0.2e1 - exp(-t0 + t) * eps * cos(omega * t0) / 0.2e1 + exp(-t0 + t) * eps * omega * sin(omega * t0) / 0.2e1 - exp(-t0 + t) / 0.2e1 - exp(t0 - t) * y0 / 0.2e1 - exp(t0 - t) * eps * cos(omega * t0) / 0.2e1 - exp(t0 - t) * eps * omega * sin(omega * t0) / 0.2e1 - exp(t0 - t) / 0.2e1 - eps * omega * omega * cos(omega * t)) * (exp(-t0 + t) - exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * eps * cos(omega * t0) + exp(-t0 + t) * eps * omega * sin(omega * t0) - exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * eps * cos(omega * t0) + exp(t0 - t) * eps * omega * sin(omega * t0) + exp(t0 - t) - 0.2e1 * eps * omega * sin(omega * t));
		aux[0][1]=exp(-t0 + t) * eps * omega * omega * cos(omega * t0) / 0.2e1 - exp(-t0 + t) * y0 / 0.2e1 + exp(-t0 + t) * eps * cos(omega * t0) / 0.2e1 + exp(-t0 + t) / 0.2e1 + exp(t0 - t) * eps * omega * omega * cos(omega * t0) / 0.2e1 + exp(t0 - t) * y0 / 0.2e1 + exp(t0 - t) * eps * cos(omega * t0) / 0.2e1 + exp(t0 - t) / 0.2e1 - (exp(-t0 + t) * y0 / 0.2e1 - exp(-t0 + t) * eps * cos(omega * t0) / 0.2e1 + exp(-t0 + t) * eps * omega * sin(omega * t0) / 0.2e1 - exp(-t0 + t) / 0.2e1 - exp(t0 - t) * y0 / 0.2e1 - exp(t0 - t) * eps * cos(omega * t0) / 0.2e1 - exp(t0 - t) * eps * omega * sin(omega * t0) / 0.2e1 - exp(t0 - t) / 0.2e1 - eps * omega * omega * cos(omega * t)) * (exp(-t0 + t) * eps * omega * omega * cos(omega * t0) - exp(-t0 + t) * y0 + exp(-t0 + t) * eps * cos(omega * t0) + exp(-t0 + t) - exp(t0 - t) * eps * omega * omega * cos(omega * t0) - exp(t0 - t) * y0 - exp(t0 - t) * eps * cos(omega * t0) - exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * eps * cos(omega * t0) + exp(-t0 + t) * eps * omega * sin(omega * t0) - exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * eps * cos(omega * t0) + exp(t0 - t) * eps * omega * sin(omega * t0) + exp(t0 - t) - 0.2e1 * eps * omega * sin(omega * t));

		aux[1][0]= -(exp(-t0 + t) - exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * eps * cos(omega * t0) + eps * omega * sin(omega * t0) * exp(-t0 + t) - exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * eps * cos(omega * t0) + eps * omega * sin(omega * t0) * exp(t0 - t) + exp(t0 - t) - 0.2e1 * eps * omega * sin(omega * t));
		aux[1][1]=-(exp(-t0 + t) * eps * omega * omega * cos(omega * t0) - exp(-t0 + t) * y0 + exp(-t0 + t) * eps * cos(omega * t0) + exp(-t0 + t) - exp(t0 - t) * eps * omega * omega * cos(omega * t0) - exp(t0 - t) * y0 - exp(t0 - t) * eps * cos(omega * t0) - exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * eps * cos(omega * t0) + eps * omega * sin(omega * t0) * exp(-t0 + t) - exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * eps * cos(omega * t0) + eps * omega * sin(omega * t0) * exp(t0 - t) + exp(t0 - t) - 0.2e1 * eps * omega * sin(omega * t));
	}
	else{
//		aux[0][0]=exp(-t0 + t) / 0.2e1 + exp(t0 - t) / 0.2e1;
//		aux[0][1]=exp(-t0 + t) * eps * omega * omega * cos(omega * t0) / 0.2e1 - exp(-t0 + t) * y0 / 0.2e1 + exp(-t0 + t) * eps * cos(omega * t0) / 0.2e1 - exp(-t0 + t) / 0.2e1 + exp(t0 - t) * eps * omega * omega * cos(omega * t0) / 0.2e1 + exp(t0 - t) * y0 / 0.2e1 + exp(t0 - t) * eps * cos(omega * t0) / 0.2e1 - exp(t0 - t) / 0.2e1;
		aux[0][0]=exp(-t0 + t) / 0.2e1 + exp(t0 - t) / 0.2e1 - (exp(-t0 + t) * y0 / 0.2e1 - exp(-t0 + t) * eps * cos(omega * t0) / 0.2e1 + exp(-t0 + t) * eps * omega * sin(omega * t0) / 0.2e1 + exp(-t0 + t) / 0.2e1 - exp(t0 - t) * y0 / 0.2e1 - exp(t0 - t) * eps * cos(omega * t0) / 0.2e1 - exp(t0 - t) * eps * omega * sin(omega * t0) / 0.2e1 + exp(t0 - t) / 0.2e1 - eps * omega * omega * cos(omega * t)) * (exp(-t0 + t) - exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * eps * cos(omega * t0) + exp(-t0 + t) * eps * omega * sin(omega * t0) + exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * eps * cos(omega * t0) + exp(t0 - t) * eps * omega * sin(omega * t0) - exp(t0 - t) - 0.2e1 * eps * omega * sin(omega * t));
		aux[0][1]= exp(-t0 + t) * eps * omega * omega * cos(omega * t0) / 0.2e1 - exp(-t0 + t) * y0 / 0.2e1 + exp(-t0 + t) * eps * cos(omega * t0) / 0.2e1 - exp(-t0 + t) / 0.2e1 + exp(t0 - t) * eps * omega * omega * cos(omega * t0) / 0.2e1 + exp(t0 - t) * y0 / 0.2e1 + exp(t0 - t) * eps * cos(omega * t0) / 0.2e1 - exp(t0 - t) / 0.2e1 - (exp(-t0 + t) * y0 / 0.2e1 - exp(-t0 + t) * eps * cos(omega * t0) / 0.2e1 + exp(-t0 + t) * eps * omega * sin(omega * t0) / 0.2e1 + exp(-t0 + t) / 0.2e1 - exp(t0 - t) * y0 / 0.2e1 - exp(t0 - t) * eps * cos(omega * t0) / 0.2e1 - exp(t0 - t) * eps * omega * sin(omega * t0) / 0.2e1 + exp(t0 - t) / 0.2e1 - eps * omega * omega * cos(omega * t)) * (exp(-t0 + t) * eps * omega * omega * cos(omega * t0) - exp(-t0 + t) * y0 + exp(-t0 + t) * eps * cos(omega * t0) - exp(-t0 + t) - exp(t0 - t) * eps * omega * omega * cos(omega * t0) - exp(t0 - t) * y0 - exp(t0 - t) * eps * cos(omega * t0) + exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * eps * cos(omega * t0) + exp(-t0 + t) * eps * omega * sin(omega * t0) + exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * eps * cos(omega * t0) + exp(t0 - t) * eps * omega * sin(omega * t0) - exp(t0 - t) - 0.2e1 * eps * omega * sin(omega * t));

		aux[1][0]=-(exp(-t0 + t) - exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * eps * cos(omega * t0) + eps * omega * sin(omega * t0) * exp(-t0 + t) + exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * eps * cos(omega * t0) + eps * omega * sin(omega * t0) * exp(t0 - t) - exp(t0 - t) - 0.2e1 * eps * omega * sin(omega * t));
		aux[1][1]=-(exp(-t0 + t) * eps * omega * omega * cos(omega * t0) - exp(-t0 + t) * y0 + exp(-t0 + t) * eps * cos(omega * t0) - exp(-t0 + t) - exp(t0 - t) * eps * omega * omega * cos(omega * t0) - exp(t0 - t) * y0 - exp(t0 - t) * eps * cos(omega * t0) + exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * eps * cos(omega * t0) + eps * omega * sin(omega * t0) * exp(-t0 + t) + exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * eps * cos(omega * t0) + eps * omega * sin(omega * t0) * exp(t0 - t) - exp(t0 - t) - 0.2e1 * eps * omega * sin(omega * t));
	}
	mmultiply(aux,DR,aux);
	mmultiply(DF,aux,DF);
//	printf("%g %g\n%g %g\n\n",DF[0][0],DF[0][1],DF[1][0],DF[1][1]);
//	exit(0);
}

DF[0][0]=DF[0][0]-1;
DF[1][1]=DF[1][1]-1;

//printf("%g %g\n%g %g\n\n",DF[0][0],DF[0][1],DF[1][0],DF[1][1]);

det=DF[0][0]*DF[1][1]-DF[1][0]*DF[0][1];

//printf("%g\n",det);

invDF[0][0]=DF[1][1]/det;
invDF[0][1]=-DF[0][1]/det;
invDF[1][0]=-DF[1][0]/det;
invDF[1][1]=DF[0][0]/det;


free(DF[0]);
free(DF[1]);
free(DF);
free(aux[0]);
free(aux[1]);
free(aux);
free(DR[0]);
free(DR[1]);
free(DR);
}


/*
void calcinvDF(double **invDF,double **impacts,int m,double eps,double r, double omega){
//Here we compute the differential analytcally for the linear case using the perturbation
//H_1=x*(cos(omega*t)+cos(k*omega*t))


double det,y0,t0,t;
double **DF,**aux,**DR;
double gamma1,gamma2,k;
int i;


gamma1=eps/(1+pow(omega,2));
gamma2=eps/(1+pow(k*omega,2));
k=3.0;

DF=(double**)malloc(2*sizeof(double*));
DF[0]=(double*)malloc(2*sizeof(double));
DF[1]=(double*)malloc(2*sizeof(double));
aux=(double**)malloc(2*sizeof(double*));
aux[0]=(double*)malloc(2*sizeof(double));
aux[1]=(double*)malloc(2*sizeof(double));
DR=(double**)malloc(2*sizeof(double*));
DR[0]=(double*)malloc(2*sizeof(double));
DR[1]=(double*)malloc(2*sizeof(double));

DF[0][0]=1;
DF[0][1]=0;
DF[1][0]=0;
DF[1][1]=1;

DR[0][0]=r;
DR[0][1]=0;
DR[1][0]=0;
DR[1][1]=1;

for (i=1;i<=2*m;i++){
	y0=impacts[i-1][0];
	t0=impacts[i-1][1];
	t=impacts[i][1];
	if (y0>0){

		aux[0][0]=exp(-t0 + t) / 0.2e1 + exp(t0 - t) / 0.2e1 - (exp(-t0 + t) * y0 / 0.2e1 - exp(-t0 + t) * gamma1 * cos(omega * t0) / 0.2e1 + exp(-t0 + t) * omega * gamma1 * sin(omega * t0) / 0.2e1 - exp(-t0 + t) * gamma2 * cos(k * omega * t0) / 0.2e1 + exp(-t0 + t) * k * omega * gamma2 * sin(k * omega * t0) / 0.2e1 - exp(-t0 + t) / 0.2e1 - exp(t0 - t) * y0 / 0.2e1 - exp(t0 - t) * gamma1 * cos(omega * t0) / 0.2e1 - exp(t0 - t) * omega * gamma1 * sin(omega * t0) / 0.2e1 - exp(t0 - t) * gamma2 * cos(k * omega * t0) / 0.2e1 - exp(t0 - t) * k * omega * gamma2 * sin(k * omega * t0) / 0.2e1 - exp(t0 - t) / 0.2e1 - gamma1 * omega * omega * cos(omega * t) - gamma2 * k * k * omega * omega * cos(k * omega * t)) * (exp(-t0 + t) - exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * gamma1 * cos(omega * t0) + exp(-t0 + t) * omega * gamma1 * sin(omega * t0) - exp(-t0 + t) * gamma2 * cos(k * omega * t0) + exp(-t0 + t) * k * omega * gamma2 * sin(k * omega * t0) - exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * gamma1 * cos(omega * t0) + exp(t0 - t) * omega * gamma1 * sin(omega * t0) + exp(t0 - t) * gamma2 * cos(k * omega * t0) + exp(t0 - t) * k * omega * gamma2 * sin(k * omega * t0) + exp(t0 - t) - 0.2e1 * gamma1 * omega * sin(omega * t) - 0.2e1 * gamma2 * k * omega * sin(k * omega * t));

		aux[0][1]=exp(-t0 + t) * omega * omega * gamma1 * cos(omega * t0) / 0.2e1 + exp(-t0 + t) * k * k * omega * omega * gamma2 * cos(k * omega * t0) / 0.2e1 - exp(-t0 + t) * y0 / 0.2e1 + exp(-t0 + t) * gamma1 * cos(omega * t0) / 0.2e1 + exp(-t0 + t) * gamma2 * cos(k * omega * t0) / 0.2e1 + exp(-t0 + t) / 0.2e1 + exp(t0 - t) * omega * omega * gamma1 * cos(omega * t0) / 0.2e1 + exp(t0 - t) * k * k * omega * omega * gamma2 * cos(k * omega * t0) / 0.2e1 + exp(t0 - t) * y0 / 0.2e1 + exp(t0 - t) * gamma1 * cos(omega * t0) / 0.2e1 + exp(t0 - t) * gamma2 * cos(k * omega * t0) / 0.2e1 + exp(t0 - t) / 0.2e1 - (exp(-t0 + t) * y0 / 0.2e1 - exp(-t0 + t) * gamma1 * cos(omega * t0) / 0.2e1 + exp(-t0 + t) * omega * gamma1 * sin(omega * t0) / 0.2e1 - exp(-t0 + t) * gamma2 * cos(k * omega * t0) / 0.2e1 + exp(-t0 + t) * k * omega * gamma2 * sin(k * omega * t0) / 0.2e1 - exp(-t0 + t) / 0.2e1 - exp(t0 - t) * y0 / 0.2e1 - exp(t0 - t) * gamma1 * cos(omega * t0) / 0.2e1 - exp(t0 - t) * omega * gamma1 * sin(omega * t0) / 0.2e1 - exp(t0 - t) * gamma2 * cos(k * omega * t0) / 0.2e1 - exp(t0 - t) * k * omega * gamma2 * sin(k * omega * t0) / 0.2e1 - exp(t0 - t) / 0.2e1 - gamma1 * omega * omega * cos(omega * t) - gamma2 * k * k * omega * omega * cos(k * omega * t)) * (exp(-t0 + t) * omega * omega * gamma1 * cos(omega * t0) + exp(-t0 + t) * k * k * omega * omega * gamma2 * cos(k * omega * t0) - exp(-t0 + t) * y0 + exp(-t0 + t) * gamma1 * cos(omega * t0) + exp(-t0 + t) * gamma2 * cos(k * omega * t0) + exp(-t0 + t) - exp(t0 - t) * omega * omega * gamma1 * cos(omega * t0) - exp(t0 - t) * k * k * omega * omega * gamma2 * cos(k * omega * t0) - exp(t0 - t) * y0 - exp(t0 - t) * gamma1 * cos(omega * t0) - exp(t0 - t) * gamma2 * cos(k * omega * t0) - exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * gamma1 * cos(omega * t0) + exp(-t0 + t) * omega * gamma1 * sin(omega * t0) - exp(-t0 + t) * gamma2 * cos(k * omega * t0) + exp(-t0 + t) * k * omega * gamma2 * sin(k * omega * t0) - exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * gamma1 * cos(omega * t0) + exp(t0 - t) * omega * gamma1 * sin(omega * t0) + exp(t0 - t) * gamma2 * cos(k * omega * t0) + exp(t0 - t) * k * omega * gamma2 * sin(k * omega * t0) + exp(t0 - t) - 0.2e1 * gamma1 * omega * sin(omega * t) - 0.2e1 * gamma2 * k * omega * sin(k * omega * t));

		aux[1][0]=-(exp(-t0 + t) - exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * gamma1 * cos(omega * t0) + exp(-t0 + t) * omega * gamma1 * sin(omega * t0) - exp(-t0 + t) * gamma2 * cos(k * omega * t0) + exp(-t0 + t) * k * omega * gamma2 * sin(k * omega * t0) - exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * gamma1 * cos(omega * t0) + exp(t0 - t) * omega * gamma1 * sin(omega * t0) + exp(t0 - t) * gamma2 * cos(k * omega * t0) + exp(t0 - t) * k * omega * gamma2 * sin(k * omega * t0) + exp(t0 - t) - 0.2e1 * gamma1 * omega * sin(omega * t) - 0.2e1 * gamma2 * k * omega * sin(k * omega * t));
		
		aux[1][1]=-(exp(-t0 + t) * omega * omega * gamma1 * cos(omega * t0) + exp(-t0 + t) * k * k * omega * omega * gamma2 * cos(k * omega * t0) - exp(-t0 + t) * y0 + exp(-t0 + t) * gamma1 * cos(omega * t0) + exp(-t0 + t) * gamma2 * cos(k * omega * t0) + exp(-t0 + t) - exp(t0 - t) * omega * omega * gamma1 * cos(omega * t0) - exp(t0 - t) * k * k * omega * omega * gamma2 * cos(k * omega * t0) - exp(t0 - t) * y0 - exp(t0 - t) * gamma1 * cos(omega * t0) - exp(t0 - t) * gamma2 * cos(k * omega * t0) - exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * gamma1 * cos(omega * t0) + exp(-t0 + t) * omega * gamma1 * sin(omega * t0) - exp(-t0 + t) * gamma2 * cos(k * omega * t0) + exp(-t0 + t) * k * omega * gamma2 * sin(k * omega * t0) - exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * gamma1 * cos(omega * t0) + exp(t0 - t) * omega * gamma1 * sin(omega * t0) + exp(t0 - t) * gamma2 * cos(k * omega * t0) + exp(t0 - t) * k * omega * gamma2 * sin(k * omega * t0) + exp(t0 - t) - 0.2e1 * gamma1 * omega * sin(omega * t) - 0.2e1 * gamma2 * k * omega * sin(k * omega * t));

	}
	else{
		aux[0][0]=exp(-t0 + t) / 0.2e1 + exp(t0 - t) / 0.2e1 - (exp(-t0 + t) * y0 / 0.2e1 - exp(-t0 + t) * gamma1 * cos(omega * t0) / 0.2e1 + exp(-t0 + t) * omega * gamma1 * sin(omega * t0) / 0.2e1 - exp(-t0 + t) * gamma2 * cos(k * omega * t0) / 0.2e1 + exp(-t0 + t) * k * omega * gamma2 * sin(k * omega * t0) / 0.2e1 + exp(-t0 + t) / 0.2e1 - exp(t0 - t) * y0 / 0.2e1 - exp(t0 - t) * gamma1 * cos(omega * t0) / 0.2e1 - exp(t0 - t) * omega * gamma1 * sin(omega * t0) / 0.2e1 - exp(t0 - t) * gamma2 * cos(k * omega * t0) / 0.2e1 - exp(t0 - t) * k * omega * gamma2 * sin(k * omega * t0) / 0.2e1 + exp(t0 - t) / 0.2e1 - gamma1 * omega * omega * cos(omega * t) - gamma2 * k * k * omega * omega * cos(k * omega * t)) * (exp(-t0 + t) - exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * gamma1 * cos(omega * t0) + exp(-t0 + t) * omega * gamma1 * sin(omega * t0) - exp(-t0 + t) * gamma2 * cos(k * omega * t0) + exp(-t0 + t) * k * omega * gamma2 * sin(k * omega * t0) + exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * gamma1 * cos(omega * t0) + exp(t0 - t) * omega * gamma1 * sin(omega * t0) + exp(t0 - t) * gamma2 * cos(k * omega * t0) + exp(t0 - t) * k * omega * gamma2 * sin(k * omega * t0) - exp(t0 - t) - 0.2e1 * gamma1 * omega * sin(omega * t) - 0.2e1 * gamma2 * k * omega * sin(k * omega * t));

		aux[0][1]=exp(-t0 + t) * omega * omega * gamma1 * cos(omega * t0) / 0.2e1 + exp(-t0 + t) * k * k * omega * omega * gamma2 * cos(k * omega * t0) / 0.2e1 - exp(-t0 + t) * y0 / 0.2e1 + exp(-t0 + t) * gamma1 * cos(omega * t0) / 0.2e1 + exp(-t0 + t) * gamma2 * cos(k * omega * t0) / 0.2e1 - exp(-t0 + t) / 0.2e1 + exp(t0 - t) * omega * omega * gamma1 * cos(omega * t0) / 0.2e1 + exp(t0 - t) * k * k * omega * omega * gamma2 * cos(k * omega * t0) / 0.2e1 + exp(t0 - t) * y0 / 0.2e1 + exp(t0 - t) * gamma1 * cos(omega * t0) / 0.2e1 + exp(t0 - t) * gamma2 * cos(k * omega * t0) / 0.2e1 - exp(t0 - t) / 0.2e1 - (exp(-t0 + t) * y0 / 0.2e1 - exp(-t0 + t) * gamma1 * cos(omega * t0) / 0.2e1 + exp(-t0 + t) * omega * gamma1 * sin(omega * t0) / 0.2e1 - exp(-t0 + t) * gamma2 * cos(k * omega * t0) / 0.2e1 + exp(-t0 + t) * k * omega * gamma2 * sin(k * omega * t0) / 0.2e1 + exp(-t0 + t) / 0.2e1 - exp(t0 - t) * y0 / 0.2e1 - exp(t0 - t) * gamma1 * cos(omega * t0) / 0.2e1 - exp(t0 - t) * omega * gamma1 * sin(omega * t0) / 0.2e1 - exp(t0 - t) * gamma2 * cos(k * omega * t0) / 0.2e1 - exp(t0 - t) * k * omega * gamma2 * sin(k * omega * t0) / 0.2e1 + exp(t0 - t) / 0.2e1 - gamma1 * omega * omega * cos(omega * t) - gamma2 * k * k * omega * omega * cos(k * omega * t)) * (exp(-t0 + t) * omega * omega * gamma1 * cos(omega * t0) + exp(-t0 + t) * k * k * omega * omega * gamma2 * cos(k * omega * t0) - exp(-t0 + t) * y0 + exp(-t0 + t) * gamma1 * cos(omega * t0) + exp(-t0 + t) * gamma2 * cos(k * omega * t0) - exp(-t0 + t) - exp(t0 - t) * omega * omega * gamma1 * cos(omega * t0) - exp(t0 - t) * k * k * omega * omega * gamma2 * cos(k * omega * t0) - exp(t0 - t) * y0 - exp(t0 - t) * gamma1 * cos(omega * t0) - exp(t0 - t) * gamma2 * cos(k * omega * t0) + exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * gamma1 * cos(omega * t0) + exp(-t0 + t) * omega * gamma1 * sin(omega * t0) - exp(-t0 + t) * gamma2 * cos(k * omega * t0) + exp(-t0 + t) * k * omega * gamma2 * sin(k * omega * t0) + exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * gamma1 * cos(omega * t0) + exp(t0 - t) * omega * gamma1 * sin(omega * t0) + exp(t0 - t) * gamma2 * cos(k * omega * t0) + exp(t0 - t) * k * omega * gamma2 * sin(k * omega * t0) - exp(t0 - t) - 0.2e1 * gamma1 * omega * sin(omega * t) - 0.2e1 * gamma2 * k * omega * sin(k * omega * t));

		aux[1][0]=-(exp(-t0 + t) - exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * gamma1 * cos(omega * t0) + exp(-t0 + t) * omega * gamma1 * sin(omega * t0) - exp(-t0 + t) * gamma2 * cos(k * omega * t0) + exp(-t0 + t) * k * omega * gamma2 * sin(k * omega * t0) + exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * gamma1 * cos(omega * t0) + exp(t0 - t) * omega * gamma1 * sin(omega * t0) + exp(t0 - t) * gamma2 * cos(k * omega * t0) + exp(t0 - t) * k * omega * gamma2 * sin(k * omega * t0) - exp(t0 - t) - 0.2e1 * gamma1 * omega * sin(omega * t) - 0.2e1 * gamma2 * k * omega * sin(k * omega * t));

		aux[1][1]=-(exp(-t0 + t) * omega * omega * gamma1 * cos(omega * t0) + exp(-t0 + t) * k * k * omega * omega * gamma2 * cos(k * omega * t0) - exp(-t0 + t) * y0 + exp(-t0 + t) * gamma1 * cos(omega * t0) + exp(-t0 + t) * gamma2 * cos(k * omega * t0) - exp(-t0 + t) - exp(t0 - t) * omega * omega * gamma1 * cos(omega * t0) - exp(t0 - t) * k * k * omega * omega * gamma2 * cos(k * omega * t0) - exp(t0 - t) * y0 - exp(t0 - t) * gamma1 * cos(omega * t0) - exp(t0 - t) * gamma2 * cos(k * omega * t0) + exp(t0 - t)) / (exp(-t0 + t) * y0 - exp(-t0 + t) * gamma1 * cos(omega * t0) + exp(-t0 + t) * omega * gamma1 * sin(omega * t0) - exp(-t0 + t) * gamma2 * cos(k * omega * t0) + exp(-t0 + t) * k * omega * gamma2 * sin(k * omega * t0) + exp(-t0 + t) + exp(t0 - t) * y0 + exp(t0 - t) * gamma1 * cos(omega * t0) + exp(t0 - t) * omega * gamma1 * sin(omega * t0) + exp(t0 - t) * gamma2 * cos(k * omega * t0) + exp(t0 - t) * k * omega * gamma2 * sin(k * omega * t0) - exp(t0 - t) - 0.2e1 * gamma1 * omega * sin(omega * t) - 0.2e1 * gamma2 * k * omega * sin(k * omega * t));

	}
	mmultiply(aux,DR,aux);
	mmultiply(DF,aux,DF);
//	printf("%g %g\n%g %g\n\n",DF[0][0],DF[0][1],DF[1][0],DF[1][1]);
//	exit(0);
}

DF[0][0]=DF[0][0]-1;
DF[1][1]=DF[1][1]-1;

//printf("%g %g\n%g %g\n\n",DF[0][0],DF[0][1],DF[1][0],DF[1][1]);

det=DF[0][0]*DF[1][1]-DF[1][0]*DF[0][1];

//printf("%g\n",det);

invDF[0][0]=DF[1][1]/det;
invDF[0][1]=-DF[0][1]/det;
invDF[1][0]=-DF[1][0]/det;
invDF[1][1]=DF[0][0]/det;


free(DF[0]);
free(DF[1]);
free(DF);
free(aux[0]);
free(aux[1]);
free(aux);
free(DR[0]);
free(DR[1]);
free(DR);
}
*/

void mmultiply(double **A,double **B,double **C){
//A=B*C
double **aux;

aux=(double**)malloc(2*sizeof(double*));
aux[0]=(double*)malloc(2*sizeof(double));
aux[1]=(double*)malloc(2*sizeof(double));

aux[0][0]=B[0][0]*C[0][0]+B[0][1]*C[1][0];
aux[0][1]=B[0][0]*C[0][1]+B[0][1]*C[1][1];
aux[1][0]=B[1][0]*C[0][0]+B[1][1]*C[1][0];
aux[1][1]=B[1][0]*C[0][1]+B[1][1]*C[1][1];

A[0][0]=aux[0][0];
A[0][1]=aux[0][1];
A[1][0]=aux[1][0];
A[1][1]=aux[1][1];

free(aux[0]);
free(aux[1]);
free(aux);

}

