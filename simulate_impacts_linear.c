# include <stdlib.h>
#include <stdio.h>
#include <math.h>

//This is the same program simulat_impacts_linear but using r<1

void eval_fun(double *x,double t,double y0,double t0, double piece,double eps,double omega);

//int main (){
int main (int argc, char *argv[]){

FILE *fimpacts,*fout,*fin;
double **impacts,*x,err,T,r;
double y0,t0,tf,eps,omega,h,piece;
double delta,ratio,te;
int i,j,n,m,npoints;
static double pi=3.14159265358979323846;

if (argc==1){
	fimpacts=fopen("impacts.dat","r");
	fin=fopen("system.dat","r");
	fout=fopen("orbit.dat","w");
}
else if (argc==4){
	fimpacts=fopen(argv[1],"r");
	fin=fopen(argv[2],"r");
	fout=fopen(argv[3],"w");
}
else{
	printf("Wrong number of input files.\n Please provide impacts, system and orbit\n");
	exit(0);
}


//fimpacts=fopen("impacts.dat","r");
//fout=fopen("orbit.dat","w");
//fin=fopen("system.dat","r");

impacts=(double**)malloc(50*sizeof(double*));
for (i=0;i<50;i++)
	{
	impacts[i]=(double*)malloc(2*sizeof(double));
	}
x=(double*)malloc(2*sizeof(double));


//fscanf(fin,"%i %i %lf %lf %lf %lf %lf",&n,&m,&y0,&t0,&eps,&r,&omega);

fscanf(fin,"%i %i %lf %lf %lf %lf %lf",&n,&m,&y0,&t0,&delta,&ratio,&omega);
//We choose \te=1
te=1;
eps=te*delta;
//r=pow(1-ratio*te*delta,1.0/(4.0*m));
r=1-ratio*te*delta;

for (i=0;i<=2*m;i++){
	fscanf(fimpacts,"%lf %lf",&impacts[i][0],&impacts[i][1]);
}

npoints=10000;
piece=1;
T=2*pi/omega;

for (i=0;i<2*m;i++){
	y0=impacts[i][0];
	t0=impacts[i][1];
	tf=impacts[i+1][1];
	h=(tf-t0)/(npoints+1);
	for (j=0;j<=npoints+1;j++){
		eval_fun(x,t0+j*h,y0,t0,piece,eps,omega);
		fprintf(fout,"%20.20f %20.20f %g %i\n",x[0],x[1],fmod(t0+h*j,T),i);
	}
	piece=piece*(-1);
}

printf("Errors at the jointing points:\n\n");
for (i=0;i<2*m;i++){
	eval_fun(x,impacts[i+1][1],impacts[i][0],impacts[i][1],pow(-1,i),eps,omega);
	err=r*x[1]-impacts[i+1][0];
	printf("Error in x:%g, error in y: %g\n",x[0],err);
}
printf("Error at the final point:\n");
printf("In y0: %g\n",r*x[1]-impacts[0][0]);
printf("In t0: %g\n",impacts[2*m][1]-impacts[0][1]-n*2*pi/omega);

free(x);
free(impacts);
fclose(fin);
fclose(fimpacts);
fclose(fout);
}

void eval_fun(double *x,double t,double y0,double t0, double piece,double eps,double omega){
double C1p,C1m,C2p,C2m;


//Perturbation H_1=xcos(omega*t)
///*
//Careful! The perturbation parameter in the solution that we are using here
//is the original one scaled by omega^2+1!!!
eps=eps/(pow(omega,2)+1);

C1p=(y0 - eps*cos(omega*t0) + eps*omega*sin(omega*t0) - 1)/2*exp(-t0);
C1m=(y0 - eps*cos(omega*t0) + eps*omega*sin(omega*t0) + 1)/2*exp(-t0);
C2p=(-y0 -eps*cos(omega*t0) - eps*omega*sin(omega*t0) - 1)/2*exp(t0);
C2m=(-y0 -eps*cos(omega*t0) - eps*omega*sin(omega*t0) + 1)/2*exp(t0);

if (piece>0){
	x[0]=C1p*exp(t)+C2p*exp(-t)+eps*cos(omega*t)+1;
	x[1]=C1p*exp(t)-C2p*exp(-t)-eps*omega*sin(omega*t);

	}
else {
	x[0]=C1m*exp(t)+C2m*exp(-t)+eps*cos(omega*t)-1;
	x[1]=C1m*exp(t)-C2m*exp(-t)-eps*omega*sin(omega*t);
	}
//*/

/*
//Perturbation H_1=x*(cos(omega*t)+cos(3*omega*t)):
double gamma1,gamma2,k;

k=3.0;
gamma1=eps/(1+pow(omega,2));
gamma2=eps/(1+pow(k*omega,2));
C1p=(y0-gamma1*cos(omega*t0)+omega*gamma1*sin(omega*t0)-gamma2*cos(k*omega*t0)+gamma2*k*omega*sin(k*omega*t0)-1)*exp(-t0)/2;
C1m=(y0-gamma1*cos(omega*t0)+omega*gamma1*sin(omega*t0)-gamma2*cos(k*omega*t0)+gamma2*k*omega*sin(k*omega*t0)+1)*exp(-t0)/2;
C2p=(-y0-gamma1*cos(omega*t0)-omega*gamma1*sin(omega*t0)-gamma2*cos(k*omega*t0)-gamma2*k*omega*sin(k*omega*t0)-1)*exp(t0)/2;
C2m=(-y0-gamma1*cos(omega*t0)-omega*gamma1*sin(omega*t0)-gamma2*cos(k*omega*t0)-gamma2*k*omega*sin(k*omega*t0)+1)*exp(t0)/2;
if (piece>0){
	x[0]=C1p*exp(t)+C2p*exp(-t)+1+gamma1*cos(omega*t)+gamma2*cos(k*omega*t);
	x[1]=C1p*exp(t)-C2p*exp(-t)-gamma1*omega*sin(omega*t)-gamma2*k*omega*sin(k*omega*t);
	}
else {
	x[0]=C1m*exp(t)+C2m*exp(-t)-1+gamma1*cos(omega*t)+gamma2*cos(k*omega*t);
	x[1]=C1m*exp(t)-C2m*exp(-t)-gamma1*omega*sin(omega*t)-gamma2*k*omega*sin(k*omega*t);
	}
*/

}

