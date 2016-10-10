

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

/*
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
*/


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
