# include <stdlib.h>
#include <stdio.h>
#include <math.h>


///This function computes the initial seed (\widehat{t}_0) to initial find_delta_max.sh.
///It reads from the file argv[1] n,m,by,ratio,omega and computes \wt.
///It writes everything again in the file argv[2].

int main (int argc, char *argv[])
{
double wt,bt0,delta,ratio,omega,by0;
int n,m;
FILE *fout, *fin;

fin=fopen(argv[1],"r");
fout=fopen(argv[2],"w");

fscanf(fin,"%i %i %lf %lf %lf %lf %lf",&n,&m,&by0,&bt0,&delta,&ratio,&omega);

wt=1.0/omega*acos(-(1+pow(omega,2))/2*pow(by0,2)*ratio);

fprintf(fout,"%i %i %lf %lf %lf %lf %lf\n",n,m,by0,wt,delta,ratio,omega);


fclose(fin);
fclose(fout);

}
