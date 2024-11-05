#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>

// Initial conditions
#define delt	25*4
#define t0 	-3*delt
#define tM 	500

//Paramenters and constants
#define chi0 	0.1
#define Er 	4.2
#define hbar	658.5
#define Emax	300
#define Del0	30
#define dt	2
#define pi	3.14159265359
#define T2	100

#define N 	 51
#define delE	 Emax/N
#define step	 (int) (tM - t0)/dt + dt

//double complex dist[step] for generating the system of ode;
double complex fp[step][2];
double complex fpn[step][2];

//Prefactor of parameters
double preEner = sqrt(Er)/pi*delE;
double preOme = sqrt(pi)/(2.*delt)*chi0;

//Subroutine write data to file
void savetofile(double t[], double a[][2])
{
	
	// contacting with file
	FILE *file;
	// open file data with name data.txt
	file = fopen("data.txt","w+");
	for(int i=1;i<step;i++)
	{	

		//printing data to file in order of t - distribution func - polarization func
		fprintf(file,"%.15f \t %.15f \t %.15f \n",t[i],a[i][0],a[i][1]);
	}
	fclose(file);
}


//Defining the g(n,l) subroutine
double g(int n,int l)
{	
	double preg = 1./sqrt(n*delE);
	double fac  = fabs( ( sqrt(n) + sqrt(l) )/( sqrt(n) - sqrt(l) ) ) ;
	double res  = preg*log(fac);
	return res;

}

//Defining the renormalize energy E(n) subroutine
double complex Ener(double complex ge)
{
	double complex res = preEner*ge;
	return res;


}
//Defining the Rabi frequency Ome(n,t) subroutine
double complex Ome( double t, double complex gome)
{
	double complex pref1 = preOme*exp(-t*t/(delt*delt));
	double complex pref2 = preEner/hbar*gome;
	double complex res = pref1 + pref2;
	return res;
	
}

//Defining system of ode of sbe
void func(double time, double complex fpn[][2],double complex fp[][2])
{
	double complex gf, ge;
	for(int n=1;n<N;n++)
	{
		gf = 0 + I*0; ge = 0 + I*0;
		for (int j=1;j<N;j++)
		{
			if (fabs(sqrt(n) - sqrt(j)) >1e-3)
			{
			gf = gf +  g(n,j)*2.*fpn[j][0];
			ge = ge +  g(n,j)*fpn[j][1];
		  }
		}

		
		fp[n][0] = -2.*cimag(Ome(time,ge)*conj(fpn[n][1]));
		
		fp[n][1] = (-I/hbar)*(n*delE - Del0 - Ener(gf))*fpn[n][1] + I*(1. - 2.*fpn[n][0])*Ome(time,ge) - fpn[n][1]/T2;
	}
}
//Defining runga kutta 4th order;
void rk4(double time,double complex fpn[][2])
{
	double complex k1[step][2], k2[step][2] , k3[step][2], k4[step][2] ; 
	double complex term[step][2];
	double complex res[step][2], fres[step][2];
	double t[step];
	double pref  =  delE*sqrt(delE);
	double result[step][2];
	double finalres[2] = {0.,0.};
	t[1] = time;
	for(int i=1;i<step;i++)
	{
		t[i+1] = t[i] + dt;
		fres[i][0] = 0 + 0*I;
		fres[i][1] = 0 + 0*I;
	}

	for(int i=1 ;i<step;i++)
	{
		func(t[i],fpn,res);

		for(int j=1;j<N;j++)
		{
			k1[j][0] = dt*res[j][0];
			k1[j][1] = dt*res[j][1];
			res[j][0] = 0. + 0.*I;
			res[j][1] = 0. + 0.*I;
		}
		for(int j=1;j<N;j++)
		{
			term[j][0] = fpn[j][0] + k1[j][0]/2.;
			term[j][1] = fpn[j][1] + k1[j][1]/2.;
		}
		
		func(t[i] + dt/2,term,res);
		
		for(int j=1;j<N;j++)
		{
			k2[j][0] = dt*res[j][0];
			k2[j][1] = dt*res[j][1];
			res[j][0] = 0. + 0.*I;
			res[j][1] = 0. + 0.*I;
		}
		for(int j=1;j<N;j++)
		{
			term[j][0] = fpn[j][0] + k2[j][0]/2.;
			term[j][1] = fpn[j][1] + k2[j][1]/2.;
		}
		
		func(t[i] + dt/2,term,res);
		
		for(int j=1;j<N;j++)
		{
			k3[j][0] = dt*res[j][0];
			k3[j][1] = dt*res[j][1];
			res[j][0] = 0. + 0.*I;
			res[j][1] = 0. + 0.*I;
		}
		for(int j=1;j<N;j++)
		{
			term[j][0] = fpn[j][0] + k3[j][0];
			term[j][1] = fpn[j][1] + k3[j][1];
		}
		
		func(t[i] + dt,term,res);
		
		for(int j=1;j<N;j++)
		{
			k4[j][0] = dt*res[j][0];
			k4[j][1] = dt*res[j][1];
		}

		for(int j=1;j<N;j++)
		{
			fpn[j][0] = fpn[j][0] + 1./6*(k1[j][0] + 2.*k2[j][0] + 2.*k3[j][0] + k4[j][0]);
			fpn[j][1] = fpn[j][1] + 1./6*(k1[j][1] + 2.*k2[j][1] + 2.*k3[j][1] + k4[j][1]);	
		}

		for (int j=1;j<N;j++)
		{
			finalres[0] += sqrt(j)*fpn[j][0];
			finalres[1] += sqrt(j)*fpn[j][1];
		}
			result[i][0] = pref*creal(finalres[0]);
			result[i][1] = cabs(pref*finalres[1]);
	}
		savetofile(t,result);
	
}
int main()
{
	double complex fp[step][2],res[step][2] ;
	#Creating initial values for fp array
	for(int i=1;i<step;i++)
	{
		fp[i][0] = 0. + 0.*I;
		fp[i][1] = 0. + 0.*I;
	}
	#Print time subroutine
	clock_t begin = clock();
  	rk4(t0,fp);
  	clock_t end   = clock();
  	double time = (double) (end - begin )/CLOCKS_PER_SEC;
  	printf("process time %.2f \n",time);
	return 0;
}
