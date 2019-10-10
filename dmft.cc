/*
 * dmft.cc
 * 
 * Author: oabdurazakov
 *
 */

#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <armadillo>

	using namespace std;
	using namespace arma;
	typedef complex<double> cdouble;

double t = 0.25; //Hopping amplitude
double T = 0.01; //Temperature
double U = 4.8; //e-e interaction
int N = 500;//Half grid size
cdouble I(0,1);//Complex i

double band(double, double);//Band structure
double bessj0(double);//Bessel function of first kind and order zero
double f(double, double, double);//Integrant of the DOS
double trapez(double, double, double);// Trapezoid method to integrate
double dos(double);// 2D DOS 
double scdos(double);//Semicircular DOS
double fermi(double, double);//Fermi distribution function
void prepare_grid(vector<double> &);//Prepares the grid
void prepare_ferm(vector<double>, vector<double> &);//Prepares the grid
void prepare_dens(vector<double>, vector<double> &);//Prepares the grid
void prepare_g0(vector<double>, vector<double>, vector<cdouble> &);//Prepares the impurity GF 
void calc_apn(vector<double> ,vector<cdouble>, vector<double>, vector<double> &, vector<double> &);//Calculates A+ and A-
void calc_ps(vector<double>, vector<double>, vector<double>, vector<double> &, vector<double> &, vector< vector<int> >);//Calculates P1 and P2
void calc_isigma(vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double> &, vector< vector<int> >);//Caculalates imaginary SE
void calc_rsigma(vector<double>, vector<double>, vector<double> &);// Calculates real SE
void create_index(vector< vector<int> > &);//Creates difference index
cdouble htransform(cdouble, vector<double>, vector<double>);//Does Hilbert transform

main(){
	vector<double> w(2*N+1,0.0);//Declare the grid
	vector<double> ferm(2*N+1,0.0);//Prepare the Fermi Func on the grid
	vector<double> dens(2*N+1,0.0);//Declare the dynamic DOS
	vector< vector<int> > index(2*N+1, vector<int>(2*N+1));//Difference index
	vector<double> ap(2*N+1,0.0);//A+ 
	vector<double> an(2*N+1,0.0);//A- 
	vector<double> p1(2*N+1,0.0);//P1 
	vector<double> p2(2*N+1,0.0);//P2 
	vector<double> isigma(2*N+1,0.0);//Declare the imaginary Sigma 
	vector<double> rsigma(2*N+1,0.0);//Declare the real Sigma
	vector<cdouble> sigma(2*N+1,0.0);//Declare the Sigma
	vector<cdouble> sigma_new(2*N+1,0.0);//Declare the new Sigma
	vector<cdouble> g0(2*N+1,0.0);//Declare the bare impurity GF 
	vector<cdouble> g(2*N+1,0.0);//Declare the impurity GF 
	prepare_grid(w);//Prepare the grid
	prepare_ferm(w,ferm);//Prepare fermi func on the grid
	prepare_dens(w,dens);//Prepare the DOS on the grid
	prepare_g0(w,dens,g0);//Prepare the impurity GF on the grid
	create_index(index);//Create difference index
	double change;// For convergence
	double big_change = 1.0e5;//For convergence
	double small_change = 1.0e-5;//For convergence
	double mix = 0.20;//For mixing the new solution with the old
	int itermax = 500;//Maximum number of iterations
	
			
	
	//Do the self-consistency
	for(int iter = 0; iter < itermax; iter++)
	{
		//Calculate A+ and A-
		calc_apn(w, g0, ferm, ap, an);
		//for(int i = 0; i < 2*N+1; i++){cout << w[i] << "  " << ap[i] << " " << an[i]<< endl;}
		//Calculate P1 and P2
		calc_ps(w, ap, an, p1, p2, index);	
		//for(int i = 0; i < 2*N+1; i++){cout << w[i] << "  " << p1[i] << " " << p2[i]<< endl;}
		//Calculate the imaginary part of the self-energy
		calc_isigma(w, p1, p2, ap, an, isigma, index);
		//Calculate the real part of the  self-energy 
		calc_rsigma(w, isigma, rsigma);
		//for(int i = 0; i < 2*N+1; i++){cout << w[i] << "  " << rsigma[i] << " " << isigma[i]<< endl;}
		//Calc the complex self-energy
		for(int i = 0; i < 2*N+1; i++)
		{
			sigma_new[i] = rsigma[i] + I*isigma[i];
		}
		//for(int i = 0; i < 2*N+1; i++){cout << w[i] << "  " << sigma_new[i].real() << " " << sigma_new[i].imag() << endl;}
		//Check the convergence
		double temp = 0.0;
		vec a(2*N+1);
		vec x(2*N+1);
		for(int i = 0; i < 2*N+1; i++)
		{
			double dw = w[1] - w[0];
			temp += abs(sigma_new[i] - sigma[i])*dw;
			a(i) = abs(sigma_new[i] - sigma[i]);
			x(i) = w[i];
			//cout << i << " "<< sigma_new[i].real()<< " " << sigma_new[i].imag() << " " << sigma[i].real() << " " << sigma[i].imag()<<endl;
		}
		change = temp;
		//mat haha = trapz(x, a);
		//change = as_scalar(haha);
		//cout << change << endl;
		//If convergence reached it is done!
		if(change < small_change) break;
		if(change > big_change) mix /= 1.2;
		big_change = change;
		//Slowly mix the new solution with the old one
		for(int i = 0; i < 2*N+1; i++)
		{
			sigma[i] = sigma[i]*(1.0-mix) + sigma_new[i]*mix;
		}
		//Display the convergence	
		cout << iter << " " << change << " " << mix << endl;
		for(int i = 0; i < 2*N+1; i++){
			g[i] = htransform(w[i] - sigma[i], w, dens);
		}
		//Calculate the new impurity GF
		for(int i = 0; i < 2*N+1; i++){
			g0[i] = 1.0/(1.0/g[i] + sigma[i]);
			//if(iter==50)cout << w[i] << " "<<  -g0[i].imag() <<endl;
		}
		//Enforce the causality if not hold
		for(int i = 0; i < 2*N+1; i++)
		{
			if(g0[i].imag()>0) g0[i] = g0[i].real() - 1.0e-10*I;
		}
		//Pring out the DOS at every 10s iteration step
		if(iter%10==0)
		{
			char filename[256];
			sprintf(filename, "dos/%03d.dat", iter/10);
			ofstream outfile(filename);
			for(int i = 0; i < 2*N+1; i++)
		{
				outfile << w[i] << " " << -g[i].imag() << endl;
			}
		}	
	}
	//Print out the initial DOS	
	ofstream outfile1("dos1.dat");
	for(double x = -4; x <= 4; x +=0.008)
	{
		outfile1 << x << " " << dos(x) << endl;
	}
}

double band(double kx, double ky)
{
	return -2*t*(cos(kx) + cos(ky));
}

double dos(double w)
{
	return (abs(w) > 4.0*t)?0:trapez(w, 0.0, 1000.0)/M_PI;
}

double fermi(double w, double T)
{
//	if(w>4*t) return 0.0;
//	if(w<-4*t) return 1.0;
	return 1.0/(exp(w/T) + 1.0);

}

double trapez(double w, double xmin, double xmax)
{
        int N = 10000;
        double x = xmin;
        double dx = (xmax - xmin)/N;
        double result = 0;
        for(int i = 0; i < N; i++)
	{
                result += (i!=(0|N-1))?f(x,t,w):f(x,t,w)/2.;
                x += dx;
        }
        return result*dx;
}
       
double f(double x, double t, double w)
{
	return cos(w*x)*pow(bessj0(2*t*x),2);
}

double bessj0( double x )
{
   double ax,z;
   double xx,y,ans,ans1,ans2;

   if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
         +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
      ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
         +y*(59272.64853+y*(267.8532712+y*1.0))));
      ans=ans1/ans2;
   } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-0.785398164;
      ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
         +y*(-0.2073370639e-5+y*0.2093887211e-6)));
      ans2 = -0.1562499995e-1+y*(0.1430488765e-3
         +y*(-0.6911147651e-5+y*(0.7621095161e-6
         -y*0.934935152e-7)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
   }
   return ans;
}

double scdos(double w)
{
	return (abs(w)>4*t)?0:sqrt(16*t*t - w*w)/(8*M_PI*t*t);
}

cdouble htransform(cdouble z, vector<double> w, vector<double> dens)
{
	double d0,a,b;
	d0 = (abs(z.imag() > 0.1))?0.0:dos(z.real());
	cx_vec temp(2*N+1);
	cx_vec x(2*N+1);
	for(int i = 0; i < 2*N+1; i++)
	{
		temp(i) = (dens[i] - d0)/(z - w[i]);
		x(i)= w[i];
	}
	cx_mat r = trapz(x, temp) + d0*(log(z - w[0])-log(z - w[2*N]));
	a = as_scalar(real(r));
	b = as_scalar(imag(r));
	cdouble c = a + b*I; 
	return c;
}

void prepare_grid(vector<double> &w)
{
	double wmax = 4.0;
	double dw = wmax/(N);
	for(int i = 0; i < 2*N+1; i++)
	{
		w[i] = -wmax + i*dw;
	}

}

void prepare_ferm(vector<double> w, vector<double>& ferm)
{
	for(int i = 0; i < 2*N+1; i++)
	{
		ferm[i] = fermi(w[i], T);
	}
}
void prepare_dens(vector<double> w, vector<double>& dens)
{
	for(int i = 0; i < 2*N+1; i++)
	{
		dens[i] = dos(w[i]);
	}
}

void prepare_g0(vector<double> w, vector<double> dens, vector<cdouble>& g0)
{
	for(int i = 0; i < 2*N+1; i++)
	{
		g0[i] = htransform(w[i] + 0.001*I, w, dens);
	}
}

void calc_apn(vector<double> w, vector<cdouble> g0, vector<double> ferm, vector<double> &ap, vector<double> &an)
{
	for(int i = 0; i < 2*N+1; i++)
	{
		ap[i] = -g0[i].imag()*ferm[i]/M_PI;
		an[i] = -g0[i].imag()*(1.0 - ferm[i])/M_PI;
	}
}

void calc_ps(vector<double> w, vector<double>ap, vector<double>an, vector<double> &p1, vector<double> &p2, vector< vector<int> > index)
{
	double dw = w[1] - w[0];
	//int l;
	for(int i = 0; i < 2*N+1; i++)
	{
		double temp1 = 0.0, temp2 = 0.0;
		for(int j = 0; j < 2*N +1; j++)
		{
			int l = index[j][i];
			if(l>=0)
			{
				temp1 += an[j]*ap[l];
				temp2 += ap[j]*an[l];
			}
		}
		p1[i] = temp1*dw;	
		p2[i] = temp2*dw;
	//	cout << p1[i]<< " " << p2[i]<< endl;	
	}
}

void calc_isigma(vector<double> w, vector<double> p1, vector<double>p2, vector<double>ap, vector<double>an, vector<double> &isigma, vector< vector<int> > index)
{
	double dw = w[1]-w[0];
	//int l;
	for(int i = 0; i < 2*N+1; i++)
	{
		double temp1 = 0.0, temp2 = 0.0;
		for(int j = 0; j < 2*N +1; j++)
		{
			int l = index[j][i];
			if(l>=0)
			{
				temp1 += ap[j]*p1[l];
				temp2 += an[j]*p2[l];
			}
		}
		isigma[i] = -U*U*(temp1 + temp2)*dw*M_PI;
	//	cout << w[i]<< " "<< isigma[i]<< endl;
	}
}

void calc_rsigma(vector<double> w, vector<double>isigma, vector<double> &rsigma)
{
	double dw = w[1] - w[0];
	vector<double> extra(2*N+1,0);
	vector<double> del(2*N+1,0);
	for(int i = 1; i < 2*N; i++)
		extra[i] = log(w[2*N] - w[i])/(w[i] - w[0]);
	
	del[0] = (isigma[1] - isigma[0])/dw;
	del[2*N] = (isigma[2*N] - isigma[2*N-1])/dw;
	for(int i = 1; i < 2*N; i++)
		 del[i] = (isigma[i+1] - isigma[i-1])/(2*dw);
	for(int i = 0; i < 2*N+1; i++)
	{
		double temp = 0.0;
		for(int j = 0; j < 2*N+1; j++)
			temp += (i!=j)?(isigma[j] - isigma[i])/(w[j] - w[i])*dw:del[i]*dw;
		
		rsigma[i] = (temp + isigma[i]*extra[i])/M_PI;
		//cout << w[i] << " " << rsigma[i] << endl;
	}
}

void create_index(vector< vector<int> > &index)
{
	int k;	
	for(int i = 0; i < 2*N+1; i++)
	{
		for(int j = 0; j < 2*N+1; j++)
		{
				k = j - i + N;
				if((k < 0)||(k >= 2*N + 1)) k = 2*N;
				index[j][i] = k;
		}
	}
}
