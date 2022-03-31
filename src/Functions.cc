#include "interface/Functions.h"



double funcAlpha(double* x, double* par)
{
  double xx = x[0];
  double val = 0.;
  
  for(int i = 0; i < 20; ++i)
  {
    double tau = pow(10,i-5);
    double alpha0 = ( alpha_p0 + alpha_p1*(i-5) + alpha_p2*pow(i-5,2) + alpha_p3*pow(i-5,3) ) / 20.;
    val += alpha0*exp(-xx/tau);
  }
  
  return val;
}



double timeScale(double T_a, double T_r)
{
  T_a = T_a + 273.15;
  T_r = T_r + 273.15;
  return exp( (-1./kB) * ((1.1692-0.00049*T_a*T_a/(T_a+655.))/T_a-(1.1692-0.00049*T_r*T_r/(T_r+655.))/T_r));
}



double myfunc_DCR(double* x, double* par)
{
  double xx = x[0];

  /*
  double x0 = par[4];
  
  double A = par[0];
  double C = par[1];
  double D = par[2];
  double E = par[3];
  double B = D*E*exp(E*x0) - 2.*C*x0;
  double F = A + B*x0 + C*x0*x0 - D*exp(E*x0);
  
  if( xx < x0 ) return A + B*xx + C*xx*xx;
  else          return D*exp(E*xx) + F;
  */

  double turnOff = 0.5 * ( 1. + TMath::Erf(-1.*(xx-par[0])/par[1]) );
  double turnOn  = 0.5 * ( 1. + TMath::Erf(+1.*(xx-par[0])/par[1]) );
  return exp(par[2]+par[3]*xx) + turnOff*exp(par[4]+par[5]*xx) + turnOn*exp(par[6]+par[7]*xx);
}



double myfunc_amp(double* x, double* par)
{
  double xx = x[0];
  double x0 = par[0];

  if( xx < x0 ) return (par[2]*xx + par[1]);
  else          return (par[3]*log(xx) + par[2]*x0 - par[3]*log(x0));
}
