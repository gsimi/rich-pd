#include "pmtwavef.h"

//===============================================
// Fit function for signal wave
//===============================================
Double_t sigwavef(Double_t* x, Double_t* par) {
  /*
    par[0]    ->    start time
    par[1]    ->    rising time
    par[2]    ->    falling time
    par[3]    ->    pedestal
    par[4]    ->    normalization
    par[5]    ->    bin width
   */
  double t=*x;
  double t0 = par[0];
  double tau1 = par[1];
  double tau2 = par[2];
  double pedestal = par[3];
  double norm = par[4];
  double binw = par[5];
  return t<t0 ? pedestal :pedestal + norm*binw*(1+tau1/tau2)/tau2*
	   (1-exp(-(t-t0)/tau1))*exp(-(t-t0)/tau2);
} 

//===============================================
// Fit function for signal wave (multi-signals)
//===============================================
Double_t multi_sigwavef(Double_t* x, Double_t* par) {
  const int nFun = 4;
  double t=*x;
  double pedestal = par[3*nFun]; 
  double norm = par[3*nFun + 1]; 
  double binw = par[3*nFun + 2]; 

  double fun = 0.;
  for (int i = 0; i < nFun; i++) {
    double t0 = par[i];
    double tau1 = par[nFun + i];
    double tau2 = par[2*nFun + i];
    double tmp = t<t0 ? pedestal : pedestal +
           norm*binw*(1+tau1/tau2)/tau2*(1-exp(-(t-t0)/tau1))*exp(-(t-t0)/tau2);
    fun += tmp;
  }
  
  return fun;
}

//===============================================
// Fit function for trigger wave
//===============================================
// Stefano (just Fermi)
Double_t trgwavef(Double_t* x, Double_t* par) {
  /*
    par[0]    ->    start time
    par[1]    ->    rising time
    par[2]    ->    pedestal
    par[3]    ->    normalization
    par[4]    ->    bin width
   */
  double t=*x;
  double t0 = par[0];
  double tau = par[1];
  double pedestal = par[2];
  double norm = par[3];
  double binw = par[4];

  return pedestal + norm*binw/(1+exp(-(t-t0)/tau));
}
