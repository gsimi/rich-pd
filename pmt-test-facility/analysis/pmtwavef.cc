#include "pmtwavef.h"
Double_t pmtwavef(Double_t* x, Double_t* par) {
  /*
    par[0]    ->    start time
    par[1]    ->    rising time
    par[2]    ->    falling time
    par[3]    ->    pedestal
    par[4]    ->    normalization
    par[5]    ->    bin width
   */
  double t=*x;
  double t0=par[0];
  double taur=par[1];
  double tauf=par[2];
  double pedestal=par[3];
  double norm=par[4];
  double binw=par[5];
  return t<t0 ? pedestal :pedestal + norm*binw*tauf/(1+taur/tauf)*
	   (1-exp(-(t-t0)/taur))*exp(-(t-t0)/tauf);
} 

