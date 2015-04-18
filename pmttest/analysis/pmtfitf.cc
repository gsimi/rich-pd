#include "TMath.h"
#include "TF1.h"
#include <iostream>

#include "pmtfitf.h"

/* 
   Sum of Gaus functions with poisson probability 
*/

Double_t pmtpdf_gaus(Double_t* x, Double_t* par) {
  /*===========================================================
         par[0]  ->   Poisson average N
         par[1]  ->   1 p.e. ADC separation (gain)
         par[2]  ->   gain RMS contribution
         par[3]  ->   N=0 ADC offset
         par[4]  ->   sigma noise fabs(vecxp.at(1)-vecxp.at(0))
         par[5]  ->   Absolute count normalisation
	 par[6]  ->   Gain of first dynode in number of photoelectrons
	 par[7]  ->   Poisson average 1st dynode
	 par[8]  ->   Probability of interaction on first dynode
	 par[9]  ->   Eff
	 par[10] ->   bin width, needed for normalization

  =============================================================*/
  const int maxN=6;
  const int minN=0; //min and max number of p.e.

  Double_t ci = 0,ci2=0;
  const Double_t gausfac  = 0.3989422804;

  Double_t npe = par[0];
  Double_t gain = par[1];
  Double_t grms = par[2];
  Double_t ped = par[3];
  Double_t noise = par[4];
  Double_t norm = par[5];//total normalisation coef
  //  Double_t pct  = par[6]; //unused
  Double_t g1 = par[6]; 
  Double_t npe1 = par[7];

  //1.3 is the theoretical increase of the relative fluctuation
  //due to the 2.3:1.2:1:...:1:0.5 voltage divider setup
  Double_t grms1=grms/g1*1.3;

  Double_t frac = par[8];// proba. 1st dyn.
  Double_t eff  = par[9];
  Double_t bw=par[10];

  Double_t probarr[maxN];;
  Double_t probarr1[maxN];;
  bool debug=false;

  if (debug)  printf ("%f %f %f %f %f %f\n",npe,gain,grms,ped,noise,norm);
  
  for (int i=0; i<maxN; i++) {
    if (i==0) {
      ci = exp(-npe);
      ci2 = exp(-npe1);
    } else {
      ci = ci*npe/((float)i);
      ci2 = ci2*npe1/((float)i);
    }
    
    if (debug)     printf("    ci = %f ",ci);
    if (i>=minN) {
      /* 
   	 probarr[i] is the recursive probability  of generating 
   	 i photoelectrons, taking into account the primary p.e. and
   	 those generated by cross talk
  
       */
      probarr[i]=ci;
      probarr1[i]=ci2;
    }//if (i>=nmin)
    else {probarr[i]=0;probarr1[i]=0;}
  }//for (i=0; i<maxN)
  
  /*
  double p=0.25; double N=npe/p; 
  for (int k=0; k<maxN; k++) {
    if (k<=N) 
      probarr[k]=TMath::Binomial(N,k)*pow(p,k)*pow(1-p,N-k);
    else
      probarr[k]=0; 
  }
  */  
  
  Double_t valK=0;
  Double_t val1=0;

  for (int i=minN; i<maxN; i++) {
    Double_t sigma = sqrt(noise*noise+float(i)*grms*grms);
    Double_t sigma2 = sqrt(noise*noise+float(i)*grms1*grms1);
    
    Double_t xm = ped + float(i)*gain;
    Double_t xm2 = ped + float(i)*gain/g1;

    //    if (debug)     printf (" sigma = %f, noise = %f, xm = %f, grms = %f\n",sigma,noise,xm,grms);
    
	/* 
	   val is the distribution of the signal amplitudes
	   obtained as the sum of gaussian functions with mean proportional to 
	   the number of p.e. and probability computed above
	*/
      valK += probarr[i] * exp(-0.5*((*x-xm)*(*x-xm))/sigma/sigma)*gausfac/sigma ;
      val1 += probarr1[i] * exp(-0.5*((*x-xm2)*(*x-xm2))/sigma2/sigma2)*gausfac/sigma2 ;	 
  }
  return   norm*eff*((1-frac)*valK+frac*val1)*bw;
}


Double_t pmtpdf_gaus2(Double_t* x, Double_t* par) {
  /*===========================================================
         par[0]  ->   Poisson average N
         par[1]  ->   1 p.e. ADC separation
         par[2]  ->   gain RMS contribution
         par[3]  ->   N=0 ADC offset
         par[4]  ->   sigma noise fabs(vecxp.at(1)-vecxp.at(0))
         par[5]  ->   Absolute count normalisation
	 par[6]  ->   Cross talk probability (ref. http://www.gap-optique.unige.ch/Publications/PDF/SiPM-OPticsExpress.pdf)
	 par[7]  ->   Poisson average 1st dynode
	 par[8]  ->   gain RMS 1st dynode
	 par[9]  ->   1 p.e 1st dynode ADC separation
	 par[10] ->   Eff
	 par[11] ->   bin width, needed for normalization

  =============================================================*/
  int maxN=10, minN=0; //min and max number of p.e.

  Double_t ci = 0;//,ci2=0;
  const Double_t gausfac  = 0.3989422804;

  Double_t npe = par[0];
  Double_t gain = par[1];
  Double_t grms = par[2];
  Double_t ped = par[3];
  Double_t noise = par[4];
  Double_t ln_norm = par[5];
  Double_t pct  = par[6];
  Double_t P = par[7];
  Double_t f = par[8];
  Double_t eff  = par[9];
  Double_t bw=par[10];

  //  Double_t standdev = par[11];
  //  Double_t frac = par[10];
  Double_t probarr[maxN];;
  // Double_t probarr2[maxN];;
  bool debug=false;

  if (debug)  printf ("%f %f %f %f %f %f\n",npe,gain,grms,ped,noise,ln_norm);

  for (int i=0; i<maxN; i++) {
    if (i==0) {
	    ci = exp(-npe);
	    // ci2 = exp(-npe1);
    } else {
      ci = ci*npe/((float)i);
      // ci2 = ci2*npe2/((float)i);
    }

    if (debug)     printf("    ci = %f ",ci);
    if (i>=minN) {
      /* 
	 probarr[i] is the recursive probability  of generating 
	 i photoelectrons, taking into account the primary p.e. and
	 those generated by cross talk
	 
      */
     double pi=ci*pow(1-pct,i);
     // double pi2=ci2*pow(1-pct,i);
     for (int k=1;2*k<i+1;k++){
	pi+=probarr[i-k]*pow(1-pct,i-2*k)*pow(pct,k)*TMath::Binomial(i-k,k);
	//pi2+=probarr2[i-k]*pow(1-pct,i-2*k)*pow(pct,k)*TMath::Binomial(i-k,k);
      }
      probarr[i]=pi;
      //probarr2[i]=pi2;
    }//if (i>=nmin)
    else {probarr[i]=0;}
  }//for (i=0; i<maxN)
  
  Double_t val=0;
  
  for (int i=minN; i<maxN; i++) {
    Double_t sigma = sqrt(noise*noise+float(i)*grms*grms);
     Double_t sigma2 = sqrt(noise*noise+f*grms*grms);
    
    Double_t xm = ped + float(i)*gain;
    // Double_t xm2 = ped + float(i)*gain2;

    if (debug)     printf (" sigma = %f, noise = %f, xm = %f, grms = %f\n",sigma,noise,xm,grms);
     if (i==0){
       val += eff * probarr[i] *((1-P)* exp(-0.5*((*x-xm)*(*x-xm))/sigma/sigma)*gausfac/sigma + P*(gausfac/(sigma2))*exp(-0.5*((*x-gain*f)*(*x-gain*f))/(sigma2*sigma2)));
       } else{
	/* 
	   val is the distribution of the signal amplitudes
	   obtained as the sum of gaussian functions with mean proportional to 
	   the number of p.e. and probability computed above
	*/
    val += eff * probarr[i] * exp(-0.5*((*x-xm)*(*x-xm))/sigma/sigma)*gausfac/sigma ;
				 } 
  }
  
  return val*exp(ln_norm)*bw;
}



/* 
   Probability density distribution for the 
   response of a pmt with equal gain dynodes.
   See H.H.Tan eq. 66 for response to a poisson distributed
   number of photoelectrons. To simplify the calculations we 
   express the constants in terms of the average gain G
    A=G
    B=1/2 *(G-1)/pow(G-1,1/ndynodes)
   using the expression for the gain of a single dynode g
   and the response normalized to the gain t

    g=G/(G-1)*pow(G-1,1./ndynodes) 

    t=y/G
   we obtain the constants in eq 66 in terms of g and t

    A/B = 2g
    sqrt(A/y)/B = 2g/(G*sqrt(t)) 
    y/B = t*2g
    2 sqrt(Ay)/B = sqrt(t)*4g

   so that eq 66 has the following expression

    1/G*{exp(-mu*(1-exp(-2g)))*delta(t) + 2g/(sqrt(t)) * Sum_n [ Poisson(n,mu) * exp (-2g(n+t)) * I1(4*g*sqrt(n*t)) ] }

   with mu=alpha_i 

   The term with delta represents the output when no photoelectrons are generated with the 
   hypothesys that there is no noise. In real cases there is noise due to thermal and shot noise, therefore
   the delta function should be replaced by a gaussian representing the noise distribution.
   The rms of the noise gaussian in the variable t is noise/G

    1/G*{exp(-mu*(1-exp(-2g)))*Gaus(t,0,noise/G,kTRUE) + 2g/(sqrt(t)) * Sum_n [ Poisson(n,mu) * exp (-2g(n+t)) * I1(4*g*sqrt(n*t)) ] }

*/

Double_t pmtpdf_bessel(Double_t* x, Double_t* par){
  // *x is the signal in e

  double mean_pe=par[0]; //poisson average of number of pe emitted by the cathode
  double G=par[1]; // in electrons
  double ndynodes=par[2];
  double ped=par[3]; //pedestal  = ADC offset
  double noise=par[4];//sigma noise 
  double norm=par[5];//absolute count normalization
  if (mean_pe<1e-12) return 0;

  /* g is the gain of one dynode g=G/(G-1)*pow(G-1,1./N) 
     In the limit G>>1 then g=pow(G,1./N) */
  double g=G/(G-1)*pow(G-1,1./double(ndynodes));  //gain of one dynode
  double t=(*x-ped)/G;  //raio between signal height and average total gain

  double sigpdf=0;
  int nmax=mean_pe+3*sqrt(mean_pe)+5; //limit for n to reduce computational time

  if (t>1e-12){
    for (int n=1; n<nmax; n++){//for n=0 the expression is zero
      double tmpval=TMath::Poisson(n,mean_pe);//alpha^n exp(-alpha)/n!
      double bessel_arg=4*g*sqrt(n*t); // = 2 sqrt(nAy)/B
      if (bessel_arg<709)
	tmpval*=exp(-2*g*(n+t))*sqrt(n)*TMath::BesselI1(bessel_arg);//exp(-y/B-nA/B)sqrt(n)*I1(2 sqrt(nAy)/B)
      else{
	//approximate expression of bessel function for large values of the argument
	tmpval*=exp(-2*g*(n+t-2*sqrt(n*t)))*sqrt(n/(2.*TMath::Pi()*bessel_arg));
      }
      sigpdf+=tmpval;
    }

    sigpdf *= 2*g/sqrt(t);
  }
  
  // now add contribution of noise corresponding to no photoelectrons  produced by the photocathode
  sigpdf/=G;
  double p0=exp(-mean_pe*( 1- exp(-2*g) )); // = exp(alpha_i*(exp(-A/B)-1))
  double noisepdf = p0*TMath::Gaus(*x,ped,noise,kTRUE); 
  double val = norm*(noisepdf+sigpdf); //this gives a function nomalized to norm

  return val;
}

/* 
   This is the pmt response pdf in the case where the gain of the first dynode is different from 
   the others. This is modeled as a PMT with N-1 dynodes, average gain G/g1 
   and an poisson distributed number of photoelectrons with average g1.
 */

Double_t pmtpdf_besselg1(Double_t* x, Double_t* par){
  //*x is the signal in Me
  double mean_pe=par[0]; //mean number of pe at the photocathode
  double G=par[1]; //overall gain of the pmt in Me
  double ndynodes=par[2]; //number of dynodes
  double ped=par[3]; //pedestal  N= ADC offset
  double noise=par[4];//sigma noise 
  double norm=par[5];//absolute count normalization
  double g1=par[6]; //gain of the fisrt dynode in number of photoelectrons
  double mean_pe1=par[7]; //poisson average of photons converting on dynode 1
  double frac=par[8]; //probability to convert on dynode 1  
  Double_t eff  = par[9];
  Double_t bw=par[10];  //bin width, needed for normalization

  double val=0;

  //Assume photoelectrons are poisson distributed with average mean_pe
  //and add all contributions

  // 1) add noise contribution
  val +=norm*TMath::Poisson(0,mean_pe)*TMath::Gaus(*x,ped,noise,kTRUE);

  // 2) add photoelectrons contribution assuming an
  //N-1 dynode amplification with g1 average photoelectrons
  //and overall gain G/g1
  Double_t x1[1], par1[6];
  par1[1]=G/g1; //gain of the remaining n-1 dynodes
  par1[2]=ndynodes-1;
  par1[3]=ped;
  par1[4]=noise;
  par1[5]=norm;
  
  //now consider that the photoelectrons emitted by the photocathode are
  //poisson distributed
  int nmax=mean_pe+3*sqrt(mean_pe)+5;
  for (int n1=1;n1<nmax;n1++){
    par1[0]=g1; //Poisson average of # of photoelectrons in the N-1 stages PMT
    *x1=(*x)/n1;//Response to a single photoelectron. x=x1*n1 is the response to n1 photoelectrons
    val+=TMath::Poisson(n1,mean_pe)*pmtpdf_bessel(x1,par1)/n1;

    /* this seems to give better fits, but I think it is wrong!! 

    par1[0]=g1*n1; //average #of electrons after the first dynode multiplication  
    *x1=*x;
    val+=TMath::Poisson(n1,mean_pe)*pmtpdf_bessel(x1,par1);  
    */ 
}

  // 3) add contribution of conversions on the first dynode
  par1[0]=mean_pe1;
  val=(1-frac)*val+frac*pmtpdf_bessel(x,par1);
  val*=eff*bw;
 
  return val;
}

/* This is the pmt response measured in ADC counts. Perhaps the conversion factor
   should be one of the parameters set by the user.
 */
Double_t pmtpdf_besselg1adc(Double_t* x, Double_t* par){
  //*x is the signal in ADC
  Double_t ADCoverEle=50./1e6;//estimate of conversion from ADC to electrons from circuit simulation
  //transform everything in electrons 
  Double_t parele[11];

  parele[0]=par[0]; //mean number of pe at the photocathode
  parele[1]=par[1]/ADCoverEle; //overall gain of the pmt in ADC
  parele[2]=par[2]; //number of dynodes
  parele[3]=par[3]/ADCoverEle; //pedestal  N= ADC offset
  parele[4]=par[4]/ADCoverEle;//sigma noise 
  parele[5]=par[5];//absolute count normalization
  parele[6]=par[6]; //gain of the fisrt dynode
  parele[7]=par[7]; //poisson average of photons converting on dynode 1
  parele[8]=par[8]; //probability to convert on dynode 1  
  parele[9]=par[9];
  parele[10]=par[10]/ADCoverEle;  //bin width, needed for normalization
  
  Double_t xele[1];
  *xele=*x/ADCoverEle;
  
  return pmtpdf_besselg1(xele,parele);
}

/* 
   Function to plot the pmtpdf_besselg1adc to test its behavior,
   useful to check the normalization of the function.
   Note that the pmtpdf_xxx have some approximations
   needed for speeding up calculations. This could 
   lead to bad normalization.
 */

TF1* 
testpmtpdf(double mean_pe=0.5, double G=50, int ndynodes=12, 
	   double ped=0, double noise=5, double norm=1e6, 
	   double g1=5.8, double mean_pe1=0.5, double frac=2.5/100, 
	   double eff=1., double bw=1.){

  double xmax=(3+mean_pe+5*sqrt(mean_pe+pow(noise/G,2)))*G;
  cout<<"xmax = "<<xmax<<endl;
  double xmin=-5*noise;

  TF1* fg=new TF1("pmtpdf_besselg1adc",pmtpdf_besselg1adc,xmin,xmax,11); 
  fg->SetParNames( "mean_pe","G","ndynodes","ped","noise","norm","g1", "mean_pe1", "frac", "eff", "bw");
  fg->SetParameters(mean_pe,  G,  ndynodes,  ped,  noise,  norm , g1,   mean_pe1,   frac,   eff,   bw); 

  /* plotting pmtpdf_bessel instead of pmtpdf_besselg1adc
  double ADCoverEle=50./1e6;
  xmin/=ADCoverEle; xmax/=ADCoverEle;
  TF1* fg=new TF1("pmtpdf_bessel",pmtpdf_bessel,xmin,xmax,6); 
  fg->SetParNames( "mean_pe","G","ndynodes","ped","noise","norm");//,"g1", "mean_pe1", "frac", "eff", "bw");
  fg->SetParameters(mean_pe,  G/ADCoverEle,  ndynodes,  ped/ADCoverEle,  noise/ADCoverEle,  norm );//, g1,   mean_pe1,   frac,   eff,   bw); 
  */

  fg->SetNpx(1000); 
  fg->SetLineColor(kBlue); 
  fg->Print();
  cout<<"drawing "<<endl;
  fg->Draw();
  cout<<"computing integral :";   double integral=fg->Integral(xmin,xmax); 
  cout<< "  integral = "<<integral<<endl;
  cout<< "  difference = "<< (norm-integral)/norm*100<<"% difference" <<endl;
  cout<<"computing noise fraction "; cout<<  fg->Integral(xmin,ped+3*noise)/integral*100<<" %"<<endl;
  return fg;
}
