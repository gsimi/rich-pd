 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * Simple Polya PDF
  * author: Kyle Cranmer <cranmer@cern.ch>
  *                                                                           * 
  *****************************************************************************/ 

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Polya pdf
// END_HTML
//

#include <iostream> 

#include "myRooPolya.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 

#include "RooRandom.h"
#include "RooMath.h"
#include "TMath.h"
#include "Math/ProbFuncMathCore.h"

using namespace std;

ClassImp(myRooPolya) 



//_____________________________________________________________________________
myRooPolya::myRooPolya(const char *name, const char *title, 
		       RooAbsReal& _x,
		       RooAbsReal& _x0,
		       RooAbsReal& _sigma,
		       RooAbsReal& _mean,
		       Bool_t noRounding) :
  RooAbsPdf(name,title), 
  x("x","x",this,_x),
  x0("x0","x0",this,_x0),
  sigma("sigma","sigma",this,_sigma),
  mean("mean","mean",this,_mean),
  _noRounding(noRounding),
  _protectNegative(false)
{ 
  // Constructor  
} 



//_____________________________________________________________________________
 myRooPolya::myRooPolya(const myRooPolya& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   x0("x0",this,other.x0),
   sigma("sigma",this,other.sigma),
   mean("mean",this,other.mean),
   _noRounding(other._noRounding),
   _protectNegative(other._protectNegative)
{ 
   // Copy constructor
} 




//_____________________________________________________________________________
Double_t myRooPolya::evaluate() const 
{ 
  // Implementation in terms of the TMath Polya function

  Double_t k = _noRounding ? (x-x0)/sigma : floor((x-x0)/sigma);  
  if(_protectNegative && mean<0) 
    return 1e-3;
  return TMath::Polya(k,mean) ;
} 





//_____________________________________________________________________________
Int_t myRooPolya::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
{
  if (matchArgs(allVars,analVars,x)) return 1 ;
  return 0 ;
}



//_____________________________________________________________________________
Double_t myRooPolya::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  assert(code==1) ;

  if(_protectNegative && mean<0) 
    return exp(-2*mean); // make it fall quickly

  // Implement integral over x as summation. Add special handling in case
  // range boundaries are not on integer values of x
  Double_t xmin = (x.min(rangeName) - x0)/sigma;
  Double_t xmax = (x.max(rangeName) - x0)/sigma;

  // Protect against negative lower boundaries
  if (xmin<0) xmin=0 ;

  Int_t ixmin = Int_t (xmin) ;
  Int_t ixmax = Int_t (xmax)+1 ;

  Double_t fracLoBin = 1-(xmin-ixmin) ;
  Double_t fracHiBin = 1-(ixmax-xmax) ;

  if (!x.hasMax()) {
    if (xmin<1e-6) {
      return 1 ;
    } else {
      
      // Return 1 minus integral from 0 to x.min() 

      if(ixmin == 0){ // first bin
	return TMath::Polya(0, mean)*(xmin-0);
      }      
      Double_t sum(0) ;
      sum += TMath::Polya(0,mean)*fracLoBin ;
      sum+= ROOT::Math::poisson_cdf(ixmin-2, mean) - ROOT::Math::poisson_cdf(0,mean) ;
      sum += TMath::Polya(ixmin-1,mean)*fracHiBin ;
      return 1-sum ;
    }
  }
  
  if(ixmin == ixmax-1){ // first bin
    return TMath::Polya(ixmin, mean)*(xmax-xmin);
  }  

  Double_t sum(0) ;
  sum += TMath::Polya(ixmin,mean)*fracLoBin ;
  if (RooNumber::isInfinite(xmax)){
    sum+= 1.-ROOT::Math::poisson_cdf(ixmin,mean) ;
  }  else {
    sum+= ROOT::Math::poisson_cdf(ixmax-2, mean) - ROOT::Math::poisson_cdf(ixmin,mean) ;
    sum += TMath::Polya(ixmax-1,mean)*fracHiBin ;
  }
  
  return sum ;

}







//_____________________________________________________________________________
Int_t myRooPolya::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  // Advertise internal generator in x

  if (matchArgs(directVars,generateVars,x)) return 1 ;  
  return 0 ;
}



//_____________________________________________________________________________
void myRooPolya::generateEvent(Int_t code)
{
  // Implement internal generator using TRandom::Polya 

  assert(code==1) ;
  Double_t xgen ;
  while(1) {    
    xgen = RooRandom::randomGenerator()->Polya(mean);
    if (xgen<=x.max() && xgen>=x.min()) {
      x = xgen ;
      break;
    }
  }
  return;
}


