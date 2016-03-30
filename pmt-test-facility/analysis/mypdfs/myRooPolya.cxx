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

#include "mypdfs/myRooPolya.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 

#include "RooRandom.h"
#include "RooMath.h"
#include "TMath.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/PdfFuncMathCore.h"

using namespace std;

ClassImp(myRooPolya) 



//_____________________________________________________________________________
myRooPolya::myRooPolya(const char *name, const char *title, 
		       RooAbsReal& _x,
		       RooAbsReal& _x0,
		       RooAbsReal& _mean,
		       RooAbsReal& _b) :
  RooAbsPdf(name,title), 
  x("x","x",this,_x),
  x0("x0","x0",this,_x0),
  mean("mean","mean",this,_mean),
  b("b","b",this,_b),
  _protectNegative(false)
{ 
  // Constructor  
} 



//_____________________________________________________________________________
 myRooPolya::myRooPolya(const myRooPolya& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   x0("x0",this,other.x0),
   mean("mean",this,other.mean),
   b("b",this,other.b),
   _protectNegative(other._protectNegative)
{ 
   // Copy constructor
} 




//_____________________________________________________________________________
Double_t myRooPolya::evaluate() const 
{ 
  // Implementation in terms of the TMath Polya function

  if(_protectNegative && mean<0) 
    return 1e-3;
  return ROOT::Math::negative_binomial_pdf(x-x0,1/(1+mean*b),b) ;
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
  Double_t xmin = (x.min(rangeName) - x0);
  Double_t xmax = (x.max(rangeName) - x0);

  // Protect against negative lower boundaries
  if (xmin<0) xmin=0 ;
  
  double val=0;
  if (!x.hasMax() || RooNumber::isInfinite(xmax)) {
    if (xmin<1e-6) {
      return 1 ;
    } else {
      val = 1-ROOT::Math::negative_binomial_cdf_c(x-x0,1./(1+mean*b),1./b);
    }
  } else{
    val = ROOT::Math::negative_binomial_cdf(x-x0,1./(1+mean*b),1./b);
  }
  return val ;

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
    //xgen = RooRandom::randomGenerator()->Polya(mean);
    xgen=0; cout<<"generator not implemented"<<endl;
    if (xgen<=x.max() && xgen>=x.min()) {
      x = xgen ;
      break;
    }
  }
  return;
}


