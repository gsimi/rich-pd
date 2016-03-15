 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * Simple Poisson PDF
  * author: Kyle Cranmer <cranmer@cern.ch>
  *                                                                           * 
  *****************************************************************************/ 

#ifndef myROOPOISSON
#define myROOPOISSON

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class myRooPoisson : public RooAbsPdf {
public:
  myRooPoisson() { _noRounding = kFALSE ; } ;
  myRooPoisson(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _x0, RooAbsReal& _sigma, RooAbsReal& _mean, Bool_t noRounding=kFALSE);
  myRooPoisson(const myRooPoisson& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new myRooPoisson(*this,newname); }
  inline virtual ~myRooPoisson() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  void generateEvent(Int_t code);
  
  void setNoRounding(bool flag = kTRUE){_noRounding = flag;}
  void protectNegativeMean(bool flag = kTRUE){_protectNegative = flag;}

protected:

  RooRealProxy x ;
  RooRealProxy x0 ;
  RooRealProxy sigma ;
  RooRealProxy mean ;
  Bool_t  _noRounding ;
  Bool_t  _protectNegative ;
  
  Double_t evaluate() const ;
  Double_t evaluate(Double_t k) const;
  

private:

  ClassDef(myRooPoisson,3) // A Poisson PDF
};
 
#endif
