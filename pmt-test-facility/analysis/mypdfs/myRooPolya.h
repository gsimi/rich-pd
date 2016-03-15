 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * Simple Polya PDF
  * author: Kyle Cranmer <cranmer@cern.ch>
  *                                                                           * 
  *****************************************************************************/ 

#ifndef myROOPOLYA
#define myROOPOLYA

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class myRooPolya : public RooAbsPdf {
public:
  myRooPolya() { _noRounding = kFALSE ; } ;
  myRooPolya(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _x0, RooAbsReal& _mean, RooAbsReal& _b);
  myRooPolya(const myRooPolya& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new myRooPolya(*this,newname); }
  inline virtual ~myRooPolya() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  void generateEvent(Int_t code);
  
  void protectNegativeMean(bool flag = kTRUE){_protectNegative = flag;}

protected:

  RooRealProxy x ;
  RooRealProxy x0 ;
  RooRealProxy mean ;
  RooRealProxy b ;
  Bool_t  _protectNegative ;
  
  Double_t evaluate() const ;
  Double_t evaluate(Double_t k) const;
  

private:

  ClassDef(myRooPolya,3) // A Polya PDF
};
 
#endif
