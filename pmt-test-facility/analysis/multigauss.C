#include <iostream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TPaveStats.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TSpectrum.h"

#include "TROOT.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "TCanvas.h"
#include "RooAddPdf.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooConstVar.h"
#include "RooProduct.h"
#include "RooFormulaVar.h"
#include "RooFitResult.h"
#include "RooCustomizer.h"
#include "mypdfs/myRooGamma.h"
#include "mypdfs/MyExp.h"
#include "mypdfs/myRooPoisson.h"
#include "mypdfs/myRooPolya.h"
#include "TROOT.h"
using namespace RooFit ;

/* TODO: io farei una bella classe... poi fai te!
 *
 * */

#define DIM 8
#define NPX DIM*DIM

// Fwd declaration
void dofit( TString, TString, Int_t );
RooFitResult* getFit(TH1F*, TString, Int_t );
void loadHamamatsuGain( TString serialNo );
void loadChannelMap( );

// Hamamatsu gains
Float_t HamamatsuGain[NPX]; 
//SerialNoMap[0]["9726E473"] = "/lhcbdata/rich/data/FA0011.txt"; 
//SerialNoMap[1]["9726E473"] = "/lhcbdata/rich/data/ZN0971.txt"; 
//SerialNoMap[0]["C464BF09"] = "/lhcbdata/rich/data/FA0002.txt"; 
//SerialNoMap[1]["C464BF09"] = "/lhcbdata/rich/data/FA0005.txt"; 

// Pixel id to MAROC channel map
Int_t pixel2ChannelMap[2][NPX];
Int_t pixel2ChipMap[2][NPX];

// Output branches
Int_t m_npx = NPX;
Int_t m_status[NPX] = {0};
Float_t m_chi2[NPX] = {0.};
Float_t m_npe[NPX] = {0.}, m_dnpe[NPX] = {0.};   
Float_t m_npe_ct[NPX] = {0.}, m_dnpe_ct[NPX] = {0.};   
Float_t m_mean_n[NPX] = {0.}, m_dmean_n[NPX] = {0.};   
Float_t m_sigma_n[NPX] = {0.}, m_dsigma_n[NPX] = {0.};   
Float_t m_mean_1pe[NPX] = {0.}, m_dmean_1pe[NPX] = {0.};   
Float_t m_sigma_1pe[NPX] = {0.}, m_dsigma_1pe[NPX] = {0.};   
Float_t m_mean_1pe_ct[NPX] = {0.}, m_dmean_1pe_ct[NPX] = {0.};   
Float_t m_sigma_1pe_ct[NPX] = {0.}, m_dsigma_1pe_ct[NPX] = {0.};   


//======================================================
// Multigauss main routine
//======================================================
void multigauss( TString dataFile, 
                 TString lab, 
                 TString serialNo,
                 Int_t pxID )
{
  // usage:
  // source roo_setup.csh
  // .x multigauss.C
  // runs on default pixel 32, if available
  // .x multigauss.C(12)
  // runs on pixel 12, if available
  // .x multigauss.C(-1)
  // runs on all the available pixels

  //Load Libs: TODO make a compiled dynamically loaded library
  /* 
     gROOT->ProcessLine(".L mypdfs/myRooGamma.cxx+");
     gROOT->ProcessLine(".L mypdfs/myRooPoisson.cxx+");
     gROOT->ProcessLine(".L mypdfs/myRooPolya.cxx+");
  */


  // Get Hamamatsu gains
  loadHamamatsuGain( serialNo );

  // Load channel map
  loadChannelMap( );

  // Create output file and tree 
  TFile *outFile = new TFile("multigauss.root", "recreate"); 
  TTree *tree = new TTree("tree", "");

  // Set branches
  tree->Branch("npx", &m_npx, "npx/I");

  tree->Branch("HamGain", &HamamatsuGain, "HamGain[npx]/F");

  tree->Branch("status", &m_status, "status[npx]/I");
  tree->Branch("chi2", &m_chi2, "chi2[npx]/F");

  tree->Branch("npe", &m_npe, "npe[npx]/F");
  tree->Branch("dnpe", &m_dnpe, "dnpe[npx]/F");

  tree->Branch("npe_ct", &m_npe_ct, "npe_ct[npx]/F");
  tree->Branch("dnpe_ct", &m_dnpe_ct, "dnpe_ct[npx]/F");

  tree->Branch("mean_n", &m_mean_n, "mean_n[npx]/F");
  tree->Branch("dmean_n", &m_dmean_n, "dmean_n[npx]/F");

  tree->Branch("sigma_n", &m_sigma_n, "sigma_n[npx]/F");
  tree->Branch("dsigma_n", &m_dsigma_n, "dsigma_n[npx]/F");

  tree->Branch("mean_1pe", &m_mean_1pe, "mean_1pe[npx]/F");
  tree->Branch("dmean_1pe", &m_dmean_1pe, "dmean_1pe[npx]/F");

  tree->Branch("sigma_1pe", &m_sigma_1pe, "sigma_1pe[npx]/F");
  tree->Branch("dsigma_1pe", &m_dsigma_1pe, "dsigma_1pe[npx]/F");

  tree->Branch("mean_1pe_ct", &m_mean_1pe_ct, "mean_1pe_ct[npx]/F");
  tree->Branch("dmean_1pe_ct", &m_dmean_1pe_ct, "dmean_1pe_ct[npx]/F");

  tree->Branch("sigma_1pe_ct", &m_sigma_1pe_ct, "sigma_1pe_ct[npx]/F");
  tree->Branch("dsigma_1pe_ct", &m_dsigma_1pe_ct, "dsigma_1pe_ct[npx]/F");

  // Fit only one px 
  if (pxID != -1 ) {
    dofit( dataFile, lab, pxID );
  }

  // Fit all pixels
  else {
    for (Int_t px = 1; px <= NPX; px++) {  
      dofit( dataFile, lab, px ); 
    }
  }

  // Close output file
  outFile->cd();
  tree->Fill();
  tree->Write();
  outFile->Close();
}


//======================================================
// Perform the fit
//======================================================
void dofit( TString dataFile, 
            TString lab, 
            Int_t pxID )
{
  gROOT->cd();

  TH1F *hh;
  Int_t dump_pdf_output = 1;
  Int_t nEvents = 0;
  

  //===================================
  // Prepare dataset
  //===================================
  // Padova data 
  if (lab == "Padova") {

    TChain ch("tree");
    ch.Add( dataFile );
    Float_t integral; ch.SetBranchAddress("integral", &integral);
    Float_t tauFall; ch.SetBranchAddress("tauFall", &tauFall);
    nEvents = ch.GetEntries();
    hh = new TH1F("hh"," ",150,-50.,800.);

    for (Int_t i = 0; i < nEvents; i++) {
      ch.GetEvent(i); 
      //if (tauFall < 125.0) continue; // Suppress cross-talk
      hh->Fill(integral);
    }
  }

  // MAROC data
  else if (lab == "MAROC") {

    TChain ch("mdfTree");
    ch.Add(dataFile);

    // ADC for pixels
    UShort_t B01_PMT1_ADC[64]; ch.SetBranchAddress("B01_PMT1_ADC",&B01_PMT1_ADC[0]);
    UShort_t B01_PMT2_ADC[64]; ch.SetBranchAddress("B01_PMT2_ADC",&B01_PMT2_ADC[0]);

    // ADC for MAROC channel
    UShort_t B01_MAR1_ADC[64]; ch.SetBranchAddress("B01_MAR1_ADC",&B01_MAR1_ADC[0]);
    UShort_t B01_MAR2_ADC[64]; ch.SetBranchAddress("B01_MAR2_ADC",&B01_MAR2_ADC[0]);

    nEvents = ch.GetEntries(); 
    hh = new TH1F("hh"," ", 210, -10., 200.);
    hh_ct = new TH1F("hh_ct"," ", 210, -10., 200.);

    // Loop over events
    for (Int_t i = 0; i < nEvents; i++) {
      ch.GetEvent(i); 
      // Standard
      hh->Fill( float(B01_PMT1_ADC[pxID-1]) ); // Read only PMT1...

      // Suppress MaPMT cross-talk
      //bool crosstalk = false;
      //for (int ip = -1; ip < 2; ip++) {
      //  for (int jp = -1; jp < 2; jp++) {
      //    int nbID = pxID + ip*1 + jp*8; // neighbour px
      //    if (nbID == pxID) continue;
      //    //cout << ip << " " << jp << " pxID " << pxID << " nbID " << nbID << endl;
      //    if (nbID > 0 && nbID <= 64 && B01_PMT1_ADC[nbID-1] > 30) crosstalk = true;
      //  }
      //}
      //if (crosstalk == false) hh->Fill( float(B01_PMT1_ADC[pxID-1]) );
      //else hh2->Fill( float(B01_PMT1_ADC[pxID-1]) );

      // Suppress MAROC cross-talk
      Int_t chan = pixel2ChannelMap[0][pxID-1];
      Int_t chip = pixel2ChipMap[0][pxID-1];
     
      UShort_t* ADC_count; // retrieve corresponding chip
      if (chip == 0) ADC_count = B01_MAR1_ADC;
      else ADC_count = B01_MAR2_ADC;
    
      hh_ct->Fill( float(ADC_count[chan]) );

      //if ((chan != NPX - 1 && ADC_count[chan + 1] > 20) ||
      //    (chan != 0 && ADC_count[chan - 1] > 20)) {
      //  hh_ct->Fill( float(B01_PMT1_ADC[pxID-1]) );  
      //}
    }//evts

    // Plot cross-talk spectrum
    //TCanvas* c_ct = new TCanvas("c_ct", "cross-talk spectra", 1000, 500);
    //c_ct->SetLogy();
    //hh_ct->Draw("histo");
  }

  // Edinburgh data
  else if (lab == "Edinburgh") {
    // dataFile="../SpectraFitting_Edinburgh/inputs/gainMap/FA0002_bl.root"; // high gain ?
    // dataFile="../SpectraFitting_Edinburgh/inputs/gainMap/ZN0971_bl.root"; // low gain
    //dataFile="../SpectraFitting_Edinburgh/inputs/gainMap/FA0006_bl.root"; // low gain

    TChain ch("QDCReadoutTree");
    ch.Add(dataFile);
    Int_t Value[32]; ch.SetBranchAddress("Value", &Value[0]);
    nEvents = ch.GetEntries();
    hh = new TH1F("hh"," ",500,0.,1000.);

    for (Int_t i = 0; i < nEvents; i++){
      ch.GetEvent(i); // cout <<  Value[pxID-1] << endl;
      hh->Fill( float(Value[pxID-1]) );
    }
  }

  // Not known!
  else {
    cout << "Setup/lab not known! Exit!" << endl;
    exit(0);
  }

  // Check if file exist!
  if (gSystem->AccessPathName(dataFile)) {
    cout << "File " << dataFile << " does not exist! Exit!" << endl;
    exit(0);
  }

  cout << "Analyzing data file: " << dataFile << "..." << endl;
  cout << "Number of events: " << nEvents << endl;

  RooFitResult* result =  getFit(hh,lab, pxID);
}
  
RooFitResult* getFit( TH1F *hh,
		      TString lab,
		      int pxID)
{

  Int_t dump_pdf_output = 1;
  Double_t lowe = 0.0;  Double_t cent = 0.0; Double_t uppe = 0.0;
  Double_t roofit_x_min = 0.;
  Double_t roofit_x_max = 0.;
  if (lab == "Padova"){
    roofit_x_min = -50.; roofit_x_max = 1050.;
  }   else if (lab == "MAROC") {
    roofit_x_min = -50.; roofit_x_max = 1050.;
  }  else if (lab == "Edinburgh") {
    roofit_x_min = 0.; roofit_x_max = 1000.;
  }

  hh->SetMinimum(1.);
  Double_t hh_y_max = hh->GetMaximum(); // cout << hh_y_max << endl;
  Double_t hh_x_max = hh->GetXaxis()->GetXmax(); // cout << hh_x_max << endl;


  //==========================
  // Pattern recognition
  //==========================
  Double_t x_peaks[2]={0.,200.}; // initial values
  TSpectrum *s = new TSpectrum(4);
  Int_t npeaks = s->Search(hh,4," ",5e-4); // used to configure the Search 
  cout << "Found " << npeaks << " candidate peaks !" << endl;
  Double_t *xpeaks = s->GetPositionX();
  for (Int_t p=0;p<npeaks;p++) {
    x_peaks[p] = static_cast<Double_t>(xpeaks[p]);
    Int_t bin = hh->GetXaxis()->FindBin(x_peaks[p]);
    Double_t yp = hh->GetBinContent(bin);
    cout << bin << " " << x_peaks[p] << " " << yp << endl;
  }
  if( npeaks==1 ){
    cout << "uninitialized peak:" << endl;
    cout << x_peaks[1] << endl; 
  }

  //***valley with TSpectrum***
  TH1F *hhswap = (TH1F*)hh->Clone();

     Double_t BinContent=0.;
  for(Int_t n_bins=0;n_bins<hh->GetNbinsX();n_bins++){

    if( hh->GetBinCenter(n_bins+1)<x_peaks[0] ){
      BinContent=0.;
    }else{
      BinContent=(hh_y_max-hh->GetBinContent(n_bins+1));
      hhswap->SetBinContent((n_bins+1),BinContent);
    }

  }
  Double_t x_swap[2]={20.,1000.}; // initial values
  TSpectrum *swap = new TSpectrum(2);
  npeaks = swap->Search(hhswap,2," ",1e-6); // used to configure the Search
  cout << "Found " << npeaks << " <<swapped>> candidate peaks !" << endl;
  Double_t *xswap = swap->GetPositionX();
  for (Int_t ps=0;ps<2;ps++) {
    x_swap[ps] = static_cast<Double_t>(xswap[ps]);
    Int_t bin = hhswap->GetXaxis()->FindBin(x_swap[ps]);
    Double_t yp = hhswap->GetBinContent(bin);
    cout << bin << " " << x_swap[ps] << " " << yp << endl;
  }
  
  

  TCanvas *cInput = new TCanvas("cInput","Pattern reco.",600,600); cInput->Divide(1,2);
  cInput->cd(1); hh->Draw();
  cInput->cd(2); hhswap->Draw();


  //================================
  // Declare fit variable and data
  //================================
  RooRealVar x("x", "x", roofit_x_min, roofit_x_max);
  RooDataHist dh("dh", "dh", x, Import(*hh));

  
  //===================
  // Noise pdf params
  //===================
  cent = x_peaks[0]; lowe = -5.*TMath::Abs(x_peaks[0]); uppe = +5.*TMath::Abs(x_peaks[0]);
  RooRealVar mean_n("mean_n", "mean of noise gaussians", cent, lowe, uppe); // MAROC
  RooRealVar sigma_n("sigma_n", "width of noise gaussians", 1.0, 0.5, 10.0); 
  //RooRealVar mean_n("mean_n", "mean of noise gaussians", 0.0, -20.0, +20.0); // Padova
  //RooRealVar sigma_n("sigma_n", "width of noise gaussians", 10.0, 5.0, 30.); 


  //========================
  // Multi-gaussian params
  //========================
  // mean number of p.e. per pixel
  RooRealVar npe("npe", "mean number of p.e.", 0.20, 1e-3, 1.0);

  // mean and error of 1 p.e. component for signal
  cent = x_peaks[1]; lowe = 0.5*cent; uppe = 1.5*cent;
  RooRealVar mean_1pe("mean_1pe", "mean of 1pe signal gaussian", cent, lowe, uppe);
  RooRealVar sigma_1pe("sigma_1pe", "width of 1pe signal gaussian", 10., 5.0, 20.0); // MAROC
  //  RooFormulaVar sigma_1pe("sigma_1pe", "mean_1pe*0.527", mean_1pe); //0.527=sqrt(g/g1/(g-1)*(1-1./g+1./g2)) fix its value from pmtsim.C a

  // mean and error of 1 p.e. component for cross talk
  RooRealVar mean_1pe_ct("mean_1pe_ct", "mean of 1pe cross-talk gaussian", 5., 0.0, 10.0); // MAROC
  //RooRealVar sigma_1pe_ct("sigma_1pe_ct", "width of 1pe cross-talk gaussian", 1, 0.0, 10.); // fit prefers small values...
  RooFormulaVar sigma_1pe_ct("sigma_1pe_ct", "mean_1pe_ct*0.4", mean_1pe_ct); //0.527=sqrt(g/g1/(g-1)*(1-1./g+1./g2)) fix its value from pmtsim.C and also from theorethical calculation using gain = 1e6 
  // mean number of neighbouring pixels
  // RooRealVar npx("npx", "mean number of neighbouring pixels", 1.0, 1e-3, 10.0); 
  // mean number of cross talk signals per pixels
  RooRealVar npe_ct("npe_ct", "mean number of xtalk signals per pixels", 0.1, 1e-3, 10.0); 


  // Fraction of the extra component for MAROC spectra
  RooRealVar *frac_extra = new RooRealVar("frac_extra", "Fraction of the extra gaussian/tot", 0.0); // fixed
  
  if (lab == "MAROC") {
    frac_extra->setVal(0.1);
    frac_extra->setRange(0.,1.);
    frac_extra->setConstant(false);
  }
  //fraciton of events in the tail of the gaus function
  RooRealVar *frac_tail = new RooRealVar("frac_tail", "Fraction of tails of gaussian", 0.0,0.0,1.0); 
  RooRealVar *sigma_tail = new RooRealVar("sigma_tail","sigma of the tail gaussian for signal model",60,5,200);

  //==============
  // CONFIGURATION
  //==============
  bool useextragamma=false;
  bool usetail=false;
  bool multistagefit=false;
  enum imodel {gaussmodel, poissonmodel};
  int signal_model=gaussmodel;

  //==================
  // Build pdfs
  //==================
  const int NGauss = 5;
  const int NGaussct = 3;
  RooArgList *pdfs = new RooArgList;
  RooArgList *coeffs = new RooArgList;
  TString name, expr, title;

  switch(signal_model){
  case gaussmodel:
  for (int ns = 0; ns < NGauss + 1 ; ns++) {
    for (int nct = 0; nct < NGaussct + 1; nct++) {
      name.Form("mean_G%i_%i", ns, nct); // mean
      expr.Form("mean_n + %i*mean_1pe + %i*mean_1pe_ct", ns, nct);
      RooFormulaVar *mean = new RooFormulaVar(name, expr, RooArgList(mean_n, mean_1pe, mean_1pe_ct));

      name.Form("sigma_G%i_%i", ns, nct); // sigma
      expr.Form("sqrt(sigma_n*sigma_n + %i*sigma_1pe*sigma_1pe + %i*sigma_1pe_ct*sigma_1pe_ct)", ns, nct);
      RooFormulaVar *sigma = new RooFormulaVar(name, expr, RooArgList(sigma_n, sigma_1pe, sigma_1pe_ct));
      

	name.Form("G%i_%i", ns, nct); // gaussian
	TString title; title.Form("Gaussian component %i-%i", ns, nct);
	RooAbsPdf *gauss = new RooGaussian(name, title, x, *mean, *sigma); 

	name.Form("T%i_%i", ns, nct); // tail gaussian
        title.Form("Gaussian tail component %i-%i", ns, nct);
	RooAbsPdf *tail = new RooGaussian(name, title, x, *mean, *sigma_tail); 
      
	name.Form("S%i_%i", ns, nct); // gauss + tail
	title.Form("signal component %i-%i", ns, nct);
	RooAbsPdf *signal;
	if (usetail)
	  signal=new RooAddPdf(name,title,RooArgList(*gauss,*tail),RooArgList(*frac_tail));
	else{
	  signal=gauss;
	  signal->SetName(name);
	  signal->SetTitle(title);
	}

	name.Form("N_G%i_%i", ns, nct); // coefficient
	expr.Form("TMath::Poisson(%i, npe)*TMath::Poisson(%i, npe_ct)", ns, nct);
	RooFormulaVar *coeff = new RooFormulaVar(name, expr, RooArgList( npe, npe_ct));
      
	pdfs->add( *signal );
	coeffs->add( *coeff );
      }
    }
    break;
  case poissonmodel:
    
    for (int ns = 0; ns < NGauss + 1 ; ns++) {
      for (int nct = 0; nct < NGaussct + 1; nct++) {
	RooRealVar *g1  = new RooRealVar("g1", "mg1",  6,   1., 10.);

	name.Form("sigma_G%i_%i", ns, nct); // sigma
	expr.Form("sqrt(sigma_n*sigma_n + %i*sigma_1pe*sigma_1pe + %i*sigma_1pe_ct*sigma_1pe_ct )", ns, nct);
	RooFormulaVar *sigma = new RooFormulaVar(name, expr, RooArgList(sigma_n, sigma_1pe, sigma_1pe_ct));

	name.Form("S%i_%i", ns, nct); // poisson
	TString title; title.Form("Poisson component %i-%i", ns, nct);
	myRooPoisson *signal = new myRooPoisson(name,title,x,mean_n,*sigma,*g1,true);

	name.Form("N_G%i_%i", ns, nct); // coefficient
	expr.Form("TMath::Poisson(%i, npe)*TMath::Poisson(%i, npe_ct)", ns, nct);
	RooFormulaVar *coeff = new RooFormulaVar(name, expr, RooArgList( npe, npe_ct));

	pdfs->add( *signal );
	coeffs->add( *coeff );
      }
    }
    
    break;
  }
  if (useextragamma){
    // Add an extra GAMMA to fit MAROC data
    if (lab == "MAROC") {

    //cent = 0.5*(x_peaks[0] + x_peaks[1]); lowe = 0.5*cent; uppe = 1.1*cent;  // better using pos. last ct gaussian

  //RooFormulaVar *mean = new RooFormulaVar("mean_Gextra", "mean_n + 3*mean_1pe_ct + mean_extra", RooArgList(mean_n, mean_1pe_ct, *mean_extra));
  //RooFormulaVar *sigma = new RooFormulaVar("sigma_Gextra", "sqrt(sigma_n*sigma_n + sigma_extra*sigma_extra)", RooArgList(sigma_n, *sigma_extra));
  //extra = new RooGaussian("extra", "Extra gaussian component", mean, sigma)

  //  // Step function
  //  //leftEdge = new RooRealVar("leftEdge", "Left edge", 20., 15., 30.); 
  //  //width = new RooRealVar("width", "Width", 5., 0., 20.); 
  //  //rightEdge = new RooFormulaVar("rightEdge", "(leftEdge + width)", RooArgList(*leftEdge, *width)); 
  //  //stepFun = new RooGenericPdf("stepFun", "step function PDF", "(@0 >= @1) && (@0 < @2)", RooArgList(x, *leftEdge, *rightEdge));

  //  // Try with expo.
  //  //dxped = new RooFormulaVar("dxped", "x - mean_n", RooArgList(x, mean_n));
  //  //sl = new RooRealVar("sl", "sl", 1./cent, 1e-1/x_peaks[1], 10/x_peaks[1]); 
  //  //extra = new RooGenericPdf("extra", "Extra", "(dxped > 0.) * sl * exp(-sl*x)", RooArgList(x, *dxped, *sl));

  //  //cent = 0.5*(x_peaks[0] + x_peaks[1]); 
  //  //sl = new RooRealVar("sl", "exponential coefficient", 1./cent, 1e-1/x_peaks[1],10/x_peaks[1] ); 
  //  //dxped = new RooFormulaVar("dxped", "x - mean_n", RooArgList(x, mean_n));
  //  //extra = new RooGenericPdf("extra", "Extra ep component", "theta*sl*TMath::Exp(-sl*dxped)",RooArgSet(theta, *dxped, *sl));

  //  //cent = 0.5*(x_peaks[0] + x_peaks[1]); 
  //  //RooRealVar alpha("alpha", "exponential coefficient", 1./cent, 1e-1/x_peaks[1],10/x_peaks[1] ); 
  //  //RooFormulaVar xped("xped","x-mean_n",RooArgList(x,mean_n));
  //  //extra = new RooGenericPdf("extra","Extra ep component","theta*alpha*TMath::Exp(-alpha*xped)",RooArgSet(theta,xped,alpha));
  //  
  //  // Try with convolution
  //  //csigma = new RooRealVar("csigma", "width of gaussian", 0.1); 
  //  //cgauss = new RooGaussian("cgauss", "Gaussian", x, RooConst(0.), *csigma); 
  //  //extra = new RooFFTConvPdf("extra", "Extra step (X) gaussian", x, *stepFun, *cgauss) ; 

  //  gammapar = new RooRealVar("gammapar", "gamma", 10, 0, 100);
  //  betapar = new RooRealVar("betapar", "beta", 1, 0, 10);
  //  extra = new myRooGamma("extra", "gammapdf", x, *gammapar, *betapar, mean_n);

  //  pdfs->add ( *extra );
  //  coeffs->add( *frac_extra );
      RooRealVar *gammapar = new RooRealVar("gamma","gamma", 1.25, 1., 100.);
      RooRealVar *betapar  = new RooRealVar("beta", "beta",  40,   0., 100.);
      myRooGamma* extra = new myRooGamma("extra", "gammapdf", x, *gammapar, *betapar, mean_n);
      
      pdfs->add ( *extra );
      //      coeffs->add(*frac_extra);
    }
  }
  // build signal model as sum of pdfs
  RooAddPdf model("model", "Model for signal + cross-talk + noise", *pdfs, *coeffs);


  //========
  // Fit!
  //========
  RooPlot* frame = x.frame( Title(" ") );
  dh.plotOn( frame );
  // x.setRange("R1",0,22);
  // model.fitTo(dh,Range("R1,R2")) ;
  model.fitTo(dh) ;
  if (multistagefit){
    sigma_n.setConstant(true);
    mean_n.setConstant(true);
    x.setRange("R2",30,200);
    model.fitTo(dh,Range("R2")) ;
  }

  //==================
  // Plot components
  //==================
  model.plotOn(frame, "", LineColor(2)); 
  for (int ns = 0; ns < NGauss + 1; ns++) {
    for (int nct = 0; nct < NGaussct + 1; nct++) {
      name.Form("S%i_%i", ns , nct ); 
      model.plotOn(frame, Components(name), LineColor(ns == 0 ? 1 : ns + 3), LineStyle(nct + 1));
    }   
  }

  if (lab == "MAROC") {
    model.plotOn(frame, Components("extra"), LineColor(kOrange), LineStyle(kSolid));
  }

  // get results 
  RooFitResult* fitres = model.fitTo(dh, RooFit::Save(true), Strategy(2));
  fitres->Print(); 
  fitres->floatParsFinal().Print("s") ;
  TCanvas* cSpec_zoom = new TCanvas("Spec_zoom","spectra zoom",1000,500) ;
  cSpec_zoom->SetLogy();
  frame->SetMinimum(1.);
  frame->Draw();
  if( dump_pdf_output==1 ){
    if( pxID==1 ){ cSpec_zoom->SaveAs("dump.pdf["); }
    if( pxID>=1 && pxID<=64 ){ cSpec_zoom->SaveAs("dump.pdf"); }
    if( pxID==64 ){ cSpec_zoom->SaveAs("dump.pdf]"); }
  }
  cout << "chi2 = " << frame->chiSquare() << endl;
  cout << "status = " << fitres->status() << endl; 

  TLatex caption;
  caption.SetTextSize(0.04);
  const char *cap_str = "#color[3]{green: 1 p.e.}, #color[4]{blue: 2 p.e.}, #color[5]{yellow: 3 p.e.}, #color[6]{pink: 4 p.e.}";
  caption.DrawLatex(0.1*hh_x_max,0.1*hh_y_max,cap_str);

  TCanvas* cSpec = new TCanvas("Spec","spectra",1000,500);
  cSpec->SetLogy();
  hh->Draw("histo");
  frame->Draw("same");


  //==============
  // Fill NTuple 
  //==============
  m_status[pxID -1] = fitres->status();
  m_chi2[pxID -1] = frame->chiSquare();

  m_npe[pxID -1] = npe.getVal();
  m_dnpe[pxID -1] = npe.getPropagatedError(*fitres);

  m_npe_ct[pxID -1] = npe_ct.getVal();
  m_dnpe_ct[pxID -1] = npe_ct.getPropagatedError(*fitres);
  
  m_mean_n[pxID -1] = mean_n.getVal();
  m_dmean_n[pxID -1] = mean_n.getPropagatedError(*fitres);
  
  m_sigma_n[pxID -1] = sigma_n.getVal();
  m_dsigma_n[pxID -1] = sigma_n.getPropagatedError(*fitres);

  m_mean_1pe[pxID -1] = mean_1pe.getVal();
  m_dmean_1pe[pxID -1] = mean_1pe.getPropagatedError(*fitres);

  m_sigma_1pe[pxID -1] = sigma_1pe.getVal();
  m_dsigma_1pe[pxID -1] = sigma_1pe.getPropagatedError(*fitres);

  m_mean_1pe_ct[pxID -1] = mean_1pe_ct.getVal();
  m_dmean_1pe_ct[pxID -1] = mean_1pe_ct.getPropagatedError(*fitres);

  m_sigma_1pe_ct[pxID -1] = sigma_1pe_ct.getVal();
  m_dsigma_1pe_ct[pxID -1] = sigma_1pe_ct.getPropagatedError(*fitres);

  return fitres;
}


//========================================
// Load Hamamatsu gains 
//========================================
void loadHamamatsuGain( TString serialNo )
{
  ifstream file( "/lhcbdata/rich/data/" + serialNo + ".txt", ios::in );
  string line;
  float val;

  std::cout << "Load Hamamatsu gains" << std::endl;

  for (size_t row = 0; row < DIM; row++) {
    std::getline( file, line );
    std::istringstream iss( line );

    for (size_t col = 0; col < DIM; col++) {
      iss >> val;
      HamamatsuGain[DIM*row + col] = val;
    }
  }
}


//========================================
// Load PMT pixels to MAROC channels map
//========================================
void loadChannelMap() 
{
  ifstream file( "MAROC3v2_channelMap.txt", ios::in );
  string line;
  Int_t chan, maroc, searray, px, pmt;

  std::cout << "Load channel map" << std::endl;

  // Skip first line w/ column description
  std::getline( file, line );

  while (std::getline( file, line )) {
    std::istringstream iss( line );

    iss >> chan >> maroc >> searray >> px >> pmt;
    pixel2ChannelMap[pmt-1][px-1] = chan; // chain: 0-63
    pixel2ChipMap[pmt-1][px-1] = maroc; // maroc chip: 0-1
    
    //std::cout << line << std::endl;
    //std::cout << ">>>>> " << pmt << " " << px << " " << chan << " " 
    //<< pixel2ChannelMap[pmt-1][px-1] << std::endl;
  }
}

