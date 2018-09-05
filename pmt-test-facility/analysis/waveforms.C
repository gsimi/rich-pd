//1;3409;0c
#include "TGraph.h"
#include "TGraphErrors.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "pmtwavef.cc"
#include "TMinuit.h"
using namespace::std;




//===================================================
// Channel class
//===================================================
class channel
{
public:
  // Number of records in the channel 
  const static int rlength = 1024;

  // Constructor and destructor 
  channel(const char* fData, 
          const char* fPed, 
          int nPed = 0);
  channel() {};
  ~channel() {};

  // Set channel according to input
  void setChannel(const char* fData, 
                  const char *fPed, 
                  int nPed = 0);

  // Load pedestals 
  void loadPedestals(const char* fPed, int nPed = 0);
  
  // Load next event 
  void loadNextEvent();
 
  // Load an event 
  void loadEvent(int index = 0);
  
  // Plot functions
  TGraph* plot(int index);
  TGraph* plot();
  TGraph* plotNext();
  TGraph* plotRaw();
  TGraph* plotPedestal();
  
  // Various ways to compute estimate height
  double integral (int sigmin = 200);
  double integral2 (int bkgmax = 200, int sigmax = rlength);
  double maxval (int bkgmax = 200);
  double minval (int nchannels = 4);
  double smoothmax (int nchannels = 4);
   
  // Fit wave using a user-defined function
  bool fit(TF1 *fun, int min = 0, int max = rlength);

  // Returns time
  float time(bool negative = false, int bkgmax = 0, int sigmmax = rlength);
  
  // Returns end-of-file 
  bool eof() { return waves.eof(); }

  // Returns if it's a binary file 
  bool isbinary(const char *fname);
  
  // Read a single float from stream 
  void readfloat(ifstream &is, float& val, bool binary);

  // Returns pedestals
  float* Pedestals() { return pedestals; }
  //float Pedestals (unsigned idx) { return pedestals[idx]; }
  
  // Returns raw values
  float* Rawval() { return rawval; }
  
  // Returns calibrated values
  float* Calibval() { return calibval; }
  
  // Returns times
  float* T() { return t; }


private:
  // Private functions
  void calibrate();

  // Data members
  float pedestals[rlength];
  float rawval[rlength];
  float calibval[rlength];
  float t[rlength];
  const char* fname;
  bool calibrated;
  ifstream waves;
  bool binaryfile;
};


//===================================
// Constructor
//===================================
channel::channel(const char* fData, const char *fPed, int nPed)
{
  setChannel(fData, fPed, nPed);
}


//===================================
// Set channel 
//===================================
void channel::setChannel(const char* fData, const char *fPed, int nPed)
{
  fname = fData;
  calibrated = false;

  // Reset arrays
  for (int i=0; i < rlength; i++) {
    pedestals[i] = rawval[i] = calibval[i] = 0.0;
    t[i] = i;
  }
 
  // Load pedestal
  cout << "loading pedestals" << endl;
  loadPedestals(fPed, nPed);

  // Open input data file 
  if (isbinary(fname)) {
    waves.open(fname,ios::binary);
    binaryfile = true;
  }
  else {
    waves.open(fname);
    binaryfile = false;
  }
}

//===================================
// Returns true if inputfile is binary 
//===================================
bool channel::isbinary(const char *f){
  string fs(f);
  bool bin = false;  
  if (fs.find(".dat",fs.size()-4) != string::npos) bin=true;    
  return bin;
}


//===================================
// Read a float from data
//===================================
void channel::readfloat(ifstream &is, float &pVal, bool binary){
  if (binary){
    static char ch[sizeof(float)];//ch is shared between calls of readfloat
    static int size = sizeof(float) ;//size is shared between calls of readlfoat
    is.read(ch, size);
    memcpy(&pVal,ch,size);
  }
  else 
    is>>pVal;
}


//===================================
// Load next event 
//===================================
void channel::loadNextEvent(){
  if (!waves.good()){cout<<"error opening file"<<endl; return;}

  for (int i=0;i<rlength;i++){
    readfloat(waves,rawval[i], binaryfile);
    //debug
    //    cout<<" record" << i<<endl;
    if(waves.eof()) {
      if (i!=0){
	cout<<"error loading next event: corrupted data, eof before end of record ";
	cout<<" at position "<<i <<" in the record"<<endl;
      }
      break;
    }
  }

  calibrate();
  calibrated=true;
}


//===================================
// Load event reading channel #index
//===================================
void channel::loadEvent(int index){
  ifstream tmpwaves;
  if (binaryfile)
    tmpwaves.open(fname,ios::binary);
  else
    tmpwaves.open(fname);    

  if (!tmpwaves.good()){cout<<"error opening file"<<endl; return;}
  //skip to the right position in the file
  float dummy;
  char ch[sizeof(float)];
  if (binaryfile) 
    for (int i=0;i<rlength*index;i++) 
      tmpwaves.read(ch, sizeof(float));
  else 
    for (int i=0;i<rlength*index;i++) 
      tmpwaves>>dummy;

  //load a waveform
  for (int i=0;i<rlength;i++){
    if(tmpwaves.eof()) continue;
    readfloat(tmpwaves,rawval[i],binaryfile);
  }

  //apply pedestal subtraction
  calibrate();
  calibrated = true;
  tmpwaves.close();
}


//===================================
// Plot a specific event 
//===================================
TGraph* channel::plot(int index) {
  loadEvent(index);
  return plot();
}


//===================================
// Plot calibration values 
//===================================
TGraph* channel::plot() {
  TGraph *g = new TGraph(rlength, t, calibval);
  g->Draw("APL");
  return g;
}


//===================================
// Plot next event
//===================================
TGraph* channel::plotNext() {
  loadNextEvent();
  return plot();
}


//===================================
// Plot raw values
//===================================
TGraph* channel::plotRaw(){
  TGraph *g=new TGraph(rlength,t,rawval);
  g->Draw("APL");
  return g;
}


//===================================
// Plot pedestal
//===================================
TGraph* channel::plotPedestal(){
  TGraph *g=new TGraph(rlength,t,pedestals);
  g->Draw("APL");
  return g;
}


//===================================
// Load pedestal
//===================================
void channel::loadPedestals(const char* pedestalfilename, int npedestals){
  for(int i=0;i<rlength;i++){pedestals[i]=0;}
  ifstream pedestalsfile;
  bool binarypedestals=isbinary(pedestalfilename);
  if (binarypedestals)
    pedestalsfile.open(pedestalfilename,ios::binary);    
  else
    pedestalsfile.open(pedestalfilename);
  float dummy;
  int nwaves=0;
  int iline=0;
  while(!pedestalsfile.eof() && (nwaves<npedestals || npedestals<=0)  ) {
    for (int i=0;i<rlength;i++){
      readfloat(pedestalsfile,dummy,binarypedestals);
      
      if(pedestalsfile.eof()) {
	if (i!=0){
	  cout<<"error loading pedestals: corrupted data, eof before end of record ";
	  cout<<" at position "<<i <<" in the record "<<" line "<<iline<<endl;
	}
	break;
      }
      iline++;
      if (dummy>pow(2,8*sizeof(float)) || dummy<0){
	cout<<"pedestals: error loading wave "<<nwaves<<", "<<i<<"-th value corrupted: "<<dummy<<endl;
	continue;
      }
      pedestals[i]+=dummy;
      
    }//end wave
    nwaves++;
    if (iline%rlength){cout<<"error: wave "<<nwaves<<" has "<<iline<<" records"<<endl;}
    iline=0;
    
  }  

  for(int i=0;i<rlength;i++){pedestals[i]/=nwaves;}
}


//===================================
// Subtract pedestal from raw data
//===================================
void channel::calibrate(){
  for(int i=0;i<rlength;i++)
    calibval[i]=rawval[i]-pedestals[i];
  calibrated=true;
}


//===================================
// Returns integral signal starting from sigmin
//===================================
double channel::integral(int sigmin){
  if (!calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  double sum=0;
  for (int i=sigmin;i<rlength;i++) sum+=calibval[i];
  return sum;
}


//===================================
// Other method of integration
//===================================
double channel::integral2(int bkgmax, int sigmax) 
{
  if (!calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  double sumbkg = 0, sumsig = 0;
  int bkgmin = 0;
  int sigmin = bkgmax;//attention hardcoded numbers
  for (int i = bkgmin; i < bkgmax; i++) sumbkg += calibval[i];
  for (int i = sigmin; i < sigmax; i++) sumsig += calibval[i];
  sumsig -= sumbkg * (sigmax-sigmin) / (bkgmax-bkgmin);
  //this division is dangerous becuase the 
  //result depends on the integration window even if the signal is localized
  //completely inside the window
  //sumsig = sumsig / (sigmax-sigmin);

  return sumsig;
}


//===================================
//===================================
double channel::maxval(int bkgmax){
  if (!calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  double bkg=0,max=-1e3;
  int sigmin=bkgmax; int sigmax=rlength;//attention hardcoded numbers
  for (int i=0;i<bkgmax;i++) bkg+=calibval[i];
  bkg=bkg/bkgmax;
  for (int i=sigmin;i<sigmax;i++) {if (calibval[i]>max) max=calibval[i];}
  /* note : names are terrible: fix them */
  return max-bkg;
}


//===================================
//===================================
double channel::smoothmax(int nchannels){
  if (!calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  double bkg=0,max=-1e3;
  int bkgmax=350; int sigmin=bkgmax; int sigmax=rlength;//attention hardcoded numbers
  for (int i=0;i<bkgmax;i++) bkg+=calibval[i];
  bkg=bkg/bkgmax;
  for (int i=sigmin;i<sigmax-nchannels;i++) {
    double average=0;
    for (int j=0;j<nchannels;j++){average+=calibval[i+j];}
    average=average/nchannels;
    if (average>max) max=average;}
  /* note : names are terrible: fix them */
  return max-bkg;
}


//===================================
//===================================
double channel::minval(int nchannels){
  if (!calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }

  double min=4096;
  int sigmin=0; int sigmax=rlength;//attention hardcoded numbers
  for (int i=sigmin;i<sigmax;i++) {
    double average=0;
    for (int j=0;j<nchannels;j++){average+=calibval[i+j];}
    average=average/nchannels;
    if (average<min) min=average;}
  return -1*min;
  
}


//===========================================
// Fit wave with a user-defined function
//===========================================
bool channel::fit(TF1 *fun, int min, int max)
{
  //fitGraph.Set(rlength);
  //for (int i = 0; i < rlength; i++) { 
  //   fitGraph.SetPoint(i, t[i], calibval[i]);
  //   fitGraph.SetPointError(i, 0.5, 15.0);
  //}
  
  TGraphErrors fitGraph(rlength, t, calibval);
  fitGraph.Fit(fun, "Q", "", min, max);
  bool ok = (gMinuit->fCstatu == "CONVERGED ");
  return ok;  
}


//===================================
// Compute time 
//===================================
float channel::time(bool negative, int bkgmax, int sigmax ){
   if (!calibrated){
     cout<<"error: integral on non calibrated channel"<<endl;
     return 0;
   }
   int bkgmin=0;
   int sigmin=bkgmax;

   // Compute mean value of the pedestal
   double pedestal=0;
   int sgn=negative?-1:+1;
   for (int i=bkgmin;i<bkgmax;i++) pedestal+=sgn*calibval[i];
   if (abs(bkgmax-bkgmin)>0) pedestal=pedestal/(bkgmax-bkgmin);
   //   cout<<"pedestal "<<pedestal<<endl;//debug

   // Find maximum of the signal
   double max=-1e4;
   int nchannels=3;
   for (int i=sigmin;i<sigmax-nchannels;i++) {
     //mooving average
    double average=0;
    for (int j=0;j<nchannels;j++){average+=sgn*calibval[i+j];}
    average=average/nchannels;
    if (average>max) max=average;
   }

   float minfrac=0.3; float maxfrac=0.7;
   //find time for which signal is >40% and 60% of maximum
   //   cout<<"max "<<max<<endl;
   float tmin=0,tmax=0; int iprev=sigmin;
   for (int i=sigmin;i<sigmax;i++) {
     if ((sgn*calibval[i]-pedestal)>minfrac*(max-pedestal)) {
       tmin=(t[i-1]+t[i])/2; break;
       iprev=i;
     }
   }
   for (int i=iprev;i<sigmax;i++) {
     if ((sgn*calibval[i]-pedestal)>maxfrac*(max-pedestal)) {
       tmax=(t[i-1]+t[i])/2; 
       break;
     }
   }
   float t0=tmin-(tmax-tmin)*minfrac/(maxfrac-minfrac);
   //cout<<"t0 "<<t0<<endl;//debug
   return t0;
 }
//=============================================================================
// End of channel class
//=============================================================================


// Define algorithms for pulse height computation
enum algorithm{integral, integral2, maxval, smoothmax, minval};
float getvalue(channel* ch, algorithm alg){
  float value=-999;
  switch (alg){
  case integral:
    value = ch->integral();
    break;
  case integral2:
    value = ch->integral2();
    break;
  case maxval:
    value = ch->maxval();
    break;
  case smoothmax:
    value = ch->smoothmax();
    break;
  case minval:
    value = ch->minval();
    break;
  //case fit: // obsolete...
  //  ch->fit();
  //  if (gMinuit->fCstatu == "CONVERGED ")
  //    value=ch->fitted_ph()/100; //normalize to fit the histogram limits
  //  else
  //    value=0;
  //  break;
  }
  return value;
}


//=============================================================================
// ChannelWithTrigger class
//=============================================================================
// - 20160112: Stefano 
//   fitSignal: now it fits from 0-600 (instead up to rlength to better fit bipolar signal)           
// - 20160115: Stefano 
//   fitSignal: now tau_r range is [1, 30] and tau_f [30, 500] 
class channelWithTrigger
{
public:
  // Constructor and destructor 
  channelWithTrigger(string fDataSig,
                     string fDataTrg, 
                     string fPedSig,
                     string fPedTrg,
                     int nPed = 0);
  channelWithTrigger() {}
  ~channelWithTrigger(); 
  
  // Read data and write tuple/histos in the outputfile
  void read(int maxEvents = -1);

  // Load next evnt for both channels
  void loadNextEvent() { 
    m_sigChannel.loadNextEvent();  
    m_trgChannel.loadNextEvent(); 
  }

  // Returns signal channel
  channel* sigChannel() { return &m_sigChannel; }

  // Returns trigger channel
  channel* trgChannel() { return &m_trgChannel; }

  // Returns signal function
  TF1* sigFun() { return m_sigFun; }
  
  // Returns trigger function
  TF1* trgFun() { return m_trgFun; }
  
  // Fit signal
  //bool fitSignal(int min = 0, int max = 600);
  bool fitSignal(int min = 25, int max = channel::rlength);

  // Fit trigger
  //bool fitTrigger(int min = 100, int max = 200); // laser default
  bool fitTrigger(int min = 25, int max = 65); //diode default

  // Plot channels and go to next event
  void plotNext(string ch = "signal");

private:
  channel m_trgChannel; // trigger ch.
  channel m_sigChannel; // signal ch. 
  TFile *m_outFile;     // output file
  TTree *m_tree;        // output tree
  TF1 *m_sigFun;
  TF1 *m_trgFun;
};

//======================================================
// Constructor
//======================================================
channelWithTrigger::channelWithTrigger(string fDataSig,
                                       string fDataTrg, 
                                       string fPedSig, 
                                       string fPedTrg, 
                                       int nPed)
{
  //m_sigFun = new TF1("sigFun", multi_sigwavef, 0, channel::rlength, 6);
  m_sigFun = new TF1("sigFun", sigwavef, 0, channel::rlength, 6);
  m_trgFun = new TF1("trgFun", trgwavef, 0, channel::rlength, 5);

  m_sigChannel.setChannel(fDataSig.c_str(), fPedSig.c_str(), nPed);
  m_trgChannel.setChannel(fDataTrg.c_str(), fPedTrg.c_str(), nPed);
}

//======================================================
// Destructor
//======================================================
channelWithTrigger::~channelWithTrigger() {}

//======================================================
// Fit signal (multi-waves)
//======================================================
//bool channelWithTrigger::fitSignal(int min, int max) 
//{
//  m_sigFun->SetRange(min, max);
//  //m_sigFun->SetParNames("t0", "tau_r", "tau_f", "ped", "norm", "bw");
// 
//  const int nFun = 4;
//
//  for (int i = 0; i < nFun; i++) {
//    m_sigFun->SetParLimits(i, 0, channel::rlength); //t0
//    m_sigFun->SetParLimits(nFun + i, 1., 30.); //tr
//    m_sigFun->SetParLimits(2*nFun + i, 30., 500.); //tf  
//
//    m_sigFun->SetParameter(i, 150 + 40*i);
//    m_sigFun->SetParameter(nFun + i, 5.);
//    m_sigFun->SetParameter(2*nFun + i, 100.);
//  }
// 
//  m_sigFun->SetParameter(3*nFun, 0);
//  m_sigFun->SetParameter(3*nFun + 1, 2e4);
//  m_sigFun->FixParameter(3*nFun + 2, 1); //bw
//
//
//  return m_sigChannel.fit(m_sigFun, min, max);
//}

bool channelWithTrigger::fitSignal(int min, int max) 
{
  m_sigFun->SetRange(min, max);
  m_sigFun->SetParNames("t0", "tau_r", "tau_f", "ped", "norm", "bw");
 
  m_sigFun->SetParLimits(0, min, max); //t0
  m_sigFun->SetParLimits(1, 0.25, 30.); //tr
  m_sigFun->SetParLimits(2, 10., 500.); //tf  
  m_sigFun->FixParameter(5, 1); //bw

  m_sigFun->SetParameters(min, 5, 100, 0, 2e4, 1); //(t0, tr, tf, ped, norm, bw)
  m_sigChannel.fit(m_sigFun, min, max);
  return m_sigChannel.fit(m_sigFun, min, m_sigFun->GetParameter("t0") + m_sigFun->GetParameter("tau_f"));
}

//======================================================
// Fit trigger
//======================================================
bool channelWithTrigger::fitTrigger(int min, int max) 
{
  m_trgFun->SetRange(min, max);
  m_trgFun->SetNpx(1000);
  
  // New settings
  m_trgFun->SetParNames("t0", "tau", "ped", "norm", "bw");
  m_trgFun->SetParLimits(0, min, max); //t0
  m_trgFun->SetParLimits(1, 0., 50.); //tau
  m_trgFun->SetParLimits(2, -50, 50.); //ped
  m_trgFun->SetParLimits(3, 1e2, 1.e4); //norm

  m_trgFun->FixParameter(4, 1); //bw

  m_trgFun->SetParameters(min, 1., 0., 1.4e3, 1); //(t0, tau, ped, norm, bw)

  // Old settings
  //m_trgFun->SetParNames("t0", "tau", "tstop", "ped", "norm", "bw");
  //m_trgFun->SetParLimits(0, 0, channel::rlength); //t0
  //m_trgFun->SetParLimits(1, 0., 50.); //tau
  //m_trgFun->SetParLimits(2, 0., 1023.); //tstop
 
  //m_trgFun->FixParameter(2, 1023); //tstop 
  //m_trgFun->FixParameter(5, 1); //bw
 
  //m_trgFun->SetParameters(100, 1, 1023, -20, 2e3, 1); //(t0, tau, tstop, ped, norm, bw)

  return m_trgChannel.fit(m_trgFun, min, max);
}

//======================================================
// Plot a channel and go to next event 
//======================================================
void channelWithTrigger::plotNext(string what)
{
  loadNextEvent(); 

  // Plot trigger 
  if (what == "trg" || what == "trigger" || what == "trig") {
    if (!fitTrigger()) cout << ">>> TRIGGER NOT CONVERGED!" << endl;
    m_trgFun->Print(); m_trgChannel.plot(); m_trgFun->Draw("PLx same"); 
  }
 
  // Plot signal
  else if (what == "sig" || what == "signal") {
    if (!fitSignal()) cout << ">>> SIGNAL NOT CONVERGED!" << endl;
    m_sigFun->Print(); m_sigChannel.plot(); m_sigFun->Draw("PLx same");
  }

  // Plot both channels
  else if (what == "both") {
    if (!fitTrigger()) cout << ">>> TRIGGER NOT CONVERGED!" << endl;
    if (!fitSignal()) cout << ">>> SIGNAL NOT CONVERGED!" << endl;
    m_trgFun->Print(); m_sigFun->Print(); 

    TGraph *t = m_trgChannel.plot(); TGraph *s = m_sigChannel.plot();
    t->Draw("APL"); m_trgFun->Draw("PLx same"); 
    s->SetLineColor(kBlue); m_sigFun->SetLineColor(kBlue); 
    s->Draw("PLx same"); m_sigFun->Draw("PLx same");
  }

  else { cout << "PlotNext has not such an option!" << endl; }
} 

//======================================================
// Read data and write histograms in the outputfile
//======================================================
void channelWithTrigger::read(int maxEvents)
{
  int nEvents = 0;
  float sigTime = 0.0, trgTime = 0.0, pulseHeight = 0.0;
  float integral = 0.0;
  float tauRise = 0.0, tauFall = 0.0, chi2 = 0.0;
  float tau = 0.0;

  if (maxEvents == -1) 
    maxEvents = 999999999;

  // Output file and tree
  m_outFile = new TFile("out.root", "recreate");
  m_tree = new TTree("tree", ""); 

  // Declare histograms
  //TCanvas *canv = new TCanvas("canvas", "");
  TH1F *hTrgTime = new TH1F("hTrgTime", "Trigger Time", 1024, 0, 1023);
  TH1F *hSigTime = new TH1F("hSigTime", "Signal_Time", 1024, 0, 1023);
  TH1F *hTimeDiff = new TH1F("hTimeDiff", "Time Difference", 200, -20., 100.);
  TH1F *hPh = new TH1F("hPh", "Pulse Height", 850, -50., 800.);
  TH1F *hChi2 = new TH1F("hChi2", "Fit chi2/nDoF", 200, 0.0, 200.0);
  TH1F *hTauRise = new TH1F("hTauRise", "Fit tau(rise)", 200, 0.0, 100.0);
  TH1F *hTauFall = new TH1F("hTauFall", "Fit tau(fall)", 200, 0.0, 500.0);
  
  // Declare branches 
  m_tree->Branch("trgTime", &trgTime);
  m_tree->Branch("sigTime", &sigTime);
  m_tree->Branch("pulseHeight", &pulseHeight);
  m_tree->Branch("integral", &integral);
  m_tree->Branch("chi2", &chi2);
  m_tree->Branch("tauRise", &tauRise);
  m_tree->Branch("tauFall", &tauFall);
  m_tree->Branch("tau", &tau);

  
  // Loop over events
  while (!(m_trgChannel.eof()) && (nEvents < maxEvents)) {
    // Load event
    loadNextEvent();
    nEvents++;
    if (nEvents%10000 == 0) cout << "." << flush;

    // Signal integral
    integral = m_sigChannel.integral2(200, channel::rlength) / 100.; //(bkgmax, sigmax)

    // Fit trigger
    trgTime = tau = 9999999.;
    if (fitTrigger()) {
      trgTime = m_trgFun->GetParameter("t0"); 
      tau = m_trgFun->GetParameter("tau");

      // Debugging 
      //m_trgChannel.plot();
      //m_trgFun->Print();
      //m_trgFun->Draw("PLx same");  
      //canv->Update();
      //canv->WaitPrimitive();
    }

    // Fit signal
    sigTime = pulseHeight = tauRise = tauFall = chi2 = 9999999.; 
    if (fitSignal()) {
      sigTime = m_sigFun->GetParameter("t0");
      pulseHeight = m_sigFun->GetParameter("norm") / 100.;
      chi2 = m_sigFun->GetChisquare() / float(m_sigFun->GetNDF());
      tauRise = m_sigFun->GetParameter("tau_r");
      tauFall = m_sigFun->GetParameter("tau_f");

      // Debugging 
      //float dt = sigTime - trgTime;
      //bool cleancut = (tauRise > 2. && tauRise < 48.) && (tauFall > 55. && tauFall < 490.) && chi2 < 1.;
      //bool tcut = (dt > 400. && dt < 450.);
      //bool phcut = true; //pulseHeight < 50.;

      //if (true) {
      ////if (cleancut && tcut && phcut) {
      //  cout << "> time: " << sigTime << ", ph: " << pulseHeight 
      //       << ", dt: " << dt << ", chi2/nDoF: " << chi2 << endl;
      //  m_sigChannel.plot();
      //  m_sigFun->Print();
      //  m_sigFun->Draw("PLx same");  
      //  canv->Update();
      //  canv->WaitPrimitive();
      //}
    }

    // Fill histograms
    hSigTime->Fill(sigTime);
    hTrgTime->Fill(trgTime);
    hTimeDiff->Fill(sigTime - trgTime);
    //    hPh->Fill(pulseHeight);
    hPh->Fill(integral);
    hChi2->Fill(chi2);
    hTauRise->Fill(tauRise);
    hTauFall->Fill(tauFall);
    m_tree->Fill();
  }
  
  cout << endl;

  // Write histograms
  hTrgTime->Write();
  hSigTime->Write();
  hTimeDiff->Write();
  hPh->Write();
  hChi2->Write();
  hTauRise->Write();
  hTauFall->Write();

  m_tree->Write();
  m_outFile->Close();
}
//=============================================================================
// End of channelWithTrigger class
//=============================================================================


