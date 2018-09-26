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
#include "TSpectrum.h"
#include "TError.h" 

using namespace::std;

typedef std::pair<float, float> Pair;

//===================================================
// Channel class
//===================================================
class channel {
public:

  // Number of records in the channel 
  const static int rlength = 1024;

  // Max number of allowed peaks
  const static int max_peaks = 10;

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
  TGraph* plotRaw();
  TGraph* plotPedestal();
  TH1F* plotHistogram();
  TGraph* plotNext();
  
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
  bool eof() { return m_waves.eof(); }

  // Returns if it's a binary file 
  bool isbinary(const char *fname);
  
  // Read a single float from stream 
  void readfloat(ifstream &is, float& val, bool binary);

  // Returns pedestals
  float* pedestals() { return m_pedestals; }
  //float Pedestals (unsigned idx) { return m_pedestals[idx]; }
  
  // Returns raw values
  float* rawval() { return m_rawval; }
  
  // Returns calibrated values
  float* calibval() { return m_calibval; }
  
  // Returns times
  float* t() { return m_t; }

  vector<Pair>& peaks() { return m_peaks; }

  // Find peaks
  vector<Pair>& findPeaks(int rebin = 1, 
		          float adc_min = 0., 
		          float sigma = 2, float thr = 0.05);

private:
  // Private functions
  void calibrate();

  // Data members
  float m_pedestals[rlength];
  float m_rawval[rlength];
  float m_calibval[rlength];
  float m_t[rlength];
  const char* m_fname;
  bool m_calibrated;
  ifstream m_waves;
  bool m_binaryfile;
  vector<Pair> m_peaks;
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
  // check if files exist
  if (gSystem->AccessPathName(fData)) {
    cout << "Data file " << fData << " does not exist!" << endl;
    return ;
  }

  if (gSystem->AccessPathName(fPed)) {
    cout << "Pedestal file " << fPed << " does not exist!" << endl;
    return;
  }

  m_fname = fData;
  m_calibrated = false;

  // Reset arrays
  for (int i=0; i < rlength; i++) {
    m_pedestals[i] = m_rawval[i] = m_calibval[i] = 0.0;
    m_t[i] = i;
  }
 
  // Load pedestal
  cout << "loading pedestals" << endl;
  loadPedestals(fPed, nPed);

  // Open input data file 
  if (isbinary(m_fname)) {
    m_waves.open(m_fname,ios::binary);
    m_binaryfile = true;
  }
  else {
    m_waves.open(m_fname);
    m_binaryfile = false;
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
  if (!m_waves.good()){cout<<"error opening file"<<endl; return;}

  for (int i=0;i<rlength;i++){
    readfloat(m_waves,m_rawval[i], m_binaryfile);
    //debug
    //    cout<<" record" << i<<endl;
    if(m_waves.eof()) {
      if (i!=0){
	cout<<"error loading next event: corrupted data, eof before end of record ";
	cout<<" at position "<<i <<" in the record"<<endl;
      }
      break;
    }
  }

  calibrate();
  m_calibrated=true;
}


//===================================
// Load event reading channel #index
//===================================
void channel::loadEvent(int index){
  ifstream tmpwaves;
  if (m_binaryfile)
    tmpwaves.open(m_fname,ios::binary);
  else
    tmpwaves.open(m_fname);    

  if (!tmpwaves.good()){cout<<"error opening file"<<endl; return;}
  //skip to the right position in the file
  float dummy;
  char ch[sizeof(float)];
  if (m_binaryfile) 
    for (int i=0;i<rlength*index;i++) 
      tmpwaves.read(ch, sizeof(float));
  else 
    for (int i=0;i<rlength*index;i++) 
      tmpwaves>>dummy;

  //load a waveform
  for (int i=0;i<rlength;i++){
    if(tmpwaves.eof()) continue;
    readfloat(tmpwaves,m_rawval[i],m_binaryfile);
  }

  //apply pedestal subtraction
  calibrate();
  m_calibrated = true;
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
  TGraph *g = new TGraph(rlength, m_t, m_calibval);
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
  TGraph *g=new TGraph(rlength,m_t,m_rawval);
  g->Draw("APL");
  return g;
}


//===================================
// Plot pedestal
//===================================
TGraph* channel::plotPedestal(){
  TGraph *g=new TGraph(rlength,m_t,m_pedestals);
  g->Draw("APL");
  return g;
}


//===================================
// Make histogram
//===================================
TH1F* channel::plotHistogram() {
  TH1F *h = new TH1F("h", "", rlength, m_t[0], m_t[rlength -1]);
  for (int i = 1; i <= rlength; i++) h->SetBinContent(i, m_t[i-1], m_calibval[i-1]);
  h->Draw();
  return h;
}


//===================================
// Load pedestal
//===================================
void channel::loadPedestals(const char* pedestalfilename, int npedestals){
  for(int i=0;i<rlength;i++){m_pedestals[i]=0;}
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
      m_pedestals[i]+=dummy;
      
    }//end wave
    nwaves++;
    if (iline%rlength){cout<<"error: wave "<<nwaves<<" has "<<iline<<" records"<<endl;}
    iline=0;
    
  }  

  for(int i=0;i<rlength;i++){m_pedestals[i]/=nwaves;}
}


//===================================
// Subtract pedestal from raw data
//===================================
void channel::calibrate(){
  for(int i=0;i<rlength;i++) {
    m_calibval[i]=m_rawval[i]-m_pedestals[i];
  }
  m_calibrated=true;
}


//===================================
// Returns integral signal starting from sigmin
//===================================
double channel::integral(int sigmin){
  if (!m_calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  double sum=0;
  for (int i=sigmin;i<rlength;i++) sum+=m_calibval[i];
  return sum;
}


//===================================
// Other method of integration
//===================================
double channel::integral2(int bkgmax, int sigmax) 
{
  if (!m_calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  double sumbkg = 0, sumsig = 0;
  int bkgmin = 0;
  int sigmin = bkgmax;//attention hardcoded numbers
  for (int i = bkgmin; i < bkgmax; i++) sumbkg += m_calibval[i];
  for (int i = sigmin; i < sigmax; i++) sumsig += m_calibval[i];
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
  if (!m_calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  double bkg=0,max=-1e3;
  int sigmin=bkgmax; int sigmax=rlength;//attention hardcoded numbers
  for (int i=0;i<bkgmax;i++) bkg+=m_calibval[i];
  bkg=bkg/bkgmax;
  for (int i=sigmin;i<sigmax;i++) {if (m_calibval[i]>max) max=m_calibval[i];}
  /* note : names are terrible: fix them */
  return max-bkg;
}


//===================================
//===================================
double channel::smoothmax(int nchannels){
  if (!m_calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  double bkg=0,max=-1e3;
  int bkgmax=350; int sigmin=bkgmax; int sigmax=rlength;//attention hardcoded numbers
  for (int i=0;i<bkgmax;i++) bkg+=m_calibval[i];
  bkg=bkg/bkgmax;
  for (int i=sigmin;i<sigmax-nchannels;i++) {
    double average=0;
    for (int j=0;j<nchannels;j++){average+=m_calibval[i+j];}
    average=average/nchannels;
    if (average>max) max=average;}
  /* note : names are terrible: fix them */
  return max-bkg;
}


//===================================
//===================================
double channel::minval(int nchannels){
  if (!m_calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }

  double min=4096;
  int sigmin=0; int sigmax=rlength;//attention hardcoded numbers
  for (int i=sigmin;i<sigmax;i++) {
    double average=0;
    for (int j=0;j<nchannels;j++){average+=m_calibval[i+j];}
    average=average/nchannels;
    if (average<min) min=average;}
  return -1*min;
  
}


//===========================================
// Fit wave with a user-defined function
//===========================================
bool channel::fit(TF1 *fun, int min, int max)
{
  TGraphErrors fitGraph(rlength, m_t, m_calibval);
  fitGraph.Fit(fun, "Q", "", min, max);
  bool ok = (gMinuit->fCstatu == "CONVERGED ");
  return ok;  
}


//===================================
// Compute time 
//===================================
float channel::time(bool negative, int bkgmax, int sigmax ){
  if (!m_calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  int bkgmin=0;
  int sigmin=bkgmax;

  // Compute mean value of the pedestal
  double pedestal=0;
  int sgn=negative?-1:+1;
  for (int i=bkgmin;i<bkgmax;i++) pedestal+=sgn*m_calibval[i];
  if (abs(bkgmax-bkgmin)>0) pedestal=pedestal/(bkgmax-bkgmin);
  //   cout<<"pedestal "<<pedestal<<endl;//debug

  // Find maximum of the signal
  double max=-1e4;
  int nchannels=3;
  for (int i=sigmin;i<sigmax-nchannels;i++) {
    //mooving average
    double average=0;
    for (int j=0;j<nchannels;j++){average+=sgn*m_calibval[i+j];}
    average=average/nchannels;
    if (average>max) max=average;
  }

  float minfrac=0.3; float maxfrac=0.7;
  //find time for which signal is >40% and 60% of maximum
  //   cout<<"max "<<max<<endl;
  float tmin=0,tmax=0; int iprev=sigmin;
  for (int i=sigmin;i<sigmax;i++) {
    if ((sgn*m_calibval[i]-pedestal)>minfrac*(max-pedestal)) {
      tmin=(m_t[i-1]+m_t[i])/2; break;
      iprev=i;
    }
  }
  for (int i=iprev;i<sigmax;i++) {
    if ((sgn*m_calibval[i]-pedestal)>maxfrac*(max-pedestal)) {
      tmax=(m_t[i-1]+m_t[i])/2; 
      break;
    }
  }
  float t0=tmin-(tmax-tmin)*minfrac/(maxfrac-minfrac);
  //cout<<"t0 "<<t0<<endl;//debug
  return t0;
}

//======================================================
// Peak finder
//======================================================
vector<Pair>& channel::findPeaks(int rebin, 
                                 float adc_min,
		                 float sigma, float thr) {
  // find peaks
  TSpectrum *s = new TSpectrum(max_peaks);
  TH1F *h = plotHistogram();
  h->Rebin(rebin);
  s->Search(h, sigma, "", thr);

  // fill peaks vector
  m_peaks.clear();
  m_peaks.reserve(s->GetNPeaks());

  //const float sigTime = 175.; // ns
  //const float tauFall = 40.; // ns
  //int i_sig = 0; // assume first peak if signal peak not found...
  //for (int i = 0; i < s->GetNPeaks(); i++) { 
  //  if (std::abs(s->GetPositionX()[i] - sigTime) < 2.) { 
  //    i_sig = i;
  //    break;
  //  }
  //}
  //float sig_time = s->GetPositionX()[i_sig];
  //float sig_adc = s->GetPositionY()[i_sig] / rebin;   

  //for (int i = 0; i < s->GetNPeaks(); i++) {
  //  s->GetPositionY()[i] /= rebin;

  //  float time = s->GetPositionX()[i];
  //  float adc = s->GetPositionY()[i];
  //  if (adc - sig_adc*std::exp(-time/tauFall) < adc_min) continue;

  //  m_peaks.push_back( {s->GetPositionX()[i], s->GetPositionY()[i]} );
  //}

  // old way of selecting peaks
  for (int i = 0; i < s->GetNPeaks(); i++) {
    s->GetPositionY()[i] /= rebin;
    if (s->GetPositionY()[i] < adc_min) continue;
    m_peaks.push_back( {s->GetPositionX()[i], s->GetPositionY()[i]} );
  }

  // sort peaks by increasing time
  auto pred = [&](const Pair &a, const Pair &b) { return a.first < b.first; }; 
  std::sort(m_peaks.begin(), m_peaks.end(), pred); 

  h->SetDirectory(gROOT);
  delete h;

  return m_peaks;
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
class channelWithTrigger {
public:
  // Constructor and destructor 
  channelWithTrigger(string fDataSig,
                     string fDataTrg, 
                     string fPedSig,
                     string fPedTrg,
		     string out_filename = "out.root",
                     int nPed = 1000);
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
  bool fitSignal(int min = 100, int max = 400);
  bool fitSignal_v2(int rebin, 
		    float adc_min, 
		    float sigma, float thr);

  // Fit trigger
  bool fitTrigger(int min = 60, int max = 85); 
  //bool fitTrigger(int min = 25, int max = 65); 

  // Plot channels and go to next event
  void plotNext(string ch = "signal");

private:
  channel m_trgChannel; // trigger ch.
  channel m_sigChannel; // signal ch. 
  TFile *m_outFile;     // output file
  TTree *m_tree;        // output tree
  TF1 *m_sigFun;
  TF1 *m_sigFun_v2;
  TF1 *m_trgFun;
  string m_out_filename;
};

//======================================================
// Constructor
//======================================================
channelWithTrigger::channelWithTrigger(string fDataSig,
                                       string fDataTrg, 
                                       string fPedSig, 
                                       string fPedTrg, 
				       string out_filename,
                                       int nPed) 
{
  //m_sigFun = new TF1("sigFun", multi_sigwavef, 0, channel::rlength, 6);
  m_sigFun = new TF1("sigFun", sigwavef, 0, channel::rlength, 6);
  m_trgFun = new TF1("trgFun", trgwavef, 0, channel::rlength, 5);
  m_sigFun_v2 = new TF1("sigFun_v2", multi_sigwavef, 0, channel::rlength, 
      n_glob_par + n_peak_par*channel::max_peaks);

  m_sigChannel.setChannel(fDataSig.c_str(), fPedSig.c_str(), nPed);
  m_trgChannel.setChannel(fDataTrg.c_str(), fPedTrg.c_str(), nPed);
 
  m_out_filename = out_filename;
}

//======================================================
// Destructor
//======================================================
channelWithTrigger::~channelWithTrigger() {}

//======================================================
// Fit signal
//======================================================
bool channelWithTrigger::fitSignal(int min, int max) 
{
  m_sigFun->SetRange(min, max);
  m_sigFun->SetParNames("t0", "tau_r", "tau_f", "ped", "norm", "bw");
 
  m_sigFun->SetParLimits(0, min, max); //t0
  m_sigFun->SetParLimits(1, 0.25, 30.); //tr
  m_sigFun->SetParLimits(2, 10., 500.); //tf  
  m_sigFun->FixParameter(5, 1); //bw

  m_sigFun->SetParameters(min, 5, 100, 0, 2e4, 1); //(t0, tr, tf, ped, norm, bw)

  return m_sigChannel.fit(m_sigFun, min, max);
}

//======================================================
// Fit signal v2
//======================================================
bool channelWithTrigger::fitSignal_v2(int rebin, 
	                    	      float adc_min, 
		                      float sigma, float thr)
{
  // find peaks
  auto peaks = m_sigChannel.findPeaks(rebin, adc_min, sigma, thr);
  int n_peaks = peaks.size();
  if (!n_peaks) return true;

  // fit each peak
  float min = 0.7*peaks.front().first;
  float max = 2.*peaks.back().first;

  m_sigFun_v2->SetRange(min, max);

  m_sigFun_v2->SetParName(0, "n_peaks"); //npeaks
  m_sigFun_v2->FixParameter(0, peaks.size());
  
  m_sigFun_v2->SetParName(1, "bw"); //bw
  m_sigFun_v2->FixParameter(1, 1); 

  m_sigFun_v2->SetParName(2, "ped"); //pedestal (common to all peaks)
  m_sigFun_v2->SetParLimits(2, -10., 10.); 
  m_sigFun_v2->SetParameter(2, 0.); 

  for (int i = 0; i < n_peaks; i++) {
    float t0 = peaks[i].first;
    m_sigFun_v2->SetParName(3 + i*n_peak_par, Form("t0_%i",i)); //t0
    m_sigFun_v2->ReleaseParameter(3 + i*n_peak_par); 
    m_sigFun_v2->SetParLimits(3 + i*n_peak_par, 0.80*t0, 1.10*t0); 
    m_sigFun_v2->SetParameter(3 + i*n_peak_par, 0.95*t0); 

    m_sigFun_v2->SetParName(4 + i*n_peak_par, Form("tr_%i",i)); //tr
    m_sigFun_v2->ReleaseParameter(4 + i*n_peak_par);
    m_sigFun_v2->SetParLimits(4 + i*n_peak_par, 0.25, 30.); 
    m_sigFun_v2->SetParameter(4 + i*n_peak_par, 5.); 

    m_sigFun_v2->SetParName(5 + i*n_peak_par, Form("tf_%i",i)); //tf
    m_sigFun_v2->ReleaseParameter(5 + i*n_peak_par); 
    m_sigFun_v2->SetParLimits(5 + i*n_peak_par, 2., 100.); 
    //m_sigFun_v2->SetParLimits(5 + i*n_peak_par, 5., 500.); // old working point
    m_sigFun_v2->SetParameter(5 + i*n_peak_par, 100.); 

    m_sigFun_v2->SetParName(6 + i*n_peak_par, Form("norm_%i",i)); //norm
    m_sigFun_v2->ReleaseParameter(6 + i*n_peak_par); 
    m_sigFun_v2->SetParLimits(6 + i*n_peak_par, 0., 1e5); 
    m_sigFun_v2->SetParameter(6 + i*n_peak_par, 2e4); 
  }

  // fix remaining peaks to zero
  for (int i = n_peaks; i < channel::max_peaks; i++) {
    m_sigFun_v2->SetParName(3 + i*n_peak_par, Form("t0_%i",i)); //t0
    m_sigFun_v2->FixParameter(3 + i*n_peak_par, 0.);

    m_sigFun_v2->SetParName(4 + i*n_peak_par, Form("tr_%i",i)); //tr
    m_sigFun_v2->FixParameter(4 + i*n_peak_par, 0.); 

    m_sigFun_v2->SetParName(5 + i*n_peak_par, Form("tf_%i",i)); //tf
    m_sigFun_v2->FixParameter(5 + i*n_peak_par, 0.); 
  
    m_sigFun_v2->SetParName(6 + i*n_peak_par, Form("norm_%i",i)); //norm
    m_sigFun_v2->FixParameter(6 + i*n_peak_par, 0.); 
  }

  return m_sigChannel.fit(m_sigFun_v2, min, max);
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
  // Suppress warnings and infos
  gErrorIgnoreLevel = kError;

  // TSpectrum configuration
  const int rebin = 10;
  const float adc_min = 20.;
  const float sigma = 1.;
  const float thr = 0.05;

  int nEvents = 0;
  int status = 0;
  float trgTime = 0.0;
  float integral = 0.0;
  float chi2 = 0.0;
  float tau = 0.0;
  int n_peaks = 0;
  float x_peaks[channel::max_peaks]= {0.};
  float y_peaks[channel::max_peaks] = {0.};
  float tau_rise[channel::max_peaks]= {0.};
  float tau_fall[channel::max_peaks] = {0.};

  if (maxEvents == -1) maxEvents = 999999999;

  // Output file and tree
  m_outFile = new TFile((m_out_filename + ".root").c_str(), "recreate");
  m_tree = new TTree("tree", ""); 

  // Declare histograms
  TCanvas *canv = new TCanvas("canvas", "");
  TH1F *hTrgTime = new TH1F("hTrgTime", "Trigger Time", 1024, 0, 1023);
  TH1F *hSigTime = new TH1F("hSigTime", "Signal_Time", 1024, 0, 1023);
  TH1F *hTimeDiff = new TH1F("hTimeDiff", "Time Difference", 200, -20., 100.);
  TH1F *hPh = new TH1F("hPh", "Pulse Height", 850, -50., 800.);
  TH1F *hChi2 = new TH1F("hChi2", "Fit chi2/nDoF", 200, 0.0, 200.0);
  TH1F *hTauRise = new TH1F("hTauRise", "Fit tau(rise)", 200, 0.0, 100.0);
  TH1F *hTauFall = new TH1F("hTauFall", "Fit tau(fall)", 200, 0.0, 500.0);
  
  // Declare branches 
  m_tree->Branch("trgTime", &trgTime);
  m_tree->Branch("integral", &integral);
  m_tree->Branch("chi2", &chi2);
  m_tree->Branch("status", &status);
  m_tree->Branch("tau", &tau);
  m_tree->Branch("n_peaks", &n_peaks, "n_peaks/I");
  m_tree->Branch("tau_rise", &tau_rise, "tau_rise[n_peaks]/F");
  m_tree->Branch("tau_fall", &tau_fall, "tau_fall[n_peaks]/F");
  m_tree->Branch("x_peaks", &x_peaks, "x_peaks[n_peaks]/F");
  m_tree->Branch("y_peaks", &y_peaks, "y_peaks[n_peaks]/F");

  
  // Start!
  time_t start, end;
  time(&start);
  
  cout << "Start processing!" << endl;
 
  // Loop over events
  while (!(m_trgChannel.eof()) && (nEvents < maxEvents)) {
    // Load event
    loadNextEvent();
    nEvents++;
    if (nEvents%1000 == 0) cout << nEvents << " processed..." << endl;

    //if (nEvents != 181 &&
    //    nEvents != 263 &&
    //    nEvents != 409 &&
    //    nEvents != 467 &&
    //    nEvents != 588 &&
    //    nEvents != 647 && 
    //    nEvents != 763 && 
    //    nEvents != 898 &&
    //    nEvents != 310 && 
    //    nEvents != 374 &&
    //    nEvents != 613 &&
    //    nEvents != 626 &&
    //    nEvents != 656 && 
    //    nEvents != 782 &&
    //    nEvents != 840 &&
    //    nEvents != 929
    //    ) continue;

    //if (nEvents != 2891 &&
    //	nEvents != 3718) continue;

    // Init to default
    integral = HUGE_VAL; 
    status = chi2 = HUGE_VAL;
    trgTime = tau = HUGE_VAL;
    n_peaks = 0;

    // Signal integral
    integral = m_sigChannel.integral2(200, channel::rlength) / 100.; //(bkgmax, sigmax)

    // Fit trigger
    if (fitTrigger()) {
      trgTime = m_trgFun->GetParameter("t0"); 
      tau = m_trgFun->GetParameter("tau");
    }

    // Fit signal!
    status = fitSignal_v2(rebin, adc_min, sigma, thr); 

    auto peaks = m_sigChannel.peaks();
    //for (auto p : peaks) cout << "x: " << p.first << ", y: " << p.second << endl;

    n_peaks = peaks.size();
    if (!n_peaks) continue;

    for (int i = 0; i < n_peaks; i++) {
      float t0 = m_sigFun_v2->GetParameter( Form("t0_%i",i) );
      float tr = m_sigFun_v2->GetParameter( Form("tr_%i",i) );
      float tf = m_sigFun_v2->GetParameter( Form("tf_%i",i) );
      float norm = m_sigFun_v2->GetParameter( Form("norm_%i",i) );
      float bw = m_sigFun_v2->GetParameter( "bw" ) ;
      float ped = m_sigFun_v2->GetParameter( "ped" );

      m_sigFun->SetParameters(t0, tr, tf, ped, norm, bw);
      tau_rise[i] = tr;
      tau_fall[i] = tf;
      x_peaks[i] = t0; 
      y_peaks[i] = m_sigFun->GetMaximum(); 
    }

    chi2 = m_sigFun_v2->GetChisquare() / float(m_sigFun_v2->GetNDF());

    // debug 
    //cout << endl;
    //cout << "******* Event " << nEvents << endl;
    //
    //if (n_peaks >= 2) {
    //for (int i = 0; i < peaks.size(); i++) {
    //  if ((x_peaks[i]-x_peaks[0] < 10.) && (y_peaks[i] > 25. && y_peaks[i] < 35.)) {

    //    cout << "npeaks: " << n_peaks << endl;
    //    cout << "chi2: " << chi2 << endl;
    //    cout << "x_peaks: " << x_peaks[i] << endl;
    //    cout << "y_peaks: " << y_peaks[i] << endl;
    //    cout << "tau_rise: " << tau_rise[i] << endl;
    //    cout << "tau_fall: " << tau_fall[i] << endl;
    //    
    //    auto h = m_sigChannel.plotHistogram();
    //    auto pm = new TPolyMarker();
    //    pm->SetMarkerStyle(23); pm->SetMarkerColor(kRed); pm->SetMarkerSize(2);
    //    for (int j = 0; j < peaks.size(); j++) {
    //      pm->SetNextPoint(peaks[j].first, peaks[j].second);
    //    }

    //    auto pm2 = new TPolyMarker();
    //    pm2->SetMarkerStyle(23); pm2->SetMarkerColor(kBlue); pm2->SetMarkerSize(2);
    //    for (int j = 0; j < n_peaks; j++) {
    //      pm2->SetNextPoint(x_peaks[j], y_peaks[j]);
    //    }

    //    pm->Draw("same");
    //    pm2->Draw("same");
    //    m_sigFun_v2->Draw("same");
    //    canv->Update();
    //    canv->WaitPrimitive();

    //    break;
    //  }
    //}

    // Fill histograms
    hTrgTime->Fill(trgTime);
    hPh->Fill(integral);
    hChi2->Fill(chi2);
    m_tree->Fill();
  }
  
  // Write histograms
  hTrgTime->Write();
  hPh->Write();
  hChi2->Write();

  m_tree->Write();
  m_outFile->Close();

  // Stop!
  time(&end);
  int timing = difftime(end, start);
  cout << "Analysis took: " 
       << int(timing / 60) << " min and " << timing % 60 << " sec" << std::endl;
  cout << "Done!" << std::endl;
}
//=============================================================================
// End of channelWithTrigger class
//=============================================================================


