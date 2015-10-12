//1;3409;0c
#include "TGraph.h"
#include <iostream>
#include <fstream>
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "pmtwavef.cc"
#include "TMinuit.h"
using namespace::std;

class channel{
public:
  channel(const char* f="data/wave_1.txt", const char* p="pedestals/wave_1.txt",int npedestals=0);
  // ~channel(){if (pedestals!=0) delete pedestals; waves.close();}
  void loadPedestals(const char* pedestalFile="pedestals/wave_1.txt", int npedestals=0);
  void loadNextEvent();
  void loadEvent(int index=0);
  //functions return ownershi to the user
  TGraph* plot(int index);
  TGraph* plot();
  TGraph* plotNext();
  TGraph* plotRaw();
  TGraph* plotPedestal();
  //various ways to compute estimate height
  double integral(int sigmin=350);
  double integral2(int bkgmax=350);
  double maxval(int bkgmax=350);
  double minval(int nchannels=4);
  double smoothmax(int nchannels=4);
  TGraph* fit();
  float fitted_time(){return fitf->GetParameter("t0");}
  float fitted_ph(){return fitf->GetParameter("norm");}
  float time(bool negative=false, int bkgmax=0, int sigmmax=1024);
  bool eof(){return waves.eof();};
  bool isbinary(const char *fname);
  void readfloat(ifstream &is, float& val, bool binary);
private:
  //private functions
  void calibrate();
  //data members
  int rlength;
  float *pedestals;
  float *rawval;
  float *calibval;
  float *t;
  const char* fname;
  bool calibrated;
  ifstream waves;
  bool binaryfile;
  TF1* fitf;
};
channel::channel(const char* f, const char *p, int npedestals){
  rlength=1024;
  //initialize arrays
  pedestals=new float[rlength];
  rawval=new float[rlength];
  calibval=new float[rlength];
  t=new float[rlength];
  for (int i=0;i<rlength;i++){t[i]=i;}
  fname=f;
  calibrated=false;
  cout<<"loading pedestals"<<endl;
  loadPedestals(p,npedestals);
  if (isbinary(fname)) {
    waves.open(fname,ios::binary);
    binaryfile=true;}
  else {
    waves.open(fname);
    binaryfile=false;}
  fitf=new TF1("pmtwavef",pmtwavef,0,1024,6);
  fitf->SetParNames("t0","tau_r","tau_f","ped","norm","bw");
}

bool channel::isbinary(const char *f){
  string fs(f);
  bool bin = false;  
  if (fs.find(".dat",fs.size()-4) != string::npos) bin=true;    
  return bin;
}

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
  calibrated=true;
  tmpwaves.close();
}

TGraph* channel::plot(int index){
  loadEvent(index);
  return plot();
}

TGraph* channel::plot(){
  TGraph *g=new TGraph(rlength,t,calibval);
  g->Draw("APL");
  return g;
}

TGraph* channel::plotNext(){
  loadNextEvent();
  return plot();
}

TGraph* channel::plotRaw(){
  TGraph *g=new TGraph(rlength,t,rawval);
  g->Draw("APL");
  return g;
}

TGraph* channel::plotPedestal(){
  TGraph *g=new TGraph(rlength,t,pedestals);
  g->Draw("APL");
  return g;
}



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
      if (dummy>pow(2,12) || dummy<0){
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

void channel::calibrate(){
  for(int i=0;i<rlength;i++)
    calibval[i]=rawval[i]-pedestals[i];
  calibrated=true;
}

double channel::integral(int sigmin){
  if (!calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  double sum=0;
  for (int i=sigmin;i<rlength;i++) sum+=calibval[i];
  return sum;
}

double channel::integral2(int bkgmax){
  if (!calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  double sumbkg=0,sumsig=0;
  int bkgmin=0;
  int sigmin=bkgmax, sigmax=800;//attention hardcoded numbers
  for (int i=bkgmin;i<bkgmax;i++) sumbkg+=calibval[i];
  for (int i=sigmin;i<sigmax;i++) sumsig+=calibval[i];
  sumsig-=sumbkg*(sigmax-sigmin)/(bkgmax-bkgmin);
  //this division is dangerous becuase the 
  //result depends on the integration window even if the signal is localized
  //completely inside the window
  sumsig=sumsig/(sigmax-sigmin); 
  return sumsig;

}


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

TGraph* channel::fit(){
  TGraph *g=new TGraph(rlength,t,calibval);
  fitf->SetParameters(420,5,100,0,2e4,1);
  fitf->FixParameter(5,1);//bw
  g->Fit(fitf,"Q","",0,1024);
  return g;  
}

float channel::time(bool negative, int bkgmax, int sigmax ){
   if (!calibrated){
     cout<<"error: integral on non calibrated channel"<<endl;
     return 0;
   }
   int bkgmin=0;
   int sigmin=bkgmax;
   //compute pedestal
   double pedestal=0;
   int sgn=negative?-1:+1;
   for (int i=bkgmin;i<bkgmax;i++) pedestal+=sgn*calibval[i];
   if (abs(bkgmax-bkgmin)>0) pedestal=pedestal/(bkgmax-bkgmin);
   //   cout<<"pedestal "<<pedestal<<endl;//debug

   //find maximum
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

enum algorithm{integral, integral2,maxval, smoothmax, minval, fit};
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
  case fit:
    ch->fit();
    if (gMinuit->fCstatu == "CONVERGED ")
      value=ch->fitted_ph()/100; //normalize to fit the histogram limits
    else
      value=0;
    break;
  }
  return value;
}

TH1F* pulseheight(const char * fdata="data/wave.01_1.txt", 
		  const char* fpedestals="pedestals/pedestal.01_1.txt", algorithm alg=integral2, int tmin=0, int tmax=1024, bool negative=false){
  channel *ch=new channel(fdata,fpedestals);
  TH1F *h = new TH1F("ph","pulse_height",2100,-52,2048);
  TH1F *t= new TH1F("th","pulse_time",1024, 0,1023);
  int ievent=0;
  while (!(ch->eof())){
    ch->loadNextEvent();
    ievent++;
    if (ievent%10000==0) cout<<"."<<flush;
    float value=getvalue(ch,alg);
    float time=alg!=fit? ch->time(negative) : ch->fitted_time() ;
    if (gMinuit->fCstatu == "CONVERGED ")
      t->Fill(time);{
      if (tmin<time && time<tmax)
	h->Fill(value,1);
    }
  }
  cout<<endl;
  //  TCanvas *c1=new TCanvas("c1","c1",800,600);
  //h->Draw();
  t->Draw();
  //  c1->Update();
  return h;
}




void waveforms(){
  TH1F* h=pulseheight("data24/2/950V tune 0%.txt","data24/2/pedestal 850V.txt");
  TFile f("ph.root","update");
  h->Write("ph");
  f.Close();
}


TH1F** phsuperimposition(){
  const int nh=64;
  TH1F** harr=new TH1F*[nh];
  for (int i=0; i<nh; i++){ 
    char fn[256]; 
    sprintf(fn,"/home/lhcb/rich-pd/pmttest/workdir/data/dataset2/px%.2d/950V/ph.root",i+1); 
    //    cout<<fn<<endl; 
    TFile *f=new TFile(fn); 
    TCanvas* c1=(TCanvas*)f->FindObjectAny("c1");
    if (i==0) c1->Draw();
    TH1F* h=(TH1F*)c1->FindObject("ph"); 
    harr[i]=(TH1F*)h->Clone();
  }
  TCanvas *c=new TCanvas("canv","canv",800,600);
  c->Draw();
  for (int i=0;i<nh;i++){
    const int ngroup=5;
    harr[i]->Rebin(ngroup);
    //    harr[i]->Smooth();
    harr[i]->SetMaximum(200*ngroup);
    harr[i]->SetLineColor(i);
    harr[i]->GetXaxis()->SetRange(int(10./ngroup),int(250./ngroup));
    harr[i]->SetStats(kFALSE);
    if (harr[i]->GetMean() > 5) harr[i]->Draw(i==0?"l":"l same");
    else cout<<"excluded pix "<<i+1<<endl;
  }
  return harr;

}

#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
void stage(){
  TCanvas *c=new TCanvas("c","c",800,600);

  channel ch("data/scint-coincid123-ph4-th60mV.dat","data/scint-pedestal.dat");
  ch.plot(25);
  //c->SaveAs("stage_pulse_shape.eps");
  //c->SaveAs("stage_pulse_shape.png");

  c->SetLogy(1);

  TH1F* h = pulseheight("data/scint-coincid123-ph4-th60mV.dat","data/scint-pedestal.dat",minval, 165,300);
  TH1F* th=(TH1F*)gROOT->FindObjectAny("th");
  th->GetXaxis()->SetTitle("time (1 unit = 0.2 ns)");
  //if (th!=nullptr){
    th->Fit("gaus","","",100,400);
    th->Draw("e");
    //c->SaveAs("stage_pulse_time.eps");
    //c->SaveAs("stage_pulse_time.png");
    //}
  
  h->GetXaxis()->SetTitle("pulse height (1 unit = 0.24 mV)");
  h->Rebin(42);
  float bw=h->GetBinWidth(0);
  TF1 * f=new TF1("f","[0]*[7]*([3]*TMath::Gaus(x,[1],[2],1) + [5]*exp(-x/[4])/[4]+(1-[3]-[5])*exp(-x/[6])/[6] )",0,2000);
  f->SetParNames("norm","mu","sigma","frac","alpha1","frac1","alpha2","bw");
  f->SetParameters(h->GetEntries(),500,300,0.75,50,0.2,1000,bw);
  f->SetParLimits(0,0,1e5);
  f->SetParLimits(1,0,1024);
  f->SetParLimits(2,0,1024);
  f->SetParLimits(3,0,1);
  f->SetParLimits(4,0,100);
  f->SetParLimits(5,0,1);
  f->SetParLimits(6,100,1e4);
  f->FixParameter(7,bw);

  gStyle->SetOptFit(1);
  c->SetLogy(0);
  h->Fit(f,"l","",0,1500);
  h->Draw("e");
  //c->SaveAs("stage_pulse_height.eps");
  //c->SaveAs("stage_pulse_height.png");

}

// my version 

void stageLeo(){
  TCanvas *c=new TCanvas("c","c",800,600);
  //channel ch("data/scint4-3-custom-2907-10000-2.dat","data/scint4-3-custom-2907-10000-pedestal-2.dat");
  channel ch("data/scint4-3-custom-3007-100000-2.dat","data/scint4-3-custom-3007-100000-pedestal-1.dat");
  ch.plot(5);
  //c->SaveAs("data/stage_pulse_shape_leo.eps");
  //c->SaveAs("data/stage_pulse_shape_leo.png");

  //c->SetLogy(1);
  
  //TH1F* h = pulseheight("data/scint4-3-custom-2907-10000-2.dat","data/scint4-3-custom-2907-10000-pedestal-2.dat",minval,165,300);
  TH1F* h = pulseheight("data/scint4-3-custom-3007-100000-2.dat","data/scint4-3-custom-3007-100000-pedestal-1.dat",minval,165,300);
  TH1F* th=(TH1F*)gROOT->FindObjectAny("th");
  th->GetXaxis()->SetTitle("time (1 unit = 0.2 ns)");
 //  TF1 * g = new TF1("g","log([0]*[1]/[2])-([3]+x*[1])/[2]",0,1000);
//   g->SetParNames("N0","delta_t_bin","mu_tau","t0");
//   g->SetParameters(th->GetEntries(),1,11,0);
//   g->SetParLimits(0,0,1e05);
//   g->SetParLimits(1,0.5,1.5);
//   g->SetParLimits(2,5,15);
//   g->SetParLimits(3,0,1e-05);
// gStyle->SetOptFit(1);
//   c->SetLogy(0);

  //if (th!=nullptr){
  th->Fit("gaus","","",100,400);
  //th->Fit(g,"l","",0,1000);
  th->Draw("e");
    //c->SaveAs("data/stage_pulse_time_leo.eps");
    //c->SaveAs("data/stage_pulse_time_leo.png");
    //}
      
  h->GetXaxis()->SetTitle("pulse height (1 unit = 0.24 mV)");
  h->Rebin(42);
  float bw=h->GetBinWidth(0);
  TF1 * f=new TF1("f","[0]*[7]*([3]*TMath::Gaus(x,[1],[2],1) + [5]*exp(-x/[4])/[4]+(1-[3]-[5])*exp(-x/[6])/[6] )",0,2000);
  f->SetParNames("norm","mu","sigma","frac","alpha1","frac1","alpha2","bw");
  f->SetParameters(h->GetEntries(),400,200,0.8,50,0.15,1000,bw);
  f->SetParLimits(0,0,1e5);
  f->SetParLimits(1,0,1024);
  f->SetParLimits(2,0,1024);
  f->SetParLimits(3,0,1);
  f->SetParLimits(4,0,100);
  f->SetParLimits(5,0,1);
  f->SetParLimits(6,100,1e4);
  f->FixParameter(7,bw);

  gStyle->SetOptFit(1);
  c->SetLogy(0);
  h->Fit(f,"l","",0,1500);
  h->Draw("e");
  //c->SaveAs("data/stage_pulse_height_leo.eps");
  //c->SaveAs("data/stage_pulse_height_leo.png");
  
}
