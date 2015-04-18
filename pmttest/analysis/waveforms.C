#include "TGraph.h"
#include <iostream>
#include <fstream>
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
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
  double integral();
  double integral2();
  double maxval();
  double smoothmax(int nchannels=4);
  float time();
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
  loadPedestals(p,npedestals);
  if (isbinary(fname)) {
    waves.open(fname,ios::binary);
    binaryfile=true;}
  else {
    waves.open(fname);
    binaryfile=false;}
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

double channel::integral(){
  if (!calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  double sum=0;
  for (int i=0;i<rlength;i++) sum+=calibval[i];
  return sum;
}

double channel::integral2(){
  if (!calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  double sumbkg=0,sumsig=0;
  int bkgmin=0, bkgmax=350; 
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


double channel::maxval(){
  if (!calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  double bkg=0,max=-1e3;
  int bkgmax=350; int sigmin=bkgmax; int sigmax=rlength;//attention hardcoded numbers
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

 float channel::time(){
   if (!calibrated){
     cout<<"error: integral on non calibrated channel"<<endl;
     return 0;
   }
   int bkgmin=0, bkgmax=350; 
   int sigmin=bkgmax, sigmax=800;//attention hardcoded numbers
   //compute pedestal
   double pedestal=0;
   for (int i=bkgmin;i<bkgmax;i++) pedestal+=calibval[i];
   pedestal=pedestal/(bkgmax-bkgmin);
   cout<<"pedestal "<<pedestal<<endl;//debug
   //find maximum
   double max=-1e3;
   for (int i=sigmin;i<sigmax;i++) {if (calibval[i]>max) max=calibval[i];}
   //find time for which signal is >50% of maximum
   cout<<"max "<<max<<endl;
   float t0=0;
   for (int i=sigmin;i<sigmax;i++) {
     if ((calibval[i]-pedestal)>0.5*(max-pedestal)) {
       t0=t[i]; 
       break;
     }
   }
   cout<<"t0 "<<t0<<endl;
   return t0;
 }

TH1F* pulseheight(const char * fdata="data/wave.01_1.txt", 
		  const char* fpedestals="pedestals/pedestal.01_1.txt"){
  //  channel *ch=new channel(fdata,fpedestals,10000);
  channel *ch=new channel(fdata,fpedestals);
  TH1F *h= new TH1F("ph","pulse_height",550,-50,500);
  TH1F *t= new TH1F("t","pulse_time",550, 300,850);
  int ievent=0;
  while (!(ch->eof())){
    ch->loadNextEvent();
    ievent++;
    if (ievent%10000==0) cout<<"."<<flush;
    float integr=ch->integral2();
    //    float integr=ch->smoothmax(3);
    h->Fill(integr,1);
    if (integr>4)
      t->Fill(ch->time());
  }
  cout<<endl;
  TCanvas *c1=new TCanvas("c1","c1",800,600);
  h->Draw();
  c1->Update();
  return h;
}




void waveforms(){
  TH1F* h=pulseheight("data24/2/950V tune 0%.txt","data24/2/pedestal 850V.txt");
  TFile f("ph.root","update");
  h->Write("ph1");
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

