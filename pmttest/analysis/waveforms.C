#include "TGraph.h"
#include <iostream>
#include <fstream>
#include "TH1F.h"
#include "TFile.h"
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

  bool eof(){return waves.eof();};
private:
  //private functions
  void calibrate();
  //data members
  int rlength;
  double *pedestals;
  double *rawval;
  double *calibval;
  double *t;
  const char* fname;
  bool calibrated;
  ifstream waves;
};
channel::channel(const char* f, const char *p, int npedestals){
  rlength=1024;
  //initialize arrays
  pedestals=new double[rlength];
  rawval=new double[rlength];
  calibval=new double[rlength];
  t=new double[rlength];
  for (int i=0;i<rlength;i++){t[i]=i;}
  fname=f;
  calibrated=false;
  loadPedestals(p,npedestals);
  waves.open(fname);
}

void channel::loadNextEvent(){
  if (!waves.good()){cout<<"error opening file"<<endl; return;}
  for (int i=0;i<rlength;i++){
    waves>>rawval[i];
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
  ifstream tmpwaves(fname);
  if (!tmpwaves.good()){cout<<"error opening file"<<endl; return;}
  double dummy;
  for (int i=0;i<rlength*index;i++){tmpwaves>>dummy;}
  
  for (int i=0;i<rlength;i++){
    if(tmpwaves.eof()) continue;
    tmpwaves>>rawval[i];
  }
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
  ifstream pedestalsfile(pedestalfilename);
  double dummy;
  int nwaves=0;
  int iline=0;
  while(!pedestalsfile.eof() && (nwaves<npedestals || npedestals<=0)  ) {
    for (int i=0;i<rlength;i++){
      pedestalsfile>>dummy;
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
  int bkgmax=350; int sigmin=bkgmax; int sigmax=800;
  for (int i=0;i<bkgmax;i++) sumbkg+=calibval[i];
  for (int i=sigmin;i<sigmax;i++) sumsig+=calibval[i];
  sumsig-=sumbkg*(sigmax-sigmin)/bkgmax;
  sumsig=sumsig/(sigmax-sigmin);
  return sumsig;

}


double channel::maxval(){
  if (!calibrated){
    cout<<"error: integral on non calibrated channel"<<endl;
    return 0;
  }
  double bkg=0,max=0;
  int bkgmax=350; int sigmin=bkgmax; int sigmax=rlength;
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
  double bkg=0,max=0;
  int bkgmax=350; int sigmin=bkgmax; int sigmax=rlength;
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





TH1F* pulseheight(const char * fdata="data/wave.01_1.txt", 
		  const char* fpedestals="pedestals/pedestal.01_1.txt"){
  channel *ch=new channel(fdata,fpedestals,10000);
  //TH1F *h= new TH1F("ph","pulse_height",100,-100,30000);
  TH1F *h= new TH1F("ph","pulse_height",250,0,200);
  int ievent=0;
  while (!(ch->eof())){
    ch->loadNextEvent();
    ievent++;
    if (ievent%10000==0) cout<<"."<<flush;
    double integr=ch->integral2();
    h->Fill(integr,1);
  }
  cout<<endl;
  h->Draw();
  return h;
}




void waveforms(){
  TH1F* h=pulseheight("data24/2/950V tune 0%.txt","data24/2/pedestal 850V.txt");
  TFile f("ph.root","update");
  h->Write("ph1");
  f.Close();
}
