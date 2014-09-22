/*
  Gabriele Simi, University of Padova, 2011

  A macro to extract a spectrum from a thscan of a SiPM,
  compute initial parameters for the fit function
  and fit.
  
  depends on spectrfitf.cc and thscan.C

*/


#include <algorithm>
#include <fstream> 
#include "waveforms.C"
#include "TCanvas.h"
using namespace std;

TF1*
fitscan(char* fname="", 
	double fmin=0, double fmax=1){

  // fmin,fmax are the fit limits expressed as fraction of the histogram range

  bool debug=true;

  //get ph distribution
  TFile* histf=new TFile(fname);
  //  TH1F *h = (TH1F*)gDirectory->Get("ph");
  TCanvas* c1= (TCanvas*)gDirectory->Get("c1");
  TH1F *h = (TH1F*)c1->FindObject("ph");
  return fitscan(h,fmin,fmax);
}

TF1*
fitscan(TH1F* h, double fmin=0, double fmax=1){
  double bw=h->GetBinWidth(1);
  int dim=h->GetNbinsX();
  double xmin=h->GetXaxis()->GetXmin();
  double xmax=h->GetXaxis()->GetXmax();  
  /* note that the range has to be tuned 
     in odre not to contain regions with 
     no data, otherwise the ML fit
     converges on a local minimum */
  //  double xmin=-15,xmax=200.; 
  for (int i=0;i<dim;i++){
    double ni=h->GetBinContent(i);
    if (ni>1e-5) {
      xmin=h->GetBinCenter(i);
      break;
    }
  }
  for (int i=dim;i>0;i--){
    double ni=h->GetBinContent(i);
    if (ni>1e-5) {
      xmax=h->GetBinCenter(i);
      break;
    }
  }

  //compute the initial estimate for norm
  double norm=0; 
  int binmin=int(fmin*dim); int binmax=int(fmax*dim)-1;
  if (binmax>dim-1) {cout<<"max out of range"<<endl; return 0;}
  norm=h->Integral(binmin,binmax);
  cout<<"norm "<<norm<<endl;

  //compute the initial estimates for gain,offset,mu,noise

  //Search for peaks to compute initial values
  //note: *** should search within fit limits ***

  TSpectrum pf;
  pf.Search(h,8,"nobackground",1e-2);
  pf.Print("V");
  const int npeaks=pf.GetNPeaks();
  double gain=25,offset=0,noise=12,mu=0.5;
  if (npeaks==1){
    double x[npeaks];
    int ix[npeaks];
    x[0]=pf.GetPositionX();
    offset=x[0];
    mu=0.01;
  }
  double pxmin,pxmax;
  if (npeaks>1){
    double x[npeaks];
    int ix[npeaks];
    for (int i=0;i<npeaks;i++) x[i]=pf.GetPositionX()[i];
    sort(x,x+npeaks-1);
    pxmin=x[0];
    pxmax=x[npeaks-1];
    gain=(pxmax-pxmin)/(npeaks-1);
    offset=pxmin;
    mu=pf.GetPositionY()[1]/pf.GetPositionY()[0]*20;//approximate
    if (npeaks>2) {
      pxmin=x[1]; pxmax=x[npeaks-1];
      gain=(pxmax-pxmin)/(npeaks-2);
      offset=pxmin-gain;
      mu=pf.GetPositionY()[2]/pf.GetPositionY()[0]*20;//approximate
    }
    noise=gain/2.;
  }

  printf("npeaks = %d\n",npeaks);
  if (npeaks>1){
    printf("pxmin = %2.2f\n",pxmin);
    printf("pxmax = %2.2f\n",pxmax);
  }
  printf("gain = %2.2f\n",gain);
  printf("offset = %2.2f\n",offset);
  printf("mu = %2.2f\n",mu);
  printf("noise = %2.2f\n",noise);

  //now setup the fit function
  gROOT->ProcessLine(".L ../analysis/spectrfitf.cc+");
  int fitmodel=3; //1=fitf, 2=fitf_g2 double gaussian model, 3=pmtfit, 4=pmtfit2
  int npar; TF1* f;

  switch (fitmodel){

  case 1:
     npar=9;  
     f=new TF1("spectrfit",fitf,xmin,xmax,npar);
     f->SetParNames( "mu","gain","gNoise","offset","iNoise","ln_norm","ct");//,"eff","bw");
     f->SetParameters(mu,  gain,  noise,  offset,    noise/10,   log(norm),0.03);
     f->FixParameter(6,0);//ct fixed
     f->FixParameter(7,1);//eff 
     f->FixParameter(8,bw);//bw
     break;

  case 2:
    npar=12; 
    f=new TF1("spectrfit",fitf_g2,xmin,xmax,npar);
    f->SetParNames(  "mu","gain","gNoise","offset","iNoise","ln_norm","ct","g2p","g2off","g2sigma");//,"eff","bw");    break;
    f->SetParameters(mu,gain     ,noise  ,offset       ,noise/10     ,log(norm)  ,0.03, 0.1, gain/5,2/5*noise);
    f->FixParameter(6,0);//ct fixed
    f->SetParLimits(7,0,0.9);//g2p
    //  f->FixParameter(7,0);
    f->SetParLimits(8,0,gain);//g2off
    //  f->FixParameter(8,0);
    f->SetParLimits(9,bw/5,50*bw);//g2sigma
    //  f->FixParameter(9,bw);
    f->FixParameter(10,1);//eff 
    f->FixParameter(11,bw);//bw
    break;
  
  default;
  case 3:
    npar=11;  
    f=new TF1("spectrfit",pmtfit,xmin,xmax,npar);
    f->SetParNames( "mu","gain","gNoise","offset","iNoise","norm","fgain1","mu1",  "P1dyn");//,"eff","bw");
    f->SetParameters(mu,  gain,  noise,   offset,  noise/10,   norm, 1./(2.3*3.2), mu, 0.05);//
    f->FixParameter(6,1./(2.3*3.2));//fgain1
    f->SetParLimits(7,0.0,10*mu+1);//mu1
    f->SetParLimits(8, 0,  1);//P1dyn
    if (npeaks==1) { 
      f->FixParameter(7,0.01);
      f->FixParameter(8,0);
    }
    f->FixParameter(9,1);//eff 
    f->FixParameter(10,bw);//bw
    break;

  case 4:
    npar=11;  
    f=new TF1("spectrfit",pmtfit2,xmin,xmax,npar);
    f->SetParNames( "mu","gain","gNoise","offset","iNoise","ln_norm","ct","P","f");//,"eff","bw");
    f->SetParameters(mu,  gain,  noise,  offset,   noise/10,   log(norm),0.0);//,0.1,0.3 );
    f->FixParameter(6,0);//ct fixed
    f->SetParLimits(7,0,1);
    f->SetParLimits(8,0,1);
    f->FixParameter(9,1);//eff 
    f->FixParameter(10,bw);//bw
    break;

  }

  f->SetNpx(1000);
  f->SetParLimits(0,0,      20*mu); //mu
  f->SetParLimits(1,bw,     10*gain); //gain
  f->SetParLimits(2,bw/5,   4*gain); //gnoise
  if (npeaks==1) { 
    f->FixParameter(0,0.01);
    f->FixParameter(1,0); 
    f->FixParameter(2,bw);
  }
  f->SetParLimits(3,-100*bw,100*bw);//offset
  f->SetParLimits(4,bw/50,  2*noise);//inoise
  f->SetParLimits(5,norm-20*sqrt(norm) , norm+20*sqrt(norm));//norm
   
  cout<<"initial fit parameters:"<<endl;
  double pmin,pmax;
  for (int i=0;i<npar;i++){
    f->GetParLimits(i,pmin,pmax);
    printf("%s\t = %.3f \t [%.3f;%.3f]\n",f->GetParName(i),f->GetParameter(i),pmin,pmax);
  }

  //Fit
  cout<<"fitting spectrum in ["<<xmin<<";"<<xmax<<"] ..."<<endl;
  h->Fit(f,"EML","",xmin,xmax);
  // f->ReleaseParameter(6);
  // f->SetParLimits(6,0.01,0.5);//fgain1
  // h->Fit(f,"EML","",xmin,xmax);

  cout<<"done"<<endl;

  //Draw
  gStyle->SetOptFit(111);
  h->SetMinimum(0.1);
  cout<<"drawing"<<endl;
  TCanvas * c= new TCanvas("c","c",800,600);
  c->SetLogy(1);
  h->Draw("e");
  return f;
}


/* 
   A macro to fit a list of threshold scans in a directory.
 */

double*
mfitscan(char *dir="/data/superb/SiPM/RadHard/Feb2011LNL/"){
  vector<string> flist=GetListOfFiles(dir,"threshold-scan-");
  cout<<"found "<<flist.size()<<" scans to analyze"<<endl;
  return mfitscan(flist);
}

/* fits a list of files passed in flist
   returns an array of double.
   The fitted values for file ifile of parametr jpar is:
   par[ifile,jpar] = pararray[ifile*npars + 2*jpar]
   epar[ifile,jpar] = pararray[ifile*npars + 2*jpar + 1]
*/

double* mfitscan(char * filelist="filelist.txt", char* fitted_data="fitted_data.txt"){
  ifstream list(filelist);
  vector<string> flist;
  string s;
  //int lines = std::count(std::istreambuf_iterator<char>( file ),std::istreambuf_iterator<char>(), '\n' ); 
  for(int i=0;i<64; i++){
    list>>s;
    flist.push_back("/home/lhcb/rich-pd/pmttest/"+s);
  }
  ofstream out; out.open(fitted_data);
  const int nfiles(flist.size());
  const int npars=11; //this is hardwired, should be possible to extract it from f
  double* parmatrix = new double[nfiles*npars*2];
  for (int ifile=0;ifile<nfiles; ifile++){
    cout<<"\n"<<"file directory : "<<flist[ifile].c_str()<<endl;
    TF1* f=fitscan(flist[ifile].c_str(),5./100,1);
    if (f==0){ cout<<"Error, fit for "<<flist[ifile].c_str()<<
	" not found, skipping"<<endl;
      continue;}
    if (ifile==0){
      out<<"ifile\t"; 
      for(int j=0;j<npars;j++){
	out<<f->GetParName(j)<<"\t"<<"e"<<f->GetParName(j)<<"\t";
      }
      out<<"pxNumber"<<endl;

    }
    double *fitpar=f->GetParameters();
    double *fiterr=f->GetParErrors();
    out<<ifile<<"\t";
    for (int jpar=0;jpar<npars;jpar++){
      parmatrix[ifile*npars+2*jpar]=fitpar[jpar];
      parmatrix[ifile*npars+2*jpar+1]=fiterr[jpar];
      out<<fitpar[jpar]<<"\t"<<fiterr[jpar]<<"\t";
    }
    //    out<<darkrate(f,flist[ifile].c_str(),gate)<<"\t";
    //    out<<filetime(flist[ifile])<<"\t";
    string pxname;
    string path;
    size_t pos;
    path = flist[ifile].c_str();
    pos = path.find("px");
    pxname = path.substr(pos+2,2);
    out<<pxname<<endl;
    if (f!=0) delete f;
  }
  return parmatrix;
}


/* helper functions to list the contents of a directory */

Bool_t IsItDirectory(const char *name, char * dirfile) const
{
  // Check if name is a directory.
  Long64_t size;
  Long_t id, flags, modtime;

  gSystem->ChangeDirectory(dirfile);
  flags = id = size = modtime = 0;
  gSystem->GetPathInfo(name, &id, &size, &flags, &modtime);
  Int_t isdir = (Int_t)flags & 2;

  return isdir ? kTRUE : kFALSE;
}
 
vector<string> GetListOfFiles(char* dirname, char* match=0){
  string pwd(gSystem->pwd());
  void *dir = gSystem->OpenDirectory(dirname);
  if (!dir) return 0;

  const char *file = 0;
  vector<string> contents;
  while ((file = gSystem->GetDirEntry(dir))) {
    if (!( IsItDirectory(file,dirname)) ) {
      string filestr(file);
      if ((match==0) || (filestr.find(match)!=string::npos) ){
	string filepath(dirname);
	filepath.append(file);
	contents.push_back(filepath);
      }
    }
  }
  gSystem->FreeDirectory(dir);
  gSystem->ChangeDirectory(pwd.c_str());
  return contents;
}


TH2F* datasheet(const char* rdata="datasheet"){
  TH2F *h = new TH2F("datasheet","gain",8,0,7, 8,0,7);
  h->SetXTitle("pixel 1 to 8");
  h->SetYTitle("pixel 1 to 57");
  ifstream file(rdata);
  double  ch;
  cout<<ch<<endl;
  //file>>ch;
  for(int i=1;i<65;i++){
    if (file.eof()) break;
    int x = int(i-1)%8+1;
    int y = int((i-0.001)/8)+1;
    // cout<<" x "<<x<<" y "<<y<<" val "<<value[ipar]<<endl;
    //getline(file,ch);
    file>>ch;
    h->SetBinContent(x,y,ch);
  }
  h->Draw("colz");
  return h;
}


TH2F* uniformity(const char* rdata="test.dat", int ipar=3, int normpx=61){
  char var[10];
  //sprintf(var,"%d",ipar);

  const int npar(24);
  if(ipar>=npar){break;}
  ifstream file(rdata);
  string  ch;
  string parname[npar];
  for(int i=0; i<npar;i++){
    file>>parname[i];
  }
  const char * title=parname[ipar].c_str();
  //sprintf(var,"gain",ipar);
  TH2F *h = new TH2F("uniformity",title,8,0,7, 8,0,7);
  h->SetXTitle("pixel 1 to 8");
  h->SetYTitle("pixel 1 to 57");
  double value[npar];
  double normvalue=100;
  while(!file.eof()){
    for (int i=0;i<npar;i++) {
      file>>value[i];
      if (file.eof()) break;
    }
    if (file.eof()) break;
    float pixel=value[23];
    int x = int(pixel-1)%8+1;
    int y = int((pixel-0.001)/8)+1;
    if(int(pixel+0.0001)==normpx){
      normvalue=value[ipar];
    }
    // cout<<" x "<<x<<" y "<<y<<" val "<<value[ipar]<<endl;
    h->SetBinContent(x,y,value[ipar]);
  }
  h->Scale(100./normvalue);
  h->SetStats(0);
  h->Draw("colz");
  TPaveText *pt=new TPaveText(0.7,0.85,0.98,0.98,"brNDC");
  char  legend1[100] ;
  char c1[20];
  sprintf(c1,"%d",normpx);
  strcpy(legend1,"normpx : ");
  strcat(legend1,c1);
  pt->AddText(legend1);
  char  legend2[100] ;
  char c2[100];
  sprintf(c2,"%f",normvalue);
  strcpy(legend2,"value for normpx : ");
  strcat(legend2,c2);
  pt->AddText(legend2);
  pt->Draw();
  return h;
}
