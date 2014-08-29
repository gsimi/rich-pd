/*
  Gabriele Simi, University of Padova, 2011

  A macro to extract a spectrum from a thscan of a SiPM,
  compute initial parameters for the fit function
  and fit.
  
  depends on spectrfitf.cc and thscan.C

*/


#include <algorithm>
#include "waveforms.C"
using namespace std;

TF1*
fitscan(char* fname="", 
	double fmin=0, double fmax=1){
  // fmin,fmax are the fit limits expressed as fraction of the histogram range

  bool debug=true;

  //get ph distribution
  TFile* histf=new TFile(fname);
  TH1F *h = (TH1F*)gDirectory->Get("ph");
  double bw=h->GetBinWidth(1);
  int dim=h->GetNbinsX();
  double xmin=h->GetXaxis()->GetXmin();
  double xmax=h->GetXaxis()->GetXmax();

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
  double gain=40,offset=7,noise=2.4,mu=0.28;
  if (npeaks>1){
    double x[npeaks];
    int ix[npeaks];
    for (int i=0;i<npeaks;i++) x[i]=pf.GetPositionX()[i];
    sort(x,x+npeaks-1);
    double pxmin=x[0],pxmax=x[npeaks-1];
    gain=(pxmax-pxmin)/(npeaks-1);
    offset=pxmin;
    mu=pf.GetPositionY()[1]/pf.GetPositionY()[0];//approximate
    if (npeaks>2) {
      pxmin=x[1]; pxmax=x[npeaks-1];
      gain=(pxmax-pxmin)/(npeaks-2);
      offset=pxmin-gain;
      mu=pf.GetPositionY()[2]/pf.GetPositionY()[0];//approximate
    }
    noise=gain/10.;
  }

  printf("npeaks = %d\n",npeaks);
  printf("pxmin = %2.2f\n",pxmin);
  printf("pxmax = %2.2f\n",pxmax);
  printf("gain = %2.2f\n",gain);
  printf("offset = %2.2f\n",offset);
  printf("mu = %2.2f\n",mu);
  printf("noise = %2.2f\n",noise);

  //now setup the fit function
  gROOT->ProcessLine(".L spectrfitf.cc+");
  int fitmodel=1; //1=fitf, 2=fitf_g2 double gaussian model
  int npar; TF1* f;

  switch (fitmodel){

  case 1:
     npar=9;  
     f=new TF1("spectrfit",fitf,xmin,xmax,npar);
     f->SetParNames( "mu","gain","gNoise","offset","iNoise","ln_norm","ct");//,"eff","bw");
     f->SetParameters(mu,  gain,  noise,  offset,   noise,   log(norm),0.03);
     f->FixParameter(7,1);//eff 
     f->FixParameter(8,bw);//bw
     break;

  case 2:
    npar=12; 
    f=new TF1("spectrfit",fitf_g2,xmin,xmax,npar);
    f->SetParNames(  "mu","gain","gNoise","offset","iNoise","ln_norm","ct","g2p","g2off","g2sigma");//,"eff","bw");    break;
    f->SetParameters(mu,gain     ,noise  ,offset       ,noise     ,log(norm)  ,0.03, 0.1, gain/5,2*noise);
    f->SetParLimits(7,0,0.9);//g2p
    //  f->FixParameter(7,0);
    f->SetParLimits(8,0,gain);//g2off
    //  f->FixParameter(8,0);
    f->SetParLimits(9,bw/5,50*bw);//g2sigma
    //  f->FixParameter(9,bw);
    f->FixParameter(10,1);//eff 
    f->FixParameter(11,bw);//bw
    break;
  }
  
  f->SetParLimits(0,0,20*mu); //mu
  f->SetParLimits(1,bw,10*gain); //gain
  f->SetParLimits(2,bw/5,50*bw); //gnoise
  f->SetParLimits(3,-100*bw,100*bw);//offset
  //  f->SetParLimits(4,bw/50,50*bw);//inoise
  //  f->FixParameter(4,0);//inoise
  f->SetParLimits(5,0.5*log(norm),1.5*log(norm));//norm
  f->SetParLimits(6,0,0.9);//ct
  f->FixParameter(6,0);//ct fixed
  
  cout<<"initial fit parameters:"<<endl;
  double pmin,pmax;
  for (int i=0;i<npar;i++){
    f->GetParLimits(i,pmin,pmax);
    printf("%s\t = %.3f \t [%.3f;%.3f]\n",f->GetParName(i),f->GetParameter(i),pmin,pmax);
  }

  //Fit
  cout<<"fitting spectrum in ["<<xmin<<";"<<xmax<<"] ..."<<endl;
  h->Fit(f,"","",xmin,xmax);
  cout<<"done"<<endl;

  //Draw
  gStyle->SetOptFit(111);
  h->SetMinimum(0.1);
  h->Draw();
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

double* mfitscan(vector<string> flist, char* fitted_data="fitpar.dat"){
  ofstream out; out.open(fitted_data);
  const int nfiles(flist.size());
  const int npars=10; //this is hardwired, should be possible to extract it from f
  double* parmatrix = new double[nfiles*npars*2];
  for (int ifile=0;ifile<nfiles; ifile++){
    TF1* f=fitscan(flist[ifile].c_str(),5./100,1);
    if (f==0){ cout<<"Error, fit for "<<flist[ifile].c_str()<<
	" not found, skipping"<<endl;
      continue;}
    if (ifile==0){
      out<<"ifile"; 
      for(int j=0;j<npars;j++){
	out<<":"<<f->GetParName(j)<<":e"<<f->GetParName(j);
      }
      out<<":dr:t"<<endl;

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
    out<<flist[ifile].c_str()<<endl;
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


