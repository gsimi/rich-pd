/*
  Gabriele Simi, University of Padova, 2011

  A macro to extract a spectrum from waveforms acquired with a v1742,
  compute initial parameters for the fit function
  and fit it to the data.
  
*/


#include <algorithm>
#include <fstream> 
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPaveText.h"
#include "spectrfitf.cc"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "waveforms.C"
#include "spectrfitf.h"
#include "pmtfitf.h"

using namespace std;

/*
  function to perform the fit of the pulse height spectrum
  generated by the PMT. 
  It expects at least one peak in the distribution.
  
  fmin,fmax are the fraction of the histogram range to be used in the fit,
  they should be tuned to avoid large regions without data

  Initial values of the fit  parameters are estimated from 

    norm: integral of h between binmin and binmax
    gain=rms**2/mean
    npe=mean**2/rms**2
    offset=xpeak[0]
    gain1 is the gain of the first dynode, and is used to account 
    for conversions on the first dynode and back scattering.
    The gain1 depends on the HV with a power alpha~2/3
    gain1 = k* pow(HV*2.3/13,alpha) =k pow(HV*voltage_1/[sum_i voltage_i],alpha)
    The gain therefore depends on HV^(12*alpha)=HV^8
    gainSpread=gain/sqrt(gain1) because fluctuations at the first amplification stage dominate

  if npeaks==1 and forcesignal==false the fit is performed without the signal component
  if forcesiganl==1 the fit is performed with signal even if only one peak was found. 
  In the latter case the conversion on the first dynode is removed from the fit to simplify the 
  function.
  
  
*/
TF1*
fitscan(TH1F* h, double fmin=0, double fmax=1, double HV=950, bool forcesignal=false, int fitmodel=3){
  //1=fitf, 2=fitf_g2 double gaussian model, 3=pmtfit, 4=pmtfit2, 5=bessel function, 6=convolution of gaus with exp
  double bw=h->GetBinWidth(1);
  int dim=h->GetNbinsX();
  double xmin=h->GetXaxis()->GetXmin();
  double xmax=h->GetXaxis()->GetXmax();  
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

  //compute the initial estimates for gain,offset,npe,gainSpread

  //Search for peaks to compute initial values
  //note: *** should search within fit limits ***

  TSpectrum pf;
  pf.Search(h,4,"nobackground",1e-3); //used to configure the fit
  pf.Print("V");
  const int npeaks=pf.GetNPeaks(); 
  printf("npeaks = %d\n",npeaks);
  float xpeak[npeaks];  for (int i=0;i<npeaks;i++) xpeak[i]=pf.GetPositionX()[i];
  sort(xpeak,xpeak+npeaks-1);
  double offset=0.;
  if (npeaks>0) offset=xpeak[0];

  double rms=h->GetRMS(), mean=h->GetMean()-offset;
  double k=0.122; //estimated from fit of gain vs HV of pixel 20
  double alpha=0.75;//estimated from pixel 20, usually alpha=3/2
  double gain1 = k* pow(HV*2.3/13,alpha); //2.3 is the f1=voltage_1/voltage_3, 13is the sum_i (voltage_i/voltage_3)
  double gain=rms*rms/mean,gainSpread=gain/sqrt(gain1),npe=mean*mean/rms/rms;
  TF1 ff("ff","[0]*TMath::Gaus(x,[1],[2],1)");  
  ff.SetParameters(norm,offset,rms/5.);
  ff.SetParLimits(0,norm/2,2*norm);
  ff.FixParameter(1,offset);
  ff.SetParLimits(2,bw/5,rms);
  h->Fit(&ff,"","",offset-gain,offset+gain); 
  double inoise=ff.GetParameter(2);
  //  double inoise=2; //noise estimate in ADC counts

  
  printf("gain = %2.2f\n",gain);
  printf("offset = %2.2f\n",offset);
  printf("npe = %2.2f\n",npe);
  printf("gainSpread = %2.2f\n",gainSpread);
  printf("gain1 = %2.2f\n",gain1);

  //now setup the fit function
  int npar; TF1* f;
  gROOT->ProcessLine(".L ./analysis/spectrfitf.cc+");
  gROOT->ProcessLine(".L ./analysis/pmtfitf.cc+");

  switch (fitmodel){

  case 1:
    npar=9;  
    f=new TF1("spectrfit",fitf,xmin,xmax,npar);
    f->SetParNames( "npe","gain","gainSpread","offset","iNoise","ln_norm","ct");//,"eff","bw");
    f->SetParameters(npe,  gain,  gainSpread,  offset,    inoise,   log(norm),0.03);
    f->SetParLimits(2,bw/2,   gain); //ggainSpread
    f->FixParameter(6,0);//ct fixed
    f->FixParameter(7,1);//eff 
    f->FixParameter(8,bw);//bw
    break;

  case 2:
    npar=12; 
    f=new TF1("spectrfit",fitf_g2,xmin,xmax,npar);
    f->SetParNames(  "npe","gain","gainSpread","offset","iNoise","ln_norm","ct","g2p","g2off","g2sigma");//,"eff","bw");    break;
    f->SetParameters(npe,gain     ,gainSpread  ,offset       ,inoise     ,log(norm)  ,0.03, 0.1, gain/5,2/5*gainSpread);
    f->SetParLimits(2,bw/5,   4*gain); //ggainSpread
    f->FixParameter(6,0);//ct fixed
    f->SetParLimits(7,0,0.9);//g2p
    f->SetParLimits(8,0,gain);//g2off
    f->SetParLimits(9,bw/2,50*bw);//g2sigma
    f->FixParameter(10,1);//eff 
    f->FixParameter(11,bw);//bw
    break;
  
  default:
  case 3:
    npar=11;  
    f=new TF1("spectrfit",pmtpdf_gaus,xmin,xmax,npar);
    f->SetParNames( "npe","gain","gainSpread","offset","iNoise","norm","gain1","npe1",  "frac");//,"eff","bw");
    f->SetParameters(npe,  gain,  gainSpread,   offset,  inoise,   norm, gain1, npe, 0.05);//
    f->SetParLimits(2,bw/2,   gain); //gainSpread
    //releasing gain1 tends to give higher gain1 values than expected from calculations
    f->FixParameter(6,gain1);//gain1

    f->SetParLimits(7,0.0,4*npe+1);//npe1
    f->SetParLimits(8, 0,  1);//frac
    if (npeaks==1) { //distribution with only noise
      f->FixParameter(7,0.01);//npe1
      f->FixParameter(8,0);//frac
    }
    f->FixParameter(9,1);//eff 
    f->FixParameter(10,bw);//bw
    break;

  case 4:
    npar=11;  
    f=new TF1("spectrfit",pmtpdf_gaus2,xmin,xmax,npar);
    f->SetParNames( "npe","gain","gainSpread","offset","iNoise","ln_norm","ct","P","f");//,"eff","bw");
    f->SetParameters(npe,  gain,  gainSpread,  offset,   inoise,   log(norm),0.0);//,0.1,0.3 );
    f->SetParLimits(2,bw/2,   gain); //ggainSpread
    f->FixParameter(6,0);//ct fixed
    f->SetParLimits(7,0,1);
    f->SetParLimits(8,0,1);
    f->FixParameter(9,1);//eff 
    f->FixParameter(10,bw);//bw
    break;

  case 5:
    npar=11;
    f=new TF1("spectrfit",pmtpdf_besselg1adc,xmin,xmax,npar);
    f->SetParNames( "npe","gain","ndyn","offset","iNoise","norm","gain1","npe1",  "frac");//,"eff","bw");
    f->SetParameters(npe,  gain,     12, offset,  inoise,  norm,  gain1,  npe,   0.05);//
    f->FixParameter(2,12);//number of dynodes
    f->FixParameter(6,gain1);//gain1

    f->SetParLimits(7,0.0,4*npe+1);//npe1
    f->SetParLimits(8, 0,  1);//frac
    if (npeaks==1) { 
      f->FixParameter(7,0.01);//npe1
      f->FixParameter(8,0);//frac
    }
    f->FixParameter(9,1);//eff 
    f->FixParameter(10,bw);//bw
    break;

    //gauss * exp
  case 6:
    npar=13;
    f=new TF1("spectrfit",pmtpdf_gaus_exp,xmin,xmax,npar);
    f->SetParNames( "npe","gain","gainSpread","offset","iNoise","norm","gain1","npe1",  "frac","w",  "alpha");//,"eff","bw");
    f->SetParameters(npe,  gain,     12, offset,  inoise,  norm,  gain1,  npe,   0.05,    0.05, gain1/gain);//

    f->SetParLimits(2,bw/2,   gain); //ggainSpread

    f->FixParameter(6,gain1);//gain1

    f->SetParLimits(7,0.0,4*npe+1);//npe1
    f->SetParLimits(8, 0,  1);//frac
    if (npeaks==1) { 
      f->FixParameter(7,0.01);//npe1
      f->FixParameter(8,0);//frac
    }

    f->FixParameter(9,0.05);
    f->FixParameter(10,gain1/gain);
    f->FixParameter(11,1);//eff 
    f->FixParameter(12,bw);//bw
    break;



  }

  f->SetNpx(1000);
  f->SetParLimits(0,0,      3*npe); //npe
  f->SetParLimits(1,bw,     4*gain); //gain
  if (npeaks==1 && forcesignal==false) { 
    f->FixParameter(0,0.01);
    f->FixParameter(1,0); 
    f->FixParameter(2,bw);
  }
  f->SetParLimits(3,-100*bw,100*bw);//offset
  f->SetParLimits(4,bw/50,  gainSpread);//inoise
  f->SetParLimits(5,norm*0.1 , norm*1.3);//norm
   
  cout<<"initial fit parameters:"<<endl;
  double pmin,pmax;
  for (int i=0;i<npar;i++){
    f->GetParLimits(i,pmin,pmax);
    printf("%s\t = %.3f \t [%.3f;%.3f]\n",f->GetParName(i),f->GetParameter(i),pmin,pmax);
  }

  
  //Fit
  cout<<"fitting spectrum in ["<<xmin<<";"<<xmax<<"] ..."<<endl;
  TCanvas * c= new TCanvas("c","c",800,600);
  // c->Divide(2);
  // c->cd(1);
  // f->Draw();
  // c->cd(2);
  //  h->Fit(f,"EML","",xmin,xmax);

  h->Fit(f,"","",xmin,xmax); //simpler fit
  c->Update();
  /*
    perform  new fit with a different line color, releasing some parameters

  */
  if (fitmodel==6){
    f->ReleaseParameter(9);
    f->ReleaseParameter(10);
    f->SetParLimits(9,0,1);//w
    f->SetParLimits(10,0.1/gain,1./inoise);//alpha
  }else{
    f->SetLineColor(kBlue);
    f->ReleaseParameter(6);
    f->SetParLimits(6,0.3*gain1,4*gain1);//gain1
  }
  h->Fit(f,"EML","",xmin,xmax);
  


  cout<<"done"<<endl;

  //Draw
  gStyle->SetOptFit(1);
  h->SetMinimum(0.1);
  cout<<"drawing"<<endl;
  c->SetLogy(1);
  h->Draw("e");
  cout<<"Chi2/NDoF = "<<f->GetChisquare()<<"/"<<f->GetNDF()<<endl;
  return f;
}



double get_voltage(char *fname){
  string sfname(fname);
  size_t pos = sfname.find("px");
  size_t end = sfname.find("V/");
  size_t start = pos+5;
  string svoltage=sfname.substr(start,end-start);
  double voltage = atof(svoltage.c_str());
  if (voltage<500 || voltage > 1100){
    cout<<"invalid file name format, cannot find HV value, forcing HV=950"<<endl;
    voltage=950;
  }
  return voltage;
}

TF1*
fitscan(char* fname, double fmin=0, double fmax=1, bool forcesignal=false, int fitmodel=3){

  // fmin,fmax are the fit limits expressed as fraction of the histogram range

  //get ph distribution
  TFile *file=new TFile(fname);; if(file->IsOpen()==false) return nullptr;
  TCanvas* c1= (TCanvas*)gDirectory->Get("c1");
  TH1F *h = (TH1F*)c1->FindObject("ph");
  h->Rebin(2);
  double voltage=get_voltage(fname);
 
  return fitscan(h,fmin,fmax,voltage,forcesignal,fitmodel);
}


//Draw various contributions of pmtpdf_gaus 
TF1*
drawcontributions_gaus(TF1 *f){
  Double_t probarr[4];;
  Double_t probarr2[4];;
  TF1 *phe[4];
  TF1 *dyn[4];
  Double_t ci = 0,ci2=0;
    
  Double_t grms = f->GetParameter("gainSpread");
  Double_t grms1 = f->GetParameter("gainSpread")/f->GetParameter("gain1")*1.3;
  TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
  legend->SetTextFont(72);
  legend->SetTextSize(0.02);
  legend->AddEntry(f,"Global Fit","l");
  
  for (int i=0; i<4; i++) {
    if (i==0) {
      ci = exp(-f->GetParameter("npe"));
      ci2 = exp(-f->GetParameter("npe1"));
    } else {
      ci = probarr[i-1]*f->GetParameter("npe")/((float)i);
      ci2 = probarr2[i-1]*f->GetParameter("npe1")/((float)i);
    }
    probarr[i]=ci;
    probarr2[i]=ci2;

    
    Double_t sigma = sqrt(f->GetParameter("iNoise")*f->GetParameter("iNoise")+float(i)*grms*grms);
    Double_t sigma2 = sqrt(f->GetParameter("iNoise")*f->GetParameter("iNoise")+float(i)*grms1*grms1);   
    
    Double_t xm = f->GetParameter("offset") + float(i)*f->GetParameter("gain");
    Double_t xm2 = f->GetParameter("offset") + float(i)*f->GetParameter("gain")/f->GetParameter("gain1");

    if(i==0){
      //condtribution of pedestal
      TF1 *fnoise = (TF1*)f->Clone("fnoise");
      phe[0]=fnoise;
      f->Copy(*fnoise);
      fnoise->SetLineColor(kRed);
      fnoise->SetParameter("npe",0);
      fnoise->SetParameter("npe1",0);
      double frac=f->GetParameter("frac");
      double norm=f->GetParameter("norm");
      double prob=TMath::Poisson(0,f->GetParameter("npe"));
      double prob1=TMath::Poisson(0,f->GetParameter("npe1"));
      fnoise->SetParameter("norm",norm*((1-frac)*prob+frac*prob1));
      fnoise->Draw("same");

    }else{
    
      //condtribution of 1, 2 and 3 phe
      phe[i] = new TF1("phe","(TMath::Sign(1,x)+1)/2.*(1-[6])*[4]*[5]*[3]*(0.39894228/[0])*exp(-((x-[1])*(x-[1]))/(2*[0]*[0]))",-10,500);
      phe[i]->SetParameter(0,sigma);
      phe[i]->SetParameter(1,xm);
      phe[i]->SetParameter(3,f->GetParameter("norm"));
      phe[i]->SetParameter(4,f->GetParameter(10));//bw
      phe[i]->SetParameter(5,probarr[i]);
      phe[i]->SetParameter(6,f->GetParameter("frac"));
      phe[i]->SetNpx(1000);
      phe[i]->SetLineColor(i+3);

      phe[i]->Draw("SAME");
    
      //contribution of 1, 2 and 3 phe at 1st dynode
      dyn[i] = new TF1("dyn","(TMath::Sign(1,x)+1)/2.*[6]*[4]*[5]*[3]*(0.39894228/[0])*exp(-((x-[1])*(x-[1]))/(2*[0]*[0]))",-10,500);
      dyn[i]->SetParameter(0,sigma2);
      dyn[i]->SetParameter(1,xm2);
      dyn[i]->SetParameter(3,f->GetParameter("norm"));
      dyn[i]->SetParameter(4,f->GetParameter(10));//bw
      dyn[i]->SetParameter(5,probarr2[i]);
      dyn[i]->SetParameter(6,f->GetParameter("frac"));
      dyn[i]->SetNpx(1000);

      dyn[i]->SetLineColor(i+8);
      dyn[i]->Draw("SAME");
    }


   //legend->AddEntry(backFcn,"Background fit","l");
   //legend->AddEntry(signalFcn,"Signal fit","l");
  
  }
  legend->AddEntry(phe[0],"pedestal","l");
  legend->AddEntry(phe[1],"one phe","l");
  legend->AddEntry(phe[2],"2 phe","l");
  legend->AddEntry(phe[3],"3 phe","l");
  legend->AddEntry(dyn[1],"one phe at 1st dynode","l");
  legend->AddEntry(dyn[2],"2 phe at 1st dynode","l");


  legend->Draw();
  return f;
}



void
drawcontributions_bessel(TF1* f){
  double *par=f->GetParameters();
  TF1 *phe[4];
  TF1 *dyn[4];
  //noise
  int npar=11;
  double xmin=-50, xmax=500; //this should be extracted from the TH1F object
  phe[1]=new TF1("phe",pmtpdf_bessel,xmin,xmax,npar);
  
  phe[1]->SetParameters(par[6],par[1]/par[6],par[2]-1,par[3],par[4],par[5]);

}

/* 
   A macro to fit a list of threshold scans in a directory.
 */
/*
double*
mfitscan(char *dir="/data/superb/SiPM/RadHard/Feb2011LNL/"){
  vector<string> flist=GetListOfFiles(dir,"threshold-scan-");
  cout<<"found "<<flist.size()<<" scans to analyze"<<endl;
  return mfitscan(flist);
}
*/
/* fits a list of files passed in flist
   returns an array of double.
   The fitted values for file ifile of parametr jpar is:
   par[ifile,jpar] = pararray[ifile][2*jpar]
   epar[ifile,jpar] = pararray[ifile][2*jpar + 1]
*/

double** mfitscan(const char * listfilename="filelist.txt", const char* fitted_data="fitted_data.txt", bool forcesignal=false, int fitmodel=3 ){
  ifstream listfile(listfilename);
  vector<string> flist;
  string line;
  while(getline(listfile,line)){
    flist.push_back(line);
  }
  listfile.close();
  ofstream out; out.open(fitted_data);
  const int nfiles(flist.size());
  double** parmatrix = new double*[nfiles];
  for (int ifile=0;ifile<nfiles; ifile++){
    cout<<"\n"<<"file directory : "<<flist[ifile].c_str()<<endl;
    TF1* f=fitscan((char*)flist[ifile].c_str(),5./100,1, forcesignal, fitmodel);
    if (f==0){ cout<<"Error, fit for "<<flist[ifile].c_str()<<
	" not found, skipping"<<endl;
      continue;}
    double *fitpar=f->GetParameters();
    double *fiterr=f->GetParErrors();
    const int npars=f->GetNpar(); 
    parmatrix[ifile]=new double[npars*2];
    if (ifile==0){
      out<<"ifile:"; 
      for(int j=0;j<npars;j++){
	out<<f->GetParName(j)<<":"<<"e"<<f->GetParName(j)<<":";
      }
      out<<"pxNumber"<<":";
      out<<"voltage"<<endl;

    }
    out<<ifile<<"\t";
    for (int jpar=0;jpar<npars;jpar++){
      parmatrix[ifile][2*jpar]=fitpar[jpar];
      parmatrix[ifile][2*jpar+1]=fiterr[jpar];
      out<<fitpar[jpar]<<"\t"<<fiterr[jpar]<<"\t";
    }
    //    out<<darkrate(f,flist[ifile].c_str(),gate)<<"\t";
    //    out<<filetime(flist[ifile])<<"\t";
    string pxname;
    string volt;
    string path;
    size_t pos;
    size_t start,end;

    path = flist[ifile].c_str();
    pos = path.find("px");
    pxname = path.substr(pos+2,2);
    out<<pxname<<"\t";

    end = path.find("V/");
    start = pos+5;
    volt = path.substr(start,end-start);
    out<<volt<<endl;

    if (f!=0) delete f;
    
  }
  out.close();
  return parmatrix;
}

TTree * getTree(const char* fname="fitted_data.txt"){
  TTree *t=new TTree;
  t->ReadFile(fname);
  return t;
}
