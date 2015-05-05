/* 
   Author List:

   Gabriele Simi, University of Padova, 2011

   A class to extract the spectrum from a threshold scan.

   The threshold scan data must be in an ascii file (
   containing 4 columns in this order: 
   threshold, rate(Hz), threshold_error, rate_error(Hz)
   The last two columns can be omitted but the fit will be less precise.
   
   -The constructor reads the threshold scan data into an TGraphErrors.
   -The convert function can be used add a calibration or to correct for 
   dead time.
   -The setmodel function can be used to switch between a linear, quadratic and cubic interpolation.
   -The fit_spectr function performs the fit.
   -The draw and draw_spectr functions are used for drawing.

   Typical usage:
   thscan scan(fname);
   scan.fit_spectr(5);
   scan.draw_spectr();

 */
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TList.h"
#include "TObject.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TNamed.h"
#include "TROOT.h"
using namespace std;

enum adccorrection{nocorrection,deadtime_notupd,boh1,boh2};
enum calib{identity,caenadc,caenadc_atten,concast};
enum spectrfit{p1,p2,p3};

double myp1(double*,double*);
double myp2(double*,double*);
double myp3(double*,double*);


class thscan: public TNamed{
public:
  thscan();
  thscan(const char * fname);
  ~thscan();
  void convert(int thcalib=concast, int ratecorrection=nocorrection, 
	       double dead_time=2e-6, double gatetime=0);
  double calibrated(double adcth);
  TGraphErrors* fit_spectr(int delta=12);

  void draw(TGraphErrors* g,const char* opt, TPad* pad=0, 
	    double xmin=-0.25, double ymin=0);
  void draw(const char* opt="APl"){draw(_gthscan,opt);}
  void draw_raw(const char* opt="APl"){draw(_gthscan_raw,opt);}
  void draw_spectr(const char* opt="APl"){if (_gspectr!=0)draw(_gspectr,opt);}
  TGraphErrors* draw_ratio(thscan* t, int dim=0);
  
  TGraphErrors* get(){return  _gthscan;}
  TGraphErrors* get_raw(){return _gthscan_raw;}
  TGraphErrors* get_spectr(){return _gspectr;}
  double get_rate(double threshold){return _gthscan->Eval(threshold);}
  TCanvas * get_canv(){return _canv;}
  double bin_width(){return _bw;}

  void scale_spectr(double scale);  
  TH1F* hconvert(TGraphErrors *gr);
  void setmodel(int model){_fitmodel=model;}

  void add(thscan ts, double coeff);
  void subtract(thscan noise);
  void write_spectr(char *fname);
  //  void fit();
  ClassDef(thscan, 1);
private:
  TGraphErrors* _gthscan_raw;
  TGraphErrors* _gthscan;
  TGraphErrors* _gspectr;
  int _dim;
  int _ical, _adccorrection;
  int _fitmodel;
  double _bw;
  TCanvas *_canv;
};

thscan::thscan(){
  _gthscan_raw = 0;
  _gthscan = 0;
  _gspectr=0;
  _dim=0;
  _ical=identity;
  _adccorrection=nocorrection;
  _fitmodel=p2;
  fName="thscan";
  _canv=0;
}

thscan::thscan(const char* fname){
  //this is the raw th scan (measured rate; adc channel)

  ifstream f;
  f.open(fname);
  if (f.rdstate()!=0) {
    cout<<"Error, file "<<fname<<" not found"<<endl; 
    _gthscan_raw=0;
    _gthscan=0;
    _gspectr=0;
    _dim=0;
    return;
  }
  //find out if file is column separated
  string s;
  getline(f,s);
  f.close();
  if (s.find(";")!=string::npos) 
    _gthscan_raw = new TGraphErrors(fname,"%lg;%lg;%lg;%lg");
  else _gthscan_raw = new TGraphErrors(fname);
  _gthscan_raw->SetNameTitle(fname,fname);
  if (_gthscan_raw->IsZombie()==true){
    _gthscan_raw=0;
    _gthscan=0;
    _dim=0;
  }else{
    _gthscan = (TGraphErrors*)_gthscan_raw->Clone(fname);
    _gthscan->SetNameTitle(fname,fname);
    _dim=_gthscan->GetN();
  }
  _gspectr=0;
  _bw=0;
  _ical=identity;
  _adccorrection=nocorrection;
  _fitmodel=p2;
  fName=fname;
  fTitle=fname;
  _canv=0;
}

thscan::~thscan(){
//   if (_gthscan!=0)
//     delete _gthscan;
   if (_gthscan_raw!=0)
     delete _gthscan_raw;
   if (_gspectr!=0)
     delete _gspectr;
}

/* converts the raw th scan to corret rate for dead time 
   and use the calibration of the adc threshold (to mV) */
void
thscan::convert(int thcalib, int ratecorrection, double dead_time,
		double gatetime){

  if (_gthscan_raw==0) return;
  _ical=thcalib;
  _adccorrection=ratecorrection;

  /* determine sample for which rate is max
     (used in dead time correction func.) */

  double maxrate=0; 
  int ithmax=0;
  for (int ith=0;ith<_dim;ith++){
    /* determine the ith value with max rate
       ( used int the dead time correction) */
    double rate=_gthscan_raw->GetY()[ith];
    if (rate>maxrate) {
      maxrate=rate;
      ithmax=ith;
    }
  }
  cout<<"max rate is "<<maxrate<<" for ith="<<ithmax<<endl;
  //  if (1/maxrate>T) T=1/4./maxrate;

  double* th=new double[_dim];
  double* rate=new double[_dim];
  double* eth=new double[_dim];
  double* erate=new double[_dim];
  /* fill true rate versus pulse height histogram */
  for (int ith=0;ith<_dim;ith++){
    
    th[ith]=calibrated(_gthscan_raw->GetX()[ith]); 
    double thisrate=_gthscan_raw->GetY()[ith];
    double thiserate=_gthscan_raw->GetEY()[ith];

    // allow input file with counts instead of rate
    if (gatetime!=0){
      thisrate=thisrate/gatetime;
      thiserate=thiserate/gatetime;
    }

    if (thisrate*dead_time<1e-4){
      rate[ith]=thisrate;
    }
    else{
      switch (ratecorrection){
      default:
      case nocorrection:
	rate[ith]=thisrate;
	break;
      case deadtime_notupd:
	//     //correct for dead time
	rate[ith]=thisrate/(1-thisrate*dead_time);
	//	cout<<"ith "<<ith<<", rate "<<thisrate<<" --> "<<rate[ith]<<endl;
	break;
      case boh1:      
	//     //correct for dead time
	if (ith<ithmax){
	  rate[ith]=(1 + sqrt(1 - 4.*thisrate*dead_time))/2./dead_time;
	}else{
	  rate[ith]=(1 - sqrt(1 - 4.*thisrate*dead_time))/2./dead_time;}
	break;    
      case boh2:      
	//     //correct for dead time
	if (ith<ithmax){
	  rate[ith]= 1/dead_time*(1/(dead_time*thisrate) - ( 1 -  sqrt(pow(( 1-1/(thisrate*dead_time)),2) -2 )));
	}else{
	  rate[ith]= 1/dead_time*(1/(dead_time*thisrate) - ( 1 +  sqrt(pow(( 1-1/(thisrate*dead_time)),2) -2 )));
	}
	break;    
      }//switch(ratecorrection)
    }// if (rate*dead_time<1e-4)
      // no calibrations for errors.... this should be fixed....
    eth[ith]=_gthscan_raw->GetEX()[ith];
    //erate[ith]=_gthscan_raw->GetEY()[ith];
    erate[ith]=thiserate;
    //    cout<<"ith "<<ith<<" th "<<th[ith]<<" rate "<<thisrate<<" rate corr "<<rate[ith]<<endl;
    //    cout<<"    "<<ith<<" eth "<<eth[ith]<<" erate "<<erate[ith]<<endl;
    
  }//for (ith)


  if (_gthscan!=0) delete  _gthscan;
  //  _gthscan=new  TGraphErrors(_dim,th,rate,eth,erate);
  _gthscan=new  TGraphErrors(_dim,th,rate,eth,erate);
  delete th; delete eth; delete rate; delete erate;
  _gthscan->SetNameTitle(fName,fTitle);
  _gthscan->GetXaxis()->SetTitle("threshold [V]");
  _gthscan->GetYaxis()->SetTitle("freq [Hz]");
  //  _gthscan->GetYaxis()->SetTitleOffset(1.5);

}


/* adc to pulse height calibration */
double thscan::calibrated(double rawth){
  double attenuation=1.0;
  double attenuator_offset=0.0;
  double adc_offset=0.0;
  double th=0;
  switch (_ical){
  case caenadc:
    th=rawth-11.423;
    break;
  case caenadc_atten:
    //calibration of (attenuator + caen discriminator)
    attenuation=3.7/(2.*pow(10,30./20.));
    attenuator_offset=5.5;
    adc_offset=-11.423;
    th=(rawth + adc_offset + attenuator_offset)/attenuation;
    break;
  case concast:
    //calibration of concast prototype comparator
    // it has a 500mV hysteresis.
    // The threshold sets the lower edge of the hysteresis cylce
    // therefore the actual threshold is the external threshold plus 500mV
    //    th=rawth+0.500;
    //    th=rawth+(12-rawth)*5./105.;  //without pulldown
    th=(rawth*1.688 + 0.4167)/1.774; //Flavio
    //th=rawth*1.0483 + 0.235; Enrico
    break;
  default:
    //identity
    th=rawth;
    break;
  }    
  //  cout << "rawth "<<rawth<<";  th "<<th<<endl;
  return th;
}

/* fit th scan with moving window and accumulate derivative */
TGraphErrors* thscan::fit_spectr(int delta){
  bool debug=false;
  if (_gthscan==0) return 0;
  cout<<"fitting "<<_gthscan->GetTitle()<<endl;
  int dim_spectr=_dim-2*delta;
  if (debug) cout<<"dim "<<_dim<<"; dim_spectr = "<<dim_spectr<<endl;
  double   m[dim_spectr];  for (int i=0;i<dim_spectr;i++){m[i]=0.;}
  double  em[dim_spectr];  for (int i=0;i<dim_spectr;i++){em[i]=0.;}
  double mth[dim_spectr];  for (int i=0;i<dim_spectr;i++){mth[i]=0.;}
  double *th=_gthscan->GetX();
  double *counts=_gthscan->GetY();
  double eps=0.1*(th[_dim-1]-th[0])/double(_dim);
  double bw=(th[delta+1]-th[delta]);

  /* prepare fit function */
  string fname;
  int fnpar;
  TF1* fitf;
  double (*fptr)(double*,double*);
  switch (_fitmodel){
  case p1:
    fname="myp1";
    fnpar=3;
    fptr=myp1;
    break;
  default:
  case p2:
    fname="myp2";
    fnpar=4;
    fptr=myp2;
    break;
  case p3:
    fname="myp3";
    fnpar=5;
    fptr=myp3;
    break;
  }
  cout<<"with model "<<fname<<endl;
  
  _bw=bw;
  for (int x0=delta;x0+delta<_dim;x0++){    
    cout<<"."<<flush;
    //    if (th[x0-delta]>0){
    if (debug) 
      cout<<"fitting at x0="<<x0<<" for th="<<th[x0]<<" in interval ["<<th[x0-delta]<<":"<<th[x0+delta]<<"]; bw="<<bw<<endl;//<<" fnpar="<<fnpar<<endl;
    bool allzero=true;
    for (int k=x0-delta;k<=x0+delta;k++)
      if (counts[k]!=0) {
	allzero=false ;
      }
    if (debug){
      cout<<"     x="<<counts[x0];
      cout<<"     ex="<<_gthscan->GetEY()[x0];
      cout<<"; y="<<th[x0];
      cout<<"; ey="<<_gthscan->GetEX()[x0]<<endl;
    }
    int im=x0-delta;
    m[im]=0.;
    em[im]=0.;
    mth[im]=th[x0];
    if (allzero==false){
      fitf=new TF1(fname.c_str(),fptr,th[0],th[_dim-1],fnpar);
      //  if (debug)  fitf->Print("V");
      

      //set initial values for all params to zero
      for (int i=0;i<fnpar; i++){ fitf->SetParameter(i,0.);}
      //set initial values for mu, m and q parameters
      double dy=counts[x0+delta]-counts[x0-delta];
      double dx=th[x0+delta]-th[x0-delta];
      double m0=dy/dx;
      double mu0=th[x0];
      double q0=counts[x0];
      fitf->FixParameter(0,mu0); //fix the mean to the central bin
      fitf->SetParameter(1,m0);
      fitf->SetParLimits(1,-5*fabs(m0),5*fabs(m0));
      fitf->SetParameter(2,q0);
      fitf->SetParLimits(2,-5*fabs(q0),5*fabs(q0));
      //      if (debug) gROOT->ProcessLine(".>.spetr_fit");
      string opt="q";
      if (debug) opt="v";
      _gthscan->Fit(fitf,opt.c_str(),"",th[x0-delta]-eps,th[x0+delta]+eps);
      //      if (debug) gROOT->ProcessLine(".>");
      if (debug) cout<<"      th1="<<th[x0-delta]-eps<<" th2="<<th[x0+delta]+eps<<endl;
      if (fitf != 0){
	m[im]=-1*bw*(fitf->GetParameter(1));
	em[im]=bw*(fitf->GetParError(1));
	if (debug) fitf->Print("V");
	if (debug) cout<<"     fitted error "<<fitf->GetParError(1)<<endl;
	if (debug) cout<<"      fit: m="<<m[im]<<", em="<<em[im]<<endl;
	delete fitf;
      } else {
	cout <<"fit  for bin "<<x0<<" failed"<<endl;
      }//f != 0
    }// allzero==false
    string c;
    if (debug) {
      cout<<"any char to continue.."; cin>>c;
      if (c.compare("q")==0) break;
    }
  }//
  cout<<endl;
  TGraphErrors *spectr=new TGraphErrors(dim_spectr,mth,m,0,em);
  spectr->SetNameTitle(_gthscan->GetName(),_gthscan->GetTitle());
  if (_gspectr!=0) delete _gspectr;
  _gspectr=spectr;
  _gspectr->GetXaxis()->SetTitle("pulse height [V]");
  char tit[256];
  sprintf(tit,"freq [Hz/%2.3fV]",bw);
  _gspectr->GetYaxis()->SetTitle(tit);
  return spectr;
}

void thscan::draw(TGraphErrors* g,const char* opt, TPad* pad, 
		  double xmin, double ymin){
  string sopt(opt);
  if(sopt.find("same")==string::npos){ //if it doesn't find "same"
    if (pad==0){
      _canv=new TCanvas("canv","canv",800,600);
      _canv->cd(0);
      pad=(TPad*)(_canv->GetPad(0));
    }
    if (pad==0) cout<<"Erorr: failed to create Pad"<<endl;
   double ymax=g->GetYaxis()->GetXmax();
    double xmax=g->GetXaxis()->GetXmax();
     pad->DrawFrame(xmin,ymin,xmax,ymax);
     //    pad->SetTitle(g->GetTitle());
  }
  g->Draw(opt);
}


TGraphErrors*
thscan::draw_ratio(thscan *t, int dim){
  TGraphErrors* gts=t->get_spectr();
  if (dim==0){
    if (gts->GetN() != _gspectr->GetN()){
      cout<< "Error: graph dimensions differ"<<endl;
      return 0;
    }
    dim=gts->GetN();
  }
  double *ratio= new double[dim];
  double *th= new double[dim];
  int skip1=0,skip2=0;
  for (int k=0; k<dim; k++){
    int dim1=_gspectr->GetN();
    while ( (_gspectr->GetX()[k+skip1] < gts->GetX()[k+skip2]) && (k+skip1)<dim1 ) 
      skip1++;
    int dim2=gts->GetN();
    while ( (_gspectr->GetX()[k+skip1] > gts->GetX()[k+skip2]) && (k+skip2)<dim2 )
      skip2++;
    if (_gspectr->GetX()[k+skip1] != gts->GetX()[k+skip2]) 
      cout<<"could not find common threshold, skip1="<<skip1<<" skip2="<<skip2<<endl;
    else
    cout<<"th1 = "<<_gspectr->GetX()[k+skip1]<<" th2 = "<<gts->GetX()[k+skip2]<<endl;

    double c_num=_gspectr->GetY()[k+skip1];
    double c_denom=gts->GetY()[k+skip2];

    //    cout<<"c_num ="<<c_num<<"; c_denom="<<c_denom<<endl;
    if (c_denom!=0.){
      ratio[k]=c_num/c_denom;
    }
    else
      ratio[k]=1e10;
    th[k]=_gspectr->GetX()[k+skip1];
    //    cout<<"k="<<k<<" ratio="<<ratio[k]<<" th="<<th[k]<<endl;
  }
  TGraphErrors* gratio=new TGraphErrors(dim,th,ratio,0,0);
  char newTit[512];
  sprintf(newTit,"%s/%s",GetTitle(),gts->GetTitle());
  gratio->SetNameTitle(newTit,newTit);
  draw(gratio,"PL*");
  return gratio;
}

void
thscan::scale_spectr(double scale){
  double *y=_gspectr->GetY();
  double *ey=_gspectr->GetEY();
  for (int i=0;i<_gspectr->GetN();i++){
    y[i]*=scale;
    ey[i]*=scale;
  }
  double max=_gspectr->GetMaximum();
  _gspectr->SetMaximum(scale*max);
}

TH1F*
thscan::hconvert(TGraphErrors *gr){
  TAxis *xa=gr->GetXaxis();
  int nbinsx=xa->GetNbins()-1;
  double xmin=xa->GetXmin();
  double xmax=xa->GetXmax();
  TH1F* h=new TH1F("h","h",nbinsx,xmin,xmax);
  double *x=gr->GetX();
  h->GetXaxis()->Set(nbinsx,x);
  double *y=gr->GetY();
  double *ey=gr->GetEY();
  for (int ibin=1;ibin<nbinsx;ibin++){
    h->SetBinContent(ibin,y[ibin]);
    h->SetBinError(ibin,ey[ibin]);
  }
  return h;
}

void
thscan::add(thscan ts, double coeff){
  //this is the raw th scan (measured rate; adc channel)
  int dim_ts=ts.get_raw()->GetN();
  if (dim_ts!=_dim){
    cout<<"Error: cannot add scans with different sizes"<<endl;
    return;
  }

  double x[_dim],y[_dim],ex[_dim],ey[_dim];
  for (int i=0;i<_dim;i++){
    x[i]=_gthscan_raw->GetX()[i];
    y[i]=_gthscan_raw->GetY()[i]+coeff*ts.get_raw()->GetY()[i];
    ex[i]=_gthscan_raw->GetEX()[i];
    ey[i]=sqrt(pow(_gthscan_raw->GetEY()[i],2)+pow(coeff*ts.get_raw()->GetEY()[i],2));
  }
  
  //  if (_gthscan_raw != 0) delete _gthscan_raw; leak!!!
  _gthscan_raw = new TGraphErrors(_dim,x,y,ex,ey);
  char newTitle[512];
  sprintf(newTitle,"%s+%f%s",GetTitle(),coeff,ts.GetTitle());
  _gthscan_raw->SetNameTitle(newTitle,newTitle);
  //  if (_gthscan != 0) delete _gthscan; leak!!!
  _gthscan = (TGraphErrors*)_gthscan_raw->Clone(fName);
  _gthscan->SetNameTitle(newTitle,newTitle);
  //  if (_gthscan != 0) delete _gthscan; leak!!!
  _gspectr=0;
  _bw=0;
  _ical=identity;
  _adccorrection=nocorrection;
  //  _fitmodel=p2;
  fName=newTitle;
  fTitle=newTitle;
}

void 
thscan::subtract(thscan noise){
  add(noise,-1);
}

void
thscan::write_spectr(char* fname){
  if (_gspectr!=0){
    ofstream data(fname);
    int dim=_gspectr->GetN();
    double *rate=_gspectr->GetY();
    double *th=_gspectr->GetX();
    for (int i=0;i<dim;i++){
      data<<th[i]<<"\t"<<rate[i]<<endl;
    }
    data.close();
  }
}

/* fit functions */

// void
// thscan::fit(){
//   fit_spectr();
//   TH1F* h = hconvert(get_spectr());
//   fit_co60(h,1.2,0.4,2.5);
// }


/* fit function: second degree polynomial centered at par[3]=mu */
double myp1(double *x, double* par){
  double xx=x[0];
  double mu=par[0];
  double m=par[1];
  double q=par[2];
  return q+m*(xx-mu);
}

double myp2(double *x, double* par){
  double xx=x[0];
  double mu=par[0];
  double m=par[1];
  double q=par[2];
  double b=par[3];
  return q+m*(xx-mu)+b*(pow(xx-mu,2));
}

double myp3(double *x, double* par){
  double xx=x[0];
  double mu=par[0];
  double m=par[1];
  double q=par[2];
  double b=par[3];
  double c=par[4];
  return q+m*(xx-mu)+b*(pow(xx-mu,2)) + c*(pow(xx-mu,3));
}
