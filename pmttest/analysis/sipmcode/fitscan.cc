/*
  Gabriele Simi, University of Padova, 2011

  A macro to extract a spectrum from a thscan of a SiPM,
  compute initial parameters for the fit function
  and fit.
  
  depends on spectrfitf.cc and thscan.C

*/


#include <algorithm>
using namespace std;

void 
fitradscan(char *dir, char* date, int idet){
  char *fname=new char[256];
  sprintf(fname,"%s/threshold-scan-%s-%d.dat",dir,date,idet);
  cout<<fname<<endl;
  fitscan(fname);
}

TF1*
fitscan(char* fname="thscans/RadHardnessLNL/8_9July2010/threshold-scan-08Jul10-1742-1.dat", 
	double fmin=0, double fmax=1,double width=40e-9){

  bool debug=true;
  //Fit threshold scan
  thscan *scan=new thscan(fname);
  if (debug) cout<<"converting.."<<endl;
  scan->convert(thscan::identity, thscan::deadtime_notupd,width); 
  if (debug) cout<<"fitting.."<<endl;
  scan->fit_spectr(2); 
  if (scan->get_spectr()==0){
    cout<<"Error, fit missing for "<<fname<<" returning null"<<endl;
    return 0;
  }

  TGraphErrors *spectr = new TGraphErrors(*(scan->get_spectr()));

  //Compute fit limits
  double bw=scan->bin_width();
  int dim=spectr->GetN(); int jxmin=int(fmin*dim); int jxmax=int(fmax*dim)-1;
  if (jxmax>dim-1) {cout<<"max out of range"<<endl; return 0;}
  double xmin=spectr->GetX()[jxmin];   double xmax=spectr->GetX()[jxmax];

  //Search for peaks to compute initial values
  TH1F* h=scan->hconvert(spectr);
  TSpectrum pf;
  pf.Search(h,8,"nobackground",1e-2);
  pf.Print("V");

  //Set initial values for fit function
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
    mu=pf.GetPositionY()[1]/pf.GetPositionY()[0];//incorrect...
    if (npeaks>2) {
      pxmin=x[1]; pxmax=x[npeaks-1];
      gain=(pxmax-pxmin)/(npeaks-2);
      offset=pxmin-gain;
      mu=pf.GetPositionY()[2]/pf.GetPositionY()[0];//incorrect...      
    }
    noise=gain/10.;
  }

  double norm=0; 
  for(int i=jxmin;i<jxmax;i++){norm+=spectr->GetY()[i];}
  cout<<"norm "<<norm<<endl;
  
  int npar=12;
  TF1 *f=new TF1("spectrfit",fitf_g2,xmin,xmax,npar);
  f->SetParNames(  "mu","gain","gNoise","offset","iNoise","ln_norm","ct","g2p","g2off","g2sigma");//,"eff","bw");
  f->SetParameters(mu,gain     ,noise  ,offset       ,noise     ,log(norm)  ,0.03, 0.1, gain/5,2*noise);

  f->SetParLimits(0,0,20*mu); //mu
  f->SetParLimits(1,bw,10*gain); //gain
  f->SetParLimits(2,bw/5,50*bw); //gnoise
  f->SetParLimits(3,-100*bw,100*bw);//offset
  //  f->SetParLimits(4,bw/50,50*bw);//inoise
  //  f->FixParameter(4,0);//inoise
  f->SetParLimits(5,0.5*log(norm),1.5*log(norm));//norm
  f->SetParLimits(6,0,0.9);//ct
  //  f->FixParameter(6,0);//ct fixed
  f->SetParLimits(7,0,0.9);//g2p
  //  f->FixParameter(7,0);
  f->SetParLimits(8,0,gain);//g2off
  //  f->FixParameter(8,0);
  f->SetParLimits(9,bw/5,50*bw);//g2sigma
  //  f->FixParameter(9,bw);
  f->FixParameter(10,1);//eff 
  f->FixParameter(11,scan->bin_width());//bw

  cout<<"initial fit parameters:"<<endl;
  double pmin,pmax;
  for (int i=0;i<npar;i++){
    f->GetParLimits(i,pmin,pmax);
    printf("%s\t = %.3f \t [%.3f;%.3f]\n",f->GetParName(i),f->GetParameter(i),pmin,pmax);
  }

  //Fit
  cout<<"fitting spectrum in ["<<xmin<<";"<<xmax<<"] ..."<<endl;
  spectr->Fit(f,"","",xmin,xmax);
  cout<<"done"<<endl;

  //Draw
  spectr->SetMinimum(0.1);
  spectr->Draw("APL");
  if (scan!=0) delete scan;
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
  gROOT->ProcessLine(".L thscan.C+");
  gROOT->ProcessLine(".L spectrfitf.cc+");
  ofstream out; out.open(fitted_data);
  const int nfiles(flist.size());
  const int npars=10; //this is hardwired, should be possible to extract it from f
  double* parmatrix = new double[nfiles*npars*2];
  double gate=48e-9; 
  for (int ifile=0;ifile<nfiles; ifile++){
    TF1* f=fitscan(flist[ifile].c_str(),5./100,1,gate);
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
    out<<darkrate(f,flist[ifile].c_str(),gate)<<"\t";
    out<<filetime(flist[ifile])<<"\t";
    out<<flist[ifile].c_str()<<endl;
    if (f!=0) delete f;
  }
  return parmatrix;
}

TGraphErrors*
makegraph(int jpar, double* parmatrix, const int nfiles){
  double x[nfiles], y[nfiles], ex[nfiles], ey[nfiles];
  int npars=10;
  for (int ifile=0; ifile<nfiles; ifile++){
    x[ifile]=ifile;
    ex[ifile]=0.5;
    y[ifile]=parmatrix[ifile*npars+2*jpar];
    ey[ifile]=parmatrix[ifile*npars+2*jpar+1];
  }
  TGraphErrors *gr=new TGraphErrors(nfiles,x,y,ex,ey);
  return gr;
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


// inline bool operator<(string &a, string &b){
//   int comp= strcmp(a.c_str(),b.c_str());
//   return comp<0? 0 : 1;
// }


int
filetime(string f){
  int p=f.find("threshold-scan-");
  p+=sizeof("threshold-scan");
  int d=atoi(f.substr(p,2).c_str());
  char *month=f.substr(p+2,3).c_str();
  int y=atoi(f.substr(p+5,2).c_str());
  int h=atoi(f.substr(p+8,2).c_str());
  int m=atoi(f.substr(p+10,2).c_str());
  //compute time since 18feb11 in seconds
  int t=(d-18)*24*3600;
  t+=h*3600;
  t+=m*60;
  
  cout <<" d="<<d<<" month="<<month<<" y="<<y<<" m="<<m<<endl;
  return t;
}


double darkrate(TF1* f, char* fname, double width){
  double th=f->GetParameter("offset")+f->GetParameter("gain")/2;
  thscan scan(fname);
  scan.convert(thscan::identity, thscan::deadtime_notupd,width); 
  return scan.get_rate(th);
}
