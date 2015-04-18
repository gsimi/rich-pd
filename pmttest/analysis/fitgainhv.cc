
/*Plot the datasheet values in a 2D histogram*/

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


/*Fonction takes fitted values to plot it in a 2D histogram*/

TH2F* uniformity(const char* rdata="test.dat", int ipar=3, int normpx=61){



  const int npar(25);
  if(ipar>=npar){return 0;}
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




//fit function for the gain VS high voltage
Double_t gainfit(Double_t* x, Double_t* par) {
  /*===========================================================
   *x       -> voltage in kV
    par[0]  ->   proportional constant
    par[1]  ->   exponent paramenter
    par[2]  ->   number of dynodes

  =============================================================*/
  Double_t c = par[0];
  Double_t alpha = par[1];
  Double_t n= par[2];
  Double_t nalpha = n*alpha;
  return pow(c,n)*pow(*x,nalpha);
}


//Draw the voltage dependency of a parameter and if the parameter is the gain, fit it
TGraphErrors* voltagedependency(const char* rdata="test.dat", int ipar=3){


  Double_t x[20], y[20],ex[20], ey[20];
  int nvalues=14;
  const int npar(25);
  if(ipar>=npar){return 0;}
  ifstream file(rdata);
  string  ch;
  string parname[npar];
  for(int i=0; i<npar;i++){
    file>>parname[i];
  }
  //const char * title=parname[ipar].c_str();
  //sprintf(var,"gain",ipar);
  int j=0;
  double value[npar];
  while(!file.eof()){
    for (int i=0;i<npar;i++) {
      file>>value[i];
      if (file.eof()) break;
    }
    if (file.eof()) break;
    x[j]=value[24]/1e3; //HV in kV
    ex[j]=0;
    y[j]=value[ipar];  //value of the fitted paramter
    if(ipar<npar-1) ey[j]=value[ipar+1];
    else ey[j]=0;

    
    if (ipar == 17){ //special case for probability of signal on 1st dynode)
      double frac=value[ipar];
      double npe=value[1];
      //      double npe1=value[15];
      //      y[j]=1./((1./frac-1)*(1-TMath::Exp(-npe))/(1-TMath::Exp(-npe1)) + 1);
      //      parname[ipar]="1stdynode fraction";
      y[j]=(1.-frac)*(1-TMath::Exp(-npe));// + frac*(1-TMath::Exp(-npe1));
      parname[ipar]="Quantum efficiency";
      
    }
    
    j++;
    // cout<<" x "<<x<<" y "<<y<<" val "<<value[ipar]<<endl;
  }
  TGraphErrors *g = new TGraphErrors(nvalues,x,y,ex,ey);
  const char * title=parname[ipar].c_str();
  g->SetTitle(title);
  g->SetMarkerSize(1.5);
  g->SetMarkerStyle(5);
  g->Draw("LPA");

  if (ipar==3){
    //fit gain depency on voltage
    int nfiles=j;
    TF1* f=new TF1("gainfit",gainfit,x[0],x[nfiles-1],3);
    cout<<"x0 "<<x[0]<<" x1 "<<x[nfiles-1]<<endl;
    f->SetParameter(0,50);
    f->SetParNames("const","alpha","ndynodes");
    f->SetParLimits(0,0,10);
    f->SetParameter(1,0.75);
    f->SetParLimits(1,0,2);
    f->FixParameter(2,12); // fix number of dynodes
    g->Fit(f,"","",x[0],x[nfiles-1]);
    double c= f->GetParameter(0);
    double alpha = f->GetParameter(1); //alpha should be 2./3
    double n = f->GetParameter(2);

    const int ndynodes=12;
    double fHV[ndynodes]={2.3,1.2,1,1,1,1,1,1,1,1,1,0.5}; 
    double fnorm=0; for (int i=0;i<ndynodes;i++){fnorm+=fHV[i];}
    double fproduct=1; for (int i=0;i<ndynodes;i++){fproduct*=fHV[i]/fnorm;}
    double eoveradc=3.5e4;
    //1./1000**alpha takes into account the fact that we use HV in kV
    double k = c*pow(eoveradc,1./n)/pow(fproduct,alpha/n)/pow(1000,alpha);
    cout << " k "<<k<<endl;

    double HV=950;
    double gain1 = k* pow(HV*2.3/13,alpha);
    cout << " gain1 "<<gain1<<endl;
    double fgain2 = k* pow(HV*1.2/13,alpha);
    cout << " fgain2 "<<fgain2<<endl;
    double fgain3 = k* pow(HV/13,alpha);
    cout << " fgain3 "<<fgain3<<endl;
    double gain12 = k* pow(HV*0.5/13,alpha);
    cout << " gain12 "<<gain12<<endl;
  }
  return g;
}



double
computetheoreticalgainrms_1stoverpk(){
  double g[12];
  for (int i=0;i<12;i++){g[i]=3.2;}; 
  g[0]=6; g[1]=3.7; g[11]=1.9; 

  double n[13], en[13];

  cout<<"conversion on photocade"<<endl;
  n[0]=1;en[0]=0;
  for (int i=0;i<12;i++){n[i+1]=n[i]*g[i]; en[i+1]=sqrt(n[i+1]+pow((en[i]*g[i]),2));}

  cout<<"gain:"<<endl;
  for (int i=0;i<12;i++){cout<<g[i]<<":";} cout<<endl;
  cout<<"number of electrons:"<<endl;
  for (int i=0;i<13;i++){cout<<n[i]<<":";} cout<<endl;
  cout<<"elative fluctuation:"<<endl;
  for (int i=0;i<13;i++){cout<<en[i]/n[i]<<":";} cout<<endl;
  double rk=en[12]/n[12];
  cout<<endl;

  cout<<"conversion on first dynode"<<endl;
  cout<<"gain:"<<endl;
  n[1]=1; en[1]=0;
  for (int i=1;i<12;i++){n[i+1]=n[i]*g[i]; en[i+1]=sqrt(n[i+1]+pow((en[i]*g[i]),2));}

  for (int i=1;i<12;i++){cout<<g[i]<<":";} cout<<endl;
  cout<<"number of electrons:"<<endl;
  for (int i=1;i<13;i++){cout<<n[i]<<":";} cout<<endl;
  cout<<"elative fluctuation:"<<endl;
  for (int i=1;i<13;i++){cout<<en[i]/n[i]<<":";} cout<<endl;
  double r1=en[12]/n[12];

  return r1/rk;

}
