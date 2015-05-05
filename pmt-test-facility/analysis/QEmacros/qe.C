double integral(TTree* tuv){
  const int n=tuv->GetSelectedRows();
  double *x=tuv->GetV2();
  double *y=tuv->GetV1();
  double sum=0;
  for (int i=0;i<n-1;i++){
    sum+=fabs(x[i+1]-x[i])*(y[i+1]+y[i])/2.;
  }
  return sum;
}

void qe(){
  TTree* t = new TTree("qe","qe sba");
  t->ReadFile("qe.dat","l:qe");
  TTree* tuv = new TTree("qeuv","qe sba+uv window");
  tuv->ReadFile("qeuv.dat","l:qe");
  double hc=2*TMath::Pi()*197;

  tuv->Draw("qe:1.23779e3/l","","");
  TGraph* gqeuv=new TGraph(tuv->GetSelectedRows(),tuv->GetV2(),tuv->GetV1());
  gqeuv->SetMarkerStyle(20);
  gqeuv->SetLineColor(kBlue);  
  gqeuv->SetMarkerColor(kBlue);
  gqeuv->SetTitle("QE as a function of photon energy");
  gqeuv->GetXaxis()->SetTitle("E[ev]");
  gqeuv->GetYaxis()->SetTitle("QE[%]");
  cout<<"integral tuv "<<integral(tuv)<<endl;

  t->Draw("qe:1.23779e3/l","","same");
  TGraph* gqe=new TGraph(t->GetSelectedRows(),t->GetV2(),t->GetV1());
  gqe->SetMarkerStyle(20);
  gqe->SetMarkerColor(kRed);
  gqe->SetLineColor(kRed);
  gqe->SetTitle("QE as a function of photon energy");
  gqe->GetXaxis()->SetTitle("E[ev]");
  gqe->GetYaxis()->SetTitle("QE[%]");
  cout<<"integral t "<<integral(t)<<endl;

  gqeuv->Draw("APL");
  gqe->Draw("same,PL");



  TLegend* l = new TLegend(0.6,0.7,0.9,0.9,"","NDC");
  l->AddEntry(gqe,"Borosilicate","p");
  l->AddEntry(gqeuv,"UV Glass","p");
  l->Draw();

  float effratio=integral(tuv)/integral(t)*100;
  char s[128];
  sprintf(s,"eff ratio %2.2f %%\n",effratio);
  TText* text=new TText(5,25,s);
  text->Draw();

}
