void afterpulse(){
  TFile *_file0 = TFile::Open("10k_1Gs_pmtold000_px04_new.root");
  //
  TTree* tree=(TTree*)_file0->FindObjectAny("tree");
  TCanvas *c1= new TCanvas("canv","afterpulse analysis",800,600);
  tree->Draw("y_peaks:x_peaks>>h2(200,0,1000,100,0,500)","","colz");
  c1->SaveAs("results/px4_AvsT.png");

  TCut presig("x_peaks<175");
  TCut sig("175<x_peaks && x_peaks<195");
  TCut promptap("195< x_peaks && x_peaks<320");
  TCut lateap("x_peaks>320");

  float maxA=600;
  int nbinsA=100;
  
  TH1F *htot = new TH1F("htot","total sample",nbinsA,0,maxA);
  TH1F *hpresig = new TH1F("hpresig","total sample",nbinsA,0,maxA);
  TH1F *hsig = new TH1F("hsig","total sample",nbinsA,0,maxA);
  TH1F *hpromptap = new TH1F("hpromptap","total sample",nbinsA,0,maxA);
  TH1F *hlateap = new TH1F("hlateap","total sample",nbinsA,0,maxA);

  tree->Draw("y_peaks>>htot");
  tree->Draw("y_peaks>>hsig",sig,"same");
  hsig->SetLineColor(kRed);
  tree->Draw("y_peaks>>hpresig",presig,"same");
  hpresig->SetLineColor(kGreen);
  tree->Draw("y_peaks>>hpromptap",promptap,"same");
  hpromptap->SetLineColor(kPink);
  tree->Draw("y_peaks>>hlateap",lateap,"same");
  hlateap->SetLineColor(kBlack);

  TLegend* l=new TLegend(0.5,0.5,0.9,0.9,"PH distribution","NDC");
  l->AddEntry(htot,"total","l");
  l->AddEntry(hpresig,"before signal","l");
  l->AddEntry(hsig,"signal [175-195]","l");
  l->AddEntry(hpromptap,"afterpulse [195-320]","l");
  l->AddEntry(hlateap,"late afterpulse [320-1000]","l");

  l->Draw();

  htot->GetXaxis()->SetTitle("amplitude [ADC]");
  htot->GetYaxis()->SetTitle("entries / 6 ADC");
  c1->SaveAs("results/px4_A_lin.png");
  c1->SetLogy(1);
  c1->SaveAs("results/px4_A_log.png");

  float threshold = 100;
  float bw=maxA/nbinsA;

  float   nsig = hsig->Integral(int(threshold/bw),int(maxA/bw))  ;
  float   npre = hpresig->Integral(int(threshold/bw),int(maxA/bw))  ;
  float   npromptap = hpromptap->Integral(int(threshold/bw),int(maxA/bw))  ;
  float   nlateap = hlateap->Integral(int(threshold/bw),int(maxA/bw))  ;

  cout<<" Threshold  at 100"<<endl      <<" npre = "<<npre<<endl      <<" nsig = "<<nsig<<endl      <<" npromptap = "<<npromptap<<endl      <<" nlateap = "<<nlateap<<endl;

  threshold = 10;
  
   nsig = hsig->Integral(int(threshold/bw),int(maxA/bw))  ;
   npre = hpresig->Integral(int(threshold/bw),int(maxA/bw))  ;
   npromptap = hpromptap->Integral(int(threshold/bw),int(maxA/bw))  ;
   nlateap = hlateap->Integral(int(threshold/bw),int(maxA/bw))  ;

  cout<<" Threshold  at 10"<<endl      <<" npre = "<<npre<<endl      <<" nsig = "<<nsig<<endl      <<" npromptap = "<<npromptap<<endl      <<" nlateap = "<<nlateap<<endl;

    
}

