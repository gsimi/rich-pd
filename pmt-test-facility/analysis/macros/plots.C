void plots(TString fileName = "out.root", bool isCrossTalk = false) 
{
  // Input file 
  TFile *_file0 = TFile::Open( fileName );
  
  const unsigned nbin = 410;
  const double xmin = -52; 
  const double xmax = 768; 

  // Histograms
  TH1F *hAll = new TH1F("hAll", "all", nbin, xmin, xmax); //fit
  TH1F *h1 = new TH1F("h1", "Left sideband", nbin, xmin, xmax);
  TH1F *h2 = new TH1F("h2", "Signal peak", nbin, xmin, xmax);
  TH1F *h3 = new TH1F("h3", "Right sideband", nbin, xmin, xmax);

  TH1F *hIntAll = new TH1F("hIntAll", "all", nbin, xmin, xmax); //integral
  TH1F *hInt1 = new TH1F("hInt1", "Left sideband", nbin, xmin, xmax);
  TH1F *hInt2 = new TH1F("hInt2", "Signal peak", nbin, xmin, xmax);
  TH1F *hInt3 = new TH1F("hInt3", "Right sideband", nbin, xmin, xmax);

  TH1F *hdt = new TH1F("hdt", "dt", 200, -20.0, 100.);
  TH2F *h2d = new TH2F("h2d", "dt vs ph", 200, -20., 100., nbin, xmin, xmax);
  TH2F *hInt2d = new TH2F("hInt2d", "dt vs integral", 200, -20., 100., nbin, xmin, xmax);
  TH2F *h2d_taus = new TH2F("h2d_taus", "tau_r vs tau_f", 100, 0, 500, 100, 0, 50);

  // Cuts
  tree->SetAlias("dt", "sigTime - trgTime");
  TCut cleancut = "(tauRise > 2. && tauRise < 48.) && (tauFall > 55. && tauFall < 490.)";
  TCut tauFallCut = "tauFall < 150";  
  cleancut = "";
  
  // Define dt regions
  TCut cut1 = "dt < 60." && cleancut;
  TCut cut2 = "dt > 60. && dt < 80." && cleancut; //signal peak
  TCut cut3 = "dt > 80." && cleancut;

  if (isCrossTalk) {
    cut1 = "dt < 50." && cleancut;
    cut2 = "dt > 50. && dt < 80." && cleancut; 
    cut3 = "dt > 80." && cleancut;
  }

  // Make histograms
  tree->Draw("dt>>hdt", cleancut, "goff");
  tree->Draw("pulseHeight:dt>>h2d", cleancut, "goff");
  tree->Draw("integral:dt>>hInt2d", cleancut, "goff");
  tree->Draw("tauRise:tauFall>>h2d_taus", cut2, "goff");
  tree->Draw("tauRise>>htr(100,0,50)", cleancut);
  tree->Draw("tauFall>>htf(100,0,500)", cleancut);
  tree->Draw("chi2>>hchi2(100,0,1)", cleancut);

  tree->Draw("pulseHeight>>hAll", cleancut, "goff");
  tree->Draw("pulseHeight>>h1", cut1, "goff");
  tree->Draw("pulseHeight>>h2", cut2, "goff");
  tree->Draw("pulseHeight>>h3", cut3, "goff");

  tree->Draw("integral>>hIntAll", cleancut, "goff");
  tree->Draw("integral>>hInt1", cut1, "goff");
  tree->Draw("integral>>hInt2", cut2, "goff");
  tree->Draw("integral>>hInt3", cut3, "goff");

  // Now fit method
  hAll->SetMarkerSize(1);
  hAll->SetMarkerStyle(8);
  hAll->SetLineColor(kBlack);
  hAll->SetXTitle("Pulse height");
  h1->SetLineColor(kGreen);
  h2->SetLineColor(kRed);
  h3->SetLineColor(kBlue);

  hAll->Draw("pe");
  h1->Draw("same");
  h2->Draw("same");
  h3->Draw("same");

  c1->Update();
  c1->WaitPrimitive();

  // Now integral method
  hIntAll->SetMarkerSize(1);
  hIntAll->SetMarkerStyle(8);
  hIntAll->SetLineColor(kBlack);
  hIntAll->SetXTitle("Pulse height");
  hInt1->SetLineColor(kGreen);
  hInt2->SetLineColor(kRed);
  hInt3->SetLineColor(kBlue);

  hIntAll->Draw("pe");
  hInt1->Draw("same");
  hInt2->Draw("same");
  hInt3->Draw("same");
}



void dtplotdiode(TString fname="out.root",bool save=false){
  plots(fname);//substitute out.root with the correct name
  TCut  sigcut = "(tauRise>5) && (tauFall>170) && pulseHeight>30";
   TCut cleancut = "(tauRise > 2. && tauRise < 48.) && (tauFall > 55. && tauFall < 490.)";
   tree->Draw("dt/5>>hdtns(350,0,70)","dt<2e3 & dt> 0" && cleancut && sigcut,"colz");
   hdtns->SetXTitle("ns");
   hdtns->SetYTitle("entries/0.2ns");

   if (save) c1->SaveAs("dt_diodo_5V.png");
   c1->SetLogy(true);
   gStyle->SetOptFit(11111);
   //  hdtns->Fit("gaus","","",23,37); //for diode
   hdtns->Fit("gaus","","",11,15); //for laser
   if (save) c1->SaveAs("dt_diodo_5V-logfit.png");
   
   c1->SetLogy(false);
   tree->Draw("pulseHeight:dt/5>>h2dpht(350,0,70,100,-50,450)","dt<2e3 & dt> 0" && cleancut && sigcut,"colz");
   tree->Draw("pulseHeight:dt/5>>h2dphtprof(350,0,70,100,-50,450)","dt<2e3 & dt> 0" && cleancut && sigcut,"prof");
   
   
}
