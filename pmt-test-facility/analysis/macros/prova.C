void prova(const TString filename, const int np)
{
  const int nbins_t = 200;
  const float t_min = -100.;
  const float t_max = 1000.;

  const int nbins_adc = 100;
  const float adc_min = 0.;
  const float adc_max = 600.;

  TFile *_file0 = TFile::Open(filename);
  auto tree = (TTree*)_file0->Get("tree");
  
  //tree->SetAlias("dt", "x_peaks");
  tree->SetAlias("dt", "(x_peaks - sigTime)");

  TCanvas *c0 = new TCanvas("c0", "", 800, 600);
  TH1F *hnp = new TH1F("hnp", "", 10, 0, 10);
  tree->Draw("n_peaks >> hnp", "tauFall > 30");

  // pulse height vs time
  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  c1->SetLogz();
  TH2F *h2d = new TH2F("h2d", "", nbins_t, t_min, t_max, nbins_adc, adc_min, adc_max);
  tree->Draw("y_peaks:dt >> h2d", Form("tauFall > 30. && n_peaks <= %i", np), "colz" );

  // pulse height specrtum
  TCanvas *c2 = new TCanvas("c2", "", 800, 600);
  TH1F *hadc_sig = new TH1F("hadc_sig", "", nbins_adc, adc_min, adc_max);
  TH1F *hadc_aft = new TH1F("hadc_aft", "", nbins_adc, adc_min, adc_max);
  TH1F *hadc_tot = new TH1F("hadc_tot", "", nbins_adc, adc_min, adc_max);
  
  tree->Draw("y_peaks >> hadc_sig", Form("tauFall > 30. && n_peaks <= %i && dt > -10 &&  dt < 30", np) ); // in time
  tree->Draw("y_peaks >> hadc_aft", Form("tauFall > 30. && n_peaks <= %i && dt > 30", np) ); // after-pulse

  hadc_tot->Add(hadc_sig); 
  hadc_tot->Add(hadc_aft); 

  hadc_tot->SetLineColor(kBlack); 
  hadc_tot->SetLineWidth(2); 
  hadc_sig->SetLineColor(kRed); 
  hadc_aft->SetLineColor(kBlue);  

  hadc_tot->Draw("");
  hadc_sig->Draw("same");
  hadc_aft->Draw("same");

  // diff wrt first peak
  TCanvas *c3 = new TCanvas("c3", "", 800, 600);
  TH1F* hdt[np];
  for (int i = 1; i < np; i++) {
    hdt[i] = new TH1F(Form("hdt%i", i), "tauFall > 30. && ", 200, 0., 1000.);
    tree->Draw(Form("x_peaks[%i]-sigTime >> hdt%i",i,i), Form("tauFall > 30. && n_peaks <= %i",np), "goff");
    hdt[i]->SetLineColor(i); 
    hdt[i]->Draw((i==1 ? "" : "same"));
  }

  // print stat
  float nsig = tree->GetEntries(Form("tauFall > 30. && n_peaks <= %i && dt > -10 && dt < 30.", np));
  float nafter = tree->GetEntries(Form("tauFall > 30. && n_peaks <= %i && dt > 50.",np));
  cout << "#sig events: " << nsig << endl;
  cout << "#after-pulses: " << nafter << endl;
  cout << "#ratio: " << 100.*(nafter/nsig) << "%" << endl;
}
