void prova(const TString filename, const TCut clean_cut, const int np)
{
  const int nbins_t = 200;
  const float t_min = 100.;
  const float t_max = 1100.;

  const int nbins_adc = 100;
  const float adc_min = 0.;
  const float adc_max = 600.;

  // Open file
  TFile *_file0 = TFile::Open(filename);
  auto tree = (TTree*)_file0->Get("tree");
  
  // Define cleaning cuts
  TCut cut = clean_cut && Form("n_peaks < %i", np);

  // Aliases
  tree->SetAlias("t", "x_peaks");
  tree->SetAlias("adc", "y_peaks");

  // chi2 
  TCanvas *c00 = new TCanvas("c00", "", 800, 600);
  TH1F *hchi2 = new TH1F("hchi2", "", 50, 0, 100);
  tree->Draw("chi2 >> hchi2", cut);

  // npeaks
  TCanvas *c0 = new TCanvas("c0", "", 800, 600);
  TH1F *hnp = new TH1F("hnp", "", 10, 0, 10);
  tree->Draw("n_peaks >> hnp", cut);

  // pulse height vs time
  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  c1->SetLogz();
  TH2F *h2d = new TH2F("h2d", "", nbins_t, t_min, t_max, nbins_adc, adc_min, adc_max);
  tree->Draw("adc:t >> h2d", cut, "colz" );

  // pulse height specrtum
  TCanvas *c2 = new TCanvas("c2", "", 800, 600);
  TH1F *hadc_sig = new TH1F("hadc_sig", "", nbins_adc, adc_min, adc_max);
  TH1F *hadc_aft = new TH1F("hadc_aft", "", nbins_adc, adc_min, adc_max);
  TH1F *hadc_tot = new TH1F("hadc_tot", "", nbins_adc, adc_min, adc_max);
  
  tree->Draw("adc >> hadc_sig", cut && "t > 170 && t < 180"); // in time
  tree->Draw("adc >> hadc_aft", cut && "t > 200"); // after-pulse

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
    hdt[i] = new TH1F(Form("hdt%i", i), "", 50, 0., 900);
    tree->Draw(Form("t[%i] - t[0] >> hdt%i",i,i), cut, "goff");
    hdt[i]->SetLineColor(i); 
    hdt[i]->Draw((i==1 ? "" : "same"));
  }

  // tau rise vs tau fall
  TCanvas *c4 = new TCanvas("c4", "", 800, 600);
  c4->SetLogz();
  TH2F *htaus = new TH2F("htaus", "", 100, 0., 30., 100, 0., 100);
  tree->Draw("tau_fall:tau_rise >> htaus", cut, "colz");


  // print stat
  float ntot = tree->GetEntries();
  float nsel = tree->GetEntries(cut);
  float nsig = tree->GetEntries(cut && "t > 170 && t < 180.");
  float nafter = tree->GetEntries(cut && "t > 200.");

  cut.Print();
  cout << "# tot events: " << tree->GetEntries() << endl;
  cout << "# nsel/ntot : " << nsel/ntot << endl;
  cout << "# sig events: " << nsig << endl;
  cout << "# after-pulses: " << nafter << endl;
  cout << "# fraction of after pulses: " << 100.*(nafter/nsig) << "%" << endl;
}
