#define DIM 8
#define NPX DIM*DIM

void plotFitResults()
{
  const Float_t maxGain = 100.; // allow fitted gain less than this value

  Float_t fitGain[NPX]; 
  Float_t HamGain[NPX]; 

  // Define dataset
  vector<TString> dataset;
  dataset.push_back("multigauss_9726E473-001827_v3.root");


  //=========================
  // Loop over dataset
  //=========================
  for (size_t idat = 0; idat < dataset.size(); idat++) {
    std::cout << "> Processing file: " << dataset[idat] << std::endl;
    // Do fit
    //multigauss( runs[i], "MAROC", -1 );
  
    // Open input file
    TFile *inputFile = new TFile(dataset[idat], "read");
    TTree *tree = (TTree*)inputFile->Get("tree");
  
    // Set branch address
    tree->SetBranchAddress("mean_1pe", fitGain);
    tree->SetBranchAddress("HamGain", HamGain);
   
    // Read gains
    tree->GetEntry(0);
  
    // Find max. gain
    std::pair<Int_t, Float_t> maxHamGain(-1, -99999.);
    std::pair<Int_t, Float_t> maxFitGain(-1, -99999.);
    for (size_t ipx = 0; ipx < NPX; ipx++) {
      if (fitGain[ipx] > maxGain) fitGain[ipx] = 0.; // WARNING: remove this cut when fit is working well...
  
      if (HamGain[ipx] > maxHamGain.second) { 
        maxHamGain.first = ipx;
        maxHamGain.second = HamGain[ipx];
      }
      if (fitGain[ipx] > maxFitGain.second) { 
        maxFitGain.first = ipx;
        maxFitGain.second = fitGain[ipx];
      }
    }
  
    std::cout << "> Max gain (Hamamatsu): " << "px = " << maxHamGain.first 
              << ", gain = " << maxHamGain.second << std::endl;
    std::cout << "> Max gain (fitted): " << "px = " << maxFitGain.first 
              << ", gain = " << maxFitGain.second << std::endl;
  
    // Normalize fitted gains
    for (size_t ipx = 0; ipx < NPX; ipx++) {
      fitGain[ipx] /= maxFitGain.second / 100.;
    }
  
    // Draw 
    TCanvas *canv = new TCanvas(Form("canv_%zu", idat), "", 10, 10, 850, 650);
    canv->Divide(2, 2);
  
    // Book histograms
    TH2F *hHamGains = new TH2F("hHamGains", "Hamamatsu gains", 8, 0., 8., 8, 0., 8.);
    TH2F *hFitGains = new TH2F("hFitGains", "Fitted gains", 8, 0., 8., 8, 0., 8.);
    TH2F *hGainCorrel = new TH2F("hGainCorrel", "Gain correlation", 50, 50., 100., 50, 50., 100.);
    TH2F *hSigmaGainCorrel = new TH2F("hSigmaGainCorrel", "Sigma 1pe - Gain correlation", 50, 12., 22., 50, 5., 7.);
  
    for (size_t row = 0; row < DIM; row++) {
      for (size_t col = 0; col < DIM; col++) { 
        Int_t hamG = HamGain[DIM*row + col];
        Int_t fitG = fitGain[DIM*row + col];
  
        //cout << row << " " << col << " " << HamamatsuGain[DIM*row + col] << endl;
        hHamGains->Fill(col, DIM -1 -row, hamG);
        hFitGains->Fill(col, DIM -1 -row, fitG);
        hGainCorrel->Fill(fitG, hamG);
      }
    }
  
    canv->cd(1);
    hHamGains->SetMinimum(60.);
    hHamGains->SetMaximum(100.);
    hHamGains->Draw("colz text");
  
    canv->cd(2);
    hFitGains->SetMinimum(60.);
    hFitGains->SetMaximum(100.);
    hFitGains->Draw("colz text");
    
    canv->cd(3);
    hGainCorrel->GetXaxis()->SetTitle("Fitted gain");
    hGainCorrel->GetYaxis()->SetTitle("Hamamatsu gain");
    hGainCorrel->Draw("colz");
    hGainCorrel->Fit("pol1");
   
    canv->cd(4);
    hSigmaGainCorrel->GetXaxis()->SetTitle("Fitted sigma 1pe");
    hSigmaGainCorrel->GetYaxis()->SetTitle("sqrt(fitted gain)");
    tree->Draw("sqrt(mean_1pe):sigma_1pe >> hSigmaGainCorrel", "", "colz");
    hSigmaGainCorrel->Fit("pol1");
  
    // Save plots
    canv->SaveAs(Form("outputPlots_%zu.C", idat));
  }//dataset

}

