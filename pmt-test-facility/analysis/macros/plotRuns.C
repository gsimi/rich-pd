{
  // Define data runs 
  vector<TString> runs;
  runs.push_back("20160111_px28_diode075V.root");
  runs.push_back("20160111_px28_diode08V.root");
  runs.push_back("20160111_px28_diode085V.root");
  runs.push_back("20160111_px28_diode09V.root");
  runs.push_back("20160111_px28_diode095V.root");
  runs.push_back("20160111_px28_diode1V.root");
  runs.push_back("20160111_px28_diode105V.root");
  runs.push_back("20160111_px28_diode11V.root");
  runs.push_back("20160111_px28_diode115V.root");
  runs.push_back("20160111_px28_diode12V.root");
  runs.push_back("20160111_px28_diode125V.root");
  runs.push_back("20160111_px28_diode13V.root");

  TCut cleanCut = "(tauRise > 2. && tauRise < 28.) && (tauFall > 35.1 && tauFall < 490.) && chi2 < 100";
  TCut timeCut = "abs(sigTime - trgTime - 200) < 40";
  TCut tauFallCut = "tauFall < 70";  

  TString what;
  TCanvas *cInt = new TCanvas("cInt", "Integral");
  TCanvas *cTime = new TCanvas("cTime", "sigTime -trgTime");
  TCanvas *cTaus = new TCanvas("cTaus", "tauRise vs tauFall");
  TCanvas *cChi2 = new TCanvas("cChi2", "chi2");
  cInt->Divide(3, 4);
  cTime->Divide(3, 4);
  cTaus->Divide(3, 4);
  cChi2->Divide(3, 4);

  // Loop over dataset
  for (size_t i = 0; i < runs.size(); i++) {
    cout << "> Running run: " << runs[i] << endl;
    TFile *file = new TFile(runs[i], "", "read");
    TTree *tree = (TTree*)file->Get("tree");

    if (i < 3) tauFallCut = "tauFall < 110";
    if (i >= 3 && i < 6) tauFallCut = "tauFall < 80";
    if (i >= 6) tauFallCut = "tauFall < 60";

    // dt
    TCanvas *c2 = (TCanvas*)cTime->cd(i + 1);
    //c2->SetLogy();
    what.Form("sigTime - trgTime>>hTime%i(100, 100, 400)", i);
    tree->Draw(what, cleanCut);

    // Integral
    TCanvas *c1 = (TCanvas*)cInt->cd(i + 1);
    c1->SetLogy();
    what.Form("integral>>hInt%i(100, -200, 600)", i);
    //what.Form("pulseHeight>>hInt%i(100, -200, 600)", i);
    tree->Draw(what, !tauFallCut);

    // tauRise vs tauFall
    TCanvas *c3 = (TCanvas*)cTaus->cd(i + 1);
    what.Form("tauRise:tauFall>>hTaus%i(80, 30., 250., 80, 0., 30.", i);
    tree->Draw(what, cleanCut && timeCut, "colz");

    // chi2
    TCanvas *c4 = (TCanvas*)cChi2->cd(i + 1);
    what.Form("chi2>>hChi2%i(80, 0, 200)", i);
    tree->Draw(what, cleanCut && timeCut);
  }

}
