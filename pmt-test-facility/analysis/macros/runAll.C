void runAll() 
{
  // Non funziona se uso gSystem->LoadMacro... si blocca...
  // Quindi prima fare .L waveforms.C++ poi .x run.C

  // Dataset path
  string path = "../data/dataset_stefano/";

  // Define data runs 
  vector<string> runs;
  runs.push_back("20160111_px28_diode075V");
  runs.push_back("20160111_px28_diode08V");
  runs.push_back("20160111_px28_diode085V");
  runs.push_back("20160111_px28_diode09V");
  runs.push_back("20160111_px28_diode095V");
  runs.push_back("20160111_px28_diode1V");
  runs.push_back("20160111_px28_diode105V");
  runs.push_back("20160111_px28_diode11V");
  runs.push_back("20160111_px28_diode115V");
  runs.push_back("20160111_px28_diode12V");
  runs.push_back("20160111_px28_diode125V");
  runs.push_back("20160111_px28_diode13V");
  
  // Loop over dataset
  for (size_t i = 0; i < runs.size(); i++) {
    cout << "> Running run: " << runs[i] << endl;

    // Load analyzer
    channelWithTrigger ch( path + runs[i] + "_sig.dat", 
                           path + runs[i] + "_trg.dat", 
                           path + "20160111_px28_sigped.dat", 
                           path + "20160111_px28_trgped.dat" );
    
    // Run
    ch.read(50000);

    // Change name of the output
    gSystem->Exec(("mv out.root " + runs[i] + ".root").c_str());
  }

}
