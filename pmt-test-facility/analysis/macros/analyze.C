void analyze(const string path, 
             const string run_tag, 
	     const string ped_tag, 
	     const int nEvents = -1)
{
  // Non funziona se uso gSystem->LoadMacro... si blocca...
  // Quindi prima fare .L waveforms.C++
  cout << "> Analyzing run: " << run_tag << " with pedestal: " << ped_tag << endl;

  // Load analyzer
  channelWithTrigger ch( path + "/sig_" + run_tag + ".dat", 
			 path + "/trg_" + run_tag + ".dat", 
			 path + "/sig_ped_" + ped_tag + ".dat", 
			 path + "/trg_ped_" + ped_tag + ".dat", 
                         run_tag);
  
  // Run
  ch.read(nEvents);
}
