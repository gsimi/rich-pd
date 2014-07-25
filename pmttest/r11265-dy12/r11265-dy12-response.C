float round(float f){
  return floor(f+0.5);
}
bool comment(string s, std::fstream& dati){
  if (!strncmp(s.c_str(),"#",strlen("#"))){
    string header; getline(dati,header);
    cout<<"skipped comment: "<<header<<endl; 
    return true;
  }
  return false;
}
bool empty(string s){
  if (!strcmp(s.c_str(),"") ||
      !strcmp(s.c_str(),"\n")){
    cout << "skipping empty "<<endl<<flush;
    return true;
  }
  return false;
}


void r11265_dy12_response(char* fn="./r11265-dy12-response.dat"){
  std::fstream dati(fn,std::ios_base::in);
  string s; //string for reading a line from the file
  const int n(8); //dimension of pixel array
  float muarr[n*n],x[n*n],y[n*n]; //arrays for storing pixel index and measured signal
  float pixel[n*n];
  float mean=0, pedestal=0;//mean value of the signal for normalization
  int i=0;//arbitrary sequential index identifying a pixel
  TH2F *hdy12=new TH2F("hdy12","dynode 12 response",n,0,n,n,0,n);
  TH2I *hpmap=new TH2I("pmap","pixel map",n,0,n,n,0,n);
  //read each line of the file
  do
    {
      dati>>s; 
      //remove headers begining with #
      if (comment(s,dati)) continue;
      if (empty(s)) continue;
      if (!strncmp(s.c_str(),"pedestal",strlen("pedestal"))){
	dati>>pedestal; getline(dati,s);
	cout<<"read pedestal "<<pedestal<<endl;
	continue;
      }
      
      //    else {
      cout<<"reading value "<<s<<endl;
      //fill arrays for plotting
      x[i]=i%8;    //pixel columns along y axis
      y[i]=int(i/8);
      stringstream(s)>>muarr[i];
      muarr[i]=muarr[i]-pedestal;
      pixel[i]=i+1;
      //update mean
      mean+=muarr[i];
      //update counter
      i++;
      //    }
    } while (!dati.eof())

  //check if we read all the data
  if (i != n*n){cout<<"read unexpected number of data : "<<i<<endl; return 0;}
  mean=mean/i;  //mean normalized to 100;
  //  getchar();

  // normlize the data such that the mean=100;
  //printout what we collected (for debugging)
  for (int j=0;j<n*n;j++){
    muarr[j]=100*muarr[j]/mean;
    //    cout<<j%8<<","<<int(j/8)<<"=";
    cout<<muarr[j];
    if (j%8==7) cout<<endl;
    hdy12->SetBinContent(x[j]+1,y[j]+1,round(muarr[j]));
    hpmap->SetBinContent(x[j]+1,y[j]+1,pixel[j]);
    
  }
  hdy12->SetStats(false);
  
  TCanvas *canv= new TCanvas("canv","canv",800,600);
  canv->Divide(2,2);
  canv->cd(1);
  hdy12->Draw("colz");
  hpmap->Draw("same text");
  canv->cd(2);
  hdy12->Draw("colz");
  hdy12->Draw("same text");

  TGraph2D *dy12map = new TGraph2D(n*n,x,y,muarr);
  dy12map->SetName("d12"); 
  dy12map->SetTitle("dynode 12 response");
  TGraph2D *dy12pixel = new TGraph2D(n*n,x,y,pixel);
  dy12pixel->SetName("pixmap"); 
  dy12pixel->SetTitle("pixel numbering");

  canv->cd(3);
  dy12map->Draw("Acolz");

  canv->cd(4);
  dy12map->Draw("tri1");
  //  dy12pixel->Draw("same");
  
}
