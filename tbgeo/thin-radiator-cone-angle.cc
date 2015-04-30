{
  gROOT->Reset();

   /* 
      generation of cherenkov light from thin radiator
   */
   
  //Radiator: index of refraction, 
  //angle of inclination of the radiatior with respect to th beam
  //thickness
  float n=1.54; //this is SiO2 at 600nm
  float alpha=45./180*TMath::Pi();
  float radiator_thick=1;//[cm]

  //cherenkov angle
  float thetac=acos(1./n);

  //max exit angle of cherenkov radiation
  float maxtheta_o=sqrt(n*n-1)*sin(alpha)-cos(alpha);

  //max angle of cherenkov light with respect to beam
  float maxbeta=maxtheta_o+alpha;
  cout<<"max angle of cherenkov light with respect to beam "<<maxtheta_o<<endl;



  TVector3 beam(0,0,1);
  ///direction orthogonal to the surface of the radiator 
  TVector3 radiator_dir(0,-1*sin(alpha),cos(alpha));

  //angle of the cherenkov radiation with respect to the x axis 
  float phi=-1*TMath::Pi()/2;

  //direction of cherenkov photons inside radiator
  TVector3 inside_dir(sin(thetac)*cos(phi),
		      sin(thetac)*sin(phi),
		      cos(thetac));

  //angle of incidence of cherenkov radiation on the 
  //inner surface of the radiator
  float theta_i =inside_dir.Angle(radiator_dir);

  if (n*sin(theta_i)>1) {
    cout<<"internal reflection for phi = "<<phi<<endl;
    continue;
  }

  // angle of refracted light
  float n_air=1.;
  float theta_o=asin(n/n_air*sin(theta_i));


  TVector3 u=radiator_dir.Cross(inside_dir);
  TVector3 v=u.Cross(radiator_dir);
  TVector3 refracted_dir=radiator_dir*sin(theta_o)+v*cos(theta_o);
  float beta=refracted_dir.Angle(beam);

  cout<<"Using three vectors the angle is "<<beta*180/TMath::Pi()<<endl;
  
  c1 = new TCanvas("c1","PolyLine3D & PolyMarker3D Window",200,10,700,500);
  
  // 3D Graph
   // create a pad
  p1 = new TPad("p1","p1",0.05,0.02,0.95,0.82,46,3,1);
  p1->Draw();
  p1->cd();
  // creating a view
   view = TView::CreateView(1);
   view->SetRange(-30,30,-30,30,30,30);
  TPolyLine3D *ray=new TPolyLine3D(3);
  ray->SetPoint(0,

 
  
}
