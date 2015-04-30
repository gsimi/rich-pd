#include "TGeoManager.h"
#include "TVector3.h"
#include "TTree.h"
#include "TMath.h"
#include <iostream>
#include "TVirtualGeoTrack.h"
#include "TCanvas.h"
#include "TView3D.h"
class tbconfig{
public:
  tbconfig();
  ~tbconfig(){};
  //box size
  double lxbox;
  double lybox;
  double lzbox;
  //radiator size
  double radheight;
  double radwidth;
  double radthick;
  //radiator center position inside box
  double pxrad;
  double pyrad;
  double pzrad;
  //radiator inclination around x axis
  double incly;
  //origin of gammas
  double ox;
  double oy;
  double oz;

  TVector3 dir_track;
  double nradiator;
  double nair;
};
   

tbconfig::tbconfig(){
  //box size
  lxbox=60; 
  lybox=60; 
  lzbox=100;

  //radiator size  
  radheight=60; 
  radthick=60; 
  radwidth=100;

  //radiator center position inside box
  pxrad=0;
  //  pyrad=-lybox/4;
  pyrad=0.;
  pzrad=0;

  //radiator inclination around x axis
  incly=0;  
  
  //origin of gammas
  ox=pxrad;
  oy=pyrad;//+radthick/3;
  oz=pzrad;

  nradiator=1.47;
  //nradiator=1.3;
  nair=1.;

  dir_track.SetXYZ(0,0,1);

}


class refractor{
public:
  refractor(TGeoManager *geom);
  void propagateToBndry();
  bool refract();
  TVector3 get_dir(){return _next_dir;}
  bool debug;
private:
  TGeoManager* _geom;
  TVector3 _normal_dir;
  double _n2;
  double _n1;
  TVector3 _next_dir;
};
refractor::refractor(TGeoManager *geom){
  _geom=geom;
  debug=false;
  geom->FindNode();
  _n1=_geom->GetCurrentVolume()->GetMedium()->GetParam(0);
  
  if (debug){
    cout<<"current node is "<<_geom->GetCurrentNode()->GetName()<<endl;
    cout<<"current refraction index "<<_n1<<endl;
  }
}

//track to next boundary and cross it
void refractor::propagateToBndry(){
  if (debug){
    cout<<"propagating to next boundary"<<endl;
    const double* current_dir=_geom->GetCurrentDirection();
    TVector3 vcurrent_dir(current_dir[0],current_dir[1],current_dir[2]);
    cout<<"current dir "; vcurrent_dir.Print();
  }


  //compute distance to next boundary
  TGeoNode* nextNode=_geom->FindNextBoundary(-1000); 
  //TGeoNode* nextNode=_geom->FindNextBoundary(); 
  if (debug){
    cout<<"Found boundary, returned node is "<<nextNode->GetName()<<endl;
    cout<<"Distance to Boundary is "<<_geom->GetStep()<<" cm."<<endl;
  }

  //find normal to surface
  double* normal_dir_master=_geom->FindNormalFast(); 
  //  TVector3 vnormal_dir_master(normal_dir_master[0],normal_dir_master[1],normal_dir_master[2]);
  //  double* normal_dir=new double[3];
  //  _geom->MasterToLocalVect(normal_dir_master,normal_dir);
  _normal_dir.SetXYZ(normal_dir_master[0],normal_dir_master[1],normal_dir_master[2]);
  //_geom->LocalToMaster(normal_dir_master,normal_dir);
  if (debug) {
    cout<<"normal direction in master coordinates is"<<endl;
    //    vnormal_dir_master.Print();
    //    cout<<"normal direction in local coordinates is"<<endl;
    _normal_dir.Print();
  }

  //now actually move to the new location
  _geom->Step();
  _geom->FindNode();
  if (debug) cout<<"next node is "<<_geom->GetCurrentNode()->GetName()<<endl;       

  _n2=_geom->GetCurrentVolume()->GetMedium()->GetParam(0);
  if (debug) {
    if (_n2!=_n1){
      cout<<"refraction index has changed: normal direction  ";
      cout<<"current n "<<_n1<<endl;
      cout<<"next n "<<_n2<<endl;
    } else {
      cout<<"refraction index unchanged"<<endl;
    }
  }
}

bool refractor::refract(){

  //define u,v,normal vectors forming an orthogonal ref system  
  const double* current_dir=_geom->GetCurrentDirection();
  TVector3 vcurrent_dir(current_dir[0],current_dir[1],current_dir[2]);
  TVector3 w(_normal_dir);
  if (debug){
    cout<<"** computing refraction"<<endl;
    cout<<"current direction "; vcurrent_dir.Print();
    cout<<"normal direction "; _normal_dir.Print();
  }
  w=w.Unit();
  TVector3 v=(w.Cross(vcurrent_dir)).Unit();
  TVector3 u=(v.Cross(w)).Unit();
  if (debug){ 
    cout<<"current dir "; vcurrent_dir.Print();
    cout<<"u dir "; u.Print();
    cout<<"v dir "; v.Print();
    cout<<"w dir "; w.Print();
  }

  //angle of incidence of cherenkov radiation on the 
  //inner surface of the radiator
  float theta_i =vcurrent_dir.Angle(_normal_dir);
  if (debug) cout<<"theta_i "<<theta_i/TMath::Pi()*180<<endl;

  float theta_o;
  bool refracted;
  if (fabs(_n1/_n2*sin(theta_i))>1) {
    //internal reflection
    if (debug) cout<<"internal reflection "<<endl;
    theta_o=TMath::Pi()-theta_i;
    refracted=false;
  } else {
    //refraction
    if (debug) cout<<"refraction "<<endl;
    theta_o=asin(_n1/_n2*sin(theta_i));
    refracted=true;
  }
  if (debug) cout<<"theta_o "<<theta_o/TMath::Pi()*180<<endl;

  //  TVector3 vrefracted_dir=u*sin(theta_o)+w*cos(theta_o);
  TVector3 vrefracted_dir=u*sin(theta_o)+w*cos(theta_o);

  _next_dir=vrefracted_dir.Unit();
  _geom->SetCurrentDirection(_next_dir[0],_next_dir[1],_next_dir[2]);

  if (debug){ cout<<"refracted dir "; vrefracted_dir.Print();}

  return refracted;
}


TVector3 generate_gamma(TVector3 vtrk ,double thetac, double phi){
  TVector3 vg(vtrk);
  TVector3 u(vg.Orthogonal());
  vg.Rotate(thetac,u); 
  vg.Rotate(phi,vtrk);
  return vg;
}

void addgamma(TGeoManager* geom, TVector3 pgamma, TVector3 vgamma, bool debug=false){
     //--------------------set initial point and direction--------------------
     double dir_gamma[3]; 
     vgamma.GetXYZ(dir_gamma);
     if (debug) {cout<<"initial gamma direction "; vgamma.Print();}
     double pt_gamma[3];
     pgamma.GetXYZ(pt_gamma);
     geom->InitTrack(pt_gamma,dir_gamma);
     //geom->SetPdgName(22,"gamma");
     geom->AddTrack(1,22);
     TVirtualGeoTrack *gamma = geom->GetLastTrack();
     geom->SetCurrentTrack(gamma);
     gamma->SetPDG(22);
     //     gamma->SetName("gamma");
     gamma->SetLineWidth(0.5);
     gamma->SetMarkerStyle(21);
     gamma->SetMarkerSize(1);
}


void addpoint(TGeoManager *geom, double time, bool debug=false){
       //add point on track
  const double* current_pt=geom->GetCurrentPoint();
       if (debug){
	 cout<<"adding on current trk position "<<current_pt[0]<<" "<<current_pt[1]<<" "<<current_pt[2]<<endl;
       }
       
       geom->GetCurrentTrack()->AddPoint(current_pt[0],current_pt[1],current_pt[2],time);
}



TGeoManager * buildGeometry(tbconfig *config){


   TGeoManager *geom = new TGeoManager("geom","TB 3D GEO");

//------------------Create materials-----------------------------
   TGeoMaterial *matVacuum = new TGeoMaterial("matVacuum",0,0,0);
   TGeoMaterial *matQuartz = new TGeoMaterial("matQuartz",68.0843,30,2.2);

//------------------Create media----------------------------------
   double nair[1]={config->nair};
   TGeoMedium *Air = new TGeoMedium("Air",0,matVacuum,nair);
   double nradiator[1]={config->nradiator};
   TGeoMedium *Quartz = new TGeoMedium("Quartz",1,matQuartz,nradiator);

//------------------Create TOP volume----------------------------
   TGeoVolume *top = geom->MakeBox("top",Air,config->lxbox/2,config->lybox/2,config->lzbox/2);//box
   geom->SetTopVolume(top);
   geom->SetTopVisible(1);
   // If you want to see the boundary, please input the number, 1 instead of 0.
   // Like this, geom->SetTopVisible(1); 


//-----------------Create Radiator volume--------------------------

   TGeoVolume *radiator=geom->MakeBox("radiator",Quartz,config->radheight/2,config->radthick/2,config->radwidth/2);
   radiator->SetFillColor(kGray);
   //   geom->RandomPoints(radiator);
   //   radiator->Raytrace();

   //--------------------Rotate--------------------
   TGeoRotation *rot=new TGeoRotation("rot",0.,config->incly,0.);//rotation around y axis
   TGeoTranslation *tr=new TGeoTranslation("tr",config->pxrad,config->pyrad,config->pzrad);
   TGeoCombiTrans *combi=new TGeoCombiTrans(*tr,*rot);
   top->AddNode(radiator,1,combi);

   geom->CloseGeometry();

   return geom;
}



TTree* tbgeo(double inclination=45., double n=1.47, double thickness=0.1) 
{
  bool debug=false;
  tbconfig *config=new tbconfig;
  config->incly=inclination;
  config->nradiator=n;
  config->radthick=thickness;
  TGeoManager *geom = buildGeometry(config);
   //--------------------prepare ntuple--------------------
  TTree* t=new TTree("xy","xy position of cherenkov photons");
  //position of photon on the box
   double axbox; t->Branch("xbox",&axbox);
   double aybox; t->Branch("ybox",&aybox);
   double azbox; t->Branch("zbox",&azbox);

   double ngammas=100; // number of rays in the cherenkof cone
   //double ngammas=3; // number of rays in the cherenkof cone
   double pi=TMath::Pi();   
   /*
   double pt_trk[3]={0,config->pyrad,0};
   double dir_track[3]={0,1,0}
   geom->InitTrack(pt_trk,dir_track);
   geom->SetPdgName(121,"pi+");
   TVirtualGeoTrack *track = geom->GetLastTrack();
   track->SetPDG(121);
   track->SetLineWidth(0.5);
   //generate cherenkov tracks
   generateCherenkov(geom,track);
   //getListOfTracks
   TObjArray* tracklist=geom->GetListOfTracks();
   //loop on tracks that are gammas
   
   //propagate each track

   */

   double thetac=acos(1./config->nradiator);
   if(debug) cout<<"thetac = "<<thetac*180./TMath::Pi()<<endl;
   for (double phi=0;phi<2*pi-1e-10;phi+=2*pi/ngammas){
     cout<<endl<<phi<<"."<<endl;

     //generate a gamma
     TVector3 dir_g=generate_gamma(config->dir_track,thetac,phi);
     if (debug) {
       cout<<"track direction"<<endl; config->dir_track.Print();
       cout<<"generated gamma dir"<<endl; dir_g.Print();}

     //add a gamma to the list of tracks
     TVector3 pos_g(config->ox,config->oy,config->oz); //point of origin of light
     addgamma(geom,pos_g,dir_g);

     double time=0,dt=1;
     //add a point on the current track
     addpoint(geom,time);

     //propagate till outside
     while (!geom->IsOutside() && (time<10)) {   


       refractor r(geom);
       r.debug=debug;

       //propagate to next surface
       r.propagateToBndry();

       //add a point on the current track
       addpoint(geom,time);

       //refract or reflect direction
       bool refracted=r.refract();

       double next_dir[3];
       r.get_dir().GetXYZ(next_dir);
       
       const double *next_pt=geom->GetCurrentPoint();
       if (debug) {
	 cout<<"next direction    "<<next_dir[0]<<" "<<next_dir[1]<<" "<<next_dir[2]<<endl;
	 cout<<"next position "<<next_pt[0]<<" "<<next_pt[1]<<" "<<next_pt[2]<<endl;
       }
              
       //if reflection then go back into the original volume
       if ( !refracted ) {geom->FindNextBoundary(); geom->Step();}

  
       if (debug)     cout<<endl;
       time+=dt;
     } 

     //print info on final position on boundary box
      printf("final : %s\n", geom->GetPath());
     TGeoVolume *cvol=geom->GetCurrentVolume();
     TGeoMedium *cmed=cvol->GetMedium();

     cout<<"refraction index "<<cmed->GetParam(0)<<endl;
     const double* final_pt=geom->GetCurrentPoint();
     if (fabs(final_pt[0])>1e28 ||fabs(final_pt[1])>1e28 ||fabs(final_pt[2])>1e28) {
       debug=true;
       cout<<"*** found infinite final position"<<endl;
       continue;
     } else {
       debug=false;
     }
     cout<<"final position "<<final_pt[0]<<" "<<final_pt[1]<<" "<<final_pt[2]<<endl;
     const double* final_dir=geom->GetCurrentDirection();
     cout<<"final direction "<<final_dir[0]<<" "<<final_dir[1]<<" "<<final_dir[2]<<endl;
     //add point on track
     geom->GetCurrentTrack()->AddPoint(final_pt[0],final_pt[1],final_pt[2],time);

     axbox=final_pt[0];
     aybox=final_pt[1];
     azbox=final_pt[2];

     if (fabs(axbox)<1e28  && fabs(aybox)<1e28  && fabs(azbox)<1e28) 
       t->Fill();
     else
       cout<<"discarding invalid numbers in ttree "<<axbox<<" "<<aybox<<" "<<azbox<<endl;
   }   
   
   geom->SetVisLevel(4);
   geom->SetVisOption(0);   

   TCanvas *canv=new TCanvas("canv","canv",800,600);
   geom->FindVolumeFast("top")->Draw("ogle"); 
   TView3D *view=(TView3D*)canv->GetView3D(); 
   if (view !=0) view->ShowAxis();
   else cout<<"null view"<<endl;
   geom->DrawTracks("same");
   
   return t;
}



