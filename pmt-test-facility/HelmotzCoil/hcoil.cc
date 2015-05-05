/* 
   hcoil is a class representing an Helmotz coil
   with radius, distance from coils and number of turns.
   B field can be calculated at the center or a a point along
   the z axis.

   widehcoil takes into account also the finite width of the 
   coils
 */

#include "TMath.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include <string>
#include "TLegend.h"

class hcoil{
public:
  hcoil(float radius,int turns); 
  hcoil(float radius,int turns, float distance); 
  virtual float B_center(float current);
  virtual float B_axis(float current, float d_center);
  void draw(float zmin, float zmax, float current, const char* opt="APL");
  TGraph* getgraph(){return _Bvsz;}
  virtual ~hcoil(){}
private:
  float _r;
  int _n;
  float _d;
  TGraph* _Bvsz;
};
hcoil::hcoil(float radius, int turns){
  _r=radius;
  _n=turns;
  _d=radius;  
  _Bvsz=nullptr;
}
hcoil::hcoil(float radius, int turns, float distance){
  _r=radius;
  _n=turns;
  _d=distance;
  _Bvsz=nullptr;
}
float hcoil::B_center(float current){
  return 32*TMath::Pi()*_n*current*1e-7/(5*sqrt(5)*_r);
  //  return b(current,0);
}
float hcoil::B_axis(float current, float z_center){
  double mu0=4*TMath::Pi()*1e-7;
  return mu0*_n*current/_r*0.5*
    (1./pow( 1+pow(z_center+_d/2,2)/pow(_r,2) , 1.5)+
     1./pow( 1+pow(z_center-_d/2,2)/pow(_r,2) , 1.5));
}
void hcoil::draw(float zmin, float zmax, float current, const char* opt){
  //  TCanvas * canv = new TCanvas("canv","canv",800,600);
  //  canv->Draw();
  const int n=1000;
  float z[n],b[n];
  float zstep=(zmax-zmin)/n;
  for (int i=0;i<n;i++){
    z[i]=zmin+i*zstep+zstep/2;
    b[i]=B_axis(current,z[i])*1e4;
  }
  TGraph *gr= new TGraph(n,z,b);
  gr->GetXaxis()->SetTitle("z[m]");
  gr->GetYaxis()->SetTitle("B[Gauss]");
  gr->SetName("bvsz");
  gr->SetTitle("Helmotz coil B field versus distance from center");
  gr->Draw(opt);
  _Bvsz=gr;
  

}


class widehcoil: public hcoil{
public:
  widehcoil(float radius,int turns, float distance,float width); 
  virtual float B_center(float current);
  virtual float B_axis(float current, float d_center);
  virtual ~widehcoil(){}
private:
  float _width;
};
widehcoil::widehcoil(float radius,int turns, float distance,float width):
  hcoil(radius,turns,distance){
  _width=width;
}
float widehcoil::B_center(float current){
  return widehcoil::B_axis(current,0);
}
float widehcoil::B_axis(float current,float d_center){
  const int n=1000;
  //  float dz=_width/n;
  float dz=_d/n;
  float z=-_width/2;
  float B=hcoil::B_axis(current,z+d_center);
  //  for (int i=0;i<n;i++){
  int i=1;
  while(z<_width/2){
    z=-_width/2+dz*i;
    B+=hcoil::B_axis(current,z+d_center);
    i++;
  }
  B=B/i;
  return B;
}

				 
void comparewidth(float r=0.2, float distance=0.2){
  //coil configuration
  int turns=36;
  float current=5;

  //graph limits
  float zmin=-0.05, zmax=0.05;

  //loop with varying width
  float wmin=0.,wmax=0.1; const int n=10;
  float step=(wmax-wmin)/n;
    TLegend *l=new TLegend(0.7,0.5,0.9,0.9,"widht[m]","NDC");
  for (int i=0;i<=n;i++){
    float width=wmin+i*step;
    widehcoil hc(r,turns,distance,width);
    if (i==0) hc.draw(zmin,zmax,current);
    else hc.draw(zmin,zmax,current,"PL same");
    int color=i+1;
    hc.getgraph()->SetLineColor(color);
    hc.getgraph()->SetMarkerColor(color);
    char ci[6]; sprintf(ci,"%2.2f",width);
    l->AddEntry(hc.getgraph(),ci,"l");
  }
  l->Draw();
}

void comparedistance(float r=0.2, float width=0.05){
  //coil configuration
  int turns=36;
  float current=5;

  //graph limits
  float zmin=-0.05, zmax=0.05;

  //loop with varying distance
  float dmin=r-0.005,dmax=r+0.005; const int n=10;
  float step=(dmax-dmin)/n;
    TLegend *l=new TLegend(0.7,0.5,0.9,0.9,"distance[m]","NDC");
  for (int i=0;i<=n;i++){
    float distance=dmin+i*step;
    widehcoil hc(r,turns,distance,width);
    if (i==0) hc.draw(zmin,zmax,current);
    else hc.draw(zmin,zmax,current,"PL same");
    int color=i+1;
    hc.getgraph()->SetLineColor(color);
    hc.getgraph()->SetMarkerColor(color);
    char ci[7]; sprintf(ci,"%2.3f",distance);
    l->AddEntry(hc.getgraph(),ci,"l");
  }
  l->Draw();
}

