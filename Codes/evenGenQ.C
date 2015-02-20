// Event generation-single particle generation. Decayer Input is 
//Mass, pT, Y of mother particle and out put is pT1,pT2 and eta1,eta2 
//of decayed particles.
#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TF1.h" 
#include "TGraph.h"
#include "TRandom.h"
#include "TLine.h"
#include <TObjArray.h>
#include <TTree.h>
#include "stdio.h"
using namespace std;

class Particle{

  public:
    Particle();
    Particle(double);
    double   px, py, pz, E, m, pAbs, pt, p[4], eta; 
    void     p4(double, double, double, double);
    void     boost(Particle);
    void     twoBodyDecay(Particle&, Particle&);
    double   InvariantMass(Particle&, Particle&);
    void     print();
};


Particle::Particle(){
  px = py = pz = E = m = pAbs = pt = eta = 0.0;
  p[0] = p[1] = p[2] = p[3] = 0.0;
}

Particle::Particle(double mass){
  m = mass;
  px = py = pz = E = pAbs = pt = eta = 0.0;
  p[0] = p[1] = p[2] = p[3] = 0.0;
}

//*** Sets components of 4-momentum vector -------------------------------------

void Particle::p4(double momx, double momy, double momz, double energy){
  // components of 4-momenta
  px = p[0] = momx;
  py = p[1] = momy;
  pz = p[2] = momz;
  E  = p[3] = energy;
  // transverse momentum and the magnitude of the space momentum
  pt      = sqrt(momx*momx + momy*momy);
  pAbs    = sqrt(momx*momx + momy*momy + momz*momz);

  eta = 0.5*TMath::Log((pAbs+pz)/(pAbs-pz));
}


void Particle::print(){
  cout << "(" << p[0] <<",\t" << p[1] <<",\t"<< p[2] <<",\t"<< p[3] << ")" << endl;
}



//*** Evaluates 4-vector of decay product in the rest frame of parent. ---------
void Particle::twoBodyDecay(Particle& prod1, Particle& prod2) {

   double m1 = prod1.m;
   double m2 = prod2.m;

   // check if particle m can decay
   if( m < m1+m2 ){
     cout << "Particle::twoBodyDecay  parent mass is less than sum of products." 
          << endl;
      prod1.p4 ( 0., 0., 0., m1);
      prod2.p4 ( 0., 0., 0., m2);
      return;
   }

   // CM energies and momentum
   double e1 = (m*m + m1*m1 - m2*m2) / (2.0*m);
   double e2 = (m*m - m1*m1 + m2*m2) / (2.0*m);
   double P  = sqrt(e1*e1 - m1*m1);

   // Isotropic random angles
   //TF1 *angle = new TF1 ("angle","1", -20.0, 20.0);
   //   double phi = angle->GetRandom(0,2*TMath::Pi());
   //   double theta = acos(angle->GetRandom(-1.0,1.0));


   double phi = 2*TMath::Pi()* gRandom->Rndm();
   double theta = acos(2.0*gRandom->Rndm()-1.0);
 
   double pX = P*sin(theta)*cos(phi);
   double pY = P*sin(theta)*sin(phi);
   double pZ = P*cos(theta);

   // set the 4-momenta
   prod1.p4(  pX,  pY,  pZ, e1 );
   prod2.p4( -pX, -pY, -pZ, e2 );
}


double Particle::InvariantMass(Particle& d1, Particle& d2) {
  return sqrt( (d1.E+d2.E)*(d1.E+d2.E) - (d1.px+d2.px)*(d1.px+d2.px) -
           (d1.py+d2.py)*(d1.py+d2.py) - (d1.pz+d2.pz)* (d1.pz+d2.pz));
}



//*** Lorentz Boost  -----------------------------------------------------------
void Particle::boost(Particle parent){

  // beta and gamma values
  double betax = parent.px / parent.E;
  double betay = parent.py / parent.E;
  double betaz = parent.pz / parent.E;
  double beta2 = betax*betax + betay*betay + betaz*betaz;
  double gamma = 1.0/sqrt(1.0-beta2);
  double dot   = betax*px + betay*py + betaz*pz;
  double prod  = gamma*( gamma*dot/(1.0+gamma) + E );

  double pX = px + betax*prod;
  double pY = py + betay*prod;
  double pZ = pz + betaz*prod;
  double e  = gamma*(E + dot);

  p4( pX, pY, pZ, e );
}

void evenGenQ(){
  gStyle->SetPalette(1);

  void TwoParDecay(double, double, double, double&, double&, double&, double&, double&,double&,double&);
  //long i; 
  double Pt1,Pt2,InvMass,eta1,eta2,P1,P2;
  TH1F *MotherRap = new TH1F("MotherRap"," MotherRap",50,-5,5);
  TH1F *MotherPt = new TH1F("MotherPt","MotherPt",50,0,15);
  
  TH1F *MotherMass = new TH1F("MotherMass","MotherMass",100,0,50);
  
  
  TH1F *MotherMassEta8Pt0 = new TH1F("MotherMassEta8Pt0","MotherMassEta8Pt0",100,0,50);
  TH1F *MotherMassEta8Pt1 = new TH1F("MotherMassEta8Pt1","MotherMassEta8Pt1",100,0,50);
  TH1F *MotherMassEta8Pt2 = new TH1F("MotherMassEta8Pt2","MotherMassEta8Pt2",100,0,50);
  TH1F *MotherMassEta8Pt3 = new TH1F("MotherMassEta8Pt3","MotherMassEta8Pt3",100,0,50);
  TH1F *MotherMassEta8PtLT1 = new TH1F("MotherMassEta8PtLT1","MotherMassEta8PtLT1",100,0,50);



  TH1F *MotherMassEta24Pt0 = new TH1F("MotherMassEta24Pt0","MotherMassEta24Pt0",100,0,50);
  TH1F *MotherMassEta24Pt1 = new TH1F("MotherMassEta24Pt1","MotherMassEta24Pt1",100,0,50);
  TH1F *MotherMassEta24Pt2 = new TH1F("MotherMassEta24Pt2","MotherMassEta24Pt2",100,0,50);
  TH1F *MotherMassEta24Pt3 = new TH1F("MotherMassEta24Pt3","MotherMassEta24Pt3",100,0,50);
  TH1F *MotherMassEta24Pt4 = new TH1F("MotherMassEta24Pt4","MotherMassEta24Pt4",100,0,50);
  TH1F *MotherMassEta24Pt10 = new TH1F("MotherMassEta24Pt10","MotherMassEta24Pt10",100,0,50);
  

  TH1F *MotherMassEta4Pt0 = new TH1F("MotherMassEta4Pt0","MotherMassEta4Pt0",100,0,50);
  TH1F *MotherMassEta4Pt1 = new TH1F("MotherMassEta4Pt1","MotherMassEta4Pt1",100,0,50);
  TH1F *MotherMassEta4P4 = new TH1F("MotherMassEta4P4","MotherMassEta4P4",100,0,50);



  TH1F *Dau1Eta = new TH1F("Dau1Eta","Dau1Eta",50,-10,10);
  TH1F *Dau2Eta = new TH1F("Dau2Eta","Dau2Eta",50,-10,10);
  TH1F *Dau1Pt = new TH1F("Dau1Pt","Dau1Pt",50,0,5);
  TH1F *Dau2Pt = new TH1F("Dau2Pt","Dau2Pt",50,0,5);
  TH1F *Dau1P = new TH1F("Dau1P","Dau1P",100,0,50);
  TH1F *Dau2P = new TH1F("Dau2P","Dau2P",100,0,50);
  TH1F *MassOut = new TH1F("MassOut","MassOut",100,0,50);
  
  
  char fileName[500];
  sprintf(fileName,"thermalhist.root");
  cout<<" open file "<<fileName<<endl;
  TFile *infile=new TFile(fileName);
  TH2F *diMuonsInvMassPt = (TH2F*)infile->Get("diMuonsMassPt");
  TH1F *diMuonsRap = (TH1F*)infile->Get("diMuonsRap");
  
  TFile *outfile=new TFile("ThermalDiMuonShdwMass.root","Recreate");
  double MassIntegral=diMuonsInvMassPt->Integral("width"); 
  
  cout<<"MassIntegral: "<<MassIntegral<<endl;
  
  new TCanvas;
  diMuonsInvMassPt->ProjectionX()->Draw();
  
  new TCanvas;
  diMuonsInvMassPt->Draw("lego");

  //return;

  
  // long int NEvents = 100000000;
  long int NEvents = 10000000;
  
  for(long int i=1;i<=NEvents ;i++) {
    Double_t Rap=  diMuonsRap->GetRandom();
    Double_t Mass, Pt;
    diMuonsInvMassPt->GetRandom2(Mass,Pt);
    
    if(Mass<0.25)continue;
    
    if(i%1000000==0)cout<<i<<" mass "<<Mass<<" pT "<<Pt<< " rap "<<Rap<<endl;  
    

    TwoParDecay(Mass, Pt, Rap, Pt1, Pt2, P1, P2, eta1, eta2, InvMass);
    
    MotherPt->Fill(Pt);
    MotherRap->Fill(Rap);
    MotherMass->Fill(InvMass);
    

    //cout<<" mass in "<<Mass<<" mass out "<<InvMass<<endl;
    //cout<<Pt1<<" pt  "<<Pt2<<endl;
    
    if( (TMath::Abs(eta1)<=0.8 && TMath::Abs(eta2)<=0.8) && (Pt1>0.0 && Pt2>0.0) )MotherMassEta8Pt0->Fill(InvMass);
    if( (TMath::Abs(eta1)<=0.8 && TMath::Abs(eta2)<=0.8) && (Pt1>1.0 && Pt2>1.0) )MotherMassEta8Pt1->Fill(InvMass);
    if( (TMath::Abs(eta1)<=0.8 && TMath::Abs(eta2)<=0.8) && (Pt1>2.0 && Pt2>2.0) )MotherMassEta8Pt2->Fill(InvMass);
    if( (TMath::Abs(eta1)<=0.8 && TMath::Abs(eta2)<=0.8) && (Pt1>3.0 && Pt2>3.0) )MotherMassEta8Pt3->Fill(InvMass);
    if( (TMath::Abs(eta1)<=0.8 && TMath::Abs(eta2)<=0.8) && (Pt1<1.0 && Pt2<1.0) )MotherMassEta8PtLT1->Fill(InvMass);
    
    if( (TMath::Abs(eta1)<=2.4 && TMath::Abs(eta2)<=2.4) && (Pt1>0.0 && Pt2>0.0) )MotherMassEta24Pt0->Fill(InvMass);
    if( (TMath::Abs(eta1)<=2.4 && TMath::Abs(eta2)<=2.4) && (Pt1>1.0 && Pt2>1.0) )MotherMassEta24Pt1->Fill(InvMass);
    if( (TMath::Abs(eta1)<=2.4 && TMath::Abs(eta2)<=2.4) && (Pt1>2.0 && Pt2>2.0) )MotherMassEta24Pt2->Fill(InvMass);
    if( (TMath::Abs(eta1)<=2.4 && TMath::Abs(eta2)<=2.4) && (Pt1>3.0 && Pt2>3.0) )MotherMassEta24Pt3->Fill(InvMass);
    if( (TMath::Abs(eta1)<=2.4 && TMath::Abs(eta2)<=2.4) && (Pt1>4.0 && Pt2>4.0) )MotherMassEta24Pt4->Fill(InvMass);
    if( (TMath::Abs(eta1)<=2.4 && TMath::Abs(eta2)<=2.4) && (Pt1>10.0 && Pt2>10.0) )MotherMassEta24Pt10->Fill(InvMass);



    if( (  (eta1>2.5 && eta1<4.0)  && (eta2>2.5 && eta2<4.0)  ) && (Pt1>0.0 && Pt2>0.0) )MotherMassEta4Pt0->Fill(InvMass);
    if( (  (eta1>2.5 && eta1<4.0)  && (eta2>2.5 && eta2<4.0)  ) && (Pt1>1.0 && Pt2>1.0) )MotherMassEta4Pt1->Fill(InvMass);
    if( (  (eta1>2.5 && eta1<4.0)  && (eta2>2.5 && eta2<4.0)  ) && (P1>4.0 && P2>4.0) )MotherMassEta4P4->Fill(InvMass);


  

    Dau1Eta->Fill(eta1);
    
    Dau2Eta->Fill(eta2);
    
    Dau1Pt->Fill(Pt1);
    
    Dau2Pt->Fill(Pt2);
    

    Dau1P->Fill(P1);
    Dau2P->Fill(P2);
    MassOut->Fill(InvMass);
  }
 
  double ScaleThermal= MassIntegral/MotherMass->GetEffectiveEntries();
  ScaleThermal=ScaleThermal/0.00566;
  
  cout<<" entries in mother mass "<<MotherMass->GetEntries()<<endl;
  cout<<" effective entries in mother mass "<<MotherMass->GetEffectiveEntries()<<endl;
  

  cout<<" ScaleThermal  "<<ScaleThermal<<endl;

  
  new TCanvas;
  MotherRap->Scale(ScaleThermal);
  MotherRap->Draw();
  MotherRap->Write();
  
  new TCanvas;
  MotherPt->Scale(ScaleThermal);
  MotherPt->Draw();
  MotherPt->Write();
     
  new TCanvas;
  MassOut->Scale(ScaleThermal);
  MassOut->Draw();
  MassOut->Write();


  new TCanvas;
  MotherMass->Scale(ScaleThermal);
  MotherMass->Draw();
  MotherMass->Write();
  
  cout<<" mother mass after scale "<<MotherMass->Integral()<<endl;
  
 
  
  new TCanvas;
  if(!(MotherMassEta8Pt0->GetEntries()==0)){MotherMassEta8Pt0->Scale(ScaleThermal);}
  MotherMassEta8Pt0->Draw();
  MotherMassEta8Pt0->Write();
 
  new TCanvas;
  MotherMass->Draw();
  MotherMass->Write();
  MotherMassEta8Pt0->SetLineColor(2);
  MotherMassEta8Pt0->Draw("same");


  new TCanvas;
  if(!(MotherMassEta8Pt1->GetEntries()==0)){MotherMassEta8Pt1->Scale(ScaleThermal);}
  MotherMassEta8Pt1->Draw();
  MotherMassEta8Pt1->Write();
  
  new TCanvas;
  if(!(MotherMassEta8Pt2->GetEntries()==0)){MotherMassEta8Pt2->Scale(ScaleThermal);}
  MotherMassEta8Pt2->Draw();
  MotherMassEta8Pt2->Write();

  new TCanvas;
  if(!(MotherMassEta8Pt3->GetEntries()==0)){MotherMassEta8Pt3->Scale(ScaleThermal);}
  MotherMassEta8Pt3->Draw();
  MotherMassEta8Pt3->Write();
  
  new TCanvas;
  if(!(MotherMassEta8PtLT1->GetEntries()==0)){MotherMassEta8PtLT1->Scale(ScaleThermal);}
  MotherMassEta8PtLT1->Draw();
  MotherMassEta8PtLT1->Write();


 new TCanvas;
 if(!(MotherMassEta24Pt0->GetEntries()==0)){MotherMassEta24Pt0->Scale(ScaleThermal);}
 MotherMassEta24Pt0->Draw();
 MotherMassEta24Pt0->Write();
 
  new TCanvas;
  if(!(MotherMassEta24Pt1->GetEntries()==0)){MotherMassEta24Pt1->Scale(ScaleThermal);}
  MotherMassEta24Pt1->Draw();
  MotherMassEta24Pt1->Write();
  
  new TCanvas;
  if(!(MotherMassEta24Pt2->GetEntries()==0)){MotherMassEta24Pt2->Scale(ScaleThermal);}
  MotherMassEta24Pt2->Draw();
  MotherMassEta24Pt2->Write();

  new TCanvas;
  if(!(MotherMassEta24Pt3->GetEntries()==0)){MotherMassEta24Pt3->Scale(ScaleThermal);}
  MotherMassEta24Pt3->Draw();
  MotherMassEta24Pt3->Write();


  new TCanvas;
  if(!(MotherMassEta24Pt4->GetEntries()==0)){MotherMassEta24Pt4->Scale(MassIntegral/MotherMassEta24Pt4->GetEntries());}
  MotherMassEta24Pt4->Draw();
  MotherMassEta24Pt4->Write();


  new TCanvas;
  if(!(MotherMassEta24Pt10->GetEntries()==0)){MotherMassEta24Pt10->Scale(ScaleThermal);}
  MotherMassEta24Pt10->Draw();
  MotherMassEta24Pt10->Write();


  new TCanvas;
  if(!(MotherMassEta4Pt0->GetEntries()==0)){MotherMassEta4Pt0->Scale(ScaleThermal);}
  MotherMassEta4Pt0->Draw();
  MotherMassEta4Pt0->Write();
  
  new TCanvas;
  if(!(MotherMassEta4Pt1->GetEntries()==0)){MotherMassEta4Pt1->Scale(ScaleThermal);}
  MotherMassEta4Pt1->Draw();
  MotherMassEta4Pt1->Write();
  
  new TCanvas;
  if(!(MotherMassEta4P4->GetEntries()==0)){MotherMassEta4P4->Scale(ScaleThermal);}
  MotherMassEta4P4->Draw();
  MotherMassEta4P4->Write();


 
  new TCanvas;
  Dau1Eta->SetLineColor(2);
  Dau1Eta->GetXaxis()->SetTitle("Y1");
  Dau1Eta->GetYaxis()->SetTitle("Entries");
  Dau1Eta->Draw();
  
  new TCanvas;
  Dau2Eta->SetLineColor(2);
  Dau2Eta->GetXaxis()->SetTitle("Y2");
  Dau2Eta->GetYaxis()->SetTitle("Entries");
  Dau2Eta->Draw();
  
  new TCanvas;
  Dau1Pt->SetLineColor(2);
  Dau1Pt->GetXaxis()->SetTitle("PT1");
  Dau1Pt->GetYaxis()->SetTitle("Entries");
  Dau1Pt->Draw();
  
  new TCanvas;
  Dau2Pt->SetLineColor(2);
  Dau2Pt->GetXaxis()->SetTitle("PT2");
  Dau2Pt->GetYaxis()->SetTitle("Entries");
  Dau2Pt->Draw(); 
  
  new TCanvas;
  Dau1P->SetLineColor(2);
  Dau1P->GetXaxis()->SetTitle("P1");
  Dau1P->GetYaxis()->SetTitle("Entries");
  Dau1P->Draw();
  
  new TCanvas;
  Dau2P->SetLineColor(2);
  Dau2P->GetXaxis()->SetTitle("P2");
  Dau2P->GetYaxis()->SetTitle("Entries");
  Dau2P->Draw(); 

 
 
  outfile->Write();
  //outfile->Close();
}


void TwoParDecay(double Mass, double Pt, double Rap,  double& Pt1, double& Pt2, double& P1, double& P2, double& eta1, double& eta2, double& InvMass)

{ 
  double m1 = 0.105658369;

  Particle dimuon(Mass);
  Particle muon1(m1);
  Particle muon2(m1);

  dimuon.twoBodyDecay(muon1, muon2);

  //TF1 *PhiFunc = new TF1 ("PhiFunc","1",0,2*TMath::Pi());
  //double PHI = PhiFunc->GetRandom(0,2*TMath::Pi());
  double PHI = 2*TMath::Pi()* gRandom->Rndm(); 
  double Mt = sqrt(Pt*Pt + Mass*Mass);
  double E  = Mt *TMath::CosH(Rap);
  double Pz  = Mt *TMath::SinH(Rap);
  double Px  = Pt*cos(PHI);
  double Py  = Pt*sin(PHI);
   
  dimuon.p4(Px, Py, Pz, E);

  // boost the muons
  muon1.boost(dimuon);
  muon2.boost(dimuon);

  P1 = muon1.pAbs ;
  P2 = muon2.pAbs ;

  Pt1 = muon1.pt;
  Pt2 = muon2.pt;

  eta1 = muon1.eta;
  eta2 = muon2.eta;

  InvMass = dimuon.InvariantMass(muon1,muon2);

}
