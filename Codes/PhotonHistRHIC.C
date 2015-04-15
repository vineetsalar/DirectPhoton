// c++ classes
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <fstream>
#include <ctime>
// ROOT classes
#include "TROOT.h"
#include "TStyle.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TTree.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TF1.h" 
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include <TObjArray.h>
#include <TTree.h>
#include "stdio.h"
#include "TVector3.h"
#include <TLegend.h>
#include "TCanvas.h"
#include <TPad.h>
#include "TFile.h"
#include "TRandom.h"
#include "TLine.h"
#include "TSystem.h"


Double_t pi = TMath::Pi();
Double_t pi2 = pi*pi;
Double_t pi3 = pi*pi2;
Double_t pi4 = pi2*pi2;
Double_t hbarc = 0.197327;
Double_t hbarc2 = hbarc*hbarc;
Double_t hbarc3 = hbarc*hbarc*hbarc;
Double_t hbarc4 = hbarc*hbarc*hbarc*hbarc;

Double_t AA = 197.0;//Au
Double_t Radius = 1.2*TMath::Power(AA, 1/3.);
Double_t contime = Radius*Radius*pi;
Double_t conrate = 1.0/(137.0*137.0*2.0*pi3);

const Double_t mPi = 0.140; 
const Double_t RAu = 1.2*TMath::Power(AA, 1/3.);
const Double_t R03 = 0.95*RAu;

const Double_t tau0 = 0.3;

const Double_t SS=5.0*1.5*700;

const Double_t Nf=2.5;
const Double_t aq = (7.0*Nf/60.0 + 16.0/90.0)*pi2;
const Double_t ah = 4.5*pi2/90.0;

//Double_t aT = 0.1;
Double_t aT = 0.1;
Double_t z0=0.0; //0
Double_t vZ=1.0;     //1.0

const Double_t VTau0 = (R03+0.5*aT*tau0*tau0)*(R03+0.5*aT*tau0*tau0)*(z0+vZ*tau0)*pi;
const Double_t ss03 = SS/VTau0;


const Double_t T0=TMath::Power(SS/(4.0*aq*VTau0),1.0/3.0)*hbarc;
const Double_t TC = 0.170;
const Double_t TF = 0.140;

const Double_t nPart03 = 358;
//const Double_t nColl0 = ;

//Double_t EtaStart = -1;
//Double_t EtaEnd = 1;

//==================================== Lattice EOS ================================================//
TFile *fileEOS=new TFile("LatticeEOS_s95p-v1.2.root","R");
TGraph *grfSSVsTemp = (TGraph*)fileEOS->Get("grfSSVsTemp");
TGraph *TempVsFQGP = (TGraph*)fileEOS->Get("TempVsFQGP");
TGraph *TempVsFQGP2 = (TGraph*)fileEOS->Get("TempVsFQGP2");




Double_t Npart(int BinLow, int BinHigh);
Double_t NColl(int BinLow, int BinHigh);
Double_t CalculateTandf_LatticeEOS( Double_t ssCent, Double_t R0Cent);




// Temp and Hadronic fraction as a function of time 
Double_t tau[10000], Temp[10000], h[10000];
Int_t Ntime;
double Steptime;


//Direct Photon
//Double_t RatePhoton(Double_t RCent, Double_t Pt);
Double_t RateQGP_IntTau(Double_t Pt, Double_t RCent);
Double_t RateQGP_IntEtaf(Double_t Pt, Double_t T);
Double_t RateQGP(Double_t Pt, Double_t T, Double_t Etaf);


Double_t RateHadron_IntTau(Double_t Pt, Double_t RCent);
Double_t RateHadron_IntEtaf(Double_t Pt, Double_t T);
Double_t RateHadron(Double_t Pt, Double_t T, Double_t Etaf);

Double_t RateQCD_IntTau(Double_t Pt);
Double_t RateQCD(Double_t Pt);


//======= Data Graphs =========================//
void Draw_PHENIX_DirectPhotonRate020_Pt(TLegend *lgd);
void Draw_PHENIX_DirectPhotonRate2040_Pt(TLegend *lgd);  
void Draw_PHENIX_DirectPhotonRate4060_Pt(TLegend *lgd);  
void Draw_PHENIX_DirectPhotonRate6092_Pt(TLegend *lgd);  


void Draw_AllDataGraphs();



void PhotonHistRHIC()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  //gStyle->SetPadTopMargin(0.10);
  //gStyle->SetPadBottomMargin(0.10);
  gStyle->SetPadLeftMargin(0.19);
  //gStyle->SetPadRightMargin(0.15);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0);
  gStyle->SetHistLineWidth(2.0);
  gStyle->SetCanvasDefH(600);
  gStyle->SetCanvasDefW(600);  
  
  // ================== Checking Data Graphs ===============================//
  cout<<" ========== checking data graphs =================================="<<endl;
  //Draw_AllDataGraphs();

  // return;

  TFile *OutFile = new TFile("RHIC_DirecPhoton.root","Recreate");


  cout<<" simulating QGP evolution : "<<endl;

  // ================ dn/deta graph for making Temp as a function of nPart ==========================================//
  Double_t NPartdNdEta[11] = {358.0,331.0,298.0,256.0,217.0,183.0,152.0,124.0,103.0,83.0,65.0};
  Double_t Err_NPartdNdEta[11] = {12.0,10.0,9.0,8.0,8.0,7.0,6.0,6.0,5.0,5.0,4.0};
  
  Double_t dNdEtabyNpartby2[11] = {3.91,3.81,3.68,3.56,3.48,3.41,3.32,3.25,3.19,3.14,2.90};
  Double_t Err_dNdEtabyNpartby2[11] = {0.20,0.18,0.18,0.18,0.18,0.18,0.19,0.19,0.21,0.22,0.22}; 
  
  
  
  TGraphErrors *grdNDetaNpart = new TGraphErrors(11,NPartdNdEta,dNdEtabyNpartby2,Err_NPartdNdEta,Err_dNdEtabyNpartby2);
  grdNDetaNpart->SetLineWidth(2);
  grdNDetaNpart->SetMarkerStyle(20);
  grdNDetaNpart->GetYaxis()->SetRangeUser(2.5,5.0);
  grdNDetaNpart->GetYaxis()->SetTitle("#frac{dN}{d#eta}/(#frac{N_{Part}}{2})");
  grdNDetaNpart->GetXaxis()->SetTitle("N_{Part}");
  
  
  // Calculate Temp and Hadronic fraction as a function of time
  double TauLatt[10000], TempTauLatt[10000];
  double fQGPLatt[10000];

  Double_t R0020 = R03*TMath::Power( (Npart(0,5)/nPart03) ,0.5);
  Double_t ss020 = ss03*(grdNDetaNpart->Eval(Npart(0,5))/dNdEtabyNpartby2[0]);

  Ntime=0;
  CalculateTandf_LatticeEOS(ss020, R0020);

  //cout<<" Ntime "<<Ntime <<" ss010 "<<ss010<<"  R0010: "<<R0010<<endl;
  
  for(int i=0;i<Ntime;i++){
    TauLatt[i]=tau[i];
    TempTauLatt[i]=Temp[i];
    fQGPLatt[i]=(1.0-h[i]);
    //cout<<tau[i]<<" Lattice  "<<Temp[i]<<endl;
  }
    
  TGraph *grTempVsTauAnaC = new TGraph(Ntime,TauLatt,TempTauLatt);
  grTempVsTauAnaC->SetName("grTempVsTauLatt");
  grTempVsTauAnaC->SetTitle("grTempVsTauLatt");
  grTempVsTauAnaC->SetLineColor(1);
  grTempVsTauAnaC->SetLineWidth(2);
  grTempVsTauAnaC->SetLineStyle(9);
  grTempVsTauAnaC->GetYaxis()->SetTitle("Temperature (GeV)");
  grTempVsTauAnaC->GetXaxis()->SetTitle("#tau (fm)");
  new TCanvas;
  grTempVsTauAnaC->Draw("AL");
  

  TGraph *grFQGPVsTauAnaC = new TGraph(Ntime,TauLatt,fQGPLatt);
  grFQGPVsTauAnaC->SetName("grFQGPVsTauLatt");
  grFQGPVsTauAnaC->SetTitle("grFQGPVsTauLatt");
  grFQGPVsTauAnaC->SetLineColor(1);
  grFQGPVsTauAnaC->SetLineWidth(2);
  grFQGPVsTauAnaC->SetLineStyle(9);
  grFQGPVsTauAnaC->GetYaxis()->SetTitle("FQGP");
  grFQGPVsTauAnaC->GetXaxis()->SetTitle("#tau (fm)");

  new TCanvas;
  grFQGPVsTauAnaC->Draw("AL");



  //=============================== PQCD Direct Photon ===================================//
  
  Double_t tPt[1000]={0.0};
  Double_t DirectPhotonPQCD[1000]={0.0};

  Double_t tPtMin =0.1;
  Double_t tPtMax =12.0;
  Double_t tPtStep =0.1;
  Int_t tNPt = (tPtMax - tPtMin)/tPtStep;
  for(int i =0;i<tNPt;i++)
    {
      tPt[i]=tPtMin+i*tPtStep;
      DirectPhotonPQCD[i] = RateQCD(tPt[i]);
    }


  TGraph *grfDirectPhotonPQCD = new TGraph(tNPt,tPt,DirectPhotonPQCD);
  grfDirectPhotonPQCD->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonPQCD->GetYaxis()->SetTitle("#frac{1}{2#pi p_{T}}#frac{d^{2}#sigma}{dp_{T}dy}[(mb.GeV/c)^{-2}]");
  grfDirectPhotonPQCD->GetYaxis()->SetTitleOffset(1.8);
  grfDirectPhotonPQCD->SetLineColor(4);


  new TCanvas;
  gPad->SetLogy();
  gPad->SetTicks();
  grfDirectPhotonPQCD->Draw("AL");




  
  //Double_t R0020 = R03*TMath::Power( (Npart(0,5)/nPart03) ,0.5);
  //Double_t ss020 = ss03*(grdNDetaNpart->Eval(Npart(0,5))/dNdEtabyNpartby2[0]);
  
  //Ntime=0;
  //CalculateTandf_LatticeEOS(ss020, R0020);
  
  //0-20-40-60-92
  Double_t NCollPhenix[4]={779.0,297.0,91.0,15.0};
  Double_t Sigmapp = 42.2; //mb inelastic

  //Double_t NPartPhenix[4]={280.0,140.0,60.0,18.0};
  //Double_t NCollPhenix[4]={1.0,1.0,1.0,1.0};

  Double_t PtMin =0.1;
  Double_t PtMax =6.0;
  Double_t PtStep =0.1;
  Int_t NPt = (PtMax - PtMin)/PtStep;
  
  Double_t Pt[100]={0.0};
  Double_t DirectPhoton1[100]={0.0};
  Double_t DirectPhotonQGP1[100]={0.0};
  Double_t DirectPhotonHadron1[100]={0.0};
  Double_t DirectPhotonPQCD1[100]={0.0};

  cout<<"Pt:  "<<"         "<<"Hadron      "<<"QGP       "<<"RatePhoton"<<endl;

  cout<<"R0020  "<<R0020<<endl;
 
  for(int i =0;i<NPt;i++)
    {
      Pt[i]=PtMin+i*PtStep;
      
      DirectPhotonQGP1[i] = RateQGP_IntTau(Pt[i],R0020);
      DirectPhotonHadron1[i] = RateHadron_IntTau(Pt[i],R0020);
      DirectPhotonPQCD1[i] = NCollPhenix[0]*RateQCD(Pt[i])/Sigmapp;
      DirectPhoton1[i] =  DirectPhotonPQCD1[i] + DirectPhotonQGP1[i] + DirectPhotonHadron1[i];

      cout<<Pt[i]<<"    "<<DirectPhotonHadron1[i]<<"    "<<DirectPhotonQGP1[i]<<"   "<<DirectPhoton1[i]<<endl;
    }  


  TGraph *grfDirectPhotonQGP1 = new TGraph(NPt,Pt,DirectPhotonQGP1);
  grfDirectPhotonQGP1->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonQGP1->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonQGP1->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhotonQGP1->SetLineColor(2);


  TGraph *grfDirectPhotonHadron1 = new TGraph(NPt,Pt,DirectPhotonHadron1);
  grfDirectPhotonHadron1->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonHadron1->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonHadron1->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhotonHadron1->SetLineColor(1);



  TGraph *grfDirectPhotonPQCD1 = new TGraph(NPt,Pt,DirectPhotonPQCD1);
  grfDirectPhotonPQCD1->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonPQCD1->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonPQCD1->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhotonPQCD1->SetLineColor(4);



  TGraph *grfDirectPhoton1 = new TGraph(NPt,Pt,DirectPhoton1);
  grfDirectPhoton1->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhoton1->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhoton1->GetYaxis()->SetTitleOffset(1.4);
    
  TLegend *legd5 = new TLegend( 0.60,0.70,0.82,0.85);
  legd5->SetBorderSize(0);
  legd5->SetFillStyle(0);
  legd5->SetFillColor(0);
  legd5->SetTextSize(0.040);

  legd5->AddEntry(grfDirectPhoton1,"Sum","L");
  legd5->AddEntry(grfDirectPhotonQGP1,"QGP","L");
  legd5->AddEntry(grfDirectPhotonHadron1,"Hadron","L");
  legd5->AddEntry(grfDirectPhotonPQCD1,"pp*N_{coll}","L");

  new TCanvas;
  gPad->SetTicks();
  gPad->SetLogy();
  Draw_PHENIX_DirectPhotonRate020_Pt(legd5);
  grfDirectPhoton1->SetLineWidth(2);
  grfDirectPhoton1->SetLineColor(4);
  grfDirectPhoton1->Draw("Lsame");
  grfDirectPhotonHadron1->Draw("Lsame");
  grfDirectPhotonQGP1->Draw("Lsame");
  grfDirectPhotonPQCD1->Draw("Lsame"); 
  legd5->Draw("Lsame");



  //return;
  //====================================== 20-40% ======================================//
  //0-3-6-10-15-20-25-30-35-40-45-50-55-60-65-70-75-80-100 
  Double_t R02040 = R03*TMath::Power( (Npart(5,9)/nPart03) ,0.5);
  Double_t ss2040 = ss03*(grdNDetaNpart->Eval(Npart(5,9))/dNdEtabyNpartby2[0]);
  Ntime=0;
  CalculateTandf_LatticeEOS(ss2040, R02040);
  
  Double_t DirectPhoton2[100]={0.0};
  Double_t DirectPhotonQGP2[100]={0.0};
  Double_t DirectPhotonPQCD2[100]={0.0};
  Double_t DirectPhotonHadron2[100]={0.0};
   

  for(int i =0;i<NPt;i++)
    {
      Pt[i]=PtMin+i*PtStep;
      
      DirectPhotonQGP2[i] = RateQGP_IntTau(Pt[i],R02040);
      DirectPhotonHadron2[i] = RateHadron_IntTau(Pt[i],R02040);
      DirectPhotonPQCD2[i] = NCollPhenix[1]*RateQCD(Pt[i])/Sigmapp;

      DirectPhoton2[i] =  DirectPhotonPQCD2[i] + DirectPhotonQGP2[i] + DirectPhotonHadron2[i];
      
      cout<<Pt[i]<<"    "<<DirectPhotonHadron2[i]<<"    "<<DirectPhotonQGP2[i]<<"   "<<DirectPhoton2[i]<<endl;
    }  


  TGraph *grfDirectPhotonQGP2 = new TGraph(NPt,Pt,DirectPhotonQGP2);
  grfDirectPhotonQGP2->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonQGP2->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonQGP2->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhotonQGP2->SetLineColor(2);


  TGraph *grfDirectPhotonHadron2 = new TGraph(NPt,Pt,DirectPhotonHadron2);
  grfDirectPhotonHadron2->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonHadron2->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonHadron2->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhotonHadron2->SetLineColor(1);

  TGraph *grfDirectPhotonPQCD2 = new TGraph(NPt,Pt,DirectPhotonPQCD2);
  grfDirectPhotonPQCD2->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonPQCD2->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonPQCD2->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhotonPQCD2->SetLineColor(4);




  TGraph *grfDirectPhoton2 = new TGraph(NPt,Pt,DirectPhoton2);
  grfDirectPhoton2->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhoton2->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhoton2->GetYaxis()->SetTitleOffset(1.4);
    
  

  new TCanvas;
  gPad->SetTicks();
  gPad->SetLogy();
  Draw_PHENIX_DirectPhotonRate2040_Pt(legd5);
  grfDirectPhoton2->SetLineWidth(2);
  grfDirectPhoton2->SetLineColor(4);
  grfDirectPhoton2->Draw("Lsame");
  grfDirectPhotonHadron2->Draw("Lsame");
  grfDirectPhotonQGP2->Draw("Lsame");
  grfDirectPhotonPQCD2->Draw("Lsame");
  legd5->Draw("Lsame");

  
  //====================================== 40-60% ======================================//
  //0-3-6-10-15-20-25-30-35-40-45-50-55-60-65-70-75-80-100 
  Double_t R04060 = R03*TMath::Power( (Npart(9,13)/nPart03) ,0.5);
  Double_t ss4060 = ss03*(grdNDetaNpart->Eval(Npart(9,13))/dNdEtabyNpartby2[0]);
  Ntime=0;
  CalculateTandf_LatticeEOS(ss4060, R04060);
  

  Double_t DirectPhoton3[100]={0.0};
  Double_t DirectPhotonQGP3[100]={0.0};
  Double_t DirectPhotonHadron3[100]={0.0};
  Double_t DirectPhotonPQCD3[100]={0.0};

    

  for(int i =0;i<NPt;i++)
    {
      Pt[i]=PtMin+i*PtStep;
      
      DirectPhotonQGP3[i] = RateQGP_IntTau(Pt[i],R04060);
      DirectPhotonHadron3[i] = RateHadron_IntTau(Pt[i],R04060);
      DirectPhotonPQCD3[i] = NCollPhenix[2]*RateQCD(Pt[i])/Sigmapp;
      DirectPhoton3[i] =  DirectPhotonPQCD3[i] + DirectPhotonQGP3[i] + DirectPhotonHadron3[i];
      
      cout<<Pt[i]<<"    "<<DirectPhotonHadron3[i]<<"    "<<DirectPhotonQGP3[i]<<"   "<<DirectPhoton3[i]<<endl;
    }  


  TGraph *grfDirectPhotonQGP3 = new TGraph(NPt,Pt,DirectPhotonQGP3);
  grfDirectPhotonQGP3->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonQGP3->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonQGP3->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhotonQGP3->SetLineColor(2);


  TGraph *grfDirectPhotonHadron3 = new TGraph(NPt,Pt,DirectPhotonHadron3);
  grfDirectPhotonHadron3->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonHadron3->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonHadron3->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhotonHadron3->SetLineColor(1);

  TGraph *grfDirectPhotonPQCD3 = new TGraph(NPt,Pt,DirectPhotonPQCD3);
  grfDirectPhotonPQCD3->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonPQCD3->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonPQCD3->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhotonPQCD3->SetLineColor(4);


  TGraph *grfDirectPhoton3 = new TGraph(NPt,Pt,DirectPhoton3);
  grfDirectPhoton3->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhoton3->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhoton3->GetYaxis()->SetTitleOffset(1.4);
    
  

  new TCanvas;
  gPad->SetTicks();
  gPad->SetLogy();
  Draw_PHENIX_DirectPhotonRate4060_Pt(legd5);
  grfDirectPhoton3->SetLineWidth(2);
  grfDirectPhoton3->SetLineColor(4);
  grfDirectPhoton3->Draw("Lsame");
  grfDirectPhotonHadron3->Draw("Lsame");
  grfDirectPhotonQGP3->Draw("Lsame");
  grfDirectPhotonPQCD3->Draw("Lsame");
  legd5->Draw("Lsame");



//====================================== 60-100% ======================================//
  //0-3-6-10-15-20-25-30-35-40-45-50-55-60-65-70-75-80-100 
  Double_t R060100 = R03*TMath::Power( (Npart(13,18)/nPart03) ,0.5);
  Double_t ss60100 = ss03*(grdNDetaNpart->Eval(Npart(13,18))/dNdEtabyNpartby2[0]);
  Ntime=0;
  CalculateTandf_LatticeEOS(ss60100, R060100);
  

  Double_t DirectPhoton4[100]={0.0};
  Double_t DirectPhotonQGP4[100]={0.0};
  Double_t DirectPhotonPQCD4[100]={0.0};
  Double_t DirectPhotonHadron4[100]={0.0};
    

  for(int i =0;i<NPt;i++)
    {
      Pt[i]=PtMin+i*PtStep;
      
      DirectPhotonQGP4[i] = RateQGP_IntTau(Pt[i],R060100);
      DirectPhotonHadron4[i] = RateHadron_IntTau(Pt[i],R060100);
      DirectPhotonPQCD4[i] = NCollPhenix[3]*RateQCD(Pt[i])/Sigmapp;

      DirectPhoton4[i] =  DirectPhotonPQCD4[i] + DirectPhotonQGP4[i] + DirectPhotonHadron4[i];
      
      cout<<Pt[i]<<"    "<<DirectPhotonHadron4[i]<<"    "<<DirectPhotonQGP4[i]<<"   "<<DirectPhoton4[i]<<endl;
    }  


  TGraph *grfDirectPhotonQGP4 = new TGraph(NPt,Pt,DirectPhotonQGP4);
  grfDirectPhotonQGP4->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonQGP4->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonQGP4->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhotonQGP4->SetLineColor(2);


  TGraph *grfDirectPhotonHadron4 = new TGraph(NPt,Pt,DirectPhotonHadron4);
  grfDirectPhotonHadron4->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonHadron4->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonHadron4->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhotonHadron4->SetLineColor(1);

  TGraph *grfDirectPhotonPQCD4 = new TGraph(NPt,Pt,DirectPhotonPQCD3);
  grfDirectPhotonPQCD4->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonPQCD4->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonPQCD4->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhotonPQCD4->SetLineColor(4);


  TGraph *grfDirectPhoton4 = new TGraph(NPt,Pt,DirectPhoton4);
  grfDirectPhoton4->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhoton4->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhoton4->GetYaxis()->SetTitleOffset(1.4);
    
  

  new TCanvas;
  gPad->SetTicks();
  gPad->SetLogy();
  Draw_PHENIX_DirectPhotonRate6092_Pt(legd5);
  grfDirectPhoton4->SetLineWidth(2);
  grfDirectPhoton4->SetLineColor(4);
  grfDirectPhoton4->Draw("Lsame");
  grfDirectPhotonHadron4->Draw("Lsame");
  grfDirectPhotonQGP4->Draw("Lsame");
  grfDirectPhotonPQCD4->Draw("Lsame");
  legd5->Draw("Lsame");

  

  grTempVsTauAnaC->SetName("grTempVsTauAnaC");
  grFQGPVsTauAnaC->SetName("grFQGPVsTauAnaC");
  grfDirectPhotonPQCD->SetName("grfDirectPhotonPQCD");
  
  grfDirectPhoton1->SetName("grfDirectPhoton1");
  grfDirectPhotonHadron1->SetName("grfDirectPhotonHadron1");
  grfDirectPhotonQGP1->SetName("grfDirectPhotonQGP1");
  grfDirectPhotonPQCD1->SetName("grfDirectPhotonPQCD1"); 


  grfDirectPhoton2->SetName("grfDirectPhoton2");
  grfDirectPhotonHadron2->SetName("grfDirectPhotonHadron2");
  grfDirectPhotonQGP2->SetName("grfDirectPhotonQGP2");
  grfDirectPhotonPQCD2->SetName("grfDirectPhotonPQCD2"); 


  grfDirectPhoton3->SetName("grfDirectPhoton3");
  grfDirectPhotonHadron3->SetName("grfDirectPhotonHadron3");
  grfDirectPhotonQGP3->SetName("grfDirectPhotonQGP3");
  grfDirectPhotonPQCD3->SetName("grfDirectPhotonPQCD3"); 


  grfDirectPhoton4->SetName("grfDirectPhoton4");
  grfDirectPhotonHadron4->SetName("grfDirectPhotonHadron4");
  grfDirectPhotonQGP4->SetName("grfDirectPhotonQGP4");
  grfDirectPhotonPQCD4->SetName("grfDirectPhotonPQCD4"); 
  
  

 grTempVsTauAnaC->SetTitle("grTempVsTauAnaC");
  grFQGPVsTauAnaC->SetTitle("grFQGPVsTauAnaC");
  grfDirectPhotonPQCD->SetTitle("grfDirectPhotonPQCD");
  
  grfDirectPhoton1->SetTitle("grfDirectPhoton1");
  grfDirectPhotonHadron1->SetTitle("grfDirectPhotonHadron1");
  grfDirectPhotonQGP1->SetTitle("grfDirectPhotonQGP1");
  grfDirectPhotonPQCD1->SetTitle("grfDirectPhotonPQCD1"); 


  grfDirectPhoton2->SetTitle("grfDirectPhoton2");
  grfDirectPhotonHadron2->SetTitle("grfDirectPhotonHadron2");
  grfDirectPhotonQGP2->SetTitle("grfDirectPhotonQGP2");
  grfDirectPhotonPQCD2->SetTitle("grfDirectPhotonPQCD2"); 


  grfDirectPhoton3->SetTitle("grfDirectPhoton3");
  grfDirectPhotonHadron3->SetTitle("grfDirectPhotonHadron3");
  grfDirectPhotonQGP3->SetTitle("grfDirectPhotonQGP3");
  grfDirectPhotonPQCD3->SetTitle("grfDirectPhotonPQCD3"); 


  grfDirectPhoton4->SetTitle("grfDirectPhoton4");
  grfDirectPhotonHadron4->SetTitle("grfDirectPhotonHadron4");
  grfDirectPhotonQGP4->SetTitle("grfDirectPhotonQGP4");
  grfDirectPhotonPQCD4->SetTitle("grfDirectPhotonPQCD4"); 
  

 
  grTempVsTauAnaC->Write();
  grFQGPVsTauAnaC->Write();
  grfDirectPhotonPQCD->Write();
  grfDirectPhoton1->Write();
  grfDirectPhotonHadron1->Write();
  grfDirectPhotonQGP1->Write();
  grfDirectPhotonPQCD1->Write(); 
  grfDirectPhoton2->Write();
  grfDirectPhotonHadron2->Write();
  grfDirectPhotonQGP2->Write();
  grfDirectPhotonPQCD2->Write();
  grfDirectPhoton3->Write();
  grfDirectPhotonHadron3->Write();
  grfDirectPhotonQGP3->Write();
  grfDirectPhotonPQCD3->Write();
  grfDirectPhoton4->Write();
  grfDirectPhotonHadron4->Write();
  grfDirectPhotonQGP4->Write();
  grfDirectPhotonPQCD4->Write();


  OutFile->Close();



  return;



}



//This is the function for time vs temp


//////////

Double_t CalculateTandf_LatticeEOS(Double_t ssCent, Double_t R0Cent)
{
  Steptime = 0.1;
  Ntime=0;
  Double_t CutTemp =0.0;
  do{
    tau[Ntime] = tau0 + Steptime*Ntime;

    Double_t VTau0Cent =  (R0Cent+0.5*aT*tau0*tau0)*(R0Cent+0.5*aT*tau0*tau0)*(z0+vZ*tau0)*pi;
    Double_t VTauCent =  (R0Cent+0.5*aT*tau[Ntime]*tau[Ntime])*(R0Cent+0.5*aT*tau[Ntime]*tau[Ntime])*(z0+vZ*tau[Ntime])*pi;
    
    Double_t ssTau = (hbarc3*ssCent)*(VTau0Cent/VTauCent);
    
    Temp[Ntime]=grfSSVsTemp->Eval(ssTau);
    Double_t tt= Temp[Ntime];
    Double_t fqgp_temp = TempVsFQGP2->Eval(tt);
    h[Ntime]= (1.0-fqgp_temp);
    CutTemp = Temp[Ntime];
    
    cout<<"Lattice: Ntime "<<Ntime<<" ss "<<ssTau<<" tau "<<tau[Ntime]<<" T  "<<Temp[Ntime]<<" h  "<<h[Ntime]<<endl;
   
    Ntime++;
  }while(CutTemp > TF);
  
  return 0;

}


//0-3-6-10-15-20-25-30-35-40-45-50-55-60-65-70-75-80-100
Double_t Npart(int BinLow, int BinHigh)
{
  //Double_t NpartArray[40]={358.0,331.0,298.0,256.0,217.0,183.0,152.0,124.0,103.0,83.0,65.0};
  Double_t NpartArray[40]={358.0,331.0,298.0,256.0,217.0,183.0,152.0,124.0,103.0,83.0,65.0,50.0,39.0,29.0,21.0,
			   15.0,11.0,4.1};
  Double_t sum=0;
  for(int i=BinLow;i<BinHigh;i++){
    sum+=NpartArray[i];
  }
  Double_t NPart = sum/(BinHigh-BinLow);
  return NPart;
}


//0-3-6-10-15-20-25-30-35-40-45-50
//These are dummy values
Double_t NColl(int BinLow, int BinHigh)
{
  Double_t NCollArray[40]={358.0,331.0,298.0,256.0,217.0,183.0,152.0,124.0,103.0,83.0,65.0};
  
  Double_t sum=0;
  for(int i=BinLow;i<BinHigh;i++){
    sum+=NCollArray[i];
  }
  Double_t NColl = sum/(BinHigh-BinLow);
  return NColl;
}



//================================ Direct Photon Functions =======================================//


/*
Double_t RatePhoton(Double_t RCent, Double_t Pt)
{

  Double_t Area = pi*RCent*RCent;
  Double_t RPhoton =0.0;
  RPhoton =Area*(RateQGP_IntTau(Pt) + RateHadron_IntTau(Pt));

  return RPhoton;
}
*/


Double_t RateQGP_IntTau(Double_t Pt, Double_t RCent)
{

  Double_t Sum =0.0;
  Double_t VTau =0.0;

  for(int i =0;i<Ntime;i++)
    {

      VTau = (RCent+0.5*aT*tau[i]*tau[i])*(RCent+0.5*aT*tau[i]*tau[i])*(z0+vZ*tau[i])*pi;
      Sum = Sum + RateQGP_IntEtaf(Pt,Temp[i])*(1-h[i])*VTau; 
      
    }
  return Sum*Steptime;
}


Double_t RateQGP_IntEtaf(Double_t Pt, Double_t T)
{
  Double_t Etaf=0.0;
  Double_t EtafMin = -4.0;
  Double_t EtafMax = 4.0;
  
  Double_t EtafStep = 0.01;
  
  Int_t NEtaf = (EtafMax - EtafMin)/EtafStep;

  Double_t Sum =0.0;
  for(int i =0;i<NEtaf;i++)
    {
      Etaf = EtafMin + i*EtafStep;
      Sum =Sum + RateQGP(Pt,T,Etaf);
    }

  return Sum*EtafStep;
}



Double_t RateQGP(Double_t Pt, Double_t T, Double_t Etaf)
{
  Double_t yy =0.0; 
  Double_t E = Pt*TMath::CosH(yy-Etaf);

  Double_t Alpha =1.0/137.0;
  //Nf, TC
  Double_t AlphaS = 6.0*pi/((33.0-2.0*Nf)*TMath::Log(8.0*T/TC));

  Double_t RQGP1 =0.0; 
  Double_t Const1=(5.0*Alpha)/(18.0*pi2);
  RQGP1 = Const1*Alpha*AlphaS*TMath::Log(0.23*E/(AlphaS*T))*T*T*TMath::Exp(-(E/T));

  Double_t RQGP2 =0.0; 
  Double_t Const2=0.0219;
  RQGP2 = Const2*Alpha*AlphaS*T*T*TMath::Exp(-(E/T));

  Double_t RQGP3 =0.0; 
  Double_t Const3=0.0105;
  RQGP3 = Const3*Alpha*AlphaS*E*T*TMath::Exp(-(E/T));

  Double_t RQGP =0.0; 

  RQGP = RQGP1 + RQGP2 + RQGP3; 

  return RQGP/hbarc4;
}


// ================================= Rate Hadron ==============================//

Double_t RateHadron_IntTau(Double_t Pt, Double_t RCent)
{
  Double_t Sum =0.0;
  Double_t VTau =0.0;
  for(int i =0;i<Ntime;i++)
    {
      VTau = (RCent+0.5*aT*tau[i]*tau[i])*(RCent+0.5*aT*tau[i]*tau[i])*(z0+vZ*tau[i])*pi;
      Sum = Sum + RateHadron_IntEtaf(Pt,Temp[i])*h[i]*VTau; 

    }

  return Sum*Steptime;

}



Double_t RateHadron_IntEtaf(Double_t Pt, Double_t T)
{
  Double_t Etaf=0.0;
  Double_t EtafMin = -4.0;
  Double_t EtafMax = 4.0;
  Double_t EtafStep = 0.01;
  Int_t NEtaf = (EtafMax - EtafMin)/EtafStep;

  Double_t Sum =0.0;
  for(int i =0;i<NEtaf;i++)
    {
      Etaf = EtafMin + i*EtafStep;
      Sum =Sum + RateHadron(Pt,T,Etaf);
    }

  return Sum*EtafStep;
}



Double_t RateHadron(Double_t Pt, Double_t T, Double_t Etaf)
{
  Double_t yy =0.0; 
  Double_t E = Pt*TMath::CosH(yy-Etaf);

  Double_t RHadron =0.0; 

  Double_t Const=4.8;
  
  Double_t Exponent = 1.0/(TMath::Power((1.35*T*E),0.77));

  RHadron = Const*TMath::Power(T,2.15)*TMath::Exp(-Exponent)*TMath::Exp(-(E/T));

  return RHadron;
}




//not used
Double_t RateQCD_IntTau(Double_t Pt)
{
  Double_t Sum =0.0;
  Double_t IntEtaf = 8.0;

  //Double_t tau[10000], Temp[10000], h[10000];
  //Int_t Ntime;
  //double Steptime;

  for(int i =0;i<=Ntime;i++)
    {
      Sum = Sum + h[i]*tau[i]; 

    }

  return Sum*IntEtaf*RateQCD(Pt)*Steptime;


}

Double_t RateQCD(Double_t Pt)
{
  /*
  Double_t a = -4.1506;
  Double_t b = -1.9845;
  Double_t c = 0.0744;
  Double_t d =-0.0383;
  Double_t RQCD =0.0;
  RQCD = TMath::Exp(a+b*Pt+c*Pt*Pt+d*Pt*Pt*Pt);
  return RQCD;
  */

  Double_t a = 0.0083;
  Double_t b = 2.26;
  Double_t c = 3.45;
  Double_t PtFac =0.0;
  PtFac = (1.0 + ((Pt*Pt)/b));
  
  Double_t RQCD =0.0;
  RQCD = a*TMath::Power(PtFac,-c);
  return RQCD;

}







// ===================================== Direct Photon Graphs (PHENIX) =======================================//
void Draw_AllDataGraphs()
{
 
  TLegend *legend_ratio = new TLegend(0.1677852,0.72,0.83,0.90);
  legend_ratio->SetBorderSize(0);
  legend_ratio->SetFillStyle(0);
  legend_ratio->SetFillColor(0);
  legend_ratio->SetTextSize(0.04);


  new TCanvas;
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLogy();
  Draw_PHENIX_DirectPhotonRate020_Pt(legend_ratio);  



  new TCanvas;
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLogy();
  Draw_PHENIX_DirectPhotonRate2040_Pt(legend_ratio);  
  


  new TCanvas;
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLogy();
  Draw_PHENIX_DirectPhotonRate4060_Pt(legend_ratio);  

  new TCanvas;
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLogy();
  Draw_PHENIX_DirectPhotonRate6092_Pt(legend_ratio);  


  legend_ratio->Draw("same");
}



void Draw_PHENIX_DirectPhotonRate020_Pt(TLegend *lgd)
{

  
  Double_t Pt[11]=   {0.50,0.70,0.90,1.10,1.30,1.50,1.70,1.90,2.25,3.00,4.25}; 
  Double_t ErrPt[11]={0.0};
  
  Double_t D2NDPtDy[11]={18.384316,3.863694,1.531712,0.500764,0.248402,0.105543,0.071944,0.028950,0.007928,0.001125,0.000025}; 
  Double_t StatErrD2NDPtDy[11]={4.162295,0.908128,0.313389,0.119182,0.053970,0.025248,0.014094,0.007115,0.001929,0.000274,0.000033}; 
  Double_t SystErrD2NDPtDy[11]={4.963615,1.412782,0.499280,0.184227,0.077102,0.032901,0.016161,0.007324,0.002126,0.000242,0.000012}; 
  
  //  for(int j=0;j<11;j++){
  //cout<<Pt[j]<<"  "<<D2NDPtDy[j]<<"   "<<StatErrD2NDPtDy[j]<<"   "<<SystErrD2NDPtDy[j]<<endl;
  //  }
  

  TGraphErrors *PhotonRatePHENIXdata = new TGraphErrors(11,Pt,D2NDPtDy,ErrPt,StatErrD2NDPtDy);
  PhotonRatePHENIXdata->SetMarkerStyle(20);
  PhotonRatePHENIXdata->SetMarkerColor(kBlue+0);
  PhotonRatePHENIXdata->SetMarkerSize(1.1);
  PhotonRatePHENIXdata->GetYaxis()->SetTitleOffset(1.8);
  PhotonRatePHENIXdata->GetYaxis()->SetRangeUser(0.000001,60.0);
  TAxis *Xaxis1 = PhotonRatePHENIXdata->GetXaxis();
  Xaxis1->SetLimits(0.0,5.2);
   

  PhotonRatePHENIXdata->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  PhotonRatePHENIXdata->GetYaxis()->SetTitle("#frac{1}{2#pi p_{T}}#frac{d^{2}N}{dp_{T}dy}[(GeV/c)^{-2}]");
  
  PhotonRatePHENIXdata->Draw("AP");

    
  
  TBox *SystErrPhotonRate[11];
  Double_t PtBound[12]=   {0.4, 0.60, 0.80,1.0,1.20, 1.40, 1.60, 1.80, 2.0,  2.5, 3.5, 5.0};
  for(int j=0;j<11;j++){
    SystErrPhotonRate[j] = new TBox(PtBound[j], D2NDPtDy[j]-SystErrD2NDPtDy[j], PtBound[j+1], D2NDPtDy[j]+SystErrD2NDPtDy[j]);
  }
  
  for(int j=0;j<11;j++){
    SystErrPhotonRate[j]->SetFillStyle(9001);
    SystErrPhotonRate[j]->SetLineColor(kBlue);
    SystErrPhotonRate[j]->Draw("same"); 
  }
    
  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  tb->DrawLatex(0.22,0.20,"AuAu #sqrt{s_{NN}} = 200GeV");
  tb->DrawLatex(0.22,0.15,"0-20%");
 
  //lgd->AddEntry(PhotonRatePHENIXdata,"Inclusive Photons", "P");
  
  
  
}



void Draw_PHENIX_DirectPhotonRate2040_Pt(TLegend *lgd)
{

  
  Double_t Pt[11]=   {0.50,0.70,0.90,1.10,1.30,1.50,1.70,1.90,2.25,3.00,4.25}; 
  Double_t ErrPt[11]={0.0};
  
  Double_t D2NDPtDy[11]={5.954453,0.906468,0.396036,0.215631,0.110570,0.050032,0.019825,0.011294,0.004157,0.000554,0.000044}; 
  Double_t StatErrD2NDPtDy[11]={1.590491,0.350932,0.121966,0.050980,0.023209,0.011015,0.005521,0.003037,0.000919,0.000135,0.000018}; 
   Double_t SystErrD2NDPtDy[11]={2.116200,0.588297,0.209478,0.084104,0.036068,0.015939,0.007187,0.003575,0.001142,0.000135,0.000009}; 
   



   //for(int j=0;j<11;j++){
   //cout<<Pt[j]<<"  "<<D2NDPtDy[j]<<"   "<<StatErrD2NDPtDy[j]<<"   "<<SystErrD2NDPtDy[j]<<endl;
   //}
  

  TGraphErrors *PhotonRatePHENIXdata = new TGraphErrors(11,Pt,D2NDPtDy,ErrPt,StatErrD2NDPtDy);
  PhotonRatePHENIXdata->SetMarkerStyle(20);
  PhotonRatePHENIXdata->SetMarkerColor(kBlue+0);
  PhotonRatePHENIXdata->SetMarkerSize(1.1);
  PhotonRatePHENIXdata->GetYaxis()->SetTitleOffset(1.8);
  PhotonRatePHENIXdata->GetYaxis()->SetRangeUser(0.000001,60.0);
  TAxis *Xaxis1 = PhotonRatePHENIXdata->GetXaxis();
  Xaxis1->SetLimits(0.0,5.2);
   

  PhotonRatePHENIXdata->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  PhotonRatePHENIXdata->GetYaxis()->SetTitle("#frac{1}{2#pi p_{T}}#frac{d^{2}N}{dp_{T}dy}[(GeV/c)^{-2}]");
  
  PhotonRatePHENIXdata->Draw("AP");

    
  
  TBox *SystErrPhotonRate[11];
  Double_t PtBound[12]=   {0.4, 0.60, 0.80,1.0,1.20, 1.40, 1.60, 1.80, 2.0,  2.5, 3.5, 5.0};
  for(int j=0;j<11;j++){
    SystErrPhotonRate[j] = new TBox(PtBound[j], D2NDPtDy[j]-SystErrD2NDPtDy[j], PtBound[j+1], D2NDPtDy[j]+SystErrD2NDPtDy[j]);
  }
  
  for(int j=0;j<11;j++){
    SystErrPhotonRate[j]->SetFillStyle(9001);
    SystErrPhotonRate[j]->SetLineColor(kBlue);
    SystErrPhotonRate[j]->Draw("same"); 
  }
    
  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  tb->DrawLatex(0.22,0.20,"AuAu #sqrt{s_{NN}} = 200GeV");
  tb->DrawLatex(0.22,0.15,"20-40%");
 
  //lgd->AddEntry(PhotonRatePHENIXdata,"Inclusive Photons", "P");
  
  
  
}

void Draw_PHENIX_DirectPhotonRate4060_Pt(TLegend *lgd)
{
 
  Double_t Pt[11]=   {0.50,0.70,0.90,1.10,1.30,1.50,1.70,1.90,2.25,3.00,4.25}; 
  Double_t ErrPt[11]={0.0};
  Double_t D2NDPtDy[11]={2.064060,0.530987,0.168944,0.061616,0.041359,0.014628,0.006184,0.003108,0.000850,0.000095,0.000008}; 
  Double_t StatErrD2NDPtDy[11]={0.615983,0.141643,0.047410,0.018912,0.009241,0.004274,0.002196,0.001249,0.000359,0.000054,0.000007}; 
  Double_t SystErrD2NDPtDy[11]={0.794669,0.232527,0.080196,0.030895,0.013852,0.005990,0.002790,0.001383,0.000435,0.000052,0.000003}; 
  

  TGraphErrors *PhotonRatePHENIXdata = new TGraphErrors(11,Pt,D2NDPtDy,ErrPt,StatErrD2NDPtDy);
  PhotonRatePHENIXdata->SetMarkerStyle(20);
  PhotonRatePHENIXdata->SetMarkerColor(kBlue+0);
  PhotonRatePHENIXdata->SetMarkerSize(1.1);
  PhotonRatePHENIXdata->GetYaxis()->SetTitleOffset(1.8);
  PhotonRatePHENIXdata->GetYaxis()->SetRangeUser(0.000001,60.0);
  TAxis *Xaxis1 = PhotonRatePHENIXdata->GetXaxis();
  Xaxis1->SetLimits(0.0,5.2);
   

  PhotonRatePHENIXdata->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  PhotonRatePHENIXdata->GetYaxis()->SetTitle("#frac{1}{2#pi p_{T}}#frac{d^{2}N}{dp_{T}dy}[(GeV/c)^{-2}]");
  
  PhotonRatePHENIXdata->Draw("AP");

    
  
  TBox *SystErrPhotonRate[11];
  Double_t PtBound[12]=   {0.4, 0.60, 0.80,1.0,1.20, 1.40, 1.60, 1.80, 2.0,  2.5, 3.5, 5.0};
  for(int j=0;j<11;j++){
    SystErrPhotonRate[j] = new TBox(PtBound[j], D2NDPtDy[j]-SystErrD2NDPtDy[j], PtBound[j+1], D2NDPtDy[j]+SystErrD2NDPtDy[j]);
  }
  
  for(int j=0;j<11;j++){
    SystErrPhotonRate[j]->SetFillStyle(9001);
    SystErrPhotonRate[j]->SetLineColor(kBlue);
    SystErrPhotonRate[j]->Draw("same"); 
  }
    
  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  tb->DrawLatex(0.22,0.20,"AuAu #sqrt{s_{NN}} = 200GeV");
  tb->DrawLatex(0.22,0.15,"40-60%");
 
  //lgd->AddEntry(PhotonRatePHENIXdata,"Inclusive Photons", "P");
  
 
}


void Draw_PHENIX_DirectPhotonRate6092_Pt(TLegend *lgd)
{

  
  Double_t Pt[11]=   {0.50,0.70,0.90,1.10,1.30,1.50,1.70,1.90,2.25,3.00,4.25}; 
  Double_t ErrPt[11]={0.0};
  Double_t D2NDPtDy[11]={0.283798,0.036460,0.012077,0.008213,0.002898,0.001267,0.001107,0.000137,0.000045,0.000025,0.000007}; 
  Double_t StatErrD2NDPtDy[11]={0.126140,0.025660,0.008446,0.003597,0.001647,0.000829,0.000481,0.000243,0.000076,0.000015,0.000004}; 
  Double_t SystErrD2NDPtDy[11]={0.148048,0.039574,0.013316,0.005268,0.002167,0.000971,0.000490,0.000222,0.000074,0.000011,0.000001}; 
  

  TGraphErrors *PhotonRatePHENIXdata = new TGraphErrors(11,Pt,D2NDPtDy,ErrPt,StatErrD2NDPtDy);
  PhotonRatePHENIXdata->SetMarkerStyle(20);
  PhotonRatePHENIXdata->SetMarkerColor(kBlue+0);
  PhotonRatePHENIXdata->SetMarkerSize(1.1);
  PhotonRatePHENIXdata->GetYaxis()->SetTitleOffset(1.8);
  PhotonRatePHENIXdata->GetYaxis()->SetRangeUser(0.000001,60.0);
  TAxis *Xaxis1 = PhotonRatePHENIXdata->GetXaxis();
  Xaxis1->SetLimits(0.0,5.2);
   

  PhotonRatePHENIXdata->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  PhotonRatePHENIXdata->GetYaxis()->SetTitle("#frac{1}{2#pi p_{T}}#frac{d^{2}N}{dp_{T}dy}[(GeV/c)^{-2}]");
  
  PhotonRatePHENIXdata->Draw("AP");

    
  
  TBox *SystErrPhotonRate[11];
  Double_t PtBound[12]=   {0.4, 0.60, 0.80,1.0,1.20, 1.40, 1.60, 1.80, 2.0,  2.5, 3.5, 5.0};
  for(int j=0;j<11;j++){
    SystErrPhotonRate[j] = new TBox(PtBound[j], D2NDPtDy[j]-SystErrD2NDPtDy[j], PtBound[j+1], D2NDPtDy[j]+SystErrD2NDPtDy[j]);
  }
  
  for(int j=0;j<11;j++){
    SystErrPhotonRate[j]->SetFillStyle(9001);
    SystErrPhotonRate[j]->SetLineColor(kBlue);
    SystErrPhotonRate[j]->Draw("same"); 
  }
    
  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  tb->DrawLatex(0.50,0.85,"AuAu #sqrt{s_{NN}} = 200GeV");
  tb->DrawLatex(0.50,0.80,"60-92%");
 
  //lgd->AddEntry(PhotonRatePHENIXdata,"Inclusive Photons", "P");
  
  
  
}
