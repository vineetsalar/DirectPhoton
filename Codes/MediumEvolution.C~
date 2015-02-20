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
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <TMath.h>
#include "TVector3.h"
#include <TLegend.h>
#include "TCanvas.h"
#include <TPad.h>
#include "TStyle.h"
#include "TFile.h"
#include "TRandom.h"
#include "TLine.h"
#include "TLatex.h"
#include "TSystem.h"
// constants
Double_t pi = TMath::Pi();
Double_t pi2 = pi*pi;
Double_t pi3 = pi2*pi;
Double_t twopi3 = 8.0*pi*pi*pi;
Double_t hbarc = 0.197327;
Double_t hbarc2 = hbarc*hbarc;
Double_t hbarc3 = hbarc2*hbarc;


//Matter at extream condition
//T0 is 550 - 580 MeV and tau0 is 0.3
const Double_t mPi = 0.140; 
const Double_t RPb = 7.11;  //1.2*TMath::Power(208,1.0/3.0)
const Double_t R05 = 0.92*RPb;
const Double_t tau0 = 0.3;

//const Double_t SS=3.6*1.5*1600;
const Double_t SS=5.0*1.5*1600;
const Double_t Nf=2.5;
const Double_t aq = (7.0*Nf/60.0 + 16.0/90.0)*pi2;
const Double_t ah = 4.5*pi2/90.0;

Double_t aT = 0.1;
//Double_t aT = 0.1;
//Double_t z0=1.8*tau0; //0
//Double_t vZ=1.4;     //1.0
//for longitudnal
//Double_t aT = 0.0;   // 0
Double_t z0=0.0; //0
Double_t vZ=1.0;     //1.0

const Double_t VTau0 = (R05+0.5*aT*tau0*tau0)*(R05+0.5*aT*tau0*tau0)*(z0+vZ*tau0)*pi;
const Double_t ss05 = SS/VTau0;

const Double_t T0=TMath::Power(SS/(4.0*aq*VTau0),1.0/3.0)*hbarc;
const Double_t TC = 0.170;
const Double_t TF = 0.140;

//const Double_t tauf = pow(T0/TC, 3.)*tau0;
const Double_t nPart0 = 393;
const Double_t nColl0 = 1747;

int NTau; double stepTau;
double Tau[10000], TempTau[10000], fQGP[10000];



//======================== Function definations =============================//
// Dynamics 
Double_t Npart(int BinLow, int BinHigh);
Double_t NColl(int BinLow, int BinHigh);
Double_t CalculateTandf(int Flag, Double_t T0Cent, Double_t R0Cent, Double_t tauCCent, Double_t tauhCent, Double_t taufCent,Double_t ssCent);
Double_t CalculateTandf_Classic(int Flag, Double_t T0Cent, Double_t R0Cent, Double_t tauCCent, Double_t tauhCent, Double_t taufCent);
Double_t TauToTemp_All(Int_t Flag, Double_t T0Cent, Double_t R0Cent, Double_t Tau);
Double_t FQGP(Int_t Flag, Double_t R0Cent, Double_t TauCCent, Double_t TauHCent, Double_t Tau);
Double_t Tau0ToTauC(Int_t Flag, Double_t T0Cent, Double_t R0Cent);
Double_t TauCToTauh(Int_t Flag, Double_t TauC, Double_t R0Cent);
Double_t TauhToTauf(Int_t Flag, Double_t Tauh, Double_t R0Cent);
Double_t CalculateTandf_LatticeEOS( Double_t ssCent, Double_t R0Cent);

//==================================== Lattice EOS ================================================//
TFile *fileEOS=new TFile("LatticeEOS_s95p-v1.2.root","R");
TGraph *grfSSVsTemp = (TGraph*)fileEOS->Get("grfSSVsTemp");
TGraph *TempVsFQGP = (TGraph*)fileEOS->Get("TempVsFQGP");
TGraph *TempVsFQGP2 = (TGraph*)fileEOS->Get("TempVsFQGP2");



void MediumEvolution()
{


//==================================//
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptTitle(0);
  //gStyle->SetOptStat("nmr");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetTextSize(0.04);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.048,"xyz");
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadRightMargin(0.065);
  gStyle->SetPadLeftMargin(0.12);
  //gStyle->SetTitleXOffset(1.15);
  //gStyle->SetTitleYOffset(1.2);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetCanvasDefH(600);
  gStyle->SetCanvasDefW(600);
  gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetErrorX(0);   
  gStyle->SetEndErrorSize(0);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.0);
  gROOT->ForceStyle();
  


 TLatex *tb = new TLatex();
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  
  int Disso=1;
  int Form =1;
  
  TLegend *lcat[100];
  for(int i =0;i<90;i++)
    {
      lcat[i]= new TLegend( 0.55,0.59,0.84,0.78);;
      lcat[i]->SetBorderSize(0);
      lcat[i]->SetFillStyle(0);
      lcat[i]->SetFillColor(0);
      lcat[i]->SetTextSize(0.040);
    }




  cout<<" simulating QGP evolution : "<<endl;

 // ================ dn/deta graph for making Temp as a function of nPart ==========================================//
  Double_t NPartdNdEta[10] = {382.8,329.7,260.5,186.4,128.9,85.0,52.8,30.0,15.8};
  Double_t Err_NPartdNdEta[10] = {3.1,4.6,4.4,3.9,3.3,2.6,2.0,1.3,0.6};
  
  Double_t dNdEtabyNpartby2[10] = {8.4,7.9,7.4,7.0,6.6,6.1,5.7,5.1,4.4};
  Double_t Err_dNdEtabyNpartby2[10] = {0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.4};
  
  TGraphErrors *grdNDetaNpart = new TGraphErrors(9,NPartdNdEta,dNdEtabyNpartby2,Err_NPartdNdEta,Err_dNdEtabyNpartby2);
  grdNDetaNpart->SetLineWidth(2);
  grdNDetaNpart->SetMarkerStyle(20);
  grdNDetaNpart->GetYaxis()->SetRangeUser(3.5,9.0);
  grdNDetaNpart->GetYaxis()->SetTitle("#frac{dN}{d#eta}/(#frac{N_{Part}}{2})");
  grdNDetaNpart->GetXaxis()->SetTitle("N_{Part}");
  
  //new TCanvas;
  //gPad->SetTicks();
  //gPad->SetBottomMargin(0.14);
  //grdNDetaNpart->Draw("AP");
  //gPad->SaveAs("dNdEtaVsNpart.png");  
  //gPad->SaveAs("dNdEtaVsNpart.eps");     
  
  TF1 *funDnDetaNPart = new TF1("funDnDetaNPart","[0]*TMath::Power(x,[1])/(0.5*x)",0.0,500);
  funDnDetaNPart->SetParameters(1.317,1.190);
  funDnDetaNPart->SetLineColor(4);
  
  TF1 *funTwoComp = new TF1("funTwoComp","[0]*( [1]*x + (1-[1])*([2] + [3]*x +[4]*x*x)) /(0.5*x)",0.0,500);
  funTwoComp->SetLineColor(2);
  funTwoComp->SetParameter(0,2.441);
  funTwoComp->SetParameter(1,0.788);
  funTwoComp->SetParameter(2,-13.4708);
  funTwoComp->SetParameter(3,1.69143);
  funTwoComp->SetParameter(4,0.00709679);

  lcat[3] = new TLegend( 0.37,0.22,0.87,0.40);
  lcat[3]->SetBorderSize(0);
  lcat[3]->SetFillStyle(0);
  lcat[3]->SetFillColor(0);
  lcat[3]->SetTextSize(0.040);
  lcat[3]->AddEntry(funDnDetaNPart,"#alpha N_{Part}^{#beta}","L");
  lcat[3]->AddEntry(funTwoComp,"fN_{Part} + (1-f)N_{Coll}","L");

  new TCanvas;
  grdNDetaNpart->Fit("funDnDetaNPart","M","",10.0,400.0);
  grdNDetaNpart->Fit("funTwoComp","M","",10.0,400.0);

  grdNDetaNpart->Draw("AP");
  funDnDetaNPart->Draw("same");
  funTwoComp->Draw("same");
  lcat[3]->Draw("same");
  

  // return;

  Double_t TempNpart[400];
  Double_t TempNpart3[400];
  Double_t NpartVal[400];
  
  for(int i=0;i<=80;i++)
    {
      NpartVal[i]=i*5;
      TempNpart[i]=T0*TMath::Power(funDnDetaNPart->Eval(NpartVal[i])/dNdEtabyNpartby2[0],1.0/3); //ALICE Func
      TempNpart3[i]=T0*TMath::Power(funTwoComp->Eval(NpartVal[i])/dNdEtabyNpartby2[0],1.0/3); //Two Comp model
    }

  TGraph *grTempVsNpart = new TGraph(80, NpartVal,TempNpart);
  grTempVsNpart->SetLineWidth(2);
  grTempVsNpart->SetMarkerStyle(20);
  grTempVsNpart->SetMarkerColor(4);
  
  grTempVsNpart->GetXaxis()->SetTitle("N_{Part}");
  grTempVsNpart->GetYaxis()->SetTitle("T (GeV)");
  grTempVsNpart->GetYaxis()->SetRangeUser(0.1,1.0);
  
  TGraph *grTempVsNpart3 = new TGraph(80, NpartVal,TempNpart3);
  grTempVsNpart3->SetLineWidth(2);
  grTempVsNpart3->SetMarkerStyle(20);
  grTempVsNpart3->SetMarkerColor(2);


  new TCanvas;
  gPad->SetTicks();
  grTempVsNpart->Draw("AP");
  grTempVsNpart3->Draw("Psame");
  lcat[3]->Draw("same");
  gPad->SaveAs("NpartVsTemp.png");  
  gPad->SaveAs("NpartVsTemp.eps");  
  
  
  //================================= Temp Vs Tau Graphs ====================================//
  
  Double_t T0Cen = T0*TMath::Power((funTwoComp->Eval(Npart(0,2)))/dNdEtabyNpartby2[0],1.0/3.0);
  Double_t ssCen = ss05*(funTwoComp->Eval(Npart(0,2))/dNdEtabyNpartby2[0]);  
 

  double TauL[10000], TempTauL[10000], fQGPL[10000];
  int NTauL;
  Double_t TauC1=Tau0ToTauC(1,T0Cen,R05);
  Double_t TauH1 = TauCToTauh(1,TauC1,R05);
  Double_t TauF1 =TauhToTauf(1,TauH1,R05);
  CalculateTandf(1, T0Cen, R05, TauC1, TauH1,TauF1,ssCen);
  
  NTauL= NTau;
  cout<<" Long " << T0Cen<<"  "<<R05<<"  "<<TauC1<<"  "<<TauH1<<"   "<<TauF1<<endl;
  for(int i=0;i<NTauL;i++){
    TauL[i]=Tau[i];
    TempTauL[i] = TempTau[i];
    fQGPL[i]=fQGP[i];
  }
 
  double TauC[10000], TempTauC[10000], fQGPC[10000];
  int NTauC;  

  Double_t TauC2=Tau0ToTauC(2,T0Cen,R05);
  Double_t TauH2 = TauCToTauh(2,TauC2,R05);
  Double_t TauF2=TauhToTauf(2, TauH2, R05);
  CalculateTandf(2,T0Cen,R05,TauC2,TauH2,TauF2,ssCen);
  NTauC= NTau;
  for(int i=0;i<NTau;i++){
    TauC[i]=Tau[i];
    TempTauC[i] = TempTau[i];
    fQGPC[i]=fQGP[i];
  }

  cout<<" Cylind  " << T0Cen<<"  "<<R05<<"  "<<TauC2<<"  "<<TauH2<<"  "<<TauF2<<endl;


  double TauLatt[10000], TempTauLatt[10000];
  double fQGPLatt[10000];
  

  cout<<" Calculating Lattice equation of state "<<endl;
  NTau=0;
  CalculateTandf(3,T0Cen,R05,TauC2,TauH2,TauF2,ssCen);

 
  cout<<" NTau "<<NTau<<" ssCen "<<ssCen<<endl;
  
  for(int i=0;i<NTau;i++){
    TauLatt[i]=Tau[i];
    TempTauLatt[i]=TempTau[i];
    fQGPLatt[i]=fQGP[i];
    cout<<Tau[i]<<" Lattice  "<<TempTau[i]<<endl;
  }
    
 

  
  TGraph *grTempVsTauAnaC = new TGraph(NTau,TauLatt,TempTauLatt);
  grTempVsTauAnaC->SetName("grTempVsTauLatt");
  grTempVsTauAnaC->SetTitle("grTempVsTauLatt");
  grTempVsTauAnaC->SetLineColor(1);
  grTempVsTauAnaC->SetLineWidth(2);
  grTempVsTauAnaC->SetLineStyle(9);
  grTempVsTauAnaC->GetYaxis()->SetTitle("Temperature (GeV)");
  grTempVsTauAnaC->GetXaxis()->SetTitle("#tau (fm)");


  new TCanvas;
  grTempVsTauAnaC->Draw("AL");
        //grTempVsTauAnaC
  TGraph *grFQGPVsTauAnaC = new TGraph(NTau,TauLatt,fQGPLatt);
  grFQGPVsTauAnaC->SetName("grFQGPVsTauLatt");
  grFQGPVsTauAnaC->SetTitle("grFQGPVsTauLatt");
  grFQGPVsTauAnaC->SetLineColor(1);
  grFQGPVsTauAnaC->SetLineWidth(2);
  grFQGPVsTauAnaC->SetLineStyle(9);
  grFQGPVsTauAnaC->GetYaxis()->SetTitle("FQGP");
  grFQGPVsTauAnaC->GetXaxis()->SetTitle("#tau (fm)");

  new TCanvas;
  grFQGPVsTauAnaC->Draw("AL");
  //grTempVsTauAnaC->Draw("AL");





  TGraph *grTempVsTauAna = new TGraph(NTauL,TauL,TempTauL);
  grTempVsTauAna->SetName("grTempVsTauL");
  grTempVsTauAna->SetTitle("grTempVsTauL");
  grTempVsTauAna->SetLineColor(2);
  grTempVsTauAna->SetLineWidth(2);
  grTempVsTauAna->SetLineStyle(1);
  grTempVsTauAna->GetYaxis()->SetTitle("Temperature (GeV)");
  grTempVsTauAna->GetXaxis()->SetTitle("#tau (fm)");
  
  TGraph *grTempVsTauAnaB = new TGraph(NTauC,TauC,TempTauC);
  grTempVsTauAnaB->SetName("grTempVsTauC");
  grTempVsTauAnaB->SetTitle("grTempVsTauC");
  grTempVsTauAnaB->SetLineColor(4);
  grTempVsTauAnaB->SetLineWidth(2);
  grTempVsTauAnaB->SetLineStyle(2);

  TLegend *legd5 = new TLegend( 0.15,0.79,0.87,0.93);
  legd5->SetBorderSize(0);
  legd5->SetFillStyle(0);
  legd5->SetFillColor(0);
  legd5->SetTextSize(0.040);
  legd5->AddEntry(grTempVsTauAna,"FO-EOS, Longitudinal expansion","L");
  legd5->AddEntry(grTempVsTauAnaB,"FO-EOS, Cylindrical expansion","L");
  legd5->AddEntry(grTempVsTauAnaC,"Lattice EOS, Cylindrical expansion","L");
   

  TLine *l1 = new TLine(0.0,0.16,10.0,0.16);
  l1->SetLineWidth(2);
  l1->SetLineColor(2);

  new TCanvas;  
  gPad->SetTicks(1);
  TAxis *Xaxis3 = grTempVsTauAna->GetXaxis();
  Xaxis3->SetLimits(0.0,20.0);
  grTempVsTauAna->GetYaxis()->SetRangeUser(0.0,0.7);
  grTempVsTauAna->Draw("AL");
  grTempVsTauAnaB->Draw("Lsame");
  grTempVsTauAnaC->Draw("Lsame");

  tb->DrawLatex(0.16,0.18,"(a)");
  legd5->Draw("Lsame");
  gPad->SaveAs("Fig1a_TauVsTemp.eps");
  gPad->SaveAs("Fig1a_TauVsTemp.png");
  gPad->SaveAs("Fig1a_TauVsTemp.pdf");

  TGraph *grFQGPVsTau = new TGraph(NTauL,TauL,fQGPL);
  grFQGPVsTau->SetName("grFQGPVsTauL");
  grFQGPVsTau->SetTitle("grFQGPVsTauL");
  grFQGPVsTau->SetLineColor(2);
  grFQGPVsTau->SetLineWidth(2);
  grFQGPVsTau->SetLineStyle(1);
  grFQGPVsTau->GetYaxis()->SetTitle("QGP fraction");
  grFQGPVsTau->GetXaxis()->SetTitle("#tau (fm)");

  TGraph *grFQGPVsTau2 = new TGraph(NTauC,TauC,fQGPC);
  grFQGPVsTau2->SetName("grFQGPVsTauC");
  grFQGPVsTau2->SetTitle("grFQGPVsTauC");
  grFQGPVsTau2->SetLineColor(4);
  grFQGPVsTau2->SetLineWidth(2);
  grFQGPVsTau2->SetLineStyle(2);


  TLegend *legdA = new TLegend( 0.15,0.79,0.87,0.93);
  legdA->SetBorderSize(0);
  legdA->SetFillStyle(0);
  legdA->SetFillColor(0);
  legdA->SetTextSize(0.040);
  legdA->AddEntry(grTempVsTauAna,"FO-EOS, Longitudinal expansion","L");
  legdA->AddEntry(grTempVsTauAnaB,"FO-EOS, Cylindrical expansion","L");
  legdA->AddEntry(grTempVsTauAnaC,"Lattice EOS, Cylindrical expansion","L");
  
  new TCanvas;  
  gPad->SetTicks(1);
  grFQGPVsTau->GetYaxis()->SetRangeUser(0.0,1.4);
  TAxis *Xaxis4 = grFQGPVsTau->GetXaxis();
  Xaxis4->SetLimits(0.,20.0);
  grFQGPVsTau->Draw("AL");
  
  grFQGPVsTau2->Draw("Lsame");
  grFQGPVsTauAnaC->Draw("Lsame");
  
  legdA->Draw("Lsame");
  tb->DrawLatex(0.16,0.18,"(b)");
  
  gPad->SaveAs("Fig1b_TauVsFQGP.eps");
  gPad->SaveAs("Fig1b_TauVsFQGP.png");
  gPad->SaveAs("TFig1b_auVsFQGP.pdf");
  
  //return 
  


}


//////////////////////////////////////////////////////////////////////////////////////////////////
Double_t Npart(int BinLow, int BinHigh)
{
  Double_t NpartArray[40]={393.622,368.96,342.32,316.49,293.49,271.98,249.65,230.53,212.28,194.50,178.54,
			 163.25,149.05,135.92,123.28,111.67,100.79,90.71,80.93,72.60,64.15,56.61,49.95,
			 43.39,37.83,32.70,27.86,23.79,20.20,16.85,14.04,11.60,9.55,7.72,6.44,4.96,4.22,
			 3.50,3.17,2.79};
  Double_t sum=0;
  for(int i=BinLow;i<BinHigh;i++){
    sum+=NpartArray[i];
  }
  Double_t NPart = sum/(BinHigh-BinLow);
  return NPart;
}

Double_t NColl(int BinLow, int BinHigh)
{
  Double_t NCollArray[40]={1747.8600, 1567.5300, 1388.3900, 1231.7700, 1098.2000, 980.4390, 861.6090, 766.0420, 676.5150, 593.4730,
			   521.9120, 456.5420, 398.5460, 346.6470, 299.3050, 258.3440, 221.2160, 188.6770, 158.9860, 134.7000,
			   112.5470, 93.4537, 77.9314, 63.5031, 52.0469, 42.3542, 33.9204, 27.3163, 21.8028, 17.2037,
			   13.5881, 10.6538, 8.3555, 6.4089, 5.1334, 3.7322, 3.0663, 2.4193, 2.1190, 1.7695};
  
  Double_t sum=0;
  for(int i=BinLow;i<BinHigh;i++){
    sum+=NCollArray[i];
  }
  Double_t NColl = sum/(BinHigh-BinLow);
  return NColl;
}







Double_t CalculateTandf_LatticeEOS(Double_t ssCent, Double_t R0Cent)
{
  stepTau = 0.1;
  NTau=0;
  Double_t CutTemp =0.0;
  do{
    Tau[NTau] = tau0 + stepTau*NTau;
    Double_t VTau0Cent =  (R0Cent+0.5*aT*tau0*tau0)*(R0Cent+0.5*aT*tau0*tau0)*(z0+vZ*tau0)*pi;
    Double_t VTauCent =  (R0Cent+0.5*aT*Tau[NTau]*Tau[NTau])*(R0Cent+0.5*aT*Tau[NTau]*Tau[NTau])*(z0+vZ*Tau[NTau])*pi;
    Double_t ssTau = (hbarc3*ssCent)*(VTau0Cent/VTauCent);
    TempTau[NTau]=grfSSVsTemp->Eval(ssTau);
    Double_t tt= TempTau[NTau];
    //fQGP[NTau]=TempVsFQGP->Eval(tt);
    fQGP[NTau]=TempVsFQGP2->Eval(tt);
    CutTemp = TempTau[NTau];
    cout<<"Lattice: NTau "<<NTau<<" ss "<<ssTau<<"  "<<Tau[NTau]<<"   "<<TempTau[NTau]<<"  "<<fQGP[NTau]<<endl;
    NTau++;
  }while(CutTemp > TF);
  
  return 0;

}


Double_t CalculateTandf(int Flag, Double_t T0Cent, Double_t R0Cent, Double_t tauCCent, Double_t tauhCent, Double_t taufCent, Double_t ssCent)
{
  if(Flag ==1){CalculateTandf_Classic(1,T0Cent,R0Cent,tauCCent,tauhCent,taufCent);}
  if(Flag ==2){CalculateTandf_Classic(2,T0Cent,R0Cent,tauCCent,tauhCent,taufCent);}
  if(Flag ==3){CalculateTandf_LatticeEOS(ssCent,R0Cent);}

  return 0;

}





Double_t CalculateTandf_Classic(int Flag, Double_t T0Cent, Double_t R0Cent, Double_t tauCCent, Double_t tauhCent, Double_t taufCent)
{

  stepTau = 0.1;
  NTau = (int)((taufCent-tau0)/stepTau);

  for (int i= 0; i<=NTau; i++) {
    Tau[i] = tau0 + stepTau*i;
    TempTau[i] = TauToTemp_All(Flag,T0Cent,R0Cent,Tau[i]);
    fQGP[i] = FQGP(Flag,R0Cent,tauCCent, tauhCent,Tau[i]);
  }
  return 0;
}





Double_t FQGP(Int_t Flag, Double_t R0Cent, Double_t TauCCent, Double_t TauHCent, Double_t Tau1)
{
  Double_t fqgp=0.0;
  Double_t Nomi=0.0;
  Double_t Deno=0.0;
  
  if(Tau1 <= TauCCent){fqgp =1.0; return fqgp;}
  if(Tau1 > TauHCent){fqgp =0.0; return fqgp;}

  if(Flag ==1){
    Nomi=aq*(Tau1-TauCCent);
    Deno=(aq-ah)*Tau1;
    fqgp = 1 - (Nomi/Deno); 
  }
 
  if(Flag ==2) { 
    Double_t Rtau = R0Cent+0.5*aT*Tau1*Tau1;
    Double_t RtauC = R0Cent+0.5*aT*TauCCent*TauCCent;
    Nomi=aq*TauCCent*RtauC*RtauC - ah*Tau1*Rtau*Rtau;
    Deno= (aq - ah)*Tau1*Rtau*Rtau;
    fqgp = Nomi/Deno;
  }
  return fqgp;
}


Double_t TauToTemp_All(Int_t Flag, Double_t T0Cent, Double_t R0Cent, Double_t Tau1) 
{
  Double_t TauC = 0.0;
  Double_t TauH = 0.0;
  Double_t TauF = 0.0;
  
  Double_t Temp = 0;
  if( Flag ==1) 
    {
      TauC = Tau0ToTauC(1,T0Cent,R0Cent);
      TauH = TauCToTauh(1,TauC, R0Cent);
      TauF = TauhToTauf(1,TauH, R0Cent);

      if(Tau1 >= tau0 && Tau1 <= TauC) {Temp = pow(T0Cent*T0Cent*T0Cent*tau0/Tau1, 1.0/3.0);}
      if(Tau1 > TauC && Tau1 <= TauH ) {Temp = TC;}
      if(Tau1 > TauH && Tau1 <= TauF ) {Temp = TC*pow(TauH/Tau1,1.0/3.0);}
    }
  
  if(Flag ==2) {
    TauC = Tau0ToTauC(2,T0Cent,R0Cent);
    TauH = TauCToTauh(2,TauC, R0Cent);
    TauF = TauhToTauf(2,TauH, R0Cent);
    
    if(Tau1 >= tau0 && Tau1 <= TauC) { 
      Double_t Num = T0Cent*T0Cent*T0Cent*(R0Cent+0.5*aT*tau0*tau0)*(R0Cent+0.5*aT*tau0*tau0)*(z0+vZ*tau0);   
      Double_t Deno = (R0Cent+0.5*aT*Tau1*Tau1)*(R0Cent+0.5*aT*Tau1*Tau1)*(z0+vZ*Tau1);
      Double_t Sum = Num/Deno;
      Temp = pow(Sum, 1.0/3.0);
    }
    if(Tau1 > TauC && Tau1 <= TauH ){Temp = TC;}
    if(Tau1 > TauH && Tau1 <= TauF ){
      Double_t Num=TauH*(R0Cent+0.5*aT*TauH*TauH)*(R0Cent+0.5*aT*TauH*TauH);
      Double_t Deno=Tau1*(R0Cent+0.5*aT*Tau1*Tau1)*(R0Cent+0.5*aT*Tau1*Tau1);
      Double_t Sum=Num/Deno;
      Temp = TC*pow((Sum),1.0/3.0);   
    }
  }
  
  return Temp; 
}

Double_t Tau0ToTauC(Int_t Flag, Double_t T0Cent, Double_t R0Cent) 
{
  Double_t TauCValLong=0;
  Double_t TauCValCylnd = 0;
  TauCValLong = T0Cent*T0Cent*T0Cent*tau0/(TC*TC*TC);
  
  if(Flag ==2) {
    Double_t TauCG = TauCValLong;
    Double_t delta;
    Int_t iters=0;
    do{
      Double_t sV0 =T0Cent*T0Cent*T0Cent*tau0*(R0Cent+0.5*aT*tau0*tau0)*(R0Cent+0.5*aT*tau0*tau0);
      Double_t Fx = TC*TC*TC*TauCG*(R0Cent+0.5*aT*TauCG*TauCG)*(R0Cent+0.5*aT*TauCG*TauCG)-sV0;     
      Double_t FxPrime= TC*TC*TC*2*(R0Cent*TauCG*TauCG+0.5*aT*TauCG*TauCG*TauCG*TauCG)*(2.0*R0Cent*TauCG + 2.0*aT*TauCG*TauCG*TauCG);
      TauCValCylnd = TauCG - 5.0*(Fx/FxPrime);
      delta=fabs(TauCValCylnd-TauCG);
      TauCG=TauCValCylnd;
      iters++;
      //cout<<delta<<" delta "<<TauCValCylnd<<endl;
    }while(delta > 0.001);
   
    //Double_t sV0 =T0Cent*T0Cent*T0Cent*tau0*(R0Cent+0.5*aT*tau0*tau0)*(R0Cent+0.5*aT*tau0*tau0);
    //Double_t Fx2 = TC*TC*TC*TauCG*(R0Cent+0.5*aT*TauCG*TauCG)*(R0Cent+0.5*aT*TauCG*TauCG)-sV0;     
    //cout<<" TauCG "<<TauCG<< "  Fx2 "<<Fx2<<" iter: "<<iters<<endl;
  }
  
  Double_t TauCVal=0;
  if(Flag==1) {TauCVal=TauCValLong;}
  if(Flag==2) {TauCVal=TauCValCylnd;}
  
  return TauCVal;
}


Double_t TauCToTauh(Int_t Flag, Double_t TauC, Double_t R0Cent)
{
  Double_t TauHLong =aq*TauC/ah;
  Double_t TauHCylnd=0;

  if(Flag ==2){
    Double_t TauHG=TauHLong;
    Double_t delta;
    Int_t iters=0;
    do{
      Double_t C = aq*TauC*(R0Cent+0.5*aT*TauC*TauC)*(R0Cent+0.5*aT*TauC*TauC);
      Double_t Fx = 0.25*ah*aT*aT*pow(TauHG,5) + R0Cent*ah*aT*pow(TauHG,3) + ah*R0Cent*R0Cent*TauHG - C;
      Double_t FPrimex = 1.25*ah*aT*aT*pow(TauHG,4) + 3.0*R0Cent*ah*aT*pow(TauHG,2) + ah*R0Cent*R0Cent;
      TauHCylnd = TauHG  - (Fx/FPrimex);
      delta=fabs(TauHCylnd-TauHG);
      TauHG=TauHCylnd;
      iters++;
    }while(delta > 0.001);
    
    //Double_t C =aq*TauC*(R0Cent+0.5*aT*TauC*TauC)*(R0Cent+0.5*aT*TauC*TauC);
    //Double_t Fx2 =  0.25*ah*aT*aT*pow(TauHG,5) + R0Cent*ah*aT*pow(TauHG,3) + ah*R0Cent*R0Cent*TauHG - C;
    // cout<<TauHLong << "  "<<TauHG<<" Fx2   "<<Fx2<<"   "<<iters<<endl;
  }
  Double_t TauHVal=0;
  
  if(Flag ==1){TauHVal=TauHLong;}
  if(Flag ==2){TauHVal=TauHCylnd;}

  return TauHVal;
}

Double_t TauhToTauf(Int_t Flag, Double_t Tauh, Double_t R0Cent)
{
  Double_t TauFCylnd=0;
  Double_t TauFLong=0;
  TauFLong = TC*TC*TC*Tauh/(TF*TF*TF);
  if(Flag ==2){
    Double_t TauFG=TauFLong ;
    Double_t delta;
    Int_t iters=0;
    do{
      Double_t C = TC*TC*TC*Tauh*(R0Cent+0.5*aT*Tauh*Tauh)*(R0Cent+0.5*aT*Tauh*Tauh);
      Double_t Fx =TF*TF*TF*TauFG*(R0Cent+0.5*aT*TauFG*TauFG)*(R0Cent+0.5*aT*TauFG*TauFG)-C;
      Double_t FPrimex = TF*TF*TF*(R0Cent+0.5*aT*TauFG*TauFG)*(R0Cent+0.5*aT*TauFG*TauFG)+2.0*TF*TF*TF*aT*TauFG*TauFG*(R0Cent+0.5*aT*TauFG*TauFG);
      TauFCylnd = TauFG  - (Fx/FPrimex);
      delta=fabs(TauFCylnd-TauFG);
      TauFG=TauFCylnd;
      iters++;
    }while(delta > 0.001);
  }
  Double_t TauFVal=0;
  if(Flag==1){TauFVal=TauFLong;}
  if(Flag==2){TauFVal=TauFCylnd;}
  return TauFVal;
}
