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
#include "TGraphAsymmErrors.h"
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

Double_t Radius = 1.2*TMath::Power(208.0, 1/3.);
Double_t contime = Radius*Radius*pi;
Double_t conrate = 1.0/(137.0*137.0*2.0*pi3);



/*
Double_t aq = 37.0*pi2/90.0;
Double_t EDNDY=1.5 * 1600;


//from ArXive 0705.1591
Double_t tau_i = 0.1;
//Double_t tau_i = 0.03;
Double_t Temp_i =  hbarc*TMath::Power((3.6*EDNDY /(tau_i*contime*4*aq)),1.0/3.0);
Double_t Temp_C = 0.160;
Double_t Temp_f = 0.130;
*/


//Matter at extream condition
//T0 is 550 - 580 MeV and tau0 is 0.3
const Double_t mPi = 0.140; 
const Double_t RPb = 7.11;  //1.2*TMath::Power(208,1.0/3.0)
const Double_t R05 = 0.96*RPb;
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

Double_t EtaStart = -1;
Double_t EtaEnd = 1;

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
Double_t RateQGP_IntTau(Double_t Pt, Double_t RCent);
Double_t RateQGP_IntEtaf(Double_t Pt, Double_t T);
Double_t RateQGP(Double_t Pt, Double_t T, Double_t Etaf);
Double_t RateHadron_IntTau(Double_t Pt, Double_t RCent);
Double_t RateHadron_IntEtaf(Double_t Pt, Double_t T);
Double_t RateHadron(Double_t Pt, Double_t T, Double_t Etaf);
Double_t RateQCD_IntTau(Double_t Pt);
Double_t RateQCD(Double_t Pt);


//======= Data Graphs =========================//
void Draw_AllDataGraphs();
void Draw_ALICE_DirectPhotonRate040_Pt(TLegend *lgd);

void PhotonHistLHC()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  //gStyle->SetPadTopMargin(0.10);
  //gStyle->SetPadBottomMargin(0.10);
  gStyle->SetPadLeftMargin(0.15);
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

  //return;


  TFile *OutFile = new TFile("LHC_DirecPhoton.root","Recreate");

  
  cout<<" simulating QGP evolution : "<<endl;



  //============ QCD rates from ALICE Paper =================================//
  
  //From ALICE Paper 
  
  Double_t APt[30]={2.14,   2.38,   2.59,   2.86,   3.15,   3.48,   3.81,   4.16,   4.51,   4.86,   5.25,   5.62,   6.07,   
		    6.58,   7.00,   7.47,   7.90,   8.35,   8.83,   9.34,   9.71,   10.24,   10.74,   11.27,   11.83,   
		    12.30,   12.59,   13.08,   13.55,   13.94};
  
  Double_t D2NByPtDPtDy[30]={0.009903112,   0.006549429,   0.003952514,   0.002736074,   0.001650603,   0.001091140,   0.000721304,   
			     0.000499134,   0.000345396,   0.000228305,   0.000150882,   0.000114419,   0.000079141,   0.000049934,   
			     0.000039631,   0.000027410,   0.000021754,   0.000016491,   0.000011405,   0.000008643,   0.000007520,   
			     0.000005444,   0.000004733,   0.000003426,   0.000002845,   0.000002258,   0.000001793,   0.000001559,   
			     0.000001295,   0.000001180};

  
  TGraph *Grf_QCDRate = new TGraph(30,APt,D2NByPtDPtDy);
  Grf_QCDRate->SetLineWidth(2);
  Grf_QCDRate->SetMarkerStyle(20);
  //Grf_QCDRate->GetYaxis()->SetRangeUser(0.0000001,1.0);
  Grf_QCDRate->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  Grf_QCDRate->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  Grf_QCDRate->GetYaxis()->SetTitleOffset(1.5);
  
  new TCanvas;
  gPad->SetTicks();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.19);
  Grf_QCDRate->Draw("AP");

  //return;



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

  TLegend *lcat = new TLegend( 0.37,0.22,0.87,0.40);
  lcat->SetBorderSize(0);
  lcat->SetFillStyle(0);
  lcat->SetFillColor(0);
  lcat->SetTextSize(0.040);
  
  lcat->AddEntry(funDnDetaNPart,"#alpha N_{Part}^{#beta}","L");
  lcat->AddEntry(funTwoComp,"fN_{Part} + (1-f)N_{Coll}","L");

  new TCanvas;
  grdNDetaNpart->Fit("funDnDetaNPart","M","",10.0,400.0);
  grdNDetaNpart->Fit("funTwoComp","M","",10.0,400.0);

  grdNDetaNpart->Draw("AP");
  funDnDetaNPart->Draw("same");
  funTwoComp->Draw("same");
  lcat->Draw("same");



  // Calculate Temp and Hadronic fraction as a function of time
  double TauLatt[10000], TempTauLatt[10000];
  double fQGPLatt[10000];


  //CalTime();
  
  //CalculateTandf_LatticeEOS(Double_t ssCent, Double_t R0Cent);
  
  Double_t R0MB = R05*TMath::Power(Npart(0,40)/Npart(0,2),0.5);
  Double_t ssMB = ss05*(funTwoComp->Eval(Npart(0,40))/dNdEtabyNpartby2[0]);
  
  Ntime=0;
  CalculateTandf_LatticeEOS(ssMB, R0MB);

  cout<<" Ntime "<<Ntime <<" ssMB "<<ssMB<<"  R0MB: "<<R0MB<<endl;
  
  for(int i=0;i<Ntime;i++){
    TauLatt[i]=tau[i];
    TempTauLatt[i]=Temp[i];
    fQGPLatt[i]=(1.0-h[i]);
    cout<<tau[i]<<" Lattice  "<<Temp[i]<<endl;
  }
    
  TGraph *grTempVsTauAnaC = new TGraph(Ntime,TauLatt,TempTauLatt);
  grTempVsTauAnaC->SetName("grTempVsTauLatt");
  grTempVsTauAnaC->SetTitle("grTempVsTauLatt");
  grTempVsTauAnaC->SetLineColor(1);
  grTempVsTauAnaC->SetLineWidth(2);
  grTempVsTauAnaC->SetLineStyle(9);
  grTempVsTauAnaC->GetYaxis()->SetTitle("Temperature (GeV)");
  grTempVsTauAnaC->GetXaxis()->SetTitle("#tau (fm)");
  grTempVsTauAnaC->GetYaxis()->SetTitleOffset(1.5);
    
  new TCanvas;
  gPad->SetTicks();
  grTempVsTauAnaC->Draw("AL");
  gPad->SaveAs("LHC_TempVsTau.pdf");  
  gPad->SaveAs("LHC_TempVsTau.png");
  gPad->SaveAs("LHC_TempVsTau.eps");
  

  TGraph *grFQGPVsTauAnaC = new TGraph(Ntime,TauLatt,fQGPLatt);
  grFQGPVsTauAnaC->SetName("grFQGPVsTauLatt");
  grFQGPVsTauAnaC->SetTitle("grFQGPVsTauLatt");
  grFQGPVsTauAnaC->SetLineColor(1);
  grFQGPVsTauAnaC->SetLineWidth(2);
  grFQGPVsTauAnaC->SetLineStyle(9);
  grFQGPVsTauAnaC->GetYaxis()->SetTitle("FQGP");
  grFQGPVsTauAnaC->GetXaxis()->SetTitle("#tau (fm)");

  new TCanvas;
  gPad->SetTicks();
  grFQGPVsTauAnaC->Draw("AL");
  gPad->SaveAs("LHC_FQGPVsTau.pdf");  
  gPad->SaveAs("LHC_FQGPVsTau.png");
  gPad->SaveAs("LHC_FQGPVsTau.eps");

  Double_t R040 = R05*TMath::Power(Npart(0,16)/Npart(0,2),0.5);
  Double_t ss040 = ss05*(funTwoComp->Eval(Npart(0,16))/dNdEtabyNpartby2[0]);
  
  Ntime=0;
  CalculateTandf_LatticeEOS(ss040, R040);
  
  Double_t PtMin =0.1;
  Double_t PtMax =12.0;
  Double_t PtStep =0.5;
  Int_t NPt = (PtMax - PtMin)/PtStep;
  
  Double_t Pt[100]={0.0};

  Double_t DirectPhoton[100]={0.0};
  Double_t DirectPhotonQGP[100]={0.0};
  Double_t DirectPhotonHadron[100]={0.0};

  Double_t DirectPhotonPQCD[100]={0.0};

  cout<<"Pt:  "<<"      "<<"Hadron  "<<"QGP  "<<"RatePhoton"<<endl;
  for(int i =0;i<NPt;i++)
    {
      Pt[i]=PtMin+i*PtStep;
      
      DirectPhotonPQCD[i] = RateQCD_IntTau(Pt[i]);
      
      DirectPhotonQGP[i] = RateQGP_IntTau(Pt[i],R040);
      DirectPhotonHadron[i] = RateHadron_IntTau(Pt[i],R040);
      

      DirectPhoton[i] = DirectPhotonPQCD[i] + DirectPhotonQGP[i] + DirectPhotonHadron[i];
      
      cout<<Pt[i]<<"    "<<DirectPhotonHadron[i]<<"    "<<DirectPhotonQGP[i]<<"   "<<DirectPhoton[i]<<endl;
    }  

  TGraph *grfPhotonPQCD = new TGraph(NPt,Pt,DirectPhotonPQCD);
  grfPhotonPQCD->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfPhotonPQCD->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfPhotonPQCD->GetYaxis()->SetTitleOffset(1.4);
  grfPhotonPQCD->SetLineColor(8);



  new TCanvas;
  gPad->SetTicks(1);
  gPad->SetLogy(1);
  grfPhotonPQCD->Draw("AL");

  //  return;


  TGraph *grfDirectPhotonQGP = new TGraph(NPt,Pt,DirectPhotonQGP);
  grfDirectPhotonQGP->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonQGP->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonQGP->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhotonQGP->SetLineColor(2);


  TGraph *grfDirectPhotonHadron = new TGraph(NPt,Pt,DirectPhotonHadron);
  grfDirectPhotonHadron->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonHadron->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonHadron->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhotonHadron->SetLineColor(4);


  TGraph *grfDirectPhoton = new TGraph(NPt,Pt,DirectPhoton);
  grfDirectPhoton->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhoton->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhoton->GetYaxis()->SetTitleOffset(1.4);
  grfDirectPhoton->SetLineColor(1);


  TLegend *legd5 = new TLegend( 0.60,0.70,0.82,0.85);
  legd5->SetBorderSize(0);
  legd5->SetFillStyle(0);
  legd5->SetFillColor(0);
  legd5->SetTextSize(0.040);

  legd5->AddEntry(grfDirectPhoton,"Sum","L");
  legd5->AddEntry(grfDirectPhotonQGP,"QGP","L");
  legd5->AddEntry(grfDirectPhotonHadron,"Hadron","L");
  legd5->AddEntry(grfPhotonPQCD,"pQCD","L");

  new TCanvas;
  gPad->SetTicks();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.19);
  Draw_ALICE_DirectPhotonRate040_Pt(legd5);

  grfDirectPhoton->SetLineWidth(2);
  grfDirectPhoton->SetLineColor(1);
  grfDirectPhoton->Draw("Csame");
  grfDirectPhotonHadron->Draw("Csame");
  grfDirectPhotonQGP->Draw("Csame");
  grfPhotonPQCD->Draw("Csame");
  legd5->Draw("Lsame");

  gPad->SaveAs("LHC_DirecPhoton.pdf");
  gPad->SaveAs("LHC_DirecPhoton.eps");
  gPad->SaveAs("LHC_DirecPhoton.gif");
  gPad->SaveAs("LHC_DirecPhoton.png");

  grTempVsTauAnaC->Write(); 
  grFQGPVsTauAnaC->Write();  
  grfPhotonPQCD->Write();  
  grfDirectPhotonQGP->Write();  
  grfDirectPhotonHadron->Write();  
  grfDirectPhoton->Write();  
  

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
    
    cout<<"Lattice: Ntime "<<Ntime<<" ss "<<ssTau<<"  "<<tau[Ntime]<<"   "<<Temp[Ntime]<<"  "<<h[Ntime]<<endl;
   
    Ntime++;
  }while(CutTemp > TF);
  
  return 0;

}

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



//================================ Direct Photon Functions =======================================//

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
  for(int i =0;i<=NEtaf;i++)
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
  RQGP1 = Const2*Alpha*AlphaS*T*T*TMath::Exp(-(E/T));

  Double_t RQGP3 =0.0; 
  Double_t Const3=0.0105;
  RQGP1 = Const3*Alpha*AlphaS*E*T*TMath::Exp(-(E/T));

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
  for(int i =0;i<=NEtaf;i++)
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



Double_t RateQCD_IntTau(Double_t Pt)
{
  /*
  Double_t Sum =0.0;
  Double_t IntEtaf = 2.0;

  for(int i =0;i<=Ntime;i++)
    {
      Sum = Sum + h[i]*tau[i]; 

    }

  return Sum*IntEtaf*RateQCD(Pt)*Steptime;
  */

  Double_t xx =0.0;
  xx = RateQCD(Pt);
  return xx;

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
  //From ALICE Paper 

  Double_t APt[30]={2.14,   2.38,   2.59,   2.86,   3.15,   3.48,   3.81,   4.16,   4.51,   4.86,   5.25,   5.62,   6.07,   
		   6.58,   7.00,   7.47,   7.90,   8.35,   8.83,   9.34,   9.71,   10.24,   10.74,   11.27,   11.83,   
		   12.30,   12.59,   13.08,   13.55,   13.94};

  Double_t D2NByPtDPtDy[30]={0.009903112,   0.006549429,   0.003952514,   0.002736074,   0.001650603,   0.001091140,   0.000721304,   
			     0.000499134,   0.000345396,   0.000228305,   0.000150882,   0.000114419,   0.000079141,   0.000049934,   
			     0.000039631,   0.000027410,   0.000021754,   0.000016491,   0.000011405,   0.000008643,   0.000007520,   
			     0.000005444,   0.000004733,   0.000003426,   0.000002845,   0.000002258,   0.000001793,   0.000001559,   
			     0.000001295,   0.000001180};


  TGraph *Grf_QCDRate = new TGraph(30,APt,D2NByPtDPtDy);



  Double_t Rate =0.0;
  Rate =Grf_QCDRate->Eval(Pt,0,"");

  return Rate;

}



// ===================================== Direct Photon Graphs (ALICE) =======================================//

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
  Draw_ALICE_DirectPhotonRate040_Pt(legend_ratio);  



  legend_ratio->Draw("same");
}





void Draw_ALICE_DirectPhotonRate040_Pt(TLegend *lgd)
{
  //12.5
  //0.0000063,
  //0.0000095
  

  Double_t Pt[16]=   {0.9,   1.1,   1.3,   1.5,   1.7,   1.9,   2.1,   2.3,   2.5,   2.8,   3.2,   3.7,   4.4,   5.5,   6.9,   9.3}; 


  Double_t ErrPtH[16]={0.0};
  Double_t ErrPtL[16]={0.0};
  
  
  Double_t D2NDPtDy[16]={2.76551558933,1.10580368159,0.661689814073,0.493215247285,0.183313311224,0.0880422051384,
			 0.0526826156565,0.0303885496417,0.0117152618748,0.00841544596557,0.00301253056964,
			 0.000965961550939,0.000309505281634,0.000128003883749,3.66497554435e-05,6.49417470134e-06}; 
  
  Double_t StatErrHD2NDPtDy[16]={4.29225854121,1.65473289898,0.954488159289,0.59238882259,0.236909598955,0.12243282222,
				 0.0706260847865,0.0422587767985,0.0181838730484,0.0108765314076,0.00403835633715,0.00139340041677,
				 0.00043035391613,0.000171601448295,5.09685944456e-05,9.71793159036e-06};
  
  Double_t StatErrLD2NDPtDy[16]={1.48361294486,0.66221532887,0.475879234702,0.38169961342,0.127101825384,0.0610378802333,0.0365258631492,
				 0.0210677668804,0.00628487296794,0.00583492020574,0.0022474131348,0.000599985742994,0.000192242389131,
				 7.66469875345e-05,2.54114160376e-05,3.48431752166e-06}; 
  
  Double_t SystErrHD2NDPtDy[16]={4.96932995217,1.98723399891,1.06536392276,0.66120220947,0.274295855142,0.136662678514,
				 0.0788346600601,0.0454710743853,0.0210534356969,0.0121406652902,0.00467590686561,0.00167348440855,
				 0.00053614300058,0.000191535107406,5.48398582896e-05,1.12508620022e-05};
  
  Double_t SystErrLD2NDPtDy[16]={0.856316192661,0.295816657337,0.306575706391,0.274482569761,0.0881171249971,0.0407942381631,0.0272490860105,
				 0.017543733603,0.00390371042577,0.00404523217165,0.00167652386416,0.000372689063311,8.58761470356e-05,
				 6.15303610473e-05,1.82745438193e-05,2.41560687e-06}; 
  
 
  Double_t SystErrD2NDPtDy[16]={0.0};

  for(int j=0;j<16;j++){
    
    cout<<Pt[j]<<"  "<<D2NDPtDy[j]<<"   "<<StatErrHD2NDPtDy[j]<<"   "<<SystErrHD2NDPtDy[j]<<endl;

      
      }
  

  //gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);

  TGraphAsymmErrors *PhotonRateALICEdata = new TGraphAsymmErrors(16,Pt,D2NDPtDy,ErrPtL,ErrPtH,  StatErrLD2NDPtDy, StatErrHD2NDPtDy);
  PhotonRateALICEdata->SetMarkerStyle(20);
  PhotonRateALICEdata->SetMarkerColor(kBlack+0);
  PhotonRateALICEdata->SetMarkerSize(1.6);
  PhotonRateALICEdata->GetYaxis()->SetTitleOffset(1.8);
  PhotonRateALICEdata->GetYaxis()->SetRangeUser(0.0000001,1000.0);
  
  TAxis *Xaxis1 = PhotonRateALICEdata->GetXaxis();
  Xaxis1->SetLimits(0.0,14.0);
   

  PhotonRateALICEdata->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  PhotonRateALICEdata->GetYaxis()->SetTitle("#frac{1}{2#pi p_{T}}#frac{d^{2}N}{dp_{T}dy}[(GeV/c)^{-2}]");
  
  PhotonRateALICEdata->Draw("AP");

  
  TBox *SystErrPhotonRate[16];

  Double_t PtBound[17]=   {0.8, 1.0,   1.2,   1.4,   1.6,   1.8,   2.0,   2.2,   2.4,   2.6,   3.0,   3.5,   4.0,   5.0,   6.0,   8.0,   11.0}; 

  for(int j=0;j<16;j++){

    //SystErrD2NDPtDy[j]=(SystErrHD2NDPtDy[j]-SystErrLD2NDPtDy[j]);

    SystErrPhotonRate[j] = new TBox(PtBound[j], D2NDPtDy[j]-SystErrLD2NDPtDy[j], PtBound[j+1], D2NDPtDy[j]+SystErrHD2NDPtDy[j]);
  }
  
  for(int j=0;j<16;j++){
    SystErrPhotonRate[j]->SetFillStyle(9001);
    SystErrPhotonRate[j]->SetLineColor(kBlack);
    SystErrPhotonRate[j]->Draw("same"); 
  }
  


  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.040);
  tb->DrawLatex(0.22,0.20,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
  tb->DrawLatex(0.22,0.15,"0-40%");
 
  //lgd->AddEntry(PhotonRateALICEdata,"Inclusive Photons", "P");
  
  
  
}

