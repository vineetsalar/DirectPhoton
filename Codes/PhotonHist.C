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





Double_t EtaStart = -1;
Double_t EtaEnd = 1;



/*
//Medium Evolution Functions
int NTau; double stepTau;
double Tau[10000], TempTau[10000], fQGP[10000];
*/

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


//double StepTime;



/*
Double_t GetYdist(Double_t Y); 
Double_t InteTimeY(Double_t Mass, Double_t Y);
Double_t InteEtaPT(Double_t Mass, Double_t Y, Double_t Temp);
Double_t InteTime(Double_t Mass);
Double_t InteEtaPTY(Double_t Mass, Double_t Temp);
Double_t InteTime(Double_t Mass, Double_t PT); 
Double_t InteEtaY(Double_t Mass, Double_t PT, Double_t Temp);
*/


//Direct Photon
Double_t RatePhoton(Double_t RCent, Double_t Pt);
Double_t RateQGP_IntTau(Double_t Pt);
Double_t RateQGP_IntEtaf(Double_t Pt, Double_t T);
Double_t RateQGP(Double_t Pt, Double_t T, Double_t Etaf);

Double_t RateHadron_IntTau(Double_t Pt);
Double_t RateHadron_IntEtaf(Double_t Pt, Double_t T);
Double_t RateHadron(Double_t Pt, Double_t T, Double_t Etaf);

Double_t RateQCD_IntTau(Double_t Pt);
Double_t RateQCD(Double_t Pt);



void PhotonHist()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  //gStyle->SetPadTopMargin(0.10);
  //gStyle->SetPadBottomMargin(0.10);
  gStyle->SetPadLeftMargin(0.12);
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

  //Double_t RateQGP_IntTau(Double_t Pt)


  //Double_t R0MB = R05*TMath::Power(Npart(0,40)/Npart(0,2),0.5);

  cout<<"Pt:  "<<"      "<<"Hadron  "<<"QGP  "<<"RatePhoton"<<endl;
  for(int i =0;i<NPt;i++)
    {
      Pt[i]=PtMin+i*PtStep;
      
      DirectPhotonQGP[i] = pi*R040*R040*RateQGP_IntTau(Pt[i])/hbarc2;
      DirectPhotonHadron[i] = pi*R040*R040*RateHadron_IntTau(Pt[i])/hbarc2;

      DirectPhoton[i] = RatePhoton(R040,Pt[i])/hbarc2;
      


      cout<<Pt[i]<<"    "<<DirectPhotonHadron[i]<<"    "<<DirectPhotonQGP[i]<<"   "<<DirectPhoton[i]<<endl;
    }  


  TGraph *grfDirectPhotonQGP = new TGraph(NPt,Pt,DirectPhotonQGP);
  grfDirectPhotonQGP->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonQGP->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonQGP->GetYaxis()->SetTitleOffset(1.4);
  //grfDirectPhotonQGP->GetYaxis()->SetRangeUser(0.00001,100);
  grfDirectPhotonQGP->SetLineColor(2);


  TGraph *grfDirectPhotonHadron = new TGraph(NPt,Pt,DirectPhotonHadron);
  grfDirectPhotonHadron->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhotonHadron->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhotonHadron->GetYaxis()->SetTitleOffset(1.4);
  //grfDirectPhotonHadron->GetYaxis()->SetRangeUser(0.00001,100);
  grfDirectPhotonHadron->SetLineColor(1);




  TGraph *grfDirectPhoton = new TGraph(NPt,Pt,DirectPhoton);
  grfDirectPhoton->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  grfDirectPhoton->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T}dydp_{T}[GeV^{-2}c^{2}])");
  grfDirectPhoton->GetYaxis()->SetTitleOffset(1.4);
  //grfDirectPhoton->GetYaxis()->SetRangeUser(0.00001,100);
  
  TLegend *legd5 = new TLegend( 0.60,0.70,0.82,0.85);
  legd5->SetBorderSize(0);
  legd5->SetFillStyle(0);
  legd5->SetFillColor(0);
  legd5->SetTextSize(0.040);

  legd5->AddEntry(grfDirectPhoton,"Sum","L");
  legd5->AddEntry(grfDirectPhotonQGP,"QGP","L");
  legd5->AddEntry(grfDirectPhotonHadron,"Hadron","L");


  




  new TCanvas;
  gPad->SetTicks();
  gPad->SetLogy();
  grfDirectPhoton->SetLineWidth(2);
  grfDirectPhoton->SetLineColor(4);
  grfDirectPhoton->Draw("APL");
  grfDirectPhotonHadron->Draw("Lsame");
  grfDirectPhotonQGP->Draw("Lsame");
  legd5->Draw("Lsame");




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
    //tau[Ntime] = tau0 + stepTau*Ntime;
    tau[Ntime] = tau0 + Steptime*Ntime;

    Double_t VTau0Cent =  (R0Cent+0.5*aT*tau0*tau0)*(R0Cent+0.5*aT*tau0*tau0)*(z0+vZ*tau0)*pi;
    Double_t VTauCent =  (R0Cent+0.5*aT*tau[Ntime]*tau[Ntime])*(R0Cent+0.5*aT*tau[Ntime]*tau[Ntime])*(z0+vZ*tau[Ntime])*pi;
    
    Double_t ssTau = (hbarc3*ssCent)*(VTau0Cent/VTauCent);
    
    //Temptau[Ntime]=grfSSVsTemp->Eval(ssTau);
    Temp[Ntime]=grfSSVsTemp->Eval(ssTau);

    //Double_t tt= Temptau[Ntime];
    Double_t tt= Temp[Ntime];
    //fQGP[Ntime]=TempVsFQGP2->Eval(tt);
    Double_t fqgp_temp = TempVsFQGP2->Eval(tt);
    h[Ntime]= (1.0-fqgp_temp);
    //CutTemp = Temptau[Ntime];
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
Double_t RatePhoton(Double_t RCent, Double_t Pt)
{
  Double_t Area = pi*RCent*RCent;
  Double_t RPhoton =0.0;
  RPhoton =Area*(RateQGP_IntTau(Pt) + RateHadron_IntTau(Pt));

  return RPhoton;


}




Double_t RateQGP_IntTau(Double_t Pt)
{
  Double_t Sum =0.0;

  //Double_t tau[10000], Temp[10000], h[10000];
  //Int_t Ntime;
  //double Steptime;

  for(int i =0;i<=Ntime;i++)
    {
      Sum = Sum + RateQGP_IntEtaf(Pt,Temp[i])*(1-h[i])*tau[i]; 

    }

  return Sum*Steptime;

}


Double_t RateQGP_IntEtaf(Double_t Pt, Double_t T)
{
  Double_t Etaf=0.0;
  Double_t EtafMin = -1.0;
  Double_t EtafMax = 1.0;
  Double_t EtafStep = 0.0001;
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

  return RQGP;
}


// ================================= Rate Hadron ==============================//

Double_t RateHadron_IntTau(Double_t Pt)
{
  Double_t Sum =0.0;

  //Double_t tau[10000], Temp[10000], h[10000];
  //Int_t Ntime;
  //double Steptime;

  for(int i =0;i<=Ntime;i++)
    {
      Sum = Sum + RateHadron_IntEtaf(Pt,Temp[i])*h[i]*tau[i]; 

    }

  return Sum*Steptime*hbarc4;

}



Double_t RateHadron_IntEtaf(Double_t Pt, Double_t T)
{
  Double_t Etaf=0.0;
  Double_t EtafMin = -1.0;
  Double_t EtafMax = 1.0;
  Double_t EtafStep = 0.0001;
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
  Double_t Sum =0.0;
  Double_t IntEtaf = 2.0;

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
  Double_t a = -4.1506;
  Double_t b = -1.9845;
  Double_t c = 0.0744;
  Double_t d =-0.0383;

  Double_t RQCD =0.0;
  RQCD = TMath::Exp(a+b*Pt+c*Pt*Pt+d*Pt*Pt*Pt);
  return RQCD;


}


















/*

////////// Get dN/dy /////////////////////////////////

Double_t GetYdist(Double_t Y) {
  Double_t MassLow = 0.1;
  Double_t MassUp = 5.0;       
  Double_t Mstep = 0.1;
  Double_t Nmass = (MassUp - MassLow)/Mstep;
  double sum = 0;    
  for(Int_t i=0; i<=Nmass; i++) {
    Double_t massy = MassLow + i * Mstep;
    Double_t fun = InteTimeY(massy, Y);
    sum = sum + fun*Mstep;
  }
  return sum;
}

Double_t InteTimeY(Double_t Mass, Double_t Y) {
  // Quark and hadron form factors
  Double_t Fq = 5.0/9.0;
  Double_t Fh = 1.0/12.0;
  //Double_t Fh = 1.0;

  Double_t sum = 0.0;
  for(int i = 0; i < Ntime; i++) {
    Double_t step = tau[i+1] - tau[i];
    Double_t FF = Fq*(1-h[i]) + Fh*h[i];
    Double_t fun = FF*tau[i]*InteEtaPT(Mass, Y, Temp[i]);
    sum = sum + fun*step;
  }

  return contime*sum;
} 


Double_t InteEtaPT(Double_t Mass, Double_t Y, Double_t Temp) {
  Double_t step = 0.1;
  Int_t N = (EtaEnd - EtaStart)/step + 0.0001;
  Double_t sum = 0;

  for(int i = 0; i <= N; i++) {
    Double_t Eta = EtaStart + i*step;

    // Calculate Rate
    Double_t xx = Mass*TMath::CosH(Y-Eta)/Temp;
    if(xx > 100) xx=100;
    Double_t fun = Mass*Mass*Mass* (1.0/(xx*xx) + 1.0/xx)*exp(-xx);

    sum = sum + fun;
  }
  double val = sum*step*conrate;

  //  cout << Mass <<"   " << Y <<"   "<< val << endl;

  return val;
}



/////////////// dN/dM ///////////////

Double_t InteTime(Double_t Mass) {
  // Quark and hadron form factors
  Double_t Fq = 5.0/9.0;
  //Double_t Mrho,Mrho2,Wrho,Wrho2;
  //Mrho=1.02; Wrho=0.005; Mrho2=Mrho*Mrho; Wrho2=Wrho*Wrho;
  ///Fh=Mrho2*Mrho2/( (Mrho2 -Mass*Mass)*(Mrho2 -Mass*Mass) + Mrho2*Wrho2) /12.0;

  Double_t Fh = 1.0/12.0;
  //Double_t Fh = 1.0;

  Double_t sum = 0.0;
  for(int i = 0; i < Ntime; i++) {
    Double_t step = tau[i+1] - tau[i];
    Double_t FF = Fq*(1-h[i]) + Fh*h[i];
    Double_t fun = 0.0;

    fun = FF*tau[i]*InteEtaPTY(Mass, Temp[i]);

    sum = sum + fun*step;
  }
  return contime*sum;
} 

Double_t InteEtaPTY(Double_t Mass, Double_t Temp) {
  
  Double_t Y=0;
  Double_t step = 0.1;
  
  Int_t N = (EtaEnd - EtaStart)/step + 0.0001;
  Double_t sum = 0;

  for(int i = 0; i <= N; i++) {
    Double_t Eta = EtaStart + i*step;

    // Calculate Rate
    Double_t xx = Mass*TMath::CosH(Y-Eta)/Temp;
    if(xx > 100) xx=100;
    Double_t fun = Mass*Mass*Mass* (1.0/(xx*xx) + 1.0/xx)*exp(-xx);
    sum = sum + fun;
  }
  double val = sum*step*conrate;

  // Y Integration 
  val = 2.0*val;
  return val;
}



//////////  dN/(dMdpT) ///////////

Double_t InteTime(Double_t Mass, Double_t PT) {
  // Quark and hadron form factors
  Double_t Fq = 5.0/9.0;
  // Double_t Fh = 1.0/12.0;
  Double_t Fh = 1.0;

  Double_t sum = 0.0;
  for(int i = 0; i < Ntime; i++) {
    Double_t step = tau[i+1] - tau[i];
    Double_t FF = Fq*(1-h[i]) + Fh*h[i];
    Double_t fun = 0.0;

    fun = FF*tau[i]*InteEtaY(Mass, PT, Temp[i]);
    sum = sum + fun*step;
  }
  return contime*sum;
} 

Double_t InteEtaY(Double_t Mass, Double_t PT, Double_t Temp) {
  

  Double_t Y=0;
  Double_t step = 0.1;
  Int_t N = (EtaEnd - EtaStart)/step + 0.0001;
  Double_t sum = 0;
  
  //  Integrate over eta  
  for(int i = 0; i <= N; i++) {
    Double_t Eta = EtaStart + i*step;

    Double_t xx = TMath::Sqrt(Mass*Mass + PT*PT)*TMath::CosH(Y-Eta)/Temp;
    if(xx > 100) xx=100;
    Double_t funm = Mass*exp(-xx)*PT;
    sum = sum + funm;
  }
  
  double val = sum * conrate *step; 
  val = 2.0*val;   // Y Integral

  return val;
}

///////////////////
*/



