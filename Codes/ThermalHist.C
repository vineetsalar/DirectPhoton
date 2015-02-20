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
#include "TLine.h"
#include <TObjArray.h>
#include <TTree.h>
#include "stdio.h"

Double_t pi = TMath::Pi();
Double_t pi2 = pi*pi;
Double_t pi3 = pi*pi2;
Double_t pi4 = pi2*pi2;
Double_t hbarc = 0.197327;
Double_t hbarc4 = hbarc*hbarc*hbarc*hbarc;

Double_t Radius = 1.2*TMath::Power(208.0, 1/3.);
Double_t contime = Radius*Radius*pi;
Double_t conrate = 1.0/(137.0*137.0*2.0*pi3);
Double_t aq = 37.0*pi2/90.0;

Double_t EDNDY=1.5 * 1600;

//from ArXive 0705.1591

Double_t tau_i = 0.1;
//Double_t tau_i = 0.03;
Double_t Temp_i =  hbarc*TMath::Power((3.6*EDNDY /(tau_i*contime*4*aq)),1.0/3.0);
Double_t Temp_C = 0.160;
Double_t Temp_f = 0.130;


Double_t EtaStart = -1;
Double_t EtaEnd = 1;


// Temp and Hadronic fraction as a function of time 
Double_t tau[10000], Temp[10000], h[10000];
Int_t Ntime;

void CalTime();
Double_t GetYdist(Double_t Y); 
Double_t InteTimeY(Double_t Mass, Double_t Y);
Double_t InteEtaPT(Double_t Mass, Double_t Y, Double_t Temp);

Double_t InteTime(Double_t Mass);
Double_t InteEtaPTY(Double_t Mass, Double_t Temp);

Double_t InteTime(Double_t Mass, Double_t PT); 
Double_t InteEtaY(Double_t Mass, Double_t PT, Double_t Temp);


void ThermalHist()
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
  gStyle->SetCanvasDefH(800);
  gStyle->SetCanvasDefW(800);  
  
  
  cout<<"Temp_i "<<Temp_i<<" tau_i "<<tau_i<<endl;
  // Calculate Temp and Hadronic fraction as a function of time
  CalTime();
  //////////

  int nbins = 100;
  Double_t step = 0.1;

  TH1F *diMuonsInvMass = new TH1F("diMuonsInvMass","diMuonsInvMass", nbins,0.0,10.0);
  TH1F *diMuonsPt = new TH1F("diMuonsPt","diMuonsPt", 100,0.0,10.0);
  TH1F *diMuonsRap = new TH1F("diMuonsRap","diMuonsRap", 100,-6,6);
  TH2F *diMuonsMassPt = new TH2F("diMuonsMassPt","diMuonsMassPt", nbins,0.0,10.0, nbins, 0.0, 10.0);

  // Pt and rapidity graphs





  // Output spectra file 
  TFile *filespectra = new TFile("thermalhist.root", "recreate");

  //**********************************************************************


  ////  Calculate dN/dM
  Double_t Mass[1000], DNDM[1000];
  Double_t sumMass =0.0;

  for(int i = 1; i <= nbins; ++i){
    Mass[i] = diMuonsInvMass->GetBinCenter(i);
    DNDM[i] = InteTime(Mass[i]);

    diMuonsInvMass->SetBinContent(i,DNDM[i]);

    cout << Mass[i] <<"       " << DNDM[i] << endl;
    sumMass = sumMass + DNDM[i]*step;
  }

  cout << "Total Integral of Mass Distribution  = " << sumMass << endl;

  new TCanvas;
  gPad->SetLogy(1);
  gPad->SetTickx();
  gPad->SetTicky();
  diMuonsInvMass->SetLineColor(2);
  diMuonsInvMass->Draw();
 
  
  TGraph *grfPt = new TGraph(nbins, Mass, DNDM);
  grfPt->SetName("pt");
  grfPt->GetYaxis()->SetTitle("dN/dp_{T}(GeV/c)^{-1}");
  grfPt->GetYaxis()->SetTitleOffset(1.2);
  grfPt->GetYaxis()->CenterTitle();
  grfPt->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  //grfPt->GetYaxis()->SetRangeUser(0.0000000001,0.002);
  grfPt->GetXaxis()->CenterTitle();
  new TCanvas;
  gPad->SetLogy(1);
  gPad->SetTickx();
  gPad->SetTicky();
  grfPt->Draw("ALP");


  cout<<" integral of mass graph as pT: "<<grfPt->Integral()<<endl;





  cout << " Integral of Mass Histogram  =  "<<diMuonsInvMass->Integral("width") << endl;


  // Fill a 2D histogram of Mass and pT

  Double_t mass, apt, Yield;

  double Integ = 0.0;

  for(int i = 1; i <= nbins; ++i){
    mass = diMuonsMassPt->GetXaxis()->GetBinCenter(i);
    double sumPT =0.0;
    for(int j = 1; j <= nbins; ++j){
      
      apt = diMuonsMassPt->GetYaxis()->GetBinCenter(j);
      Yield = InteTime(mass, apt);
      diMuonsMassPt->SetBinContent(i, j, Yield);
      
      //cout<<" pt 2 d : " <<diMuonsMassPt->GetBinContent(j)<<endl;
      
      sumPT = sumPT+Yield*step;
    }
    Integ = Integ + sumPT*step;
    cout << mass <<"       " << sumPT << endl;
  }


  double gPt[100];
  double gDnDpt[100];

for(int i = 1; i <= nbins; ++i){

  gPt[i]=diMuonsMassPt->ProjectionX()->GetBinCenter(i);
  gDnDpt[i]=diMuonsMassPt->ProjectionX()->GetBinContent(i);

 }
 
/* TGraph *grfPt = new TGraph(nbins, gPt, gDnDpt);
 grfPt->SetName("pt");
 grfPt->GetYaxis()->SetTitle("dN/dp_{T}(GeV/c)^{-1}");
 grfPt->GetYaxis()->SetTitleOffset(1.2);
 grfPt->GetYaxis()->CenterTitle();
 grfPt->GetXaxis()->SetTitle("p_{T}(GeV/c)");
 grfPt->GetYaxis()->SetRangeUser(0.0000000001,0.2);
 grfPt->GetXaxis()->CenterTitle();
 new TCanvas;
 gPad->SetLogy(1);
 gPad->SetTickx();
 gPad->SetTicky();
 grfPt->Draw("ALP");*/
 
 cout << " Integral from 2D data sum = " << Integ << endl;
 cout << " Integral from 2D histo = " << diMuonsMassPt->Integral("width") << endl;
 
 
  new TCanvas;
  gPad->SetTickx();
  gPad->SetTicky();
  diMuonsMassPt->Draw("lego");
  

  /////////////////////////////////////////////////////////////////////////////////////////////////
  //////Calculate dn/dy
  cout << endl <<" dN/dY " << endl << endl;

  double yyg[100];
  double dndyg[100];
  for(Int_t i=1; i <= nbins; i++) {
    double yy = diMuonsRap->GetBinCenter(i);
    yyg[i]=yy;
    double dndy = GetYdist(yy);
    dndyg[i]=dndy;
    diMuonsRap->SetBinContent(i,dndy);
    cout << yy <<"       " << dndy << endl;
    //cout <<yyg[i]  <<" xxxx   " <<dndyg[i]<< endl;
  }
  
  TGraph *grfRap = new TGraph(nbins, yyg, dndyg);
  grfRap->SetName("Rap");
  grfRap->GetYaxis()->SetTitle("dN/dy");
  grfRap->GetYaxis()->SetTitleOffset(1.2);
  grfRap->GetYaxis()->CenterTitle();
  grfRap->GetXaxis()->SetTitle("y");
  //grfRap->GetYaxis()->SetRangeUser(0.000001,0.001);
  grfRap->GetXaxis()->CenterTitle();
  new TCanvas;
  gPad->SetLogy(1);
  gPad->SetTickx();
  gPad->SetTicky();
  grfRap->Draw("ALP");
  
  cout<<" integral of rap graph : "<< grfRap->Integral()<<endl;

  new TCanvas;
  gPad->SetTickx();
  gPad->SetTicky(); 
  diMuonsRap->SetLineColor(2);
  diMuonsRap->Draw();

  
  ///// Write the spectra 
  grfPt->Write();
  grfRap->Write();
  
  diMuonsInvMass->Write();
  //diMuonsPt->Write();
  diMuonsRap->Write();
  diMuonsMassPt->Write();

  filespectra->Close();
  return;
}


void CalTime() {
  
  // Degrees of freedom in QGP and hadron phases
  Double_t aq = 37.0*pi2/90.0;
  Double_t ah = 3.0*pi2/90.0;
  
  // times at different stages
  Double_t tau_C = TMath::Power(Temp_i/Temp_C, 3)*tau_i;
  Double_t tau_h = tau_C*aq/ah;
  Double_t tau_f = TMath::Power(Temp_C/Temp_f, 3)*tau_h;

  cout << " times " << endl;
  cout << tau_C << "    " << tau_h << "   " << tau_f << endl << endl;

  Double_t step = 0.1;

  // Quark phase
  //Double_t sumq = 0.0;
  
  Int_t Nq = (tau_C-tau_i)/step;
  for(int i = 0; i <= Nq; i++) {
    tau[i]  = tau_i + i*step;
    Temp[i] = Temp_i* TMath::Power(tau_i/tau[i], 1/3.0);
    h[i] = 0.0;
    //    cout << tau[i] << "   " << Temp[i] << "   " << h[i] << endl;
  }

  // Mixed phase
  // Double_t summ = 0.0;
 
  Int_t Nm = (tau_h-tau_C)/step;
  for(int i = 0; i <= Nm; i++) {
    int j = Nq+1+i;
    tau[j]  = tau_C + i*step;
    Temp[j] = Temp_C;
    h[j] = (aq/(aq-ah)) * (tau[j]-tau_C)/tau[j];
    //    cout << tau[j] << "   " << Temp[j] << "   " << h[j] << endl;
  }

  // hadron phase
  //Double_t sumh = 0.0;
  
  Int_t Nh = (tau_f-tau_h)/step;
  for(int i = 0; i <= Nh; i++) {
    int k = i+Nq+Nm+2;
    tau[k]  = tau_h + i*step;
    Temp[k] = Temp_C* TMath::Power(tau_h/tau[k], 1/3.0);
    h[k] = 1.0;
    //cout << tau[k] << "   " << Temp[k] << "   " << h[k] << endl;
  }
  Ntime = Nq+Nm+Nh+2;

  TGraph *grfTime = new TGraph(Ntime, tau, Temp);
 
  new TCanvas; 
  gPad->SetTickx();
  gPad->SetTicky();

  grfTime->Draw("ALP");

  return;
} 


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




