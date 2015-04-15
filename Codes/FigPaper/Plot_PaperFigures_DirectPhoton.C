#ifndef __CINT__
#endif
#include "TLatex.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TMath.h"
#include "TVector3.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include <fstream>
#include <map>
#include <string>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "TObjArray.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TColor.h"
#include "TStyle.h"
#include <iostream>
#include <iomanip>
#include "TGraphAsymmErrors.h"


Double_t Npart(int BinLow, int BinHigh);
Double_t NColl(int BinLow, int BinHigh);

void Plot_PaperFigures_DirectPhoton() 

{
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
  gStyle->SetMarkerSize(1.3);

  gROOT->ForceStyle();
  

  TLegend *legd0 = new TLegend(0.16,0.64,0.97,0.93);
  legd0->SetBorderSize(0);
  legd0->SetFillStyle(0);
  legd0->SetFillColor(0);
  legd0->SetTextSize(0.050);
  

  TLatex *tb= new TLatex;
  tb->SetNDC(); 
  tb->SetTextAlign(12);
  tb->SetTextColor(1);
  tb->SetTextSize(0.045);

  TFile *file_RHIC_DirecPhoton = new TFile("InRootFiles/RHIC_DirecPhoton.root","R");
  

  TGraph *RHIC_grTempVsTauLatt =(TGraph*)file_RHIC_DirecPhoton->Get("grTempVsTauAnaC");
  RHIC_grTempVsTauLatt->SetLineColor(1);
  RHIC_grTempVsTauLatt->SetLineStyle(1);
  RHIC_grTempVsTauLatt->SetLineWidth(2);
  RHIC_grTempVsTauLatt->GetXaxis()->SetTitleSize(0.06);
  RHIC_grTempVsTauLatt->GetYaxis()->SetTitleSize(0.06);  

  TGraph *RHIC_grFQGPVsTauLatt =(TGraph*)file_RHIC_DirecPhoton->Get("grFQGPVsTauAnaC");
  RHIC_grFQGPVsTauLatt->SetLineColor(1);
  RHIC_grFQGPVsTauLatt->SetLineStyle(1);
  RHIC_grFQGPVsTauLatt->SetLineWidth(2);
  RHIC_grFQGPVsTauLatt->GetXaxis()->SetTitleSize(0.06);
  RHIC_grFQGPVsTauLatt->GetYaxis()->SetTitleSize(0.06);
    

  TFile *file_LHC_DirecPhoton = new TFile("InRootFiles/LHC_DirecPhoton.root","R");
  

  TGraph *LHC_grTempVsTauLatt =(TGraph*)file_LHC_DirecPhoton->Get("grTempVsTauAnaC");
  LHC_grTempVsTauLatt->SetLineColor(2);
  LHC_grTempVsTauLatt->SetLineStyle(9);
  LHC_grTempVsTauLatt->SetLineWidth(2);
  LHC_grTempVsTauLatt->GetXaxis()->SetTitleSize(0.06);
  LHC_grTempVsTauLatt->GetYaxis()->SetTitleSize(0.06); 


  TGraph *LHC_grFQGPVsTauLatt =(TGraph*)file_LHC_DirecPhoton->Get("grFQGPVsTauAnaC");
  LHC_grFQGPVsTauLatt->SetLineColor(2);
  LHC_grFQGPVsTauLatt->SetLineStyle(9);
  LHC_grFQGPVsTauLatt->SetLineWidth(2);
  LHC_grFQGPVsTauLatt->GetXaxis()->SetTitleSize(0.06);
  LHC_grFQGPVsTauLatt->GetYaxis()->SetTitleSize(0.06);
    

 
  TLegend *legdA = new TLegend( 0.22,0.71,0.79,0.85);
  legdA->SetBorderSize(0);
  legdA->SetTextSize(0.040);
  legdA->SetFillStyle(0);
  legdA->SetFillColor(0);
  legdA->AddEntry(LHC_grTempVsTauLatt,"LHC","L");
  legdA->AddEntry(RHIC_grTempVsTauLatt,"RHIC","L");
  

  new TCanvas;
  gPad->SetTicks(1);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.16);
  LHC_grTempVsTauLatt->GetXaxis()->SetTitleOffset(1.1);
  LHC_grTempVsTauLatt->GetYaxis()->SetTitleOffset(1.1);

  LHC_grTempVsTauLatt->GetYaxis()->SetRangeUser(0.0,0.7);
  LHC_grTempVsTauLatt->Draw("AC");
  RHIC_grTempVsTauLatt->Draw("Csame");
  legdA->Draw("same");
  tb->DrawLatex(0.22,0.90,"Lattice EOS, Cylindrical expansion");
    
  gPad->SaveAs("plots/Fig1a_TempVsTauLatt.pdf");
  gPad->SaveAs("plots/Fig1a_TempVsTauLatt.eps");
  gPad->SaveAs("plots/Fig1a_TempVsTauLatt.gif");
  gPad->SaveAs("plots/Fig1a_TempVsTauLatt.png");








  new TCanvas;
  gPad->SetTicks(1);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.16);
  
  LHC_grFQGPVsTauLatt->GetXaxis()->SetTitleOffset(1.1);
  LHC_grFQGPVsTauLatt->GetYaxis()->SetTitleOffset(1.1);

  LHC_grFQGPVsTauLatt->GetYaxis()->SetRangeUser(0.0,1.6);
  LHC_grFQGPVsTauLatt->Draw("AL");
  RHIC_grFQGPVsTauLatt->Draw("Lsame");
  legdA->Draw("same");
  tb->DrawLatex(0.22,0.90,"Lattice EOS, Cylindrical expansion");
    
  gPad->SaveAs("plots/Fig1b_FQGPVsTauLatt.pdf");
  gPad->SaveAs("plots/Fig1b_FQGPVsTauLatt.eps");
  gPad->SaveAs("plots/Fig1b_FQGPVsTauLatt.gif");
  gPad->SaveAs("plots/Fig1b_FQGPVsTauLatt.png");







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

