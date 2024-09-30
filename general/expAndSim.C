#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLegend.h"

void expAndSim()
{
  TFile* file1 = new TFile("./simfiles/withAl_simCs137_1.root", "read");
  TFile* file2 = new TFile("./simfiles/withAl_simNa22_2.root", "read");
  // TFile* file3 = new TFile("./expfiles/new/stilbene_dividers_Plastic_ch2_Cs137_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0001.root", "read");
  TFile* file3 = new TFile("./expfiles/time/6/plastic_detectors_Cs137_ch0HV1899_ch1HV1915_ch2HV1950_ch0_9668_ch1_9806_ch2_9852_run_0_0001.root", "read");
  TFile* file4 = new TFile("./expfiles/new/stilbene_dividers_Plastic_ch2_Na22_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0001.root", "read");

  TTree* tree1 =  (TTree*) file1->Get("dEEtree");
  TTree* tree2 =  (TTree*) file2->Get("dEEtree");
  TTree* tree3 =  (TTree*) file3->Get("AnalysisxTree");
  TTree* tree4 =  (TTree*) file4->Get("AnalysisxTree");

  // histogram options
  int binSim = 1300;
  int binExp = 700;
  double xminSim = 0.;
  double xmaxSim = 1.3;
  double xminExp = 100;
  double xmaxExp = 800;

  TH1D* hist1 = new TH1D("hist1", "Cs137 Simulation", binSim, xminSim, xmaxSim);
  TH1D* hist2 = new TH1D("hist2", "Na22 Simulation", binSim, xminSim, xmaxSim);
  TH1D* hist3 = new TH1D("hist3", "Cs137 Experiment", binExp, xminExp, xmaxExp);
  TH1D* hist4 = new TH1D("hist4", "Na22 Experiment", binExp, xminExp, xmaxExp);

  TCanvas* c1 = new TCanvas("c1", "c1", 1500, 1000);
  int cols = 2;
  int rows = 2;
  c1->cd();
  c1->Divide(cols, rows);

  c1->cd(1);
  c1->cd(1)->SetGridx();
  c1->cd(1)->SetGridy();
  hist1->SetLineColor(kBlack);
  hist1->GetXaxis()->SetTitle("Energy(Mev)");
  hist1->GetYaxis()->SetTitle("Count");
  
  tree1->Draw("Scintillator>>hist1", "Scintillator>0");

  c1->cd(2);
  c1->cd(2)->SetGridx();
  c1->cd(2)->SetGridy();
  hist2->SetLineColor(kBlack);
  hist2->GetXaxis()->SetTitle("Energy(Mev)");
  hist2->GetYaxis()->SetTitle("Count");
  
  tree2->Draw("Scintillator>>hist2", "Scintillator>0");

  c1->cd(3);
  c1->cd(3)->SetGridx();
  c1->cd(3)->SetGridy();
  hist3->GetXaxis()->SetTitle("Channel");
  hist3->GetYaxis()->SetTitle("Count");
  
  tree3->Draw("NeEvent.neutAmp[0]>>hist3", "");

  TF1 *fEnergySpread = new TF1("function", "[0]*x+[1]", 250, 300);
  hist3->Fit(fEnergySpread, "", "", 274, 286);

  c1->cd(4);
  c1->cd(4)->SetGridx();
  c1->cd(4)->SetGridy();
  hist4->GetXaxis()->SetTitle("Channel");
  hist4->GetYaxis()->SetTitle("Count");
  
  tree4->Draw("NeEvent.neutAmp[0]>>hist4", "");

  TCanvas* newT = new TCanvas("newT", "Simulation and Measurement", 1500, 750);
  newT->Divide(2,1);
  newT->cd(1);
  hist1->SetLineColor(kBlue);
  tree1->Draw("Scintillator>>hist1", "Scintillator>0");
  newT->cd(2);
  hist3->SetLineColor(kRed);
  tree3->Draw("NeEvent.neutAmp[0]>>hist3", "NeEvent.neutAmp[0]>0");

}