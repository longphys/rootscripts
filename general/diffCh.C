#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLegend.h"
#include "TRandom3.h"

void diffCh()
{
  //! Files and trees
  TFile* fsimCs137 = new TFile("withAl_simCs137_new2.root", "read");
  // TFile* fsimNa22 = new TFile("withAl_simNa22.root", "read");
  TFile* fexpCs137 = new TFile("old_divider_Cs137_no_coinc_ch2HV1746_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "read");
  // TFile* fexpNa22 = new TFile("old_divider_Na22_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "read");

  TTree* tsimCs137 =  (TTree*) fsimCs137->Get("dEEtree");
  // TTree* tsimNa22 =  (TTree*) fsimNa22->Get("dEEtree");
  TTree* texpCs137 =  (TTree*) fexpCs137->Get("AnalysisxTree");
  // TTree* texpNa22 =  (TTree*) fexpNa22->Get("AnalysisxTree");

  //! Calibration
  // double a = 0.002633333333333; //first
  // double b = -0.254766666666667;

  // double a = 0.0027; //second
  // double b = -0.254766666666667;

  double a = 0.0029; //fourth
  double b = -0.33;

  //! histogram options
  int binSim = 1300;
  int binExp = 700;
  double xminSim = 0.;
  double xmaxSim = 1.3;
  double xminExp = 100;
  double xmaxExp = 800;
  double xminCal = xminExp*a+b;
  double xmaxCal = xmaxExp*a+b;

  TH1D* hsimCs137 = new TH1D("hsimCs137", "Cs137 Simulation", binExp, xminCal, xmaxCal);
  // TH1D* hsimNa22 = new TH1D("hsimNa22", "Na22 Simulation", binExp, xminCal, xmaxCal);
  TH1D* hexpCs137_0 = new TH1D("hexpCs137_0", "Cs137 Experiment", binExp, xminExp, xmaxExp);
  TH1D* hexpCs137_1 = new TH1D("hexpCs137_1", "Cs137 Experiment", binExp, xminExp, xmaxExp);
  TH1D* hexpCs137_2 = new TH1D("hexpCs137_2", "Cs137 Experiment", binExp, xminExp, xmaxExp);
  // TH1D* hexpNa22 = new TH1D("hexpNa22", "Na22 Experiment", binExp, xminExp, xmaxExp);

  hsimCs137->GetXaxis()->SetTitle("Energy(MeV)");
  hsimCs137->GetYaxis()->SetTitle("Count");
  // hsimNa22->GetXaxis()->SetTitle("Energy(MeV)");
  // hsimNa22->GetYaxis()->SetTitle("Count");

  //! Fill experiment histograms
  UShort_t x_exp[48], y_exp[48];

  texpCs137->SetBranchAddress("NeEvent.neutAmp[48]", x_exp);
  // texpNa22->SetBranchAddress("NeEvent.neutAmp[48]", y_exp);

  Long64_t entries = texpCs137->GetEntries();

  for(int i = 0; i < entries; i++){
    texpCs137->GetEntry(i);
    hexpCs137_0->Fill(x_exp[0]);
    hexpCs137_1->Fill(x_exp[1]);
    hexpCs137_2->Fill(x_exp[2]);
  }
  
  // entries = texpNa22->GetEntries();

  // for(int i = 0; i < entries; i++){
  //   texpNa22->GetEntry(i);
  //   hexpNa22->Fill(y_exp[0]);
  // }

  //! Differentiate
  TCanvas *cDiff = new TCanvas("cDiff", "cDiff", 1600, 800);
  cDiff->cd();
  cDiff->Divide(3, 1);

  double binContent, binContentPlus, binLength, newBinContent;

  // hcalNa22->Rebin(4);
  // hcalCs137->Rebin(4);
  // binExp = binExp/4;

  binLength = (xmaxExp - xminExp)/ binExp;

  // TH1D* hDifexpNa22 = new TH1D("hDifexpNa22", "Na22 Exp Differentiated", binExp, xminExp, xmaxExp);
  TH1D* hDifexpCs137 = new TH1D("hDifexpCs137", "Cs137 Exp Differentiated", binExp, xminExp, xmaxExp);

  // for(int i = 1; i < binExp; i++)
  // {
    // binContent = hexpNa22->GetBinContent(i);
    // binContentPlus = hexpNa22->GetBinContent(i+1);
    // newBinContent = (binContentPlus - binContent)/binLength;
    // hDifexpNa22->SetBinContent(i, newBinContent);

    // binContent = hexpCs137->GetBinContent(i);
    // binContentPlus = hexpCs137->GetBinContent(i+1);
    // newBinContent = (binContentPlus - binContent)/binLength;
    // hDifexpCs137->SetBinContent(i, newBinContent);
  // }

  // cDiff->cd(1);
  // hexpNa22->Draw();

  // cDiff->cd(3);
  // hDifexpNa22->SetLineColor(kRed);
  // hDifexpNa22->Draw();

  // cDiff->cd(2);
  // hexpCs137->Draw();

  // cDiff->cd(4);
  // hDifexpCs137->SetLineColor(kRed);
  // hDifexpCs137->Draw();

  cDiff->cd(1);
  hexpCs137_0->SetLineColor(kRed);
  hexpCs137_1->SetLineColor(kRed);
  hexpCs137_2->SetLineColor(kRed);
  hexpCs137_0->Draw();
  cDiff->cd(2);
  hexpCs137_1->Draw();
  cDiff->cd(3);
  hexpCs137_2->Draw();
}