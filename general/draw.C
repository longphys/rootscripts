#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"

void draw()
{
  //! Files and trees
  TFile* fileOpen = new TFile("save_new_passive_divider_Cs137_no_coinc_ch0HV1925_ch2HV1550_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "read");

  TTree* treeOpen =  (TTree*) fileOpen->Get("tree");

  //! histogram options
  int binExp = 700;
  double xminExp = 100;
  double xmaxExp = 800;

  TH1D* hist = new TH1D("hist", "Exp data", 2000, 0, 2000);

  //! Fill experiment histograms
  // UShort_t x_exp[48];
  // treeOpen->SetBranchAddress("NeEvent.neutAmp[48]", x_exp);

  // Long64_t entries = treeOpen->GetEntries();
  // for(int i = 0; i < entries; i++){
  //   treeOpen->GetEntry(i);
  //   treeSave->Fill();
  // }

  //! Canvas
  TCanvas *canvas = new TCanvas("canvas", "Exp data show", 700, 500);
  canvas->cd();

  hist->SetLineColor(kRed);
  treeOpen->Draw("Ch0>>hist");

  // fileOpen->Close();

  TFile* fileSave = new TFile("out.root", "recreate");

  TH1D* histSave = new TH1D();
  histSave = (TH1D*) hist->Clone();

  fileSave->Write();
  fileSave->Close(); 
}