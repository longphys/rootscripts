#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"

#include "TStopwatch.h"

void refill()
{
  auto timer = new TStopwatch();
  timer->Start();

  //! Files and trees
  TFile* fileOpen = new TFile("./expfiles/usable/stilbene_divider_Co60_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "read");
  TTree* treeOpen =  (TTree*) fileOpen->Get("AnalysisxTree");

  //! histogram options
  double xminExp = 100;
  double xmaxExp = 1500;
  int binExp = xmaxExp-xminExp;
  
  TH1D* hCh0 = new TH1D("hCh0", "Channel 0", binExp, xminExp, xmaxExp);
  TH1D* hCh1 = new TH1D("hCh1", "Channel 1", binExp, xminExp, xmaxExp);
  TH1D* hCh2 = new TH1D("hCh2", "Channel 2", binExp, xminExp, xmaxExp);

  TCanvas *canvas = new TCanvas("canvas", "Exp data show", 1500, 500);
  canvas->Divide(3, 1);
  
  canvas->cd(1);
  treeOpen->Draw("NeEvent.neutAmp[0]>>hCh0");

  canvas->cd(2);
  treeOpen->Draw("NeEvent.neutAmp[1]>>hCh1");

  canvas->cd(3);
  treeOpen->Draw("NeEvent.neutAmp[2]>>hCh2");

  //! Fill experiment histograms
  TFile* fileSave = new TFile("./expfiles/usable/save_stilbene_divider_Co60_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "recreate");

  TH1D* save_hCh0 = new TH1D();
  save_hCh0 = (TH1D*) hCh0->Clone();
  TH1D* save_hCh1 = new TH1D();
  save_hCh1 = (TH1D*) hCh1->Clone();
  TH1D* save_hCh2 = new TH1D();
  save_hCh2 = (TH1D*) hCh2->Clone();

  fileSave->Write();
  fileSave->Close();

  std::cout << "time: " << timer->RealTime() << " seconds \n";
}
