#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLegend.h"
#include "TRandom3.h"

void eRes()
{
  TFile* fCs137sim = new TFile("new_withAl_simCs137.root", "read");
  TFile* fCs137exp = new TFile("new_passive_divider_ch2_Cs137_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "read");
  
  TTree* tCs137sim =  (TTree*) fCs137sim->Get("dEEtree");
  TTree* tCs137exp =  (TTree*) fCs137exp->Get("AnalysisxTree");
  
  // histogram options
  int binSim = 1300;
  int binExp = 700;
  double xminSim = 0.;
  double xmaxSim = 1.3;
  double xminExp = 100;
  double xmaxExp = 800;

  double a = 0.002569802420751

;
  double b = -0.216819075991958

;
  // a = ;
  // b = -0.258;
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1400, 700);
  c1->cd();

// hist Cs137 sim
  TH1D* hCs137sim = new TH1D("hCs137sim", "Cs137 sim", binExp, xminExp*a+b, xmaxExp*a+b);

  double entries = tCs137sim->GetEntries();
  double x;
  tCs137sim->SetBranchAddress("Scintillator", &x);

  TRandom3* gausGen = new TRandom3(0);

  double coA = 0.003;
  double coB = 0.003;
  double coC = 0.007;

  int count = 0;

  for(int i = 0; i < entries; i++){
    tCs137sim->GetEntry(i);

    double sigma = sqrt(pow(coA,2) + pow((coB/sqrt(x)),2) + pow(coC/x,2));
    if(x <= 0)
    {
      continue;
    }
    double edep = gausGen->Gaus(x, sigma);
    count ++;
    if(edep <= 0.00001)
    {
      continue;
    }
    // std::cout << edep << '\n';
    hCs137sim->Fill(edep);
  }
  hCs137sim->Draw();

// hist Cs137 exp
  TH1D* hCs137exp = new TH1D("hCs137exp", "Cs137 exp", binExp, xminExp*a+b, xmaxExp*a+b);

  UShort_t y[48];
  tCs137exp->SetBranchAddress("NeEvent.neutAmp[48]", y);

  for(int i = 0; i < entries; i++){
    tCs137exp->GetEntry(i);
    hCs137exp->Fill(a*(y[0]+0.5) + b);
  }
  hCs137exp->SetLineColor(kRed);
  hCs137exp->Scale(0.072, "nosw2");
  hCs137exp->Draw("same");

  std::cout << "Cs137 sim entries: " << count << "\n";
  std::cout << "Cs137 exp entries: " << tCs137exp->GetEntries() << "\n";  
}
