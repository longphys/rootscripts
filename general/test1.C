#include "TRandom3.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"

#include <iostream>

void test()
{
  double x_sim;
  tsimCs137->SetBranchAddress("Scintillator", &x_sim);
  double y_sim;
  tsimNa22->SetBranchAddress("Scintillator", &y_sim);

  for(int i = 0; i < timeSig1; i++){
    TH1D* hsimCs137res = new TH1D("hsimCs137res", "Cs137 Simulation with resolution", binExp, xminCal, xmaxCal);
    TH1D* hsimNa22res = new TH1D("hsimNa22res", "Na22 Simulation with resolution", binExp, xminCal, xmaxCal);
    std::cout << "test";
    entries = tsimCs137->GetEntries();

    for (int i = 0; i < entries; i++)
    {
      tsimCs137->GetEntry(i);
      double sigma = x_sim*sqrt( pow(coA,2) + pow(coB[i]/sqrt(x_sim),2) + pow(coC[i]/x_sim,2) );
      hsimCs137res->Fill(ranGen->Gaus(x_sim,sigma));
    }
  }
}