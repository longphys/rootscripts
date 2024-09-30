#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TVirtualFFT.h"
#include "TStopwatch.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMarker.h"

void res(){
  double a,b,c;
  double e[2], sigma[2];
  e[0] = 383.651;
  sigma[0] = 24.9106;
  e[1] = 1104.65;
  sigma[1] = 111.069;
  TF1* res = new TF1("fit", "x*sqrt(pow([0],2)+pow([0],2)/pow(x,2))", 0., 1500.);
  TCanvas* canvas1 = new TCanvas("canvas1", "fit", 800, 600);
  TGraph* gr = new TGraph(2, e, sigma);
  gr->Fit(res);

  canvas1->cd();
  gr->Draw();
}