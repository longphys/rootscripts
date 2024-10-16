#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"

void draw()
{
  int first = 0.00001;
  int last;
  TGraph* gr = new TGraph();


  for (int i=0; i<last; i++){
    double energy = 
  }

  //! Canvas
  TCanvas *c = new TCanvas("c", "c", 700, 500);
  c->cd();

  gr->SetLineColor(kRed);
  gr->Draw();
}