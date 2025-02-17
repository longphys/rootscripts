#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TPaveText.h"
#include "TGraph.h"

//! pixels, 1 pixel - diameter = 3,5cm
// double distance_from_center[] = {0, 3.5, 7.0, 10.5, 14.0, 17.5, 21.0};
double size_of_pixel = 2.; //cm
//! Proton cross-talk
// 8x8: 2,6%; 8x9: 5,2%; 8x10: 12,8%
double proton_cross_talk[] = {2.6, 4.5, 5.0, 2.8, 1.9, 1.3, 1.0};
//! Neutron mis-hits neighbor
// 8x8: 3,2%; 8x9: 12,5%; 8x10: 148,0%
double neutron_mis_hits_neighbor[] = {3.1, 10.4, 57.1, 89.1, 93.3, 95.1, 94.4};
//! Neutron mis-hits central
// 8x8: 37,0%; 8x9: 18,1%; 8x10: 6,8%
double neutron_mis_hits_center[] = {91.8, 80.6, 32.9, 5.3, 2.8, 2.4, 3.6};
void read2(){
  TCanvas* c_1 = new TCanvas("c_1", "c_1", 1000, 700);
  TGraph *gr1 = new TGraph();
  TGraph *gr2 = new TGraph();
  TGraph *gr3 = new TGraph();

  gr1->SetLineColor(kRed);
  gr1->SetLineWidth(3);
  gr1->SetMarkerStyle(20); // Full circle marker
  gr1->SetMarkerSize(1.2);
  gr1->SetMarkerColor(kRed);

  gr2->SetLineColor(kBlue);
  gr2->SetLineWidth(3);
  gr2->SetMarkerStyle(21); // Square marker
  gr2->SetMarkerSize(1.2);
  gr2->SetMarkerColor(kBlue);

  gr3->SetLineColor(kBlack);
  gr3->SetLineWidth(3);
  gr3->SetMarkerStyle(22); // Triangle marker
  gr3->SetMarkerSize(1.2);
  gr3->SetMarkerColor(kBlack);

  gr1->GetXaxis()->SetTitle("Distance from center (cm)");
  gr1->GetXaxis()->SetLabelFont(42);
  gr1->GetXaxis()->SetTitleFont(52);
  gr1->GetXaxis()->SetTitleSize(0.04);
  gr1->GetXaxis()->CenterTitle(true);

	gr1->GetYaxis()->SetTitle("Percentage of detected events (%)");
  gr1->GetYaxis()->SetLabelFont(42);
  gr1->GetYaxis()->SetTitleFont(52);
  gr1->GetYaxis()->SetTitleSize(0.04);
  gr1->GetYaxis()->CenterTitle(true);

  for (int i = 0; i < sizeof(proton_cross_talk)/sizeof(proton_cross_talk[0]); i++){
    std::cout << i << "\n";
    gr1->AddPoint(2*i, proton_cross_talk[i]);
    std::cout << "distance = " << 2*i << "; cross-talk = " << proton_cross_talk[i] << "\n";
    gr2->AddPoint(2*i, neutron_mis_hits_neighbor[i]);
    std::cout << "distance = " << 2*i << "; mis_hits_neighbor = " << neutron_mis_hits_neighbor[i] << "\n";
    gr3->AddPoint(2*i, neutron_mis_hits_center[i]);
    std::cout << "distance = " << 2*i << "; mis_hits_center = " << neutron_mis_hits_center[i] << "\n";
  }

  c_1->cd();
  gr1->Draw("ALP");
  gr1->GetXaxis()->SetRangeUser(-1.0, 13.0);
  gr1->GetYaxis()->SetRangeUser(0., 100.);
  gr2->Draw("LP SAME");
  gr3->Draw("LP SAME");

  TLegend *legend = new TLegend(0.4, 0.55, 0.8, 0.85);
  legend->SetBorderSize(0);
  legend->SetLineWidth(2);
  legend->SetTextSize(0.04);
  legend->AddEntry(gr1, "Proton cross-talk", "l");
  legend->AddEntry(gr2, "#splitline{Neutron mis-hits to}{the neighboring module}", "l");
  legend->AddEntry(gr3, "#splitline{Neutron mis-hits to}{the center module}", "l");

  legend->Draw();

  TCanvas* c_2 = new TCanvas("c_2", "c_2", 1000, 700);
  c_2->cd();

  gr1->Draw("ALP");
  // gr1->GetXaxis()->SetRangeUser(-1.0, 10.0);
  // gr1->GetYaxis()->SetRangeUser(0., 200.);
  gr2->Draw("LP SAME");
  // gr3->Draw("LP SAME");

  TLegend *legend1 = new TLegend(0.4, 0.55, 0.8, 0.85);
  legend1->SetBorderSize(0);
  legend1->SetLineWidth(2);
  legend1->SetTextSize(0.04);
  // legend1->AddEntry(gr1, "Proton cross-talk", "l");
  // legend1->AddEntry(gr2, "#splitline{Neutron mis-hits to}{the neighboring module}", "l");
  // legend->AddEntry(gr3, "Neutron mis-hits \n to the center module \n", "l");

  // legend1->Draw();
}