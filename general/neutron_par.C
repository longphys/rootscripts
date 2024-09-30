#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"

void neutron_par()
{
    TFile* file = new TFile("./simfiles/neutron/out.root", "read");

    TTree* tree =  (TTree*) file->Get("dEEtree");

    // histogram options
    int bins = 100;
    double xmin = 0.;
    double xmax = 20.0;
    TH1D* hist0 = new TH1D("hist0", "Total EDep", bins, xmin, xmax);
    TH1D* hist1 = new TH1D("hist1", "Neutron EDep", bins, xmin, xmax);
    TH1D* hist2 = new TH1D("hist2", "Proton EDep", bins, xmin, xmax);
    TH1D* hist3 = new TH1D("hist3", "Gamma EDep", bins, xmin, xmax);
    TH1D* hist4 = new TH1D("hist4", "Alpha EDep", bins, xmin, xmax);
    TH1D* hist5 = new TH1D("hist5", "C12 EDep", bins, xmin, xmax);
    TH1D* hist6 = new TH1D("hist6", "Other EDep", bins, xmin, xmax);
    TH1D* hist7 = new TH1D("hist7", "Total Compare", bins, xmin, xmax);

    TCanvas* c1 = new TCanvas("c1", "c1", 1600, 1000);
    c1->cd();
    c1->Divide(3, 2);

    c1->cd(1);
    tree->Draw("NeutronEDep>>hist1", "NeutronEDep>8.&&NeutronEDep<10.");
    c1->cd(2);
    tree->Draw("ProtonEDep>>hist2", "ProtonEDep>8.&&ProtonEDep<10.");
    c1->cd(3);
    tree->Draw("GammaEDep>>hist3", "GammaEDep>8.&&GammaEDep<10.");
    c1->cd(4);
    tree->Draw("AlphaEDep>>hist4", "AlphaEDep>8.&&AlphaEDep<10.");
    c1->cd(5);
    tree->Draw("C12EDep>>hist5", "C12EDep>8.&&C12EDep<10.");
    c1->cd(6);
    tree->Draw("OtherEDep>>hist6", "OtherEDep>8.&&OtherEDep<10.");

    TCanvas* c2 = new TCanvas("c2", "c2", 1500, 500);
    c2->Divide(3,1);
    c2->cd(1);
    hist0->SetLineColor(kBlack);
    hist0->SetLineWidth(2);
    tree->Draw("Scintillator>>hist0", "Scintillator>0");

    hist1->SetLineColor(kGreen);
    hist1->SetLineWidth(2);
    tree->Draw("NeutronEDep>>hist1", "NeutronEDep>0", "same");
    
    hist2->SetLineColor(kBlue);
    hist2->SetLineWidth(2);
    tree->Draw("ProtonEDep>>hist2", "ProtonEDep>0", "same");

    hist3->SetLineColor(kRed);
    hist3->SetLineWidth(2);
    tree->Draw("GammaEDep>>hist3", "GammaEDep>0", "same");
    
    hist4->SetLineColor(kMagenta);
    hist4->SetLineWidth(2);
    tree->Draw("AlphaEDep>>hist4", "AlphaEDep>0", "same");

    hist5->SetLineColor(kYellow);
    hist5->SetLineWidth(2);
    tree->Draw("C12EDep>>hist5", "C12EDep>0", "same");

    hist6->SetLineColor(kOrange);
    hist6->SetLineWidth(2);
    tree->Draw("OtherEDep>>hist6", "OtherEDep>0", "same");

    TLegend *leg1 = new TLegend(0.75, 0.6, 0.98, 0.75);
    leg1->SetHeader("Particles", "C");
    leg1->SetBorderSize(2);
    leg1->AddEntry(hist0, "Total EDep", "l");
    leg1->AddEntry(hist1, "Neutron EDep", "l");
    leg1->AddEntry(hist2, "Proton EDep", "l");
    leg1->AddEntry(hist3, "Gamma EDep", "l");
    leg1->AddEntry(hist4, "Alpha EDep", "l");
    leg1->AddEntry(hist5, "C12 EDep", "l");
    leg1->AddEntry(hist6, "Other EDep", "l");
    leg1->Draw();
    
    c2->cd(2);
    hist7->SetLineColor(kRed);
    hist7->SetLineWidth(2);
    tree->Draw("NeutronEDep+ProtonEDep+GammaEDep+AlphaEDep+C12EDep+OtherEDep>>hist7", "NeutronEDep+ProtonEDep+GammaEDep+AlphaEDep+C12EDep+OtherEDep>0");
    c2->cd(3);
    tree->Draw("Scintillator>>hist0", "Scintillator>0");

    TLegend *leg2 = new TLegend(0.75, 0.6, 0.98, 0.75);
    leg2->SetHeader("Particles", "C");
    leg2->SetBorderSize(2);
    leg2->AddEntry(hist0, "Total EDep", "l");
    leg2->AddEntry(hist7, "Compare total EDep", "l");

    // file->Close();
}