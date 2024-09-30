#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"

void multiComp()
{
    TFile* file = new TFile("./simfiles/withoutAl_simCs137_1.root", "read");

    TTree* tree =  (TTree*) file->Get("dEEtree");

    // histogram options
    int bins = 1000;
    double xmin = 0.;
    double xmax = 1.0;
    TH1D* hist0 = new TH1D("hist0", "EDep", bins, xmin, xmax);
    TH1D* hist1 = new TH1D("hist1", "Compton", bins, xmin, xmax);
    TH1D* hist2 = new TH1D("hist2", "1 Compton", bins, xmin, xmax);
    TH1D* hist3 = new TH1D("hist3", "2 Comptons", bins, xmin, xmax);
    TH1D* hist4 = new TH1D("hist4", "3 Comptons", bins, xmin, xmax);

    TCanvas* c1 = new TCanvas("c1", "c1", 1600, 400);
    int cols = 4;
    int rows = 1;
    c1->cd();
    c1->Divide(cols, rows);

    c1->cd(1);
    hist0->SetLineColor(kMagenta);
    tree->Draw("Scintillator>>hist0", "ComptonEnergyDeposition>0");

    hist1->SetLineColor(kBlack);
    tree->Draw("ComptonEnergyDeposition>>hist1","ComptonEnergyDeposition>0", "same");
    
    hist2->SetLineColor(kBlue);
    // hist1->GetYaxis()->SetRangeUser(0, 3000);   
    
    tree->Draw("ComptonEnergyDeposition[0]>>hist2","ComptonEnergyDeposition[0]>0", "same");
    hist3->SetLineColor(kRed);
    // hist2->GetYaxis()->SetRangeUser(0, 3000);
    tree->Draw("ComptonEnergyDeposition[1]>>hist3", "ComptonEnergyDeposition[1]>0", "same");
    hist4->SetLineColor(kGreen);
    // hist3->GetYaxis()->SetRangeUser(0, 100);
    tree->Draw("ComptonEnergyDeposition[2]>>hist4", "ComptonEnergyDeposition[2]>0", "same");

    TLegend *leg1 = new TLegend(0.75, 0.6, 0.98, 0.75);
    leg1->SetHeader("Multiple Comptons", "C");
    leg1->SetBorderSize(2);
    leg1->AddEntry(hist0, "Total", "l");
    leg1->AddEntry(hist1, "Total", "l");
    leg1->AddEntry(hist2, "1 Compton", "l");
    leg1->AddEntry(hist3, "2 Comptons", "l");
    leg1->AddEntry(hist4, "3 Comptons", "l");
    leg1->Draw();

    c1->cd(2);
    tree->Draw("ComptonEnergyDeposition[0]>>hist2","ComptonEnergyDeposition[0]>0");
    c1->cd(3);
    tree->Draw("ComptonEnergyDeposition[1]>>hist3", "ComptonEnergyDeposition[1]>0");
    c1->cd(4);
    tree->Draw("ComptonEnergyDeposition[2]>>hist4", "ComptonEnergyDeposition[2]>0");
    // file->Close();
}