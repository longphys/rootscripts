#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"

void backScatter()
{
    TFile* file = new TFile("out.root", "read");

    TTree* tree =  (TTree*) file->Get("dEEtree");

    // histogram options
    int bins = 400;
    double xmin = 0.;
    double xmax = 0.7;
    //TH1D* hist1 = new TH1D("hist1", "Compton", bins, xmin, xmax);
    TH1D* hist1 = new TH1D("hist1", "Compton", bins, xmin, xmax);
    TH1D* hist2 = new TH1D("hist2", "BackScatteredComptons", bins, xmin, xmax);

    TCanvas* c1 = new TCanvas("c1", "c1", 1600, 400);
    int cols = 2;
    int rows = 1;
    c1->cd();
    c1->Divide(cols, rows);

    c1->cd(1);
    hist1->SetLineColor(kBlack);
    tree->Draw("ComptonEnergyDeposition>>hist1","ComptonEnergyDeposition>0");
       
    hist2->SetLineColor(kPink);
    tree->Draw("BackScatteredComptonEnergyDeposition>>hist2", "BackScatteredComptonEnergyDeposition>0", "same");

    c1->cd(2);
    tree->Draw("BackScatteredComptonEnergyDeposition>>hist2", "BackScatteredComptonEnergyDeposition>0");
    
    // file->Close();
}