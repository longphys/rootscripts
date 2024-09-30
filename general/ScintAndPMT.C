#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"

void ScintAndPMT()
{
    TFile* file = new TFile("out.root", "read");

    TTree* tree =  (TTree*) file->Get("dEEtree");

    // histogram options
    int bins = 400;
    double xmin = 0.;
    double xmax = 0.7;
    //TH1D* hist1 = new TH1D("hist1", "Compton", bins, xmin, xmax);
    TH1D* hist1 = new TH1D("hist1", "Scint", bins, xmin, xmax);
    TH1D* hist2 = new TH1D("hist2", "ScintL", bins, xmin, xmax);
    TH1D* hist3 = new TH1D("hist3", "ScintR", bins, xmin, xmax);
    TH1D* hist4 = new TH1D("hist4", "PMT", bins, xmin, xmax);
    TH1D* hist5 = new TH1D("hist5", "PMT_L", bins, xmin, xmax);
    TH1D* hist6 = new TH1D("hist6", "PMT_R", bins, xmin, xmax);
    TH1D* hist7 = new TH1D("hist7", "BackScatteredComptonEnergyDeposition", bins, xmin, xmax);

    TCanvas* c1 = new TCanvas("c1", "c1", 1500, 800);
    int cols = 3;
    int rows = 2;
    c1->cd();
    c1->Divide(cols, rows);

    c1->cd(1);
    hist1->SetLineColor(kBlack);
    tree->Draw("Scintillator>>hist1","Scintillator>0");
    hist7->SetLineColor(kGreen);
    tree->Draw("BackScatteredComptonEnergyDeposition>>hist7", "BackScatteredComptonEnergyDeposition>0", "same");
       
    c1->cd(2);
    hist2->SetLineColor(kBlack);
    tree->Draw("ScintillatorL>>hist2","ScintillatorL>0");

    c1->cd(3);
    hist3->SetLineColor(kBlack);
    tree->Draw("ScintillatorR>>hist3","ScintillatorR>0");

    c1->cd(4);
    hist4->SetLineColor(kRed);
    tree->Draw("PMT>>hist4","PMT>0");

    c1->cd(5);
    hist5->SetLineColor(kRed);
    tree->Draw("PMTL>>hist5","PMTL>0");

    c1->cd(6);
    hist6->SetLineColor(kRed);
    tree->Draw("PMTR>>hist6","PMTR>0");
    
    // file->Close();
}