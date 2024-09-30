#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"


void plotTree()
{
    TCanvas* c1 = new TCanvas();

    TFile* input = new TFile("out.root", "read");

    TTree* tree = (TTree*)input->Get("tree");

    double x, y;

    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);

    int entries = tree->GetEntries();

    std::cout << entries << std::endl;

    TH1F *hist = new TH1F("hist", "Histogram", 20, 0, 20);

    for (int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);

        std::cout << x << "" << y << std::endl;

        hist->Fill(x);
    }

    hist->Draw();

    input->Close();
}