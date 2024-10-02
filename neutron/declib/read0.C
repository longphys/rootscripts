#include "DecLibTest.hh"
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TStopwatch.h>
#include <sstream>
#include <fstream>
#include <iostream>

int64_t EventNr;
vector <double> dE;
vector <double> dTime;
vector <unsigned int> channel;

void read0()
{
  std::vector <double> *de = nullptr;

  TFile* file = new TFile("./output/test12_2_Plast_BC404_Cs137.root","read");
  TTree* tree = (TTree*)file->Get("ETree");

  tree->SetBranchAddress("EDep", &de);

  TH1D* h0 = new TH1D("h0", "h0", 1000, 0., 1000.);
  TH1D* h1 = new TH1D("h1", "h1", 1000, 0., 1000.);
  TH1D* h2 = new TH1D("h2", "h2", 1000, 0., 1000.);

  long int entries = tree->GetEntriesFast();
  for(int i = 0; i < entries ;i++){
    tree->GetEntry(i);

    int de_size = (int)de->size(); // MUST be done for every tree entry
    double *de_data = de->data();

    if(de_data[0]>0.){
      h0->Fill(de_data[0]);
    }

    if(de_data[1]>0.){
      h1->Fill(de_data[1]);
    }
  }

  TCanvas* c = new TCanvas("c", "c", 1000, 800);
  c->Divide(2,1);

  c->cd(1);
  h0->Draw();

  c->cd(2);
  h1->Draw();

}
