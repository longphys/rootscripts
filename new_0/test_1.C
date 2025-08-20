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
#include "TChain.h"
#include "TBranchElement.h"
#include "TDirectory.h"

void test_1()
{
  auto timer = new TStopwatch();
  timer->Start();

  //! Files and trees
  TFile* file = new TFile("~/data/25e04/run25_00.root", "read");
  if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file or file is corrupted.\n";
        return;
  }
        
  TTree* tree =  (TTree*) file->Get("AnalysisxTree");
  if (!tree) {
        std::cerr << "Error: TTree 'AnalysisxTree' not found in the file.\n";
        file->Close();
        return;
  }
    
  /*TChain* chain = new TChain("chain");
  chain->Add("~/25e04/run12_0001.root?#AnalysisxTree");
  chain->Add("~/25e04/run12_0001.root?#AnalysisxTree");
  chain->Add("~/25e04/run12_0001.root?#AnalysisxTree");
  */
  
  TBranch* element = (TBranch*)tree->GetBranch("NeEvent.Lea[16]");
  element->Print();
  
  
  constexpr int n_de_right = 16;
  constexpr int n_de_left = 32;
 
  constexpr int n_e_right = 16;
  constexpr int n_e_left = 32;
  
  double thresh_map_right_min = 200.;
  double thresh_map_right_max = 7600.;
  double thresh_map_left_min = 300.;
  double thresh_map_left_max = 7500.;
  
  // Values: Silicon strips
  UShort_t de_amp_right_x[n_de_right];
  UShort_t de_amp_right_y[n_de_right];
  UShort_t de_time_right_x[n_de_right];
  UShort_t de_time_right_y[n_de_right];
  
  UShort_t de_amp_left_x[n_de_left];
  UShort_t de_amp_left_y[n_de_left];
  UShort_t de_time_left_x[n_de_left];
  UShort_t de_time_left_y[n_de_left];
  
  // Values: LYSO (Left) and CsI(Tl) (Right)
  UShort_t e_amp_right[n_e_right];
  UShort_t e_time_right[n_e_right];
  
  UShort_t e_amp_left[n_e_left];
  UShort_t e_time_left[n_e_left];
  
  // Trees: Silicon strips
  tree->SetBranchAddress("NeEvent.Rxa[16]", de_amp_right_x);
  tree->SetBranchAddress("NeEvent.Rya[16]", de_amp_right_y);
  tree->SetBranchAddress("NeEvent.Rxt[16]", de_time_right_x);
  tree->SetBranchAddress("NeEvent.Ryt[16]", de_time_right_y);

  tree->SetBranchAddress("NeEvent.Lxa[32]", de_amp_left_x);
  tree->SetBranchAddress("NeEvent.Lya[32]", de_amp_left_y);
  tree->SetBranchAddress("NeEvent.Lxt[32]", de_time_left_x);
  tree->SetBranchAddress("NeEvent.Lyt[32]", de_time_left_y);
  
  // Values: LYSO (Left) and CsI(Tl) (Right)
  tree->SetBranchAddress("NeEvent.Rea[16]", e_amp_right);
  tree->SetBranchAddress("NeEvent.Ret[16]", e_time_right);

  tree->SetBranchAddress("NeEvent.Lea[16]", e_amp_left);
  tree->SetBranchAddress("NeEvent.Let[16]", e_time_left);
  
  /*Double_t de_amp_right_x[n_de_right];
  Double_t de_amp_right_y[n_de_right];
  tree->SetBranchAddress("NeEvent.Rxc[16]", de_amp_right_x);
  tree->SetBranchAddress("NeEvent.Ryc[16]", de_amp_right_y);*/
  
  // Entries
  Long64_t entries_all = tree->GetEntries();
  std::cout << "\nNumber of all entries: " << entries_all << "\n";
  // Long64_t entries_used = entries_all;
  Long64_t entries_used = 50000;
  std::cout << "\nNumber of entries used: " << entries_used << "\n";
  
  // Histograms: Silicon strips
  TH1D* hist_de_amp_right_x[n_de_right];
  TH1D* hist_de_amp_right_y[n_de_right];
  TH1D* hist_de_time_right_x[n_de_right];
  TH1D* hist_de_time_right_y[n_de_right];
  
  TH1D* hist_de_amp_left_x[n_de_left];
  TH1D* hist_de_amp_left_y[n_de_left];
  TH1D* hist_de_time_left_x[n_de_left];
  TH1D* hist_de_time_left_y[n_de_left];
  
  char* name = new char[20];
    
  for (int i = 0; i < n_de_right; i++)
  {
    sprintf(name,"hist_de_amp_right_x_%d",i);
    hist_de_amp_right_x[i] = new TH1D(name, name , 1000, 0, 10000);
    sprintf(name,"hist_de_amp_right_y_%d",i);
    hist_de_amp_right_y[i] = new TH1D(name, name , 1000, 0, 10000);
    sprintf(name,"hist_de_time_right_x_%d",i);
    hist_de_time_right_x[i] = new TH1D(name, name , 1000, 0, 10000);
    sprintf(name,"hist_de_time_right_y_%d",i);
    hist_de_time_right_y[i] = new TH1D(name, name , 1000, 0, 10000);
  }
  
  for (int i = 0; i < n_de_left; i++)
  {
    sprintf(name,"hist_de_amp_left_x_%d",i);
    hist_de_amp_left_x[i] = new TH1D(name, name , 1000, 0, 10000);
    sprintf(name,"hist_de_amp_left_y_%d",i);
    hist_de_amp_left_y[i] = new TH1D(name, name , 1000, 0, 10000);
    sprintf(name,"hist_de_time_left_x_%d",i);
    hist_de_time_left_x[i] = new TH1D(name, name , 1000, 0, 10000);
    sprintf(name,"hist_de_time_left_y_%d",i);
    hist_de_time_left_y[i] = new TH1D(name, name , 1000, 0, 10000);
  }
  
  // Histograms: LYSO (Left) and CsI(Tl) (Right)
  TH1D* hist_e_amp_right[n_e_right];
  TH1D* hist_e_time_right[n_e_right];
  
  TH1D* hist_e_amp_left[n_e_left];
  TH1D* hist_e_time_left[n_e_left];
  
  for (int i = 0; i < n_e_right; i++)
  {
    sprintf(name,"hist_e_amp_right_%d",i);
    hist_e_amp_right[i] = new TH1D(name, name , 1000, 0, 10000);
    sprintf(name,"hist_e_time_right_%d",i);
    hist_e_time_right[i] = new TH1D(name, name , 1000, 0, 10000);
  }
  
  for (int i = 0; i < n_e_left; i++)
  {
    sprintf(name,"hist_e_amp_left_%d",i);
    hist_e_amp_left[i] = new TH1D(name, name , 1000, 0, 10000);
    sprintf(name,"hist_e_time_left_%d",i);
    hist_e_time_left[i] = new TH1D(name, name , 1000, 0, 10000);
  }
  
  // Histograms: 2D plots
  TH2D* hist_map_right = new TH2D("hist_map_right", "hist_map_right", 16, 0, 16, 16, 0, 16);
  TH2D* hist_map_left = new TH2D("hist_map_left", "hist_map_left", 32, 0, 32, 32, 0, 32);
  
  TH2D* hist_de_e_right[n_e_right];
  TH2D* hist_de_e_left[n_e_left];
  
  for (int i = 0; i < n_e_right; i++)
  {
    sprintf(name,"de_e_right_%d",i);
    hist_de_e_right[i] = new TH2D(name, name, 1000, 0., 10000., 3000, 5000., 8000.);
    hist_de_e_right[i]->GetXaxis()->SetTitle("Amp");
    hist_de_e_right[i]->GetYaxis()->SetTitle("Time");
  }
  
  for (int i = 0; i < n_e_left; i++)
  {
    sprintf(name,"de_e_left_%d",i);
    hist_de_e_left[i] = new TH2D(name, name, 1000, 0., 10000., 1000, 0., 10000.);
    hist_de_e_left[i]->GetXaxis()->SetTitle("Amp");
    hist_de_e_left[i]->GetYaxis()->SetTitle("Time");
  }
  
  for (int i_entry = 0; i_entry < entries_used; i_entry++)
  {
    tree->GetEntry(i_entry);
    if(i_entry%5000 == 0)
    {
      std::cout << "entry number: " << i_entry << "\n";
    }

    // Fill: Silicon strips
    for(int strip_x = 0; strip_x < n_de_right; strip_x++)
    {
      hist_de_amp_right_x[strip_x]->Fill(de_amp_right_x[strip_x]);
      hist_de_time_right_x[strip_x]->Fill(de_time_right_x[strip_x]);
      for(int strip_y = 0; strip_y < n_de_right; strip_y++)
      {
        hist_de_time_right_y[strip_y]->Fill(de_time_right_y[strip_y]);
        hist_de_amp_right_y[strip_y]->Fill(de_amp_right_y[strip_y]);
        if(de_amp_right_y[strip_y] > thresh_map_right_min && 
        de_amp_right_y[strip_y] < thresh_map_right_max && 
        de_amp_right_x[strip_x] > thresh_map_right_min && 
        de_amp_right_x[strip_x] < thresh_map_right_max)
        {
          hist_map_right->Fill(strip_x,strip_y);
        }
      }
    }

    for(int strip_x = 0; strip_x < n_de_left; strip_x++)
    {
      hist_de_time_left_x[strip_x]->Fill(de_time_left_x[strip_x]);
      hist_de_amp_left_x[strip_x]->Fill(de_amp_left_x[strip_x]);
      for(int strip_y = 0; strip_y < n_de_left; strip_y++)
      {
        hist_de_time_left_y[strip_y]->Fill(de_time_left_y[strip_y]);
        hist_de_amp_left_y[strip_y]->Fill(de_amp_left_y[strip_y]);
        if(de_amp_left_y[strip_y] > thresh_map_left_min && 
        de_amp_left_y[strip_y] < thresh_map_left_max && 
        de_amp_left_x[strip_x] > thresh_map_left_min && 
        de_amp_left_x[strip_x] < thresh_map_left_max)
        {
          hist_map_left->Fill(strip_x,strip_y);
        }
      }
    }

    // Fill: LYSO (Left) and CsI(Tl) (Right)
    for(int strip_right = 0; strip_right < n_e_right; strip_right++)
    {
      hist_e_amp_right[strip_right]->Fill(e_amp_right[strip_right]);
      hist_e_time_right[strip_right]->Fill(e_time_right[strip_right]);
      hist_de_e_right[strip_right]->Fill(e_amp_right[strip_right], e_time_right[strip_right]);
    }

    for(int strip_left = 0; strip_left < n_e_left; strip_left++)
    {
      hist_e_amp_left[strip_left]->Fill(e_amp_left[strip_left]);
      hist_e_time_left[strip_left]->Fill(e_time_left[strip_left]);
      hist_de_e_left[strip_left]->Fill(e_amp_left[strip_left], e_time_left[strip_left]);
    }
  }
  
  //! Canvas and draw
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  c1->cd();
  //hist_map_right->Draw("colz");
  hist_de_e_right[4]->Draw("colz");
  
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
  c2->cd();
  //hist_map_left->Draw("colz");
  hist_de_e_left[4]->Draw("colz");
  
  TCanvas *c3 = new TCanvas("c3", "c3", 800, 800);
  c3->cd();
  hist_de_amp_right_x[0]->Draw(); 
  
  //! Refill histograms
  TH1D* clone_hist_de_amp_right_x[n_de_right];
  TH1D* clone_hist_de_amp_right_y[n_de_right];
  TH1D* clone_hist_de_time_right_x[n_de_right];
  TH1D* clone_hist_de_time_right_y[n_de_right];
  
  TH1D* clone_hist_de_amp_left_x[n_de_left];
  TH1D* clone_hist_de_amp_left_y[n_de_left];
  TH1D* clone_hist_de_time_left_x[n_de_left];
  TH1D* clone_hist_de_time_left_y[n_de_left];
  
  TH1D* clone_hist_e_amp_right[n_e_right];
  TH1D* clone_hist_e_time_right[n_e_right];
  
  TH1D* clone_hist_e_amp_left[n_e_left];
  TH1D* clone_hist_e_time_left[n_e_left];
  
  for (int i = 0; i < n_de_right; i++)
  {
    clone_hist_de_amp_right_x[i] = (TH1D*)hist_de_amp_right_x[i]->Clone();
    clone_hist_de_amp_right_y[i] = (TH1D*)hist_de_amp_right_y[i]->Clone();
    clone_hist_de_time_right_x[i] = (TH1D*)hist_de_time_right_x[i]->Clone();
    clone_hist_de_time_right_y[i] = (TH1D*)hist_de_time_right_y[i]->Clone();
  }
  
  for (int i = 0; i < n_de_left; i++)
  {
    clone_hist_de_amp_left_x[i] = (TH1D*)hist_de_amp_left_x[i]->Clone();
    clone_hist_de_amp_left_y[i] = (TH1D*)hist_de_amp_left_y[i]->Clone();
    clone_hist_de_time_left_x[i] = (TH1D*)hist_de_time_left_x[i]->Clone();
    clone_hist_de_time_left_y[i] = (TH1D*)hist_de_time_left_y[i]->Clone();
  }
  
  for (int i = 0; i < n_e_right; i++)
  {
    clone_hist_e_amp_right[i] = (TH1D*)hist_e_amp_right[i]->Clone();
    clone_hist_e_time_right[i] = (TH1D*)hist_e_time_right[i]->Clone();
  }
  
  for (int i = 0; i < n_e_left; i++)
  {
    clone_hist_e_amp_left[i] = (TH1D*)hist_e_amp_left[i]->Clone();
    clone_hist_e_time_left[i] = (TH1D*)hist_e_time_left[i]->Clone();
  }
  
  TH2D* clone_hist_map_right = new TH2D();
  clone_hist_map_right = (TH2D*)hist_map_right->Clone();
  
  TH2D* clone_hist_map_left = new TH2D();
  clone_hist_map_left = (TH2D*)hist_map_left->Clone();
  
  TH2D* clone_hist_de_e_right[n_e_right];
  TH2D* clone_hist_de_e_left[n_e_left];
  
  for (int i = 0; i < n_e_right; i++)
  {
    clone_hist_de_e_right[i] = (TH2D*)hist_de_e_right[i]->Clone();
  }
  
  for (int i = 0; i < n_e_left; i++)
  {
    clone_hist_de_e_left[i] = (TH2D*)hist_de_e_left[i]->Clone();
  }
  
  TFile* file_results = new TFile("test_1_results.root", "recreate");
  
  TDirectory* dir_de_amp_right_x = file_results->mkdir("dE Right Amp X");
  TDirectory* dir_de_amp_right_y = file_results->mkdir("dE Right Amp Y");
  TDirectory* dir_de_time_right_x = file_results->mkdir("dE Right Time X");
  TDirectory* dir_de_time_right_y = file_results->mkdir("dE Right Time Y");
  
  TDirectory* dir_de_amp_left_x = file_results->mkdir("dE Left Amp X");
  TDirectory* dir_de_amp_left_y = file_results->mkdir("dE Left Amp Y");
  TDirectory* dir_de_time_left_x = file_results->mkdir("dE Left Time X");
  TDirectory* dir_de_time_left_y = file_results->mkdir("dE Left Time Y");
  
  TDirectory* dir_e_amp_right = file_results->mkdir("E Right Amp");
  TDirectory* dir_e_time_right = file_results->mkdir("E Right Time");

  TDirectory* dir_e_amp_left = file_results->mkdir("E Left Amp");
  TDirectory* dir_e_time_left = file_results->mkdir("E Left Time");

  for (int i = 0; i < n_de_right; i++)
  {
    dir_de_amp_right_x->cd();
    clone_hist_de_amp_right_x[i]->Write();
    dir_de_amp_right_y->cd();
    clone_hist_de_amp_right_y[i]->Write();
    dir_de_time_right_x->cd();
    clone_hist_de_time_right_x[i]->Write();
    dir_de_time_right_y->cd();
    clone_hist_de_time_right_y[i]->Write();    
  }
  
  for (int i = 0; i < n_de_left; i++)
  {
    dir_de_amp_left_x->cd();
    clone_hist_de_amp_left_x[i]->Write();
    dir_de_amp_left_y->cd();
    clone_hist_de_amp_left_y[i]->Write();
    dir_de_time_left_x->cd();
    clone_hist_de_time_left_x[i]->Write();
    dir_de_time_left_y->cd();
    clone_hist_de_time_left_y[i]->Write();
  }
  
  for (int i = 0; i < n_e_right; i++)
  {
    dir_e_amp_right->cd();
    clone_hist_e_amp_right[i]->Write();
    dir_e_time_right->cd();
    clone_hist_e_time_right[i]->Write();
  }
  
  for (int i = 0; i < n_e_left; i++)
  {
    dir_e_amp_left->cd();
    clone_hist_e_amp_left[i]->Write();
    dir_e_time_left->cd();
    clone_hist_e_time_left[i]->Write();
  }

  file_results->cd();
  clone_hist_map_right->Write();
  clone_hist_map_left->Write();
  
  for (int i = 0; i < n_e_right; i++)
  {
    clone_hist_de_e_right[i]->Write();
  }
  
  for (int i = 0; i < n_e_left; i++)
  {
    clone_hist_de_e_left[i]->Write();
  }

  file_results->Write();
  file_results->Close();
  
  std::cout << "\ntime: " << timer->RealTime() << " seconds\n";
}
