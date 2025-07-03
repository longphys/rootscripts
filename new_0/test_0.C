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

void test_0()
{
  auto timer = new TStopwatch();
  timer->Start();

  //! Files and trees
  TFile* file = new TFile("~/25e04/run12_0001.root", "read");
  if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file or file is corrupted.\n";
        return;
  }
        
  TTree* tree =  (TTree*) file->Get("AnalysisxTree");
  if (!tree) {
        std::cerr << "Error: TTree 'mytree' not found in the file.\n";
        file->Close();
        return;
  }
    
  /*TChain* chain = new TChain("chain");
  chain->Add("~/25e04/run12_0001.root?#AnalysisxTree");
  chain->Add("~/25e04/run12_0001.root?#AnalysisxTree");
  chain->Add("~/25e04/run12_0001.root?#AnalysisxTree");
  */
  
  TBranchElement* element = (TBranchElement*)tree->GetBranch("NeEvent.Ryt[16]");
	std::cout << "\nElement name: " << element->GetName() << "\n";
	std::cout << "Type name: " << element->GetTypeName() << "\n";
  
  int n_right = 16;
  int n_left = 32;
  
  UShort_t amp_right_x[n_right];
  UShort_t amp_right_y[n_right];
  UShort_t time_right_x[n_right];
  UShort_t time_right_y[n_right];
  
  UShort_t amp_left_x[n_left];
  UShort_t amp_left_y[n_left];
  UShort_t time_left_x[n_left];
  UShort_t time_left_y[n_left];
  
  tree->SetBranchAddress("NeEvent.Rxa[16]", amp_right_x);
  tree->SetBranchAddress("NeEvent.Rya[16]", amp_right_y);
  tree->SetBranchAddress("NeEvent.Rxt[16]", time_right_x);
  tree->SetBranchAddress("NeEvent.Ryt[16]", time_right_y);

  tree->SetBranchAddress("NeEvent.Lxa[32]", amp_left_x);
  tree->SetBranchAddress("NeEvent.Lya[32]", amp_left_y);
  tree->SetBranchAddress("NeEvent.Lxt[32]", time_left_x);
  tree->SetBranchAddress("NeEvent.Lyt[32]", time_left_y);
  
  /*Double_t amp_right_x[n_right];
  Double_t amp_right_y[n_right];
  tree->SetBranchAddress("NeEvent.Rxc[16]", amp_right_x);
  tree->SetBranchAddress("NeEvent.Ryc[16]", amp_right_y);*/
  
  Long64_t entries_all = tree->GetEntries();
  std::cout << "\nNumber of all entries: " << entries_all << "\n";
  Long64_t entries_used = entries_all;
  //Long64_t entries_used = 1000;
  std::cout << "\nNumber of entries used: " << entries_used << "\n";

  TH2D* hist_map_right = new TH2D("hist_map_right", "hist_map_right", 16, 0, 16, 16, 0, 16);
  TH2D* hist_map_left = new TH2D("hist_map_left", "hist_map_left", 32, 0, 32, 32, 0, 32);
    
  TH1D* hist_amp_right_x[n_right];
  TH1D* hist_amp_right_y[n_right];
  TH1D* hist_time_right_x[n_right];
  TH1D* hist_time_right_y[n_right];
  
  TH1D* hist_amp_left_x[n_left];
  TH1D* hist_amp_left_y[n_left];
  TH1D* hist_time_left_x[n_left];
  TH1D* hist_time_left_y[n_left];
  
  char* name = new char[20];
    
  for (int i = 0; i < n_right; i++)
  {
  	sprintf(name,"hist_amp_right_x_%d",i);
  	hist_amp_right_x[i] = new TH1D(name, name , 10000, 0, 10000);
  	sprintf(name,"hist_amp_right_y_%d",i);
  	hist_amp_right_y[i] = new TH1D(name, name , 10000, 0, 10000);
  	sprintf(name,"hist_time_right_x_%d",i);
  	hist_time_right_x[i] = new TH1D(name, name , 10000, 0, 10000);
  	sprintf(name,"hist_time_right_y_%d",i);
  	hist_time_right_y[i] = new TH1D(name, name , 10000, 0, 10000);
  }
  
  for (int i = 0; i < n_left; i++)
  {
  	sprintf(name,"hist_amp_left_x_%d",i);
  	hist_amp_left_x[i] = new TH1D(name, name , 10000, 0, 10000);
  	sprintf(name,"hist_amp_left_y_%d",i);
  	hist_amp_left_y[i] = new TH1D(name, name , 10000, 0, 10000);
  	sprintf(name,"hist_time_left_x_%d",i);
  	hist_time_left_x[i] = new TH1D(name, name , 10000, 0, 10000);
  	sprintf(name,"hist_time_left_y_%d",i);
  	hist_time_left_y[i] = new TH1D(name, name , 10000, 0, 10000);
  }
  
  for (int i_entry = 0; i_entry < entries_used; i_entry++)
  {
  	tree->GetEntry(i_entry);
  	if(i_entry%5000 == 0)
  	{
  		std::cout << "entry number: " << i_entry << "\n";
  	}
  	
	for(int strip_x = 0; strip_x < n_right; strip_x++)
	{
		hist_time_right_x[strip_x]->Fill(time_right_x[strip_x]);
		hist_amp_right_x[strip_x]->Fill(amp_right_x[strip_x]);
		for(int strip_y = 0; strip_y < n_right; strip_y++)
		{	
			hist_time_right_y[strip_y]->Fill(time_right_y[strip_y]);
			hist_amp_right_y[strip_y]->Fill(amp_right_y[strip_y]);
			if(amp_right_y[strip_y] > 0 && amp_right_x[strip_x] > 0)
			{
				hist_map_right->Fill(strip_x,strip_y);
			}
		}
	}
  	
	for(int strip_x = 0; strip_x < n_left; strip_x++)
	{
		hist_time_left_x[strip_x]->Fill(time_left_x[strip_x]);
		hist_amp_left_x[strip_x]->Fill(amp_left_x[strip_x]);
		for(int strip_y = 0; strip_y < n_left; strip_y++)
		{	
			hist_time_left_y[strip_y]->Fill(time_left_y[strip_y]);
			hist_amp_left_y[strip_y]->Fill(amp_left_y[strip_y]);
			if(amp_left_y[strip_y] > 0 && amp_left_x[strip_x] > 0)
			{
				hist_map_left->Fill(strip_x,strip_y);
			}
		}
	}
  }
  
  //! Canvas and draw
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  c1->cd();
  hist_map_right->Draw("colz");
  
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
  c2->cd();
  hist_map_left->Draw("colz");
  
  TCanvas *c3 = new TCanvas("c3", "c3", 800, 800);
  c3->cd();
  hist_amp_right_x[0]->Draw();
  
  //! Refill histograms
  TH1D* clone_hist_amp_right_x[n_right];
  TH1D* clone_hist_amp_right_y[n_right];
  TH1D* clone_hist_time_right_x[n_right];
  TH1D* clone_hist_time_right_y[n_right];
  
  TH1D* clone_hist_amp_left_x[n_left];
  TH1D* clone_hist_amp_left_y[n_left];
  TH1D* clone_hist_time_left_x[n_left];
  TH1D* clone_hist_time_left_y[n_left];
  
  for (int i = 0; i < n_right; i++)
  {
    clone_hist_amp_right_x[i] = (TH1D*)hist_amp_right_x[i]->Clone();
    clone_hist_amp_right_y[i] = (TH1D*)hist_amp_right_y[i]->Clone();
    clone_hist_time_right_x[i] = (TH1D*)hist_time_right_x[i]->Clone();
    clone_hist_time_right_y[i] = (TH1D*)hist_time_right_y[i]->Clone();
  }
  
  for (int i = 0; i < n_left; i++)
  {
    clone_hist_amp_left_x[i] = (TH1D*)hist_amp_left_x[i]->Clone();
    clone_hist_amp_left_y[i] = (TH1D*)hist_amp_left_y[i]->Clone();
    clone_hist_time_left_x[i] = (TH1D*)hist_time_left_x[i]->Clone();
    clone_hist_time_left_y[i] = (TH1D*)hist_time_left_y[i]->Clone();
  }
  
  TH2D* clone_hist_map_right = new TH2D();
  clone_hist_map_right = (TH2D*)hist_map_right->Clone();
  
  TH2D* clone_hist_map_left = new TH2D();
  clone_hist_map_left = (TH2D*)hist_map_left->Clone();
  
  TFile* file_results = new TFile("test_0_results.root", "recreate");
  
  TDirectory* dir_amp_right_x = file_results->mkdir("Right Amp X");
  TDirectory* dir_amp_right_y = file_results->mkdir("Right Amp Y");
  TDirectory* dir_time_right_x = file_results->mkdir("Right Time X");
  TDirectory* dir_time_right_y = file_results->mkdir("Right Time Y");
  
  TDirectory* dir_amp_left_x = file_results->mkdir("Left Amp X");
  TDirectory* dir_amp_left_y = file_results->mkdir("Left Amp Y");
  TDirectory* dir_time_left_x = file_results->mkdir("Left Time X");
  TDirectory* dir_time_left_y = file_results->mkdir("Left Time Y");

  for (int i = 0; i < n_right; i++)
  {
    dir_amp_right_x->cd();
    clone_hist_amp_right_x[i]->Write();
    dir_amp_right_y->cd();
    clone_hist_amp_right_y[i]->Write();
    dir_time_right_x->cd();
    clone_hist_time_right_x[i]->Write();
    dir_time_right_y->cd();
    clone_hist_time_right_y[i]->Write();    
  }
  
  for (int i = 0; i < n_left; i++)
  {
    dir_amp_left_x->cd();
    clone_hist_amp_left_x[i]->Write();
    dir_amp_left_y->cd();
    clone_hist_amp_left_y[i]->Write();
    dir_time_left_x->cd();
    clone_hist_time_left_x[i]->Write();
    dir_time_left_y->cd();
    clone_hist_time_left_y[i]->Write();
  }

  file_results->cd();
  clone_hist_map_right->Write();
  clone_hist_map_left->Write();

  file_results->Write();
  file_results->Close();
  
  std::cout << "\ntime: " << timer->RealTime() << " seconds\n";
}
