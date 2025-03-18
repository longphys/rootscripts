#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"

#include "TStopwatch.h"

void refill_copy()
{
  auto timer = new TStopwatch();
  timer->Start();

  //! Files and trees
  TFile* fileOpen = new TFile("./measurement_Cs137.root", "read");
  // TFile* fileOpen = new TFile("./measurement_Na22.root", "read");
  TTree* treeOpen =  (TTree*) fileOpen->Get("Events");

  //! Get Events
	double eventOpen;
	treeOpen->SetBranchAddress("Amplitude", &eventOpen);

  double x_min = 100.;
  double x_max = 1000.;
  int bin = x_max-x_min;
  TH1D* histogram = new TH1D("histogram", "histogram", bin, x_min, x_max);
  for(int i = 0; i<treeOpen->GetEntries(); i++){
    if(i%100000==0){std::cout << "Event: " << i << "\n";}
    treeOpen->GetEntry(i);
    histogram->Fill(eventOpen);
  }
  TCanvas* canvas = new TCanvas();
  canvas->cd();
  histogram->Draw();

  //! Fill experiment histograms
  std::ofstream fileSave("measurement_ascii_Cs137.txt"); // Open file for writing

  std::cout << "histogram->GetXaxis()->GetXmin() = " << histogram->GetXaxis()->GetXmin() << "\n";
  std::cout << "histogram->GetXaxis()->GetXmax() = " << histogram->GetXaxis()->GetXmax() << "\n";  

  for(int i = histogram->GetXaxis()->GetXmin(); i<histogram->GetXaxis()->GetXmax(); i++){
    // std::cout << "i = " << i << "; bin to be counted = " << i-histogram->GetXaxis()->GetXmin() + 1 << "\n";
    fileSave << i << " " << histogram->GetBinContent(i-histogram->GetXaxis()->GetXmin()+1) << "\n";
  }

  fileSave.close();

  std::ifstream inputFile("measurement_ascii_Cs137.txt");

  int x;
  double y;
  std::vector<int> bin_input;
  std::vector<double> content_input;
  while(inputFile >> x >> y){
    std::cout << x << " " << y << "\n";
    bin_input.push_back(x);
    content_input.push_back(y);
  }

  std::cout << "bin_input.size() = " << bin_input.size() << "\n";
  std::cout << "bin_input[0] = " << bin_input[0] << "\n";
  std::cout << "bin_input[bin_input.size()-1] = " << bin_input[bin_input.size()-1] << "\n";

  TH1D* new_histogram = new TH1D("new_histogram", "new_histogram", bin_input.size(), bin_input[0], bin_input[bin_input.size()-1]);
  for (int i=1; i<=bin_input.size(); i++){
    new_histogram->SetBinContent(i, content_input[i-1]);
  }

  TCanvas* canvas_1 = new TCanvas();
  canvas_1->cd();
  new_histogram->Draw();

  std::cout << "time: " << timer->RealTime() << " seconds \n";
}