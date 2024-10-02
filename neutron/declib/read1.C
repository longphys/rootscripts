#include "DecLibTest.hh"
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TStopwatch.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <sstream>
#include <fstream>
#include <iostream>

int64_t EventNr;
vector <double> dE;
vector <double> dTime;
vector <unsigned int> channel;
int startingChannelNum = 0;
int endingChannelNum = 29;

void read1()
{
  auto timer = new TStopwatch();
  timer->Start();

  // TFile* file = new TFile("../../../EfficiencyPlasticFLNR/declib/output/test06_2_Plast_BC404_14x_14y_out_test.root","open");
  // TFile* file = new TFile("../../../EfficiencyPlasticFLNR/declib/output/test06_2_Plast_BC404_14x_14y_out_test1.root","open");
  TFile* file = new TFile("../../../EfficiencyPlasticFLNR/declib/output/test08_2_Plast_BC404_14x_14y_Pos_x0_y0_out.root","open");
  // TFile* file = new TFile("../../../EfficiencyPlasticFLNR/declib/output/test09_2_Plast_BC404_14x_14y_Pos_x20_y-20_out.root","open");
  // TFile* file = new TFile("../../../EfficiencyPlasticFLNR/declib/output/test10_2_Plast_BC404_14x_14y_Pos_x0_y-92_out.root","open");
  // TFile* file = new TFile("../../../EfficiencyPlasticFLNR/declib/output/test11_2_Plast_BC404_14x_14y_Pos_x0_y-92_wo_upper_detector-ch1_out.root","open");
  TTree* tree = (TTree*)file->Get("ETree");

  // TH1D* E[endingChannelNum-startingChannelNum+1];
  // TH1D* Time[endingChannelNum-startingChannelNum+1];

  TCanvas* c = new TCanvas("c", "c", 1600, 800);
  c->Divide(2,1);
  
  TCanvas* cn_pl0 = new TCanvas("cn_pl0", "cn_pl0", 800, 800);
  TCanvas* cn_pl1 = new TCanvas("cn_pl1", "cn_pl1", 800, 800);

  // TFile* fileSave = new TFile("./output/test06_2_Plast_BC404_14x_14y_result.root", "recreate");
  // // TFile* fileSave = new TFile("./output/test06_2_Plast_BC404_14x_14y_result_threshold.root", "recreate");

  // for (int i = startingChannelNum; i <= endingChannelNum; i++){
  //   auto num = std::to_string(i);
  //   std::string scriptE = "E" + num;
  //   const char* scriptEChar = scriptE.c_str();
    
  //   E[i] = new TH1D();

  //   std::string drawScriptE = "EDep[" + num + "]>>E[" + num + "]";
  //   const char* drawScriptEChar = drawScriptE.c_str();

  //   std::string drawScriptECon = "EDep[" + num + "]>30";
  //   const char* drawScripEConChar = drawScriptECon.c_str();

  //   c->cd(1);
  //   tree->Draw(drawScriptEChar);
  //   // tree->Draw(drawScriptEChar,drawScripEConChar);
  //   c->Modified();
  //   c->Update();
  // }

  // for (int i = startingChannelNum; i <= endingChannelNum; i++){
  //   auto num = std::to_string(i);
  //   std::string scriptTime = "Time" + num;
  //   const char* scriptTimeChar = scriptTime.c_str();
    
  //   Time[i] = new TH1D();

  //   std::string drawScriptTime = "Time[" + num + "]>>Time[" + num + "]";
  //   const char* drawScriptTimeChar = drawScriptTime.c_str();

  //   std::string drawScriptTimeCon = "Time[" + num + "]>30";
  //   const char* drawScripTimeConChar = drawScriptTimeCon.c_str();

  //   c->cd(2);
  //   tree->Draw(drawScriptTimeChar);
  //   c->Modified();
  //   c->Update();
  // }

  // TH1D* hist0 = new TH1D();
  // tree->Draw("EDep[0]>>hist0");

  // fileSave->Write();
  // fileSave->Close();

  std::vector <double> *de = nullptr;
  std::vector <double> *ch = nullptr;
  std::vector <double> *time = nullptr;

  tree->SetBranchAddress("EDep", &de);
  tree->SetBranchAddress("Channel", &ch);
  tree->SetBranchAddress("Time", &time);

  TH2F* hmap_alpha = new TH2F("hmap_alpha", "hmap_alpha", 16, 0, 16, 16, 0, 16);
  TH2F* hmap = new TH2F("hmap", "hmap", 16, 0, 16, 16, 0, 16);
  TH2F* hmap1 = new TH2F("hmap1", "hmap1", 16, 0, 16, 16, 0, 16);

  double xmin = 0.;
  double xmax = 10000.;
  double xminTime = -10.;
  double xmaxTime = 10.;

  TH1D* heff0 = new TH1D("heff0", "heff0", 10000, xmin, xmax);
  TH1D* heff1 = new TH1D("heff1", "heff1", 10000, xmin, xmax);

  TH1D* h0 = new TH1D("h0", "h0", 10000, xmin, xmax);
  TH1D* h1 = new TH1D("h1", "h1", 10000, xmin, xmax);
  TH1D* h_cross = new TH1D("h_cross", "h_cross", 10000, xmin, xmax);

  TH1D* h_time0 = new TH1D("h_time0", "h_time0", 1000, xminTime, xmaxTime);
  TH1D* h_time1 = new TH1D("h_time1", "h_time1", 1000, xminTime, xmaxTime);
  TH1D* h_timetest = new TH1D("h_timetest", "h_timetest", 1000, xminTime, xmaxTime);

  long int entries = tree->GetEntriesFast();
  entries = (int)1.0e5;
  //std::cout << "test loop:" << std::endl;

  double threshold_plas_0 = 30.;
  double threshold_plas_1 = 30.;

  int beam_xc = 8;
  int beam_yc = 9;

  int counttest = 0;

  for (long int i = 0; i < entries; i++)
  {
    // std::cout << "Entry: " << i << "\n";
    if (i%100000 == 0)
    {
      std::cout << "event " << i << " of " << entries << std::endl;
    }

    tree->GetEntry(i);

    int de_size = (int)de->size(); // MUST be done for every tree entry
    double *de_data = de->data(); 

    int time_size = (int)time->size();
    double *time_data = time->data();

    //! Fill time data
    // if(abs(time_data[0])>0){
    //   h_time0->Fill(time_data[0]);
    // }

    // if(abs(time_data[1])>0){
    //   h_time1->Fill(time_data[1]);
    // }


    // std::cout << "de->size() = " << de_size << std::endl;
    // std::cout << "time->size() = " << time_size << std::endl;
    std::vector <int> xy_allalphas;
    std::vector <int> xy_plst0;
    std::vector <int> xy_plst1;
    std::vector <int> x_plst0;
    std::vector <int> x_plst1;
    std::vector <int> y_plst0;
    std::vector <int> y_plst1;
    std::vector <int> xy_plst1_cross;

    for (int j = 2; j < de_size; j++)
    {
      // if (de_data[j] > 1)
      // if (de_plas0 > threshold_plas_0 || de_plas1 > threshold_plas_1)
      // threshold_plas_0 = -1.0;
      if (de_data[j] > 1 ) //! Signal on alpha
      {
        xy_allalphas.push_back(j);

        //! Signal on alpha and Plastic0
        if(de_data[0] > threshold_plas_0){
          //! vector 2 elements [x_channel, y_channel]
          // xy_plst0.push_back(j);
        }
        //! Signal on alpha and Plastic1
        if(de_data[1] > threshold_plas_1){
          // xy_plst1.push_back(j);
        }
      }

//! TESTING
      // if (j < 16){ //! If lands on x-strip
      //   if (de_data[j] > 1 && (de_data[0] > threshold_plas_0))
      //   {
      //     // std::cout << "x_0 = " << j << "\n";
      //     x_plst0.push_back(j);
      //   }

      //   if (de_data[j] > 1 && (de_data[1] > threshold_plas_1))
      //   {
      //     // std::cout << "x_1 = " << j << "\n";
      //     x_plst1.push_back(j);
      //   }
      // }
      // else { //!If lands on y-strip
      //   if (de_data[j] > 1 && (de_data[0] > threshold_plas_0))
      //   {
      //     // std::cout << "y_0 = " << j << "\n";
      //     y_plst0.push_back(j);
      //   }

      //   if (de_data[j] > 1 && (de_data[1] > threshold_plas_1))
      //   {
      //     // std::cout << "y_1 = " << j << "\n";
      //     y_plst1.push_back(j);
      //   }
      // }

//! TESTING
    }

//! TESTING
    // if ((x_plst0.size() > 1 || y_plst0.size() > 1)){
    //   counttest ++;
    //   // std::cout << "Entry: " << i << "------------\n";
    //   // std::cout << "size_x0 = " << x_plst0.size() << "; size_y0 = " << y_plst0.size() << "\n";

    //   for (int i_element = 0; i_element < x_plst0.size(); i_element++){
    //     // std::cout << "x_plst0[" << i_element <<"] = " << x_plst0[i_element] << "\n";
    //   }
    //   for (int i_element = 0; i_element < y_plst0.size(); i_element++){
    //     // std::cout << "y_plst0[" << i_element <<"] = " << y_plst0[i_element] << "\n";
    //   }
    // }
    // if ((x_plst1.size() > 1) || (y_plst1.size() > 1)){
    //   counttest ++;
    //   // std::cout << "Entry: " << i << "------------\n";
    //   // std::cout << "size_x1 = " << x_plst1.size() << "; size_y1 = " << y_plst1.size() << "\n";
    //   for (int i_element = 0; i_element < x_plst1.size(); i_element++){
    //     // std::cout << "x_plst1[" << i_element <<"] = " << x_plst1[i_element] << "\n";
    //   }
    //   for (int i_element = 0; i_element < y_plst1.size(); i_element++){
    //     // std::cout << "y_plst1[" << i_element <<"] = " << y_plst1[i_element] << "\n";
    //   }
    // }
//! TESTING

    if (xy_allalphas.size() == 2) //! should only happen if vector has 2 elements
    {
      int x_alpha = xy_allalphas.at(0); //! x[2-15]->x[1-14]
      int y_alpha = xy_allalphas.at(1); //! y[16-29]->y[1-14]
      if (x_alpha < 16 && y_alpha > 15)
      {
        x_alpha -= 1;
        y_alpha -= 15;
        hmap_alpha->Fill(x_alpha, y_alpha);

        if ((x_alpha == beam_xc)&& (y_alpha == beam_yc))
        {
          heff0->Fill(de_data[0]);
          heff1->Fill(de_data[1]);
        }
      }
    }

    if (xy_plst0.size() == 2) //! should only happen if vector has 2 elements
    {
      int x = xy_plst0.at(0);
      int y = xy_plst0.at(1);
      // if (x >= 15 || y <= 16 || y >= 29)
      // {
      //   std::cout << "Event: " << i << "; first: " << x << "; second: " << y << std::endl;
      // }
      if (x < 16 && y > 15)
      {
        x -= 1;
        y -= 15;
        // hmap->Fill(x, y);

        if ((x == beam_xc)&& (y == beam_yc))
        {
          // h0->Fill(de_data[0]);
          // h1->Fill(de_data[1]);

          // if (de_data[1] < threshold_plas_1)
          // {
          //   h_cross->Fill(de_data[0]);
          // }
        }

        //! time
        // if (abs(time_data[beam_xc+1])>0 && abs(time_data[beam_yc+15])>0){
        //   h_timetest->Fill(time_data[beam_xc+1] - time_data[beam_yc+15]);
        // }
      }
    }

    if (xy_plst1.size() == 2) //! should only happen if vector has 2 elements
    {
      int x1 = xy_plst1.at(0);
      int y1 = xy_plst1.at(1);
      if (x1 < 16 && y1 > 15)
      {
        x1-= 1;
        y1 -= 15;
        // hmap1->Fill(x1, y1);
        if ((x1 == beam_xc)&& (y1 == beam_yc))
        {
          // h0->Fill(de_data[0]);
          // h1->Fill(de_data[1]);

          // if (de_data[0] < threshold_plas_0)
          // {
          //   h_cross->Fill(de_data[1]);
          // }
        }
      }
    }
  }
  cn_pl0->cd();
  // hmap->Draw("COLZ");
  hmap_alpha->Draw("COLZ");
  double n_alpha = hmap_alpha->GetBinContent(beam_xc+1, beam_yc+1);
  
  cn_pl1->cd();
  hmap1->Draw("COLZ");

  // TCanvas* cn_cross = new TCanvas("cn_cross", "cn_cross", 700, 700);
  // cn_cross->cd();
  // h_cross->Draw();

  // TCanvas* ctime = new TCanvas("ctime", "ctime", 800, 800);
  // ctime->Divide(3,1);
  // ctime->cd(1);
  // h_time0->Draw();
  // ctime->cd(2);
  // h_time1->Draw();
  // ctime->cd();
  // h_timetest->Draw();

  c->cd(1);
  c->cd(1)->SetLogy();
  // heff0->Draw();
  // heff0->GetXaxis()->SetRangeUser(threshold_plas_0, xmax);
  // double n_0 = heff0->Integral(threshold_plas_0, xmax);

  // std::cout << "n_0 ["<<beam_xc<<"]["<<beam_yc<<"]= " << n_0 << "; n_alpha["<<beam_xc<<"]["<<beam_yc<<"]= " << n_alpha << "\n";
  // std::cout << "Efficiency_0 = " << (n_0/n_alpha)*100 << "%\n";

  h0->Draw();
  h0->GetXaxis()->SetRangeUser(threshold_plas_0, xmax);
  c->cd(1)->Update();

  // TPaveStats *st_eff0 = (TPaveStats*)heff0->FindObject("stats");
  // st_eff0->SetOptStat(1111111);
  // TPaveStats *st0 = (TPaveStats*)h0->FindObject("stats");
  // st0->SetOptStat(1111111);
  c->cd(1)->Modified();

  c->cd(2);
  c->cd(2)->SetLogy();

  heff1->Draw();
  heff1->GetXaxis()->SetRangeUser(threshold_plas_1, xmax);
  double n_1 = heff1->Integral(threshold_plas_1, xmax);

  std::cout << "n_1 ["<<beam_xc<<"]["<<beam_yc<<"]= " << n_1 << "; n_alpha["<<beam_xc<<"]["<<beam_yc<<"]= " << n_alpha << "\n";
  std::cout << "Efficiency_1 = " << (n_1/n_alpha)*100 << "%\n";

  // h1->Draw();
  // h1->GetXaxis()->SetRangeUser(threshold_plas_1, xmax);
  c->cd(2)->Update();

  TPaveStats *st_eff1 = (TPaveStats*)heff1->FindObject("stats");
  st_eff1->SetOptStat(1111111);
  // TPaveStats *st1 = (TPaveStats*)h1->FindObject("stats");
  // st1->SetOptStat(1111111);
  c->cd(2)->Modified();

  std::cout << "counttest = " << counttest << " = " << (double)(counttest*100)/entries << "% of " << entries << " entries.\n";

  std::cout << "\ntime: " << timer->RealTime() << " (s)\n";
}
