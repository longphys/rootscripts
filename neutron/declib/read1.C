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

//! Calibration
double ECs137 = 0.477;
double ENa22 = 1.061;

int ChCs137[2] = {365, 283};
int ChNa22[2] = {860, 642};

double a[2], b[2];

double xminCal[2], xmaxCal[2];

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
  TCanvas* cn_alphas = new TCanvas("cn_alphas", "cn_alphas", 800, 800);

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

  TH2F* hmap_alpha = new TH2F("hmap_alpha", "hmap_alpha", 14, 1, 15, 14, 1, 15);
  TH2F* hmap_total = new TH2F("hmap_total", "hmap_total", 14, 1, 15, 14, 1, 15);
  TH2F* hmap = new TH2F("hmap", "hmap", 14, 1, 15, 14, 1, 15);
  TH2F* hmap1 = new TH2F("hmap1", "hmap1", 14, 1, 15, 14, 1, 15);

  double xmin = 0.;
  double xmax = 10000.;
  double bin = (xmax-xmin);
  // double bin = (xmax-xmin)/20.;
  double xminTime = -10.;
  double xmaxTime = 10.;

  TH1D* heff0 = new TH1D("heff0", "heff0", bin, xmin, xmax);
  TH1D* heff1 = new TH1D("heff1", "heff1", bin, xmin, xmax);

  TH1D* h0 = new TH1D("h0", "h0", bin, xmin, xmax);
  TH1D* h1 = new TH1D("h1", "h1", bin, xmin, xmax);
  TH1D* h_cross0 = new TH1D("h_cross0", "h_cross0", bin, xmin, xmax);
  TH1D* h_cross1 = new TH1D("h_cross1", "h_cross1", bin, xmin, xmax);

  TCanvas* cFit = new TCanvas("cFit", "Linear calibration", 800, 600);
  for(int i = 0; i < 2; i++){
    TH1D* hCali = new TH1D("hCali", "Calibration fit", bin, xmin, xmax);
    hCali->SetBinContent(ChCs137[i] - xmin, ECs137);
    hCali->SetBinContent(ChNa22[i] - xmax, ENa22);

    TF1* fCali = new TF1("fCali", "[0]*x + [1]", xmin, xmax);

    cFit->Divide(2,1);
    cFit->cd(1);
    std::cout << "\nENERGY CALIBRATION\n";
    hCali->Fit("fCali");

    a[i] = fCali->GetParameter(0);
    b[i] = fCali->GetParameter(1);

    a[0]=0.0011798;
    b[0]=0.0469636;
    a[1]=0.00162674;
    b[1]=0.0174457;

    xminCal[i] = xmin*a[i]+b[i];
    xmaxCal[i] = xmax*a[i]+b[i];
  }


  TH1D* h0_cal = new TH1D("h0_cal", "h0_cal", bin/10, xminCal[0], xminCal[0]);
  TH1D* h_cross0_cal = new TH1D("h_cross0_cal", "h_cross0_cal", bin/10, xminCal[1], xmaxCal[1]);
  TH1D* h_alpha_amplitude = new TH1D("h_alpha_amplitude", "h_alpha_amplitude", bin, xmin, xmax);

  // TH1D* h_time0 = new TH1D("h_time0", "h_time0", 1000, xminTime, xmaxTime);
  // TH1D* h_time1 = new TH1D("h_time1", "h_time1", 1000, xminTime, xmaxTime);
  // TH1D* h_timetest = new TH1D("h_timetest", "h_timetest", 1000, xminTime, xmaxTime);

  long int entries = tree->GetEntriesFast();
  // entries = (int)1.0e6;
  // entries = (int)5.0e7;
  //std::cout << "test loop:" << std::endl;

  // double threshold_plas_0 = 30.;
  // double threshold_plas_1 = 30.;

  double threshold_plas_0 = 87.;
  double threshold_plas_1 = 81.;

  // double threshold_plas_0 = 300.;
  // double threshold_plas_1 = 300.;

  int beam_xc = 8;
  int beam_yc = 8;

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

    if(de_data[2]>0.){
      h_alpha_amplitude->Fill(de_data[2]);
    }

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
          xy_plst0.push_back(j);
        }
        //! Signal on alpha and Plastic1
        if(de_data[1] > threshold_plas_1){
          xy_plst1.push_back(j);
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

    if (xy_allalphas.size() == 2) //! ONLY SIGNAL IN SILICON
    {
      int x_alpha = xy_allalphas.at(0); //! x[2-15]->x[1-14]
      int y_alpha = xy_allalphas.at(1); //! y[16-29]->y[1-14]
      if (x_alpha < 16 && y_alpha > 15)
      {
        x_alpha -= 1;
        y_alpha -= 15;
        hmap_alpha->Fill(x_alpha, y_alpha);

        // if ((x_alpha == beam_xc) && (y_alpha == beam_yc))
        // {
        //   heff0->Fill(de_data[0]);
        //   heff1->Fill(de_data[1]);
        // }
      }
    }

    //! PLASTIC 0
    if (xy_plst0.size() == 2) //! SIGNAL IN ALPHA AND PLASTIC 0
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
        hmap->Fill(x, y);
        hmap_total->Fill(x, y);

        if ((x == beam_xc)&& (y == beam_yc))
        {
          h0->Fill(de_data[0]);
          h0_cal->Fill(a[0]*(de_data[0]+0.5) + b[0]);
          h1->Fill(de_data[1]);

          // if (de_data[1] < threshold_plas_1)
          // {
          //   h_cross1->Fill(de_data[0]);
          // }
        }

        //! time
        // if (abs(time_data[beam_xc+1])>0 && abs(time_data[beam_yc+15])>0){
        //   h_timetest->Fill(time_data[beam_xc+1] - time_data[beam_yc+15]);
        // }
      }
    }

    //! PLASTIC 1
    if (xy_plst1.size() == 2) //! SIGNAL IN ALPHA AND PLASTIC 1
    {
      int x1 = xy_plst1.at(0);
      int y1 = xy_plst1.at(1);
      if (x1 < 16 && y1 > 15)
      {
        x1-= 1;
        y1 -= 15;
        hmap1->Fill(x1, y1);
        hmap_total->Fill(x1, y1);

        if ((x1 == beam_xc)&& (y1 == beam_yc))
        {
          // h0->Fill(de_data[0]);
          // h1->Fill(de_data[1]);

          //! If undetected in PLASTIC 0 but detected in PLASTIC 1 = PLASTIC 0'S CROSSTALK
          if (de_data[0] < threshold_plas_0)
          {
            h_cross0->Fill(de_data[1]);
            h_cross0_cal->Fill(a[1]*(de_data[1]+0.5) + b[1]);
          }
        }
      }
    }
  }
  cn_pl0->cd();
  hmap->Draw("COLZ");

  cn_alphas->cd();
  hmap_alpha->Draw("COLZ");
  double n_alpha = hmap_alpha->GetBinContent(beam_xc+1, beam_yc+1);
  
  cn_pl1->cd();
  hmap1->Draw("COLZ");

  TCanvas* cn_cross = new TCanvas("cn_cross", "cn_cross", 700, 700);
  cn_cross->Divide(2,1);
  cn_cross->cd(1);
  h_cross0->Draw();
  cn_cross->cd(2);
  h_cross1->Draw();

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
  h0->Draw();
  h0->GetXaxis()->SetRangeUser(threshold_plas_0, xmax);
  double n_0 = h0->Integral(threshold_plas_0, xmax);

  std::cout << "n_0 ["<<beam_xc<<"]["<<beam_yc<<"]= " << n_0 << "; n_alpha["<<beam_xc<<"]["<<beam_yc<<"]= " << n_alpha << "\n";
  std::cout << "Efficiency_0 = " << (n_0/n_alpha)*100 << "%\n";

  double n_cross0 = h_cross0->GetEntries();
  std::cout << "n_cross0 ["<<beam_xc<<"]["<<beam_yc<<"]= " << n_cross0 << "\n";
  std::cout << "cross-talk0 = " << (n_cross0/n_0)*100 << "%\n";

  // h0->Draw();
  // h0->GetXaxis()->SetRangeUser(threshold_plas_0, xmax);
  c->cd(1)->Update();

  // TPaveStats *st_eff0 = (TPaveStats*)heff0->FindObject("stats");
  // st_eff0->SetOptStat(1111111);
  TPaveStats *st0 = (TPaveStats*)h0->FindObject("stats");
  st0->SetOptStat(1111111);
  c->cd(1)->Modified();

  c->cd(2);
  c->cd(2)->SetLogy();

  // heff1->Draw();
  // heff1->GetXaxis()->SetRangeUser(threshold_plas_1, xmax);
  // double n_1 = heff1->Integral(threshold_plas_1, xmax);
  h1->Draw();
  h1->GetXaxis()->SetRangeUser(threshold_plas_1, xmax);
  double n_1 = h1->Integral(threshold_plas_1, xmax);

  std::cout << "n_1 ["<<beam_xc<<"]["<<beam_yc<<"]= " << n_1 << "; n_alpha["<<beam_xc<<"]["<<beam_yc<<"]= " << n_alpha << "\n";
  std::cout << "Efficiency_1 = " << (n_1/n_alpha)*100 << "%\n";

  double n_cross1 = h_cross1->GetEntries();
  std::cout << "n_cross1 ["<<beam_xc<<"]["<<beam_yc<<"]= " << n_cross1 << "\n";
  std::cout << "cross-talk1 = " << (n_cross1/n_1) << "%\n";

  // h1->Draw();
  // h1->GetXaxis()->SetRangeUser(threshold_plas_1, xmax);
  c->cd(2)->Update();

  // TPaveStats *st_eff1 = (TPaveStats*)heff1->FindObject("stats");
  // st_eff1->SetOptStat(1111111);
  TPaveStats *st1 = (TPaveStats*)h1->FindObject("stats");
  st1->SetOptStat(1111111);
  c->cd(2)->Modified();

  // std::cout << "counttest = " << counttest << " = " << (double)(counttest*100)/entries << "% of " << entries << " entries.\n";

  TCanvas* c_total1 = new TCanvas("c_total1", "c_total1", 1200, 1000);
  TCanvas* c_total2 = new TCanvas("c_total2", "c_total2", 1200, 1000);
  TCanvas* c_total3 = new TCanvas("c_total3", "c_total3", 1200, 1000);

  TCanvas* c_total = new TCanvas("c_total", "c_total", 1700, 500);
  c_total->Divide(3,1);
  // c_total->cd(1);
  // c_total->cd(1)->SetRightMargin(0.15);
  c_total1->cd();
  c_total1->SetRightMargin(0.15);
  hmap_alpha->SetTitle("Amplitude from alpha recoils");

  hmap_alpha->GetXaxis()->SetTitle("X");
  hmap_alpha->GetXaxis()->SetLabelFont(42);
  hmap_alpha->GetXaxis()->SetTitleFont(42);
  hmap_alpha->GetXaxis()->SetTitleSize(0.04);
  hmap_alpha->GetXaxis()->CenterTitle(true);

  hmap_alpha->GetYaxis()->SetTitle("Y");
  hmap_alpha->GetYaxis()->SetLabelFont(42);
  hmap_alpha->GetYaxis()->SetTitleFont(42);
  hmap_alpha->GetYaxis()->SetTitleSize(0.04);
  hmap_alpha->GetYaxis()->CenterTitle(true);
  hmap_alpha->SetStats(0);
  hmap_alpha->Draw("col z");

  // c_total->cd(2);
  // c_total->cd(2)->SetRightMargin(0.15);
  c_total2->cd();
  c_total2->SetRightMargin(0.15);
  hmap_total->SetTitle("Amplitude from silicon and plastic detectors");
  hmap_total->GetXaxis()->SetTitle("X");
  hmap_total->GetXaxis()->SetLabelFont(42);
  hmap_total->GetXaxis()->SetTitleFont(42);
  hmap_total->GetXaxis()->SetTitleSize(0.04);
  hmap_total->GetXaxis()->CenterTitle(true);

  hmap_total->GetYaxis()->SetTitle("Y");
  hmap_total->GetYaxis()->SetLabelFont(42);
  hmap_total->GetYaxis()->SetTitleFont(42);
  hmap_total->GetYaxis()->SetTitleSize(0.04);
  hmap_total->GetYaxis()->CenterTitle(true);
  hmap_total->SetStats(0);
  hmap_total->Draw("col z");

  // c_total->cd(3);
  // c_total->cd(3)->SetRightMargin(0.15);
  c_total3->cd();
  c_total3->SetRightMargin(0.15);
  hmap->SetTitle("Amplitude from silicon and the centered plastic detector");
  hmap->GetXaxis()->SetTitle("X");
  hmap->GetXaxis()->SetLabelFont(42);
  hmap->GetXaxis()->SetTitleFont(42);
  hmap->GetXaxis()->SetTitleSize(0.04);
  hmap->GetXaxis()->CenterTitle(true);

  hmap->GetYaxis()->SetTitle("Y");
  hmap->GetYaxis()->SetLabelFont(42);
  hmap->GetYaxis()->SetTitleFont(42);
  hmap->GetYaxis()->SetTitleSize(0.04);
  hmap->GetYaxis()->CenterTitle(true);
  hmap->SetStats(0);
  hmap->Draw("col z");

  TCanvas* c_crosstalk_calibrated1 = new TCanvas("c_crosstalk_calibrated1", "c_crosstalk_calibrated1", 1100, 1000);
  c_crosstalk_calibrated1->SetTicks();
  c_crosstalk_calibrated1->SetGrid();
  TCanvas* c_crosstalk_calibrated2 = new TCanvas("c_crosstalk_calibrated2", "c_crosstalk_calibrated2", 1100, 1000);
  c_crosstalk_calibrated2->SetTicks();
  c_crosstalk_calibrated2->SetGrid();

  // TCanvas* c_crosstalk_calibrated = new TCanvas("c_crosstalk_calibrated", "c_crosstalk_calibrated", 1500, 700);
  // c_crosstalk_calibrated->Divide(2,1);
  // c_crosstalk_calibrated->cd(1);
  c_crosstalk_calibrated1->cd();
  h0_cal->SetTitle("");
  h0_cal->SetLineWidth(2);
  h0_cal->SetLineColor(kBlue);
  h0_cal->SetStats(0);
  h0_cal->GetXaxis()->SetTitle("Energy (MeV)");
  h0_cal->GetXaxis()->SetLabelFont(42);
  h0_cal->GetXaxis()->SetTitleFont(42);
  h0_cal->GetXaxis()->SetTitleSize(0.04);
  h0_cal->GetXaxis()->CenterTitle(true);
  h0_cal->GetYaxis()->SetTitle("Events");
  h0_cal->GetYaxis()->SetLabelFont(42);
  h0_cal->GetYaxis()->SetTitleFont(42);
  h0_cal->GetYaxis()->SetTitleSize(0.04);
  h0_cal->GetYaxis()->CenterTitle(true);
  h0_cal->Draw();

  // c_crosstalk_calibrated->cd(2);
  c_crosstalk_calibrated2->cd();
  h_cross0_cal->SetTitle("");
  h_cross0_cal->SetLineWidth(2);
  h_cross0_cal->SetLineColor(kRed);
  h_cross0_cal->SetStats(0);
  h_cross0_cal->GetXaxis()->SetTitle("Energy (MeV)");
  h_cross0_cal->GetXaxis()->SetLabelFont(42);
  h_cross0_cal->GetXaxis()->SetTitleFont(42);
  h_cross0_cal->GetXaxis()->SetTitleSize(0.04);
  h_cross0_cal->GetXaxis()->CenterTitle(true);
  h_cross0_cal->GetYaxis()->SetTitle("Events");
  h_cross0_cal->GetYaxis()->SetLabelFont(42);
  h_cross0_cal->GetYaxis()->SetTitleFont(42);
  h_cross0_cal->GetYaxis()->SetTitleSize(0.04);
  h_cross0_cal->GetYaxis()->CenterTitle(true);
  h_cross0_cal->Draw();

  TLegend *legend_SAVE = new TLegend(0.25, 0.55, 0.75, 0.85);
  legend_SAVE->SetBorderSize(0);
  legend_SAVE->SetLineWidth(2);
  legend_SAVE->SetTextSize(0.04);
  legend_SAVE->AddEntry(h0_cal, "Light ouput of the centered module", "l");
  legend_SAVE->AddEntry(h_cross0_cal, "Cross-talk to side module", "l");

  legend_SAVE->Draw();

  TCanvas* c_alpha_amp = new TCanvas("c_alpha_amp", "c_alpha_amp", 1100, 1000);
  h_alpha_amplitude->SetTitle("");
  h_alpha_amplitude->SetLineWidth(2);
  h_alpha_amplitude->SetLineColor(kBlack);
  h_alpha_amplitude->SetStats(0);
  h_alpha_amplitude->GetXaxis()->SetTitle("Channel");
  h_alpha_amplitude->GetXaxis()->SetLabelFont(42);
  h_alpha_amplitude->GetXaxis()->SetTitleFont(42);
  h_alpha_amplitude->GetXaxis()->SetTitleSize(0.04);
  h_alpha_amplitude->GetXaxis()->CenterTitle(true);
  h_alpha_amplitude->GetYaxis()->SetTitle("Events");
  h_alpha_amplitude->GetYaxis()->SetLabelFont(42);
  h_alpha_amplitude->GetYaxis()->SetTitleFont(42);
  h_alpha_amplitude->GetYaxis()->SetTitleSize(0.04);
  h_alpha_amplitude->GetYaxis()->CenterTitle(true);
  h_alpha_amplitude->Draw();

  std::cout << "\ntime: " << timer->RealTime() << " (s)\n";
}
