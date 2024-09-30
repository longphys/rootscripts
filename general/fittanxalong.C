#include<iostream>
#include<fstream>

#include<vector>
#include"TROOT.h"
#include"TMath.h"
#include"TStyle.h"
#include"TAxis.h"
#include"TH1.h"
#include"TF1.h"
#include"TCanvas.h"
#include"TGraph.h"
#include"TLegend.h"
#include"TPaveText.h"

    double fGaus1(double*x,double*par){ //first gaus peak: red
    return par[0]*TMath::Exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]));
    }

    double fGaus2(double*x,double*par){ //second gaus peak: black
    return par[0]*TMath::Exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]));
    }

    double fPol(double*x, double*par){ //polynomial scatter background
    return par[0] + x[0]*par[1] + x[0]*x[0]*par[2] + x[0]*x[0]*x[0]*par[3]+pow(x[0],4)*par[4];
    }

    double fTotal(double*x,double*par){ //total curve = gaus1 + gaus2 + pol
    return fGaus1(x,par) + fGaus2(x,&par[3]) + fPol(x,&par[6]); //&par[i] is using memory point i of the par array
    }
int fittanxalong(char*fin)
{        
    TGraph *gr = new TGraph();

    gr->SetMarkerStyle(kDot);
    gr->SetMarkerColor(kRed);
    gr->SetTitle("Scattering Energy Spectrum");
    gr->GetXaxis()->SetTitle("Energy (keV)");
	gr->GetYaxis()->SetTitle("Count");
    
    fstream file;
    file.open(fin, ios::in); //file in

    while(1)
    {
        double x, y;
        file >> x >> y; //input from file to x y
        gr->SetPoint(gr->GetN(), x, y); //GetN get the maximum number of datapoint
        if(file.eof()) break;
    }

    file.close();

    gr->GetYaxis()->SetRangeUser(0, 1400); //setting the maximum Y axis to 1400

    double fitMin = 100.;
    double fitMax = 330.; //fit range
    int nparam = 11; //parameter arrays for fit
  
    double A1= 100.; //amplitude, mean and standard deviation of the gaus peak
    double mean1 = 225.0;
    double sigma1 = 15.;
    double A2= 100.;
    double mean2 = 155.0;
    double sigma2 = 20.;
    double p0=22.; //poly parameters
    double p1=-8.;
    double p2=1.;
    double p3=0.;
    double p4=0.;

    double params[] ={A1,mean1,sigma1,A2,mean2,sigma2,p0,p1,p2,p3,p4}; 
    double pars[11];
    TF1*fG1 = new TF1("fG1",fGaus1,fitMin,fitMax,3); //fit for each single function
    fG1->SetLineColor(kRed);
    TF1*fG2 = new TF1("fG2",fGaus2,fitMin,fitMax,3);
    fG2->SetLineColor(kBlack);
    TF1*fBg = new TF1("fBg",fPol,fitMin,fitMax,5);
    fBg->SetLineColor(kGreen);

    TF1*func = new TF1("func",fTotal,fitMin,fitMax,nparam); //fit for the total function
    func->SetParameters(params); //setting paramaters from earlier
    func->SetParLimits(0,1.,700);
    func->SetParLimits(1,200,250.);
    func->SetParLimits(2,10,50); //these are fitting limit for the the first 6 parameteters, which are from the 2 gaus peaks
    func->SetParLimits(3,1.,200);
    func->SetParLimits(4,125.,170.);
    func->SetParLimits(5,10,50);
    func->SetLineWidth(1); //customize fit line width, color and style
    func->SetLineColor(kBlue);
    func->SetLineStyle(1);
    
    TCanvas *c1 = new TCanvas("canvas", "Scattering Energy Spectrum Canvas", 1280, 720);

    c1->SetTickx();
    c1->SetTicky();
    c1->SetGridx();
    c1->SetGridy();

    gr->Draw("AP");
    gr->Fit("func","R");

    double chi2 = func->GetChisquare();
    double ndf = func->GetNDF();

    std::cout << chi2/ndf << std::endl; // calculating reduced-chisquare

    //TLegend *leg = new TLegend(0.3, 0.7, 0.7, 0.7);
    //leg->SetBorderSize(0);
    //leg->AddEntry(gr, "Measured Data", "f");
    //leg->AddEntry(func, "Fit Function", "l");
    //leg->Draw();

    TPaveText *fitInfo = new TPaveText(0.7, 0.7, 0.83, 0.83, "NDC");
    fitInfo->SetFillColor(0);
    fitInfo->SetTextAlign(12); //12=left align
    fitInfo->SetTextFont(42); //font 42 according to default
    fitInfo->SetTextSize(0.03); //% of canvas size

    fitInfo->AddText(Form("r-Chi2: %.2f", chi2/ndf));
    //fitInfo->AddText(Form("NDF: %.0f", ndf));
    fitInfo->Draw();

    func->GetParameters(pars);
    fG1->SetParameters(pars);
    fG2->SetParameters(&pars[3]);
    fBg->SetParameters(&pars[6]);
    fG1->Draw("same");
    fG2->Draw("same");
    fBg->Draw("same");

    //c1->SaveAs("output.eps");
    return 0;
}