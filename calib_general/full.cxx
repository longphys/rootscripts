#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLatex.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

vector<double> readCSV(const string& filename) {
    vector<double> data;
    ifstream infile(filename);
    double value;
    while (infile >> value) data.push_back(value);
    infile.close();
    return data;
}

double CAL_A = 0.;   
double CAL_B = 0.; 

void full() {

    const int NPEAKS = 4;

    cout << "\n=== BPX61 Alpha Detector Analysis ===\n" << endl;

    // --------------------------------------------------
    // Load data
    vector<double> data =
        readCSV("file.csv");

    // Raw spectrum
    TH1F* h = new TH1F("h", "^{226}Ra Alpha Spectrum;Channel;Counts", data.size(), 0, data.size()
    );

    for (int i = 0; i < data.size(); i++)
        h->SetBinContent(i + 1, data[i]);

    // Approximate peak positions and known energies
/*    double start_pos[NPEAKS] = {1066, 1236, 1298, 1456, 1978};
    double known_E[NPEAKS]   = {4.784, 5.304, 5.490, 6.002, 7.687}; // MeV
*/
    double start_pos[NPEAKS] = {1066, 1298, 1456, 1978};
    double known_E[NPEAKS]   = {4.784, 5.490, 6.002, 7.687}; 
    double ch[NPEAKS];
    double fwhm_ch[NPEAKS];

    TCanvas* c1 = new TCanvas("c1", "BPX61 Analysis", 900, 800);
    c1->Divide(2,2);
    gStyle->SetOptStat(0);

    // 1) Gaussian fits in channel 
    c1->cd(1);
    h->SetLineWidth(2);
    h->Draw("hist");

    for (int i = 0; i < NPEAKS; i++) {

        TF1* g = new TF1(Form("g%d", i), "gaus",
                         start_pos[i] - 30, start_pos[i] + 30);

        g->SetParameters(h->GetMaximum(), start_pos[i], 20);
        h->Fit(g, "RQ");

        ch[i]      = g->GetParameter(1);
        fwhm_ch[i] = 2.355 * g->GetParameter(2);

        g->SetLineColor(kRed);
        g->Draw("same");
    }

    // 2) Energy calibration
    c1->cd(2);

    TGraph* gcal = new TGraph(NPEAKS, ch, known_E);
    gcal->SetMarkerStyle(20);
    gcal->Draw("AP");

    TF1* fcal = new TF1("fcal", "[0] + [1]*x", 900, 2100);
    gcal->Fit(fcal, "R");

    CAL_B = fcal->GetParameter(0);
    CAL_A = fcal->GetParameter(1);

    TLatex tex;
    tex.SetNDC();
    tex.SetTextSize(0.04);
    tex.DrawLatex(0.15, 0.85, Form("E = %.6f ch + %.4f", CAL_A, CAL_B));

    // 3) Calibrated spectrum
    double Emin = CAL_B;
    double Emax = CAL_A * data.size() + CAL_B;

    TH1F* hE = new TH1F(
        "hE", "^{226}Ra Alpha Spectrum (Calibrated);Energy [MeV];Counts",
        data.size(), Emin, Emax
    );

    for (int i = 0; i < data.size(); i++)
        // hE->Fill(CAL_A * i + CAL_B, data[i]);
        hE->SetBinContent(i + 1, data[i]);

    TCanvas* c2 = new TCanvas("c2", "Calibrated Spectrum", 900, 600);
    hE->SetLineWidth(2);
    hE->Draw("hist");

    double E[NPEAKS];
    double fwhm_E[NPEAKS];
    double fwhm_E_err[NPEAKS];
    double rel_E[NPEAKS];
    double rel_E_err[NPEAKS];

    for (int i = 0; i < NPEAKS; i++) {

        double E0 = CAL_A * start_pos[i] + CAL_B;

        TF1* gE = new TF1(Form("gE%d", i), "gaus",E0 - 0.065, E0 + 0.065);

        gE->SetParameters(hE->GetMaximum(), E0, 0.03);
        hE->Fit(gE, "RQ+");

        double sigma     = gE->GetParameter(2);
        double sigma_err = gE->GetParError(2);

        E[i] = gE->GetParameter(1);

        fwhm_E[i]     = 2.355*sigma;
        fwhm_E_err[i] = 2.355*sigma_err;

        rel_E[i]     = 100*fwhm_E[i]/E[i];
        rel_E_err[i] = rel_E[i]*(fwhm_E_err[i]/fwhm_E[i]);

        cout << fixed << setprecision(4)
             << "Peak " << i+1
             << ": E = " << E[i]
             << " MeV, sigma = " << sigma
             << " ± " << sigma_err << " MeV"
             << endl;

        gE->SetLineColor(kRed);
        gE->Draw("same");
    }

    // 4) Resolution
    c1->cd(3);
    TGraphErrors* gFWHM = new TGraphErrors(NPEAKS, E, fwhm_E, nullptr, fwhm_E_err);
    gFWHM->SetTitle("BPX61 Resolution;Energy [MeV];FWHM [MeV]");
    gFWHM->SetMarkerStyle(20);
    gFWHM->Draw("AP");

    c1->cd(4);
    TGraphErrors* gRel = new TGraphErrors(NPEAKS, E, rel_E, nullptr, rel_E_err);
    gRel->SetTitle("BPX61 Relative Resolution;Energy [MeV];#DeltaE/E [%]");
    gRel->SetMarkerStyle(20);
    gRel->Draw("AP");

    c1->Update();

    cout << "\n=== Analysis Complete ===\n" << endl;
}
