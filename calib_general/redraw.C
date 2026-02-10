#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

// double a = 0.00320272;
// double b = 1.35799;

double a = 0.00299768;
double b = 1.9309;

double e_min = a*0 + b;
double e_max = a*4096 + b;
	
double en_1 = 4.6;
double en_2 = 4.784;
double en_3 = 5.304;
double en_4 = 5.49;
double en_5 = 6.002;
double en_6 = 7.687;

void redraw() {
    std::ifstream file("g4.csv");
    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    std::vector<double> values;
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        double value;
        if (ss >> value)
        {
            // std::cout << "Read value: " << value << std::endl;
            // sleep(1); // Pause for 1 second
            values.push_back(value);
        }
    }
    file.close();

    int nbins = values.size();
    if (nbins == 0) return;

    TCanvas* c_csv = new TCanvas("c_csv", "CSV Data Plot", 1200, 600);
    // c_csv->Divide(3,1);
    // c_csv->cd(1);
    TH1D* h_csv = new TH1D("h_csv", "CSV Data Histogram",
                           nbins, e_min, e_max);

    for (int i = 0; i < nbins; ++i)
    {
        h_csv->SetBinContent(i + 1, values[i]);
    }
    h_csv->SetStats(0);
    h_csv->GetXaxis()->SetTitle("Energy (MeV)");
    h_csv->GetYaxis()->SetTitle("Counts");

    // h_csv->Draw("HIST");

    TF1* f_csv = new TF1("gaussian", "gaus", e_min, e_max);
    h_csv->Fit("gaussian", "+", "", en_2-0.1, en_2+0.1);
    double centroid_1 = f_csv->GetParameter(1);
    double centroid_1_error = f_csv->GetParError(1);
    double width_1 = f_csv->GetParameter(2)*2.355; // FWHM
    double width_1_error = f_csv->GetParError(2)*2.355;
    h_csv->Fit("gaussian", "+", "", en_4-0.05, en_4+0.15);
    double centroid_2 = f_csv->GetParameter(1);
    double centroid_2_error = f_csv->GetParError(1);
    double width_2 = f_csv->GetParameter(2)*2.355; // FWHM
    double width_2_error = f_csv->GetParError(2)*2.355;
    h_csv->Fit("gaussian", "+", "", en_5-0.07, en_5+0.13);
    double centroid_3 = f_csv->GetParameter(1);
    double centroid_3_error = f_csv->GetParError(1);
    double width_3 = f_csv->GetParameter(2)*2.355; // FWHM
    double width_3_error = f_csv->GetParError(2)*2.355;
    h_csv->Fit("gaussian", "+", "", en_6-0.1, en_6+0.1);
    double centroid_4 = f_csv->GetParameter(1);
    double centroid_4_error = f_csv->GetParError(1);
    double width_4 = f_csv->GetParameter(2)*2.355; // FWHM
    double width_4_error = f_csv->GetParError(2)*2.355;
    c_csv->Modified();
    c_csv->Update();

    std::cout << "Centroid 1: " << centroid_1 << " +- " << centroid_1_error << "\n";
    std::cout << "Centroid 2: " << centroid_2 << " +- " << centroid_2_error << "\n";
    std::cout << "Centroid 3: " << centroid_3 << " +- " << centroid_3_error << "\n";
    std::cout << "Centroid 4: " << centroid_4 << " +- " << centroid_4_error << "\n";

    std::cout << "FWHM 1: " << width_1 << " +- " << width_1_error << "\n";
    std::cout << "FWHM 2: " << width_2 << " +- " << width_2_error << "\n";
    std::cout << "FWHM 3: " << width_3 << " +- " << width_3_error << "\n";
    std::cout << "FWHM 4: " << width_4 << " +- " << width_4_error << "\n";

    TCanvas* c_2 = new TCanvas("c_2", "FWHM Plot", 1200, 600);
    // c_csv->cd(2);
    TGraphErrors* g_sigma = new TGraphErrors();
    g_sigma->SetPoint(0, centroid_1, width_1);
    g_sigma->SetPointError(0, centroid_1_error, width_1_error);
    g_sigma->SetPoint(1, centroid_2, width_2);
    g_sigma->SetPointError(1, centroid_2_error, width_2_error);
    g_sigma->SetPoint(2, centroid_3, width_3);
    g_sigma->SetPointError(2, centroid_3_error, width_3_error);
    g_sigma->SetPoint(3, centroid_4, width_4);
    g_sigma->SetPointError(3, centroid_4_error, width_4_error);
    g_sigma->SetTitle("FWHM vs Energy;Energy (MeV);Sigma (MeV)");
    g_sigma->SetMarkerStyle(21);
    g_sigma->Draw("AP");
    c_csv->Modified();
    c_csv->Update();  

    TCanvas* c_3 = new TCanvas("c_3", "FWHM (%) Plot", 1200, 600);
    c_csv->cd(3);
    TGraphErrors* g_resolution = new TGraphErrors();
    g_resolution->SetPoint(0, centroid_1, 100*width_1/centroid_1);
    // g_resolution->SetPointError(0, centroid_1_error, (width_1_error*centroid_1 - width_1*centroid_1_error)/(centroid_1*centroid_1));
    g_resolution->SetPoint(1, centroid_2, 100*width_2/centroid_2);
    // g_resolution->SetPointError(1, centroid_2_error, (width_2_error*centroid_2 - width_2*centroid_2_error)/(centroid_2*centroid_2));
    g_resolution->SetPoint(2, centroid_3, 100*width_3/centroid_3);
    // g_resolution->SetPointError(2, centroid_3_error, (width_3_error*centroid_3 - width_3*centroid_3_error)/(centroid_3*centroid_3));
    g_resolution->SetPoint(3, centroid_4, 100*width_4/centroid_4);
    // g_resolution->SetPointError(3, centroid_4_error, (width_4_error*centroid_4 - width_4*centroid_4_error)/(centroid_4*centroid_4));
    g_resolution->SetTitle("Resolution vs Energy;Energy (MeV);Resolution (%)");
    g_resolution->SetMarkerStyle(21);
    g_resolution->Draw("AP");
    c_csv->Modified();
    c_csv->Update();
}