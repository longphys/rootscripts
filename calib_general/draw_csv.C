#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

void draw_csv() {
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

    TCanvas* c_csv = new TCanvas("c_csv", "CSV Data Plot", 800, 600);
    TH1D* h_csv = new TH1D("h_csv", "CSV Data Histogram",
                           nbins, 0, nbins);

    for (int i = 0; i < nbins; ++i)
    {
        h_csv->SetBinContent(i + 1, values[i]);
    }
    h_csv->SetStats(0);
    h_csv->GetXaxis()->SetTitle("Bin index");
    h_csv->GetYaxis()->SetTitle("Counts");

    // h_csv->Draw("HIST");

    TF1* f_csv = new TF1("gaussian", "gaus", 0, nbins);
    // h_csv->Fit("gaussian", "+", "", 1040, 1090);
    h_csv->Fit("gaussian", "+", "", 1040, 1090);
    double centroid_1 = f_csv->GetParameter(1);
    double width_1 = f_csv->GetParameter(2);
    // h_csv->Fit("gaussian", "+", "", 1270, 1320);
    h_csv->Fit("gaussian", "+", "", 1270, 1320);
    double centroid_2 = f_csv->GetParameter(1);
    double width_2 = f_csv->GetParameter(2);
    // h_csv->Fit("gaussian", "+", "", 1415, 1485);
    h_csv->Fit("gaussian", "+", "", 1415, 1485);
    double centroid_3 = f_csv->GetParameter(1);
    double width_3 = f_csv->GetParameter(2);
    // h_csv->Fit("gaussian", "+", "", 1945, 2000);
    h_csv->Fit("gaussian", "+", "", 1945, 2000);
    double centroid_4 = f_csv->GetParameter(1);
    double width_4 = f_csv->GetParameter(2);
    c_csv->Modified();
    c_csv->Update();

    std::cout << "Centroid 1: " << centroid_1 << " +- " << width_1 << "\n";
    std::cout << "Centroid 2: " << centroid_2 << " +- " << width_2 << "\n";
    std::cout << "Centroid 3: " << centroid_3 << " +- " << width_3 << "\n";
    std::cout << "Centroid 4: " << centroid_4 << " +- " << width_4 << "\n";
}
