#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TVector3.h"
#include "TVirtualFFT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TLegend.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <sstream>
#include <getopt.h>
#include <cstdlib>  // For atof, atoi
#include <cstring>  // For strcmp

// Name of input file
std::string name_f_sim_1 = "";
std::string name_f_sim_2 = "";
std::string name_f_mea_1 = "";
std::string name_f_mea_2 = "";

// Default channel of the measurement file
int channel = 0;

// Default number of entries for FFT
int entries_sim_fft = 100000;
int entries_mea_fft = 100000;

// Default number of entries for gradient descent
int entries_sim_descent = 500000;
int entries_mea_descent = 500000;

// Default options for all histograms
double x_min_sim = 0.;
double x_max_sim = 1.3; // in MeV
int bin_sim = (x_max_sim - x_min_sim)*1000.;

double x_min_mea = 0.;
double x_max_mea = 1000.; // in channels
int bin_mea = (x_max_mea - x_min_mea);

// Default energies of Compton edges in MeV (Cs_137 and Na_22)
double energy_1 = 0.477;
double energy_2 = 1.061;

// Default first value for energy resolution in percentage
double res_1 = 8.;
double res_2 = 5.5;

// Default channel window to perform FFT
double x_min_fft_1 = 200;
double x_max_fft_1 = x_max_mea;
double x_min_fft_2 = 400;
double x_max_fft_2 = x_max_mea;

// Default energy window to perform gradient descent
double x_min_descent_1 = 0.35;
double x_max_descent_1 = 0.7;
double x_min_descent_2 = 0.8;
double x_max_descent_2 = 1.3;

// Default FFT cut-off growth rate for estimating Compton edge
double rate_fft_ini_1 = 0.1;
double rate_fft_ini_2 = 0.1;

// Default FFT threshold for estimating Compton edge
double thresh_fft_ini_1 = 50.;
double thresh_fft_ini_2 = 50.;

// Default FFT cut-off growth rate for performing gradient descent
double rate_fft_descent_1 = 0.1;
double rate_fft_descent_2 = 0.1;

// Default FFT threshold for performing gradient descent
double thresh_fft_descent_sim = 800.;
double thresh_fft_descent_mea = 800.;

// Default learning rate for gradient descent
double learning_rate_channel_1 = 2.;
double learning_rate_channel_2 = 2.; 
double learning_rate_res_1 = 0.05;
double learning_rate_res_2 = 0.05;

// Default step for calculating gradient
double step_channel_1 = 1.;
double step_channel_2 = 1.;
double step_res_1;
double step_res_2;

// Default scaling factor of histograms
double scale_sim_1 = 1.;
double scale_sim_2 = 1.;
double scale_mea_1 = 1.;
double scale_mea_2 = 1.;

bool CheckRootFile(const char* filename) {
    TFile file(filename);
    return !file.IsZombie();  // Check if the file exists
}

// Function to read the macro file and simulate command-line arguments
std::vector<std::string> ReadMacroFile(const char* filename) {
    std::vector<std::string> args;
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Cannot open macro file: " << filename << std::endl;
        exit(1);
    }
    
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string word;
        while (iss >> word) {
            args.push_back(word);  // Store each word as an argument
        }
    }
    
    return args;
}

enum OptionIDs {
    opt_name_f_sim_1 = 1000,
    opt_name_f_sim_2,
    opt_name_f_mea_1,
    opt_name_f_mea_2,
    opt_channel,
    opt_entries_sim_fft,
    opt_entries_mea_fft,
    opt_entries_sim_descent,
    opt_entries_mea_descent,
    opt_x_min_sim,
    opt_x_max_sim,
    opt_x_min_mea,
    opt_x_max_mea,
    opt_energy_1,
    opt_energy_2,
    opt_res_1,
    opt_res_2,
    opt_x_min_fft_1,
    opt_x_max_fft_1,
    opt_x_min_fft_2,
    opt_x_max_fft_2,
    opt_x_min_descent_1,
    opt_x_max_descent_1,
    opt_x_min_descent_2,
    opt_x_max_descent_2,
    opt_rate_fft_ini_1,
    opt_rate_fft_ini_2,
    opt_thresh_fft_ini_1,
    opt_thresh_fft_ini_2,
    opt_rate_fft_descent_1,
    opt_rate_fft_descent_2,
    opt_thresh_fft_descent_1,
    opt_thresh_fft_descent_2,
    opt_learning_rate_channel_1,
    opt_learning_rate_channel_2,
    opt_learning_rate_res_1,
    opt_learning_rate_res_2,
    opt_scale_sim_1,
    opt_scale_sim_2,
    opt_scale_mea_1,
    opt_scale_mea_2,
    opt_step_channel_1,
    opt_step_channel_2,
    opt_step_res_1,
    opt_step_res_2
};

template <typename T>
void ProcessArgument(const char* optarg, T& variable, T defaultValue, const std::string& name, const std::string& unit = "") {
    if (optarg == nullptr || optarg[0] == '-') {
        std::cerr << "Default " << name << " used: " << defaultValue << " " << unit << "\n";
        variable = defaultValue;
    } else {
        variable = atof(optarg);
        std::cerr << name << " used: " << variable << " " << unit << "\n";
    }
}

int ReadArgs(int argc, char** argv)
{
    static struct option longOptions[] = 
    {
        {"name_f_sim_1", required_argument, 0, opt_name_f_sim_1},
        {"name_f_sim_2", required_argument, 0, opt_name_f_sim_2},
        {"name_f_mea_1", required_argument, 0, opt_name_f_mea_1},
        {"name_f_mea_2", required_argument, 0, opt_name_f_mea_2},
        {"channel", required_argument, 0, opt_channel},
        {"entries_sim_fft", required_argument, 0, opt_entries_sim_fft},
        {"entries_mea_fft", required_argument, 0, opt_entries_mea_fft},
        {"entries_sim_descent", required_argument, 0, opt_entries_sim_descent},
        {"entries_mea_descent", required_argument, 0, opt_entries_mea_descent},
        {"x_min_sim", required_argument, 0, opt_x_min_sim},
        {"x_max_sim", required_argument, 0, opt_x_max_sim},
        {"x_min_mea", required_argument, 0, opt_x_min_mea},
        {"x_max_mea", required_argument, 0, opt_x_max_mea},
        {"energy_1", required_argument, 0, opt_energy_1},
        {"energy_2", required_argument, 0, opt_energy_2},
        {"res_1", required_argument, 0, opt_res_1},
        {"res_2", required_argument, 0, opt_res_2},
        {"x_min_fft_1", required_argument, 0, opt_x_min_fft_1},
        {"x_max_fft_1", required_argument, 0, opt_x_max_fft_1},
        {"x_min_fft_2", required_argument, 0, opt_x_min_fft_2},
        {"x_max_fft_2", required_argument, 0, opt_x_max_fft_2},
        {"x_min_descent_1", required_argument, 0, opt_x_min_descent_1},
        {"x_max_descent_1", required_argument, 0, opt_x_max_descent_1},
        {"x_min_descent_2", required_argument, 0, opt_x_min_descent_2},
        {"x_max_descent_2", required_argument, 0, opt_x_max_descent_2},
        {"rate_fft_ini_1", required_argument, 0, opt_rate_fft_ini_1},
        {"rate_fft_ini_2", required_argument, 0, opt_rate_fft_ini_2},
        {"thresh_fft_ini_1", required_argument, 0, opt_thresh_fft_ini_1},
        {"thresh_fft_ini_2", required_argument, 0, opt_thresh_fft_ini_2},
        {"rate_fft_descent_1", required_argument, 0, opt_rate_fft_descent_1},
        {"rate_fft_descent_2", required_argument, 0, opt_rate_fft_descent_2},
        {"thresh_fft_descent_sim", required_argument, 0, opt_thresh_fft_descent_1},
        {"thresh_fft_descent_mea", required_argument, 0, opt_thresh_fft_descent_2},
        {"learning_rate_channel_1", required_argument, 0, opt_learning_rate_channel_1},
        {"learning_rate_channel_2", required_argument, 0, opt_learning_rate_channel_2},
        {"learning_rate_res_1", required_argument, 0, opt_learning_rate_res_1},
        {"learning_rate_res_2", required_argument, 0, opt_learning_rate_res_2},
        {"step_channel_1", required_argument, 0, opt_step_channel_1},
        {"step_channel_2", required_argument, 0, opt_step_channel_2},
        {"step_res_1", required_argument, 0, opt_step_res_1},
        {"step_res_2", required_argument, 0, opt_step_res_2},
        {"scale_sim_1", required_argument, 0, opt_scale_sim_1},
        {"scale_sim_2", required_argument, 0, opt_scale_sim_2},
        {"scale_mea_1", required_argument, 0, opt_scale_mea_1},
        {"scale_mea_2", required_argument, 0, opt_scale_mea_2},
        {NULL, 0, NULL, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "", longOptions, 0)) != -1) {
        switch(opt) {
            case opt_name_f_sim_1:
                if (optarg == nullptr || optarg[0] == '-') {
                    std::cerr << "Error: --opt_name_f_sim_1 requires an argument but none was provided.\n";
                    return -1;
                }
                name_f_sim_1 = optarg;
                if (!CheckRootFile(name_f_sim_1.c_str())) {
                    std::cerr << "Error: Invalid simulated file 1: " << name_f_sim_1 << std::endl;
                    return -1;
                }
                break;

            case opt_name_f_sim_2:
                if (optarg == nullptr || optarg[0] == '-') {
                    std::cerr << "Error: --opt_name_f_sim_2 requires an argument but none was provided.\n";
                    return -1;
                }
                name_f_sim_2 = optarg;
                if (!CheckRootFile(name_f_sim_2.c_str())) {
                    std::cerr << "Error: Invalid simulated file 2: " << name_f_sim_2 << std::endl;
                    return -1;
                }
                break;
                
            case opt_name_f_mea_1:
                if (optarg == nullptr || optarg[0] == '-') {
                    std::cerr << "Error: --opt_name_f_mea_1 requires an argument but none was provided.\n";
                    return -1;
                }
                name_f_mea_1 = optarg;
                if (!CheckRootFile(name_f_mea_1.c_str())) {
                    std::cerr << "Error: Invalid measurement file 1: " << name_f_mea_1 << std::endl;
                    return -1;
                }
                break;

            case opt_name_f_mea_2:
                if (optarg == nullptr || optarg[0] == '-') {
                    std::cerr << "Error: --opt_name_f_mea_2 requires an argument but none was provided.\n";
                    return -1;
                }
                name_f_mea_2 = optarg;
                if (!CheckRootFile(name_f_mea_2.c_str())) {
                    std::cerr << "Error: Invalid measurement file 2: " << name_f_mea_2 << std::endl;
                    return -1;
                }
                break;
                
			case opt_channel:
				ProcessArgument(optarg, channel, channel, "channel");
				break;
				
			case opt_entries_sim_fft:
				ProcessArgument(optarg, entries_sim_fft, entries_sim_fft, "number of simulation entries for FFT");
				break;

			case opt_entries_mea_fft:
				ProcessArgument(optarg, entries_mea_fft, entries_mea_fft, "number of measurement entries for FFT");
				break;
				
			case opt_entries_sim_descent:
				ProcessArgument(optarg, entries_sim_descent, entries_sim_descent, "number of simulation entries for gradient descent");
				break;

			case opt_entries_mea_descent:
				ProcessArgument(optarg, entries_mea_descent, entries_mea_descent, "number of measurement entries for gradient descent");
				break;

			case opt_x_min_sim:
				ProcessArgument(optarg, x_min_sim, x_min_sim, "x_min for simulation");
        bin_sim = (x_max_sim - x_min_sim)*1000.;
				break;

			case opt_x_max_sim:
				ProcessArgument(optarg, x_max_sim, x_max_sim, "x_max for simulation");
        bin_sim = (x_max_sim - x_min_sim)*1000.;
				break;

			case opt_x_min_mea:
				ProcessArgument(optarg, x_min_mea, x_min_mea, "x_min for measurement");
        bin_mea = (x_max_mea - x_min_mea);
				break;

			case opt_x_max_mea:
				ProcessArgument(optarg, x_max_mea, x_max_mea, "x_max for measurement");
        bin_mea = (x_max_mea - x_min_mea);
				break;

			case opt_energy_1:
				ProcessArgument(optarg, energy_1, energy_1, "Energy of the first Compton edge", "(MeV)");
				break;

			case opt_energy_2:
				ProcessArgument(optarg, energy_2, energy_2, "Energy of the second Compton edge", "(MeV)");
				break;
				
			case opt_res_1:
				ProcessArgument(optarg, res_1, res_1, "Energy resolution of the first Compton edge", "(%)");
				break;
				
			case opt_res_2:
				ProcessArgument(optarg, res_2, res_2, "Energy resolution of the second Compton edge", "(%)");
				break;
                
			case opt_x_min_fft_1:
				ProcessArgument(optarg, x_min_fft_1, x_min_fft_1, "Minimum channel for FFT 1");
				break;
				
			case opt_x_max_fft_1:
				ProcessArgument(optarg, x_max_fft_1, x_max_fft_1, "Maximum channel for FFT 1");
				break;
				
			case opt_x_min_fft_2:
				ProcessArgument(optarg, x_min_fft_2, x_min_fft_2, "Minimum channel for FFT 2");
				break;
				
			case opt_x_max_fft_2:
				ProcessArgument(optarg, x_max_fft_2, x_max_fft_2, "Maximum channel for FFT 2");
				break;
				
			case opt_x_min_descent_1:
				ProcessArgument(optarg, x_min_descent_1, x_min_descent_1, "Minimum channel for gradient descent 1");
				break;
				
			case opt_x_max_descent_1:
				ProcessArgument(optarg, x_max_descent_1, x_max_descent_1, "Maximum channel for gradient descent 1");
				break;
				
			case opt_x_min_descent_2:
				ProcessArgument(optarg, x_min_descent_2, x_min_descent_2, "Minimum channel for gradient descent 2");
				break;
				
			case opt_x_max_descent_2:
				ProcessArgument(optarg, x_max_descent_2, x_max_descent_2, "Maximum channel for gradient descent 2");
				break;
				
			case opt_rate_fft_ini_1:
				ProcessArgument(optarg, rate_fft_ini_1, rate_fft_ini_1, "FFT cut-off growth rate 1 (estimation)");
				break;
				
			case opt_rate_fft_ini_2:
				ProcessArgument(optarg, rate_fft_ini_2, rate_fft_ini_2, "FFT cut-off growth rate 2 (estimation)");
				break;
				
			case opt_thresh_fft_ini_1:
				ProcessArgument(optarg, thresh_fft_ini_1, thresh_fft_ini_1, "FFT cut-off threshold 1 (estimation)");
				break;
				
			case opt_thresh_fft_ini_2:
				ProcessArgument(optarg, thresh_fft_ini_2, thresh_fft_ini_2, "FFT cut-off threshold 2 (estimation)");
				break;
				
			case opt_rate_fft_descent_1:
				ProcessArgument(optarg, rate_fft_descent_1, rate_fft_descent_1, "FFT cut-off growth rate 1 (descent)");
				break;
				
			case opt_rate_fft_descent_2:
				ProcessArgument(optarg, rate_fft_descent_2, rate_fft_descent_2, "FFT cut-off growth rate 2 (descent)");
				break;
				
			case opt_thresh_fft_descent_1:
				ProcessArgument(optarg, thresh_fft_descent_sim, thresh_fft_descent_sim, "FFT cut-off threshold 1 (descent)");
				break;
				
			case opt_thresh_fft_descent_2:
				ProcessArgument(optarg, thresh_fft_descent_mea, thresh_fft_descent_mea, "FFT cut-off threshold 2 (descent)");
				break;
				
			case opt_learning_rate_channel_1:
				ProcessArgument(optarg, learning_rate_channel_1, learning_rate_channel_1, "Gradient descent learning rate for estimating channel 1");
				break;
				
			case opt_learning_rate_channel_2:
				ProcessArgument(optarg, learning_rate_channel_2, learning_rate_channel_2, "Gradient descent learning rate for estimating channel 2");
				break;
				
			case opt_learning_rate_res_1:
				ProcessArgument(optarg, learning_rate_res_1, learning_rate_res_1, "Gradient descent learning rate for estimating energy resolution 1");
				break;
				
			case opt_learning_rate_res_2:
				ProcessArgument(optarg, learning_rate_res_2, learning_rate_res_2, "Gradient descent learning rate for estimating energy resolution 2");
				break;
				
			case opt_step_channel_1:
				ProcessArgument(optarg, step_channel_1, step_channel_1, "Step of channel 1 for calculating gradient");
				break;
				
			case opt_step_channel_2:
				ProcessArgument(optarg, step_channel_2, step_channel_2, "Step of channel 2 for calculating gradient");
				break;
				
			case opt_step_res_1:
				ProcessArgument(optarg, step_res_1, step_res_1, "Step of energy resolution 1 for calculating gradient");
				break;
				
			case opt_step_res_2:
				ProcessArgument(optarg, step_res_2, step_res_2, "Step of energy resolution 2 for calculating gradient");
				break;
				
			case opt_scale_sim_1:
				ProcessArgument(optarg, scale_sim_1, scale_sim_1, "Scaling factor for simulation histogram 1");
				break;
				
			case opt_scale_sim_2:
				ProcessArgument(optarg, scale_sim_2, scale_sim_2, "Scaling factor for simulation histogram 2");
				break;
				
			case opt_scale_mea_1:
				ProcessArgument(optarg, scale_mea_1, scale_mea_1, "Scaling factor for measurement histogram 1");
				break;
				
			case opt_scale_mea_2:
				ProcessArgument(optarg, scale_mea_2, scale_mea_2, "Scaling factor for measurement histogram 2");
				break;
				
            case '?': // Handle unknown options
                return -1;
        }
    }

    // Ensure both simulated files are provided
    if (name_f_sim_1.empty() || name_f_sim_2.empty() || name_f_mea_1.empty() || name_f_mea_2.empty()) {
        std::cerr << "Error: --opt_name_f_sim_1 --opt_name_f_sim_2 --opt_name_f_mea_1 --opt_name_f_mea_2 must be specified.\n";
        return -1;
    }
    
    return 0;
}

// Function that performs fast fourier transformation (FFT) on a histogram h_channel.
TH1D* fft(TH1D* h_channel, double para_k, double para_c, std::string name_h_channel, std::string title_h_channel)
{
	int bin = h_channel->GetNbinsX();
	int bin_even = 2*bin;
	
	double x_min = h_channel->GetXaxis()->GetXmin();
	double x_max = h_channel->GetXaxis()->GetXmax();
	
	// std::cout << "x_min =" << x_min << "\n";
	// std::cout << "x_max =" << x_max << "\n";
	
	const char* name_h_channel_char = name_h_channel.c_str();
	const char* title_h_channel_char = title_h_channel.c_str();

	TH1D* h_filtered = new TH1D(name_h_channel_char, title_h_channel_char, bin, x_min, x_max);
	TH1D* h_channel_even = new TH1D("h_channel_even", "Transformed to even function", bin_even, -(x_max - x_min), x_max - x_min);

	for(int i = 1; i <= bin; i++)
	{
		h_channel_even->SetBinContent(i, h_channel->GetBinContent(i));
		h_channel_even->SetBinContent(bin_even - i, h_channel->GetBinContent(i));
	}
	// h_channel_even->SetBinContent(bin_mea, h_channel_even->GetBinContent(bin_mea + 1));
	h_channel_even->SetBinContent(bin, h_channel_even->GetBinContent(bin + 1));

	TF1* f_logistic = new TF1("f_logistic", "1./(1.+exp([0]*(x-[1])))", 0, bin_even);
	f_logistic->SetParameter(0, para_k);
	f_logistic->SetParameter(1, para_c);
	f_logistic->SetNpx(10000);

	//! Magnitude
	TH1 *hm = nullptr;
	TVirtualFFT::SetTransform(nullptr);
	hm = h_channel_even->FFT(hm, "MAG");

	TH1D* newhm = new TH1D("newhm", "newhm", bin_even, -(x_max - x_min), x_max - x_min);
	for(int i = 1; i <= bin_even; i++)
	{
		newhm->SetBinContent(i, hm->GetBinContent(i)/bin_even);
	}

	TH1D* hLogis = new TH1D("hLogis", "Logis", bin_even, 0, bin_even);
	for(int i = 1; i <= bin; i++)
	{
		hLogis->SetBinContent(i, f_logistic->Eval(i));
		hLogis->SetBinContent(bin_even-i, f_logistic->Eval(i));
	}

	hLogis->SetLineColor(kRed);
	hLogis->Draw("same");

	//! Apply threshold
	TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
	Double_t *re_full = new Double_t[bin_even];
	Double_t *im_full = new Double_t[bin_even];

	fft->GetPointsComplex(re_full, im_full);

	for(int i = 1; i <= bin_even; i++)
	{
		re_full[i] = re_full[i]*hLogis->GetBinContent(i);
		im_full[i] = im_full[i]*hLogis->GetBinContent(i);
	}

	TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &bin_even, "C2R M K");
	fft_back->SetPointsComplex(re_full,im_full);
	fft_back->Transform();
	TH1 *hb = nullptr;
	hb = TH1::TransformHisto(fft_back,hb,"RE");

	TH1D* newhb = new TH1D("newhb", "newhb", bin_even, -(x_max - x_min), x_max - x_min);
	for(int i = 1; i <= bin_even; i++)
	{
		newhb->SetBinContent(i, hb->GetBinContent(i)/bin_even);
	}

	for(int i = 1; i <= bin; i++)
	{
		h_filtered->SetBinContent(i, newhb->GetBinContent(i));
	}

  // TCanvas* c_test = new TCanvas("c_test", "c_test", 800, 600);
  // c_test->cd();
  // h_channel->Draw();

	delete h_channel_even;
	delete hm;
	delete newhm;
	delete f_logistic;
	delete hLogis;
	delete[] re_full;
	delete[] im_full;
	delete fft_back;
	delete hb;
	delete newhb;

	return h_filtered;

	delete name_h_channel_char;
	delete title_h_channel_char;
}

// Function that returns the first derivative of a histogram h_channel.
TH1D* Diff(TH1D* h_channel, std::string name_h_channel, std::string title_h_channel)
{
	int bin = h_channel->GetNbinsX();	
	double x_min = h_channel->GetXaxis()->GetXmin();
	double x_max = h_channel->GetXaxis()->GetXmax();
	
	const char* name_h_channel_char = name_h_channel.c_str();
	const char* title_h_channel_char = title_h_channel.c_str();
	double new_bin_content;
	double bin_length = (x_max-x_min)/bin;

	TH1D* h_diff = new TH1D(name_h_channel_char, title_h_channel_char, bin, x_min, x_max);
	for(int i = 3; i <= bin - 2; i++)
	{
		new_bin_content = (-h_channel->GetBinContent(i+2) 
		+ 8*h_channel->GetBinContent(i+1) 
		- 8*h_channel->GetBinContent(i-1) 
		+ h_channel->GetBinContent(i-2))/12*bin_length;
		h_diff->SetBinContent(i, new_bin_content);
	}
	return h_diff;
	delete name_h_channel_char;
	delete title_h_channel_char;
}

// Function to perform unweighted chi2 test and return chi2/ndf between h_channel_1 and h_channel_2.
double chi2(TH1D* h_channel_1, TH1D* h_channel_2){
	double chi2_result = h_channel_1->Chi2Test(h_channel_2, "UU CHI2/NDF");
	return chi2_result;
}

// Function that performs linear calibration of energy and returns coefficients.
double fit_energy(int channel_1, int channel_2, TCanvas* c_fit, int i){
	TH1D* h_calibrate = new TH1D("h_calibrate", "Calibration fit", bin_mea, x_min_mea, x_max_mea);
	h_calibrate->SetBinContent(channel_1 - x_min_mea, energy_1);
	h_calibrate->SetBinContent(channel_2 - x_min_mea, energy_2);

	TF1* f_linear = new TF1("f_linear", "[0]*x + [1]", x_min_mea, x_max_mea);
	c_fit->cd(1);
	h_calibrate->Fit("f_linear", "Q");

	double coef_a = f_linear->GetParameter(0);
	double coef_b = f_linear->GetParameter(1);

	c_fit->Modified();
	c_fit->Update();

	delete h_calibrate;
	delete f_linear;

	if(i == 0){return coef_a;}
	else{return coef_b;}
}

// Function that performs calibration of energy resolution and returns coefficients.
double fit_energy_res(double res_sig_1, double res_sig_2, TH1D* h_1, TH1D* h_2, double x_min, double x_max, TCanvas* c_fit, int i){
	TH1D* h_res = new TH1D("h_res", "Resolution fit", bin_mea, x_min, x_max);
	h_res->SetBinContent(h_1->FindBin(energy_1), res_sig_1/100.);
	h_res->SetBinContent(h_2->FindBin(energy_2), res_sig_2/100.);

	//! Choose fit function
	TF1* f_res = new TF1("f_res", "sqrt([0]*[0] + [1]*[1]/(x*x))");
	c_fit->cd(2);
	h_res->Fit("f_res", "Q", "", energy_1-0.1, energy_2+0.1);

	//! Correspond to chosen function
	double coef_a = std::abs(f_res->GetParameter(0));
	double coef_b = 0.;
	double coef_c = std::abs(f_res->GetParameter(1));

	c_fit->Modified();
	c_fit->Update();

	delete h_res;
	delete f_res;

	if(i == 0){return coef_a;}
	else if(i == 1){return coef_b;}
	else{return coef_c;}
}

void GetTrees()
{
	std::cerr << "Sim1: " << name_f_sim_1 << "\n";
	std::cerr << "Sim2: " << name_f_sim_2 << "\n";
	std::cerr << "Mea1: " << name_f_mea_1 << "\n";
	std::cerr << "Mea2: " << name_f_mea_2 << "\n";
	
	//! Read the simulation and measurement files
	TFile* f_sim_1 = new TFile(name_f_sim_1.c_str(), "read");
	TFile* f_sim_2 = new TFile(name_f_sim_2.c_str(), "read");
	
	TFile* f_mea_1 = new TFile(name_f_mea_1.c_str(), "read");
	TFile* f_mea_2 = new TFile(name_f_mea_2.c_str(), "read");

	TCanvas* c_fit = new TCanvas("c_fit", "c_fit", 1000, 500);
	c_fit->Divide(2,1);

	//! Select trees and branches from the files
	TTree* t_sim_1 = (TTree*) f_sim_1->Get("dEEtree");
	TTree* t_mea_1 = (TTree*) f_mea_1->Get("AnalysisxTree");

	UShort_t event_mea_1[48];
	t_mea_1->SetBranchAddress("NeEvent.neutAmp[48]", event_mea_1);

	double event_sim_1;
	t_sim_1->SetBranchAddress("Scintillator", &event_sim_1);

	TTree* t_sim_2 =  (TTree*) f_sim_2->Get("dEEtree");
	TTree* t_mea_2 =  (TTree*) f_mea_2->Get("AnalysisxTree");

	UShort_t event_mea_2[48];
	t_mea_2->SetBranchAddress("NeEvent.neutAmp[48]", event_mea_2);

	double event_sim_2;
	t_sim_2->SetBranchAddress("Scintillator", &event_sim_2);
	
	//! Fill histograms for FFT to estimate Compton edge
	TH1D* h_mea_1_estimate = new TH1D("h_mea_1_estimate", "First measurement histogram", bin_mea, x_min_mea, x_max_mea);

	for(int i = 0; i < entries_mea_fft; i++){
		t_mea_1->GetEntry(i);
		h_mea_1_estimate->Fill(event_mea_1[channel]);
	}

	TH1D* h_mea_1_estimate_filtered = fft(h_mea_1_estimate, rate_fft_ini_1, thresh_fft_ini_1,
	"h_mea_1_estimate_filtered", "First measurement histogram filtered");
	
	h_mea_1_estimate_filtered->GetXaxis()->SetTitle("Channel");
	h_mea_1_estimate_filtered->GetYaxis()->SetTitle("Count");
	TH1D* h_mea_1_estimate_diff = Diff(h_mea_1_estimate_filtered, "h_mea_1_estimate_diff", "Derivative of First measurement histogram");

	TH1D* h_mea_2_estimate = new TH1D("h_mea_2_estimate", "Second measurement histogram", bin_mea, x_min_mea, x_max_mea);
	for(int i = 0; i < entries_mea_fft; i++){
		t_mea_2->GetEntry(i);
		h_mea_2_estimate->Fill(event_mea_2[channel]);
	}
	
	TH1D* h_mea_2_estimate_filtered = fft(h_mea_2_estimate, rate_fft_ini_2, thresh_fft_ini_2, "h_mea_2_estimate_filtered", "Second measurement histogram filtered");
	h_mea_2_estimate_filtered->GetXaxis()->SetTitle("Channel");
	h_mea_2_estimate_filtered->GetYaxis()->SetTitle("Count");
	TH1D* h_mea_2_estimate_diff = Diff(h_mea_2_estimate_filtered, "h_mea_2_estimate_diff", "Derivative of Second measurement histogram");
	
	//! Show derivatives on a canvas
	TCanvas* c_derivative = new TCanvas("c_derivative", "Derivative check", 1000, 500);
	c_derivative->Divide(2,2);
	c_derivative->cd(1);
	h_mea_1_estimate_filtered->Draw();
	c_derivative->cd(3);
	h_mea_1_estimate_diff->Draw();
	h_mea_1_estimate_diff->GetXaxis()->SetRangeUser(x_min_fft_1, x_max_fft_1);
	int ch1 = h_mea_1_estimate_diff->GetMinimumBin() + x_min_mea;
	std::cout << "ch1 = " << ch1 << "\n";

	c_derivative->cd(2);
	h_mea_2_estimate_filtered->Draw();
	c_derivative->cd(4);
	h_mea_2_estimate_diff->Draw();
	h_mea_2_estimate_diff->GetXaxis()->SetRangeUser(x_min_fft_2, x_max_fft_2);
	int ch2 = h_mea_2_estimate_diff->GetMinimumBin() + x_min_mea;
	std::cout << "ch2 = " << ch2 << "\n";

	c_derivative->Modified();
	c_derivative->Update();
	//! cShow2 is c_derivative

  //! Canvas to compare measurement and simulation spectrum while calibrating
	TCanvas* c_compare = new TCanvas("c_compare", "Compare Spectrum", 1800, 1000);
	c_compare->Divide(2,2);
	c_compare->cd(1)->SetLeftMargin(0.15);
	c_compare->cd(3)->SetLeftMargin(0.15);
	c_compare->cd(2)->SetLeftMargin(0.15);
	c_compare->cd(4)->SetLeftMargin(0.15);

  //! Chi2 value to estimate the goodness of fit
	double chi2_1, chi2_2;
	TGraph* gr_chi2_1 = new TGraph();
	TGraph* gr_chi2_2 = new TGraph();

  //! delta_chi2 value to check the convergence of Chi2 value
	double delta_chi2 = 1000000.;
	double delta_chi2_1, delta_chi2_2, delta_chi2_3 = 0.;

  //! threshold for acceptance for changing of delta_chi2 for all spectrum 
	double thresh_delta_chi2 = 0.3;
  //! threshold of acceptance for changing of delta_chi2 for each spectrum
	double thresh_delta_chi2_each = 0.5;
	int iteration = 1;

	while (delta_chi2>thresh_delta_chi2){
    std::cout << "\nIteration " << iteration << ":\n";
  
    //! Energy calibration coefficients
    double a = fit_energy(ch1, ch2, c_fit, 0);
    double b = fit_energy(ch1, ch2, c_fit, 1);
  
    double a_up_ch1 = fit_energy(ch1+step_channel_1, ch2, c_fit, 0);
    double b_up_ch1 = fit_energy(ch1+step_channel_1, ch2, c_fit, 1);
  
    double a_up_ch2 = fit_energy(ch1, ch2+step_channel_2, c_fit, 0);
    double b_up_ch2 = fit_energy(ch1, ch2+step_channel_2, c_fit, 1);
  
    //! Boundaries
    double x_min_cal = x_min_mea*a+b;
    double x_max_cal = x_max_mea*a+b;
  
    double x_min_cal_up_ch1 = x_min_mea*a_up_ch1+b_up_ch1;
    double x_max_cal_up_ch1 = x_max_mea*a_up_ch1+b_up_ch1;
  
    double x_min_cal_up_ch2 = x_min_mea*a_up_ch2+b_up_ch2;
    double x_max_cal_up_ch2 = x_max_mea*a_up_ch2+b_up_ch2;
  
    //! Histograms (1st file)
    TH1D* h_sim_1_res = new TH1D("h_sim_1_res", "Simulation 1 with resolution", bin_mea, x_min_cal, x_max_cal);
    TH1D* h_cal_1 = new TH1D("h_cal_1", "Measurement 1 Calibrated", bin_mea, x_min_cal, x_max_cal);
  
    TH1D* h_sim_1_res_up_ch1 = new TH1D("h_sim_1_res_up_ch1", "Simulation 1 with resolution", bin_mea, x_min_cal_up_ch1, x_max_cal_up_ch1);
    TH1D* h_cal_1_up_ch1 = new TH1D("h_cal_1_up_ch1", "Measurement 1 Calibrated", bin_mea, x_min_cal_up_ch1, x_max_cal_up_ch1);
  
    TH1D* h_sim_1_res_up_ch2 = new TH1D("h_sim_1_res_up_ch2", "Simulation 1 with resolution", bin_mea, x_min_cal_up_ch2, x_max_cal_up_ch2);
    TH1D* h_cal_1_up_ch2 = new TH1D("h_cal_1_up_ch2", "Measurement 1 Calibrated", bin_mea, x_min_cal_up_ch2, x_max_cal_up_ch2);
  
    TH1D* h_sim_1_res_up_sig1 = new TH1D("h_sim_1_res_up_sig1", "Simulation 1 with resolution", bin_mea, x_min_cal, x_max_cal);
  
    TH1D* h_sim_1_res_up_sig2 = new TH1D("h_sim_1_res_up_sig2", "Simulation 1 with resolution", bin_mea, x_min_cal, x_max_cal);
  
    //! Histograms (2nd file)
    TH1D* h_sim_2_res = new TH1D("h_sim_2_res", "Na22 Simulation with resolution", bin_mea, x_min_cal, x_max_cal);
    TH1D* h_cal_2 = new TH1D("h_cal_2", "Na22 Experiment Calibrated", bin_mea, x_min_cal, x_max_cal);
  
    TH1D* h_sim_2_res_up_ch1 = new TH1D("h_sim_2_res_up_ch1", "Na22 Simulation with resolution", bin_mea, x_min_cal_up_ch1, x_max_cal_up_ch1);
    TH1D* h_cal_2_up_ch1 = new TH1D("h_cal_2_up_ch1", "Na22 Experiment Calibrated", bin_mea, x_min_cal_up_ch1, x_max_cal_up_ch1);
  
    TH1D* h_sim_2_res_up_ch2 = new TH1D("h_sim_2_res_up_ch2", "Na22 Simulation with resolution", bin_mea, x_min_cal_up_ch2, x_max_cal_up_ch2);
    TH1D* h_cal_2_up_ch2 = new TH1D("h_cal_2_up_ch2", "Na22 Experiment Calibrated", bin_mea, x_min_cal_up_ch2, x_max_cal_up_ch2);
  
    TH1D* h_sim_2_res_up_sig1 = new TH1D("h_sim_2_res_up_sig1", "Na22 Simulation with resolution", bin_mea, x_min_cal, x_max_cal);
  
    TH1D* h_sim_2_res_up_sig2 = new TH1D("h_sim_2_res_up_sig2", "Na22 Simulation with resolution", bin_mea, x_min_cal, x_max_cal);
  
    //! Fill experiment histograms (1st file)
    for(int i = 0; i < entries_mea_descent; i++){
      t_mea_1->GetEntry(i);
      h_cal_1->Fill(a*(event_mea_1[channel]+0.5) + b);
      h_cal_1_up_ch1->Fill(a_up_ch1*(event_mea_1[channel]+0.5) + b_up_ch1);
      h_cal_1_up_ch2->Fill(a_up_ch2*(event_mea_1[channel]+0.5) + b_up_ch2);
    }
    TH1D* h_cal_1_filtered = fft(h_cal_1, 0.1, thresh_fft_descent_mea, "h_cal_1_filtered", "h_cal_1_filtered");
    TH1D* h_cal_1_up_ch1_filtered = fft(h_cal_1_up_ch1, 0.1, thresh_fft_descent_mea, "h_cal_1_up_ch1_filtered", "h_cal_1_up_ch1_filtered");
    TH1D* h_cal_1_up_ch2_filtered = fft(h_cal_1_up_ch2, 0.1, thresh_fft_descent_mea, "h_cal_1_up_ch2_filtered", "h_cal_1_up_ch2_filtered");
  
    //! Fill experiment histograms (2nd file)
    for(int i = 0; i < entries_mea_descent; i++){
      t_mea_2->GetEntry(i);
      h_cal_2->Fill(a*(event_mea_2[channel]+0.5) + b);
      h_cal_2_up_ch1->Fill(a_up_ch1*(event_mea_2[channel]+0.5) + b_up_ch1);
      h_cal_2_up_ch2->Fill(a_up_ch2*(event_mea_2[channel]+0.5) + b_up_ch2);
    }
    TH1D* h_cal_2_filtered = fft(h_cal_2, 0.1, thresh_fft_descent_mea, "h_cal_2_filtered", "h_cal_2_filtered");
    TH1D* h_cal_2_up_ch1_filtered = fft(h_cal_2_up_ch1, 0.1, thresh_fft_descent_mea, "h_cal_2_up_ch1_filtered", "h_cal_2_up_ch1_filtered");
    TH1D* h_cal_2_up_ch2_filtered = fft(h_cal_2_up_ch2, 0.1, thresh_fft_descent_mea, "h_cal_2_up_ch2_filtered", "h_cal_2_up_ch2_filtered");
  
    //! Prepare Energy Resolution coefficients
    double a_res = fit_energy_res(res_1, res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 0);
    double b_res = fit_energy_res(res_1, res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 1);
    double c_res = fit_energy_res(res_1, res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 2);
  
    double a_res_up_sig1 = fit_energy_res(res_1+step_res_1, res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 0);
    double b_res_up_sig1 = fit_energy_res(res_1+step_res_1, res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 1);
    double c_res_up_sig1 = fit_energy_res(res_1+step_res_1, res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 2);
  
    double a_res_up_sig2 = fit_energy_res(res_1, res_2+step_res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 0);
    double b_res_up_sig2 = fit_energy_res(res_1, res_2+step_res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 1);
    double c_res_up_sig2 = fit_energy_res(res_1, res_2+step_res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 2);
  
    TRandom3* ranGen = new TRandom3();
    //! Fill Simulation histograms (1st file)
    for (int i = 0; i < entries_sim_descent; i++)
    {
      t_sim_1->GetEntry(i);
      double sigma = event_sim_1*sqrt( pow(a_res,2) + pow(b_res/sqrt(event_sim_1),2) + pow(c_res/event_sim_1,2) );
      double energy = ranGen->Gaus(event_sim_1,sigma);
      h_sim_1_res->Fill(energy);
      h_sim_1_res_up_ch1->Fill(energy);
      h_sim_1_res_up_ch2->Fill(energy);
    }
    TH1D* h_sim_1_res_filtered = fft(h_sim_1_res, 0.1, thresh_fft_descent_sim, "h_sim_1_res_filtered", "h_sim_1_res_filtered");
    TH1D* h_sim_1_res_up_ch1_filtered = fft(h_sim_1_res_up_ch1, 0.1, thresh_fft_descent_sim, "h_sim_1_res_up_ch1_filtered", "h_sim_1_res_up_ch1_filtered");
    TH1D* h_sim_1_res_up_ch2_filtered = fft(h_sim_1_res_up_ch2, 0.1, thresh_fft_descent_sim, "h_sim_1_res_up_ch2_filtered", "h_sim_1_res_up_ch2_filtered");    for (int i = 0; i < entries_sim_descent; i++)
    {
      t_sim_1->GetEntry(i);
      double sigma = event_sim_1*sqrt( pow(a_res_up_sig1,2) + pow(b_res_up_sig1/sqrt(event_sim_1),2) + pow(c_res_up_sig1/event_sim_1,2) );
      h_sim_1_res_up_sig1->Fill(ranGen->Gaus(event_sim_1,sigma));
    }
    TH1D* h_sim_1_res_up_sig1_filtered = fft(h_sim_1_res_up_sig1, 0.1, thresh_fft_descent_sim, "h_sim_1_res_up_sig1_filtered", "h_sim_1_res_up_sig1_filtered");
  
    for (int i = 0; i < entries_sim_descent; i++)
    {
      t_sim_1->GetEntry(i);
      double sigma = event_sim_1*sqrt( pow(a_res_up_sig2,2) + pow(b_res_up_sig2/sqrt(event_sim_1),2) + pow((c_res_up_sig2)/event_sim_1,2) );
      h_sim_1_res_up_sig2->Fill(ranGen->Gaus(event_sim_1,sigma));
    }
    TH1D* h_sim_1_res_up_sig2_filtered = fft(h_sim_1_res_up_sig2, 0.1, thresh_fft_descent_sim, "h_sim_1_res_up_sig2_filtered", "h_sim_1_res_up_sig2_filtered");
  
    h_cal_1_filtered->GetXaxis()->SetRangeUser(x_min_descent_1, x_max_descent_1);
    h_sim_1_res_filtered->GetXaxis()->SetRangeUser(x_min_descent_1, x_max_descent_1);
  
    //! Fill Simulation histograms (2nd file)
    for (int i = 0; i < entries_sim_descent; i++)
    {
      t_sim_2->GetEntry(i);
      double sigma = event_sim_2*sqrt( pow(a_res,2) + pow(b_res/sqrt(event_sim_2),2) + pow(c_res/event_sim_2,2) );
      double energy = ranGen->Gaus(event_sim_2,sigma);
      h_sim_2_res->Fill(energy);
      h_sim_2_res_up_ch1->Fill(energy);
      h_sim_2_res_up_ch2->Fill(energy);
    }
    TH1D* h_sim_2_res_filtered = fft(h_sim_2_res, 0.1, thresh_fft_descent_sim, "h_sim_2_res_filtered", "h_sim_2_res_filtered");
    TH1D* h_sim_2_res_up_ch1_filtered = fft(h_sim_2_res_up_ch1, 0.1, thresh_fft_descent_sim, "h_sim_2_res_up_ch1_filtered", "h_sim_2_res_up_ch1_filtered");
    TH1D* h_sim_2_res_up_ch2_filtered = fft(h_sim_2_res_up_ch2, 0.1, thresh_fft_descent_sim, "h_sim_2_res_up_ch2_filtered", "h_sim_2_res_up_ch2_filtered");
  
    for (int i = 0; i < entries_sim_descent; i++)
    {
      t_sim_2->GetEntry(i);
      double sigma = event_sim_2*sqrt( pow(a_res_up_sig1,2) + pow(b_res_up_sig1/sqrt(event_sim_2),2) + pow(c_res_up_sig1/event_sim_2,2) );
      h_sim_2_res_up_sig1->Fill(ranGen->Gaus(event_sim_2,sigma));
    }
    TH1D* h_sim_2_res_up_sig1_filtered = fft(h_sim_2_res_up_sig1, 0.1, thresh_fft_descent_sim, "h_sim_2_res_up_sig1_filtered", "h_sim_2_res_up_sig1_filtered");
  
    for (int i = 0; i < entries_sim_descent; i++)
    {
      t_sim_2->GetEntry(i);
      double sigma = event_sim_2*sqrt( pow(a_res_up_sig2,2) + pow(b_res_up_sig2/sqrt(event_sim_2),2) + pow((c_res_up_sig2)/event_sim_2,2) );
      h_sim_2_res_up_sig2->Fill(ranGen->Gaus(event_sim_2,sigma));
    }
    TH1D* h_sim_2_res_up_sig2_filtered = fft(h_sim_2_res_up_sig2, 0.1, thresh_fft_descent_sim, "h_sim_2_res_up_sig2_filtered", "h_sim_2_res_up_sig2_filtered");
  
    h_cal_2_filtered->GetXaxis()->SetRangeUser(x_min_descent_2, x_max_descent_2);
    h_sim_2_res_filtered->GetXaxis()->SetRangeUser(x_min_descent_2, x_max_descent_2);
  
    //! Delta Chi2
    if(iteration == 2){delta_chi2_1 = abs(chi2(h_cal_1_filtered, h_sim_1_res_filtered) - chi2_1) 
    + abs(chi2(h_cal_2_filtered, h_sim_2_res_filtered) -chi2_2);}
    else if(iteration == 3){delta_chi2_2 = abs(chi2(h_cal_1_filtered, h_sim_1_res_filtered) - chi2_1) 
    + abs(chi2(h_cal_2_filtered, h_sim_2_res_filtered) -chi2_2);}
    else if(iteration == 4){delta_chi2_3 = abs(chi2(h_cal_1_filtered, h_sim_1_res_filtered) - chi2_1) 
    + abs(chi2(h_cal_2_filtered, h_sim_2_res_filtered) -chi2_2);
    delta_chi2 = delta_chi2_1 + delta_chi2_2 + delta_chi2_3;}
    else if (iteration > 4){
    delta_chi2_1 = delta_chi2_2;
    delta_chi2_2 = delta_chi2_3;
    delta_chi2_3 = abs(chi2(h_cal_1_filtered, h_sim_1_res_filtered) - chi2_1) 
    + abs(chi2(h_cal_2_filtered, h_sim_2_res_filtered) -chi2_2);
    delta_chi2 = delta_chi2_1 + delta_chi2_2 + delta_chi2_3;
    }
  
    chi2_1 = chi2(h_cal_1_filtered, h_sim_1_res_filtered);
    chi2_2 = chi2(h_cal_2_filtered, h_sim_2_res_filtered);
  
    gr_chi2_1->SetTitle("#chi^{2} distance by iteration for the 1^{st} file");
    c_compare->cd(1)->SetTicks();
    c_compare->cd(1)->SetGrid();
    gr_chi2_1->GetXaxis()->SetTitle("Iteration");
    gr_chi2_1->GetYaxis()->SetTitle("#chi^{2} distance");
  
    gr_chi2_1->SetStats(0);
    gr_chi2_1->SetLineColor(kBlack);
    gr_chi2_1->SetLineWidth(3);
  
    gr_chi2_1->GetXaxis()->SetLabelFont(42);
    gr_chi2_1->GetXaxis()->SetTitleFont(52);
    gr_chi2_1->GetXaxis()->SetTitleSize(0.04);
    gr_chi2_1->GetXaxis()->CenterTitle(true);
  
    gr_chi2_1->GetYaxis()->SetLabelFont(42);
    gr_chi2_1->GetYaxis()->SetTitleFont(52);
    gr_chi2_1->GetYaxis()->SetTitleSize(0.04);
    gr_chi2_1->GetYaxis()->CenterTitle(true);
  
    gr_chi2_1->AddPoint(iteration, chi2_1);
  
    gr_chi2_2->SetTitle("#chi^{2} distance by iteration for the 2^{nd} file");
    c_compare->cd(2)->SetTicks();
    c_compare->cd(2)->SetGrid();
    gr_chi2_2->GetXaxis()->SetTitle("Iteration");
    gr_chi2_2->GetYaxis()->SetTitle("#chi^{2} distance");
  
    gr_chi2_2->SetStats(0);
    gr_chi2_2->SetLineColor(kBlack);
    gr_chi2_2->SetLineWidth(3);
  
    gr_chi2_2->GetXaxis()->SetLabelFont(42);
    gr_chi2_2->GetXaxis()->SetTitleFont(52);
    gr_chi2_2->GetXaxis()->SetTitleSize(0.04);
    gr_chi2_2->GetXaxis()->CenterTitle(true);
  
    gr_chi2_2->GetYaxis()->SetLabelFont(42);
    gr_chi2_2->GetYaxis()->SetTitleFont(52);
    gr_chi2_2->GetYaxis()->SetTitleSize(0.04);
    gr_chi2_2->GetYaxis()->CenterTitle(true);
  
    gr_chi2_2->AddPoint(iteration, chi2_2);
  
    c_compare->cd(1);
    gr_chi2_1->Draw();
    c_compare->cd(1)->Modified();
    c_compare->cd(1)->Update();
  
    c_compare->cd(2);
    gr_chi2_2->Draw();
    c_compare->cd(2)->Modified();
    c_compare->cd(2)->Update();
  
    //! Chi2(1st file)
    h_cal_1_up_ch1_filtered->GetXaxis()->SetRangeUser(x_min_descent_1, x_max_descent_1);
    h_sim_1_res_up_ch1_filtered->GetXaxis()->SetRangeUser(x_min_descent_1, x_max_descent_1);
    double chi2_1_up_ch1 = chi2(h_cal_1_up_ch1_filtered, h_sim_1_res_up_ch1_filtered);
  
    h_cal_1_up_ch2_filtered->GetXaxis()->SetRangeUser(x_min_descent_1, x_max_descent_1);
    h_sim_1_res_up_ch2_filtered->GetXaxis()->SetRangeUser(x_min_descent_1, x_max_descent_1);
    double chi2_1_up_ch2 = chi2(h_cal_1_up_ch2_filtered, h_sim_1_res_up_ch2_filtered);
  
    h_sim_1_res_up_sig1_filtered->GetXaxis()->SetRangeUser(x_min_descent_1, x_max_descent_1);
    double chi2_1_up_sig1 = chi2(h_cal_1_filtered, h_sim_1_res_up_sig1_filtered);
  
    h_sim_1_res_up_sig2_filtered->GetXaxis()->SetRangeUser(x_min_descent_1, x_max_descent_1);
    double chi2_1_up_sig2 = chi2(h_cal_1_filtered, h_sim_1_res_up_sig2_filtered);
  
    //! Chi2(2nd file)
    h_cal_2_up_ch1_filtered->GetXaxis()->SetRangeUser(x_min_descent_2, x_max_descent_2);
    h_sim_2_res_up_ch1_filtered->GetXaxis()->SetRangeUser(x_min_descent_2, x_max_descent_2);
    double chi2_2_up_ch1 = chi2(h_cal_2_up_ch1_filtered, h_sim_2_res_up_ch1_filtered);
  
    h_cal_2_up_ch2_filtered->GetXaxis()->SetRangeUser(x_min_descent_2, x_max_descent_2);
    h_sim_2_res_up_ch2_filtered->GetXaxis()->SetRangeUser(x_min_descent_2, x_max_descent_2);
    double chi2_2_up_ch2 = chi2(h_cal_2_up_ch2_filtered, h_sim_2_res_up_ch2_filtered);
  
    h_sim_2_res_up_sig1_filtered->GetXaxis()->SetRangeUser(x_min_descent_2, x_max_descent_2);
    double chi2_2_up_sig1 = chi2(h_cal_2_filtered, h_sim_2_res_up_sig1_filtered);
  
    h_sim_2_res_up_sig2_filtered->GetXaxis()->SetRangeUser(x_min_descent_2, x_max_descent_2);
    double chi2_2_up_sig2 = chi2(h_cal_2_filtered, h_sim_2_res_up_sig2_filtered);
  
    //! Derivative
    double deri_chi2_up_ch1 = ((chi2_1_up_ch1+chi2_2_up_ch1)-(chi2_1+chi2_2))/step_channel_1;
    double deri_chi2_up_ch2 = ((chi2_1_up_ch2+chi2_2_up_ch2)-(chi2_1+chi2_2))/step_channel_2;
    double deri_chi2_up_sig1 = ((chi2_1_up_sig1+chi2_2_up_sig1)-(chi2_1+chi2_2))/step_res_1;
    double deri_chi2_up_sig2 = ((chi2_1_up_sig2+chi2_2_up_sig2)-(chi2_1+chi2_2))/step_res_2;
  
    // std::cout << "~~~~~~~~~~~CHECK~~~~~~~~~~~~\n";
    // std::cout << "deri_chi2_up_ch1 = " << deri_chi2_up_ch1 << "; deri_chi2_up_ch2 = " << deri_chi2_up_ch2 
    // << "\nderi_chi2_up_sig1 = " << deri_chi2_up_sig1 << "; deri_chi2_up_sig2 = " << deri_chi2_up_sig2 << "\n";
    // std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

    //! Moving along the gradient
    ch1 = ch1 - learning_rate_channel_1*deri_chi2_up_ch1;
    ch2 = ch2 - learning_rate_channel_2*deri_chi2_up_ch2;
  
    res_1 = res_1 - learning_rate_res_1*deri_chi2_up_sig1;
    res_2 = res_2 - learning_rate_res_2*deri_chi2_up_sig2;

    std::cout << "ch1 = " << ch1 << "; ch2 = " << ch2 
    << "\nSig1 = " << res_1 << "; Sig2 = " << res_2 << "\n";
    
    // std::cout << "\ndelta_Chi2_1 = " << delta_chi2_1 
    // << "\ndelta_Chi2_2 = " << delta_chi2_2
    // << "\ndelta_Chi2_3 = " << delta_chi2_3
    // << "\nTotal delta_chi2 = " << delta_chi2 << "\n";
  
    c_compare->cd(3);
    c_compare->cd(3)->SetTicks();
    h_cal_1_filtered->SetTitle("");
    h_cal_1_filtered->GetXaxis()->SetTitle("E(MeVee)");
    h_cal_1_filtered->GetYaxis()->SetTitle("Events");
  
    h_cal_1_filtered->SetStats(0);
    h_cal_1_filtered->SetLineWidth(2);
    h_sim_1_res_filtered->SetLineWidth(2);
  
    h_cal_1_filtered->GetXaxis()->SetLabelFont(42);
    h_cal_1_filtered->GetXaxis()->SetTitleFont(52);
    h_cal_1_filtered->GetXaxis()->SetTitleSize(0.04);
    h_cal_1_filtered->GetXaxis()->CenterTitle(true);
  
    h_cal_1_filtered->GetYaxis()->SetLabelFont(42);
    h_cal_1_filtered->GetYaxis()->SetTitleFont(52);
    h_cal_1_filtered->GetYaxis()->SetTitleSize(0.04);
    h_cal_1_filtered->GetYaxis()->CenterTitle(true);
  
    h_cal_1_filtered->GetXaxis()->UnZoom();
    h_sim_1_res_filtered->GetXaxis()->UnZoom();
    h_cal_1_filtered->SetLineColor(kRed);
    h_cal_1_filtered->Draw();
    h_sim_1_res_filtered->Scale(scale_sim_1, "noSW2");
    h_sim_1_res_filtered->Draw("same");
  
    TLegend *legend = new TLegend(0.4, 0.55, 0.8, 0.85);
    legend->SetBorderSize(0);
    legend->SetLineWidth(2);
    legend->SetTextSize(0.06);
    legend->AddEntry(h_cal_1_filtered, "1^{st} measurement file", "l");
    legend->AddEntry(h_sim_1_res_filtered, "1^{st} simulation file", "l");
  
    legend->Draw();
  
    c_compare->cd(4);
    c_compare->cd(4)->SetTicks();
    h_cal_2_filtered->SetTitle("");
    h_cal_2_filtered->GetXaxis()->SetTitle("E(MeVee)");
    h_cal_2_filtered->GetYaxis()->SetTitle("Events");
  
    h_cal_2_filtered->SetStats(0);
    h_cal_2_filtered->SetLineWidth(2);
    h_sim_2_res_filtered->SetLineWidth(2);
  
    h_cal_2_filtered->GetXaxis()->SetLabelFont(42);
    h_cal_2_filtered->GetXaxis()->SetTitleFont(52);
    h_cal_2_filtered->GetXaxis()->SetTitleSize(0.04);
    h_cal_2_filtered->GetXaxis()->CenterTitle(true);
  
    h_cal_2_filtered->GetYaxis()->SetLabelFont(42);
    h_cal_2_filtered->GetYaxis()->SetTitleFont(52);
    h_cal_2_filtered->GetYaxis()->SetTitleSize(0.04);
    h_cal_2_filtered->GetYaxis()->CenterTitle(true);
  
    h_cal_2_filtered->GetXaxis()->UnZoom();
    h_sim_2_res_filtered->GetXaxis()->UnZoom();
    h_cal_2_filtered->SetLineColor(kRed);
    h_cal_2_filtered->Draw();
    h_sim_2_res_filtered->Scale(scale_sim_2, "noSW2");
    h_sim_2_res_filtered->Draw("same");
  
    TLegend *legend1 = new TLegend(0.4, 0.55, 0.8, 0.85);
    legend1->SetBorderSize(0);
    legend1->SetLineWidth(10);
    legend1->SetTextSize(0.06);
    legend1->AddEntry(h_cal_2_filtered, "2^{nd} measurement file", "l");
    legend1->AddEntry(h_sim_2_res_filtered, "2^{nd} simulation file", "l");
  
    legend1->Draw();
  
    c_compare->cd(3)->Modified();
    c_compare->cd(3)->Update();
    c_compare->cd(4)->Modified();
    c_compare->cd(4)->Update();
  
    c_derivative->cd(1);    
    h_cal_1_up_ch1_filtered->SetLineColor(kRed);
    h_cal_1_up_ch1_filtered->Draw();
    h_sim_1_res_up_ch1_filtered->Scale(scale_sim_1, "noSW2");
    h_sim_1_res_up_ch1_filtered->Draw("same");
    c_derivative->cd(1)->Modified();
    c_derivative->cd(1)->Update();
  
    c_derivative->cd(2);    
    h_cal_1_up_ch2_filtered->SetLineColor(kRed);
    h_cal_1_up_ch2_filtered->Draw();
    h_sim_1_res_up_ch2_filtered->Scale(scale_sim_1, "noSW2");
    h_sim_1_res_up_ch2_filtered->Draw("same");
    c_derivative->cd(2)->Modified();
    c_derivative->cd(2)->Update();
  
    c_derivative->cd(3);
    h_cal_1_filtered->Draw();
    h_sim_1_res_up_sig1_filtered->Scale(scale_sim_1, "noSW2");
    h_sim_1_res_up_sig1_filtered->Draw("same");
    c_derivative->cd(3)->Modified();
    c_derivative->cd(3)->Update();
  
    c_derivative->cd(4);
    h_cal_1_filtered->Draw();
    h_sim_1_res_up_sig2_filtered->Scale(scale_sim_1, "noSW2");
    h_sim_1_res_up_sig2_filtered->Draw("same");
    c_derivative->cd(4)->Modified();
    c_derivative->cd(4)->Update();
  
    c_derivative->cd(5);    
    h_cal_2_up_ch1_filtered->SetLineColor(kRed);
    h_cal_2_up_ch1_filtered->Draw();
    h_sim_2_res_up_ch1_filtered->Scale(scale_sim_2, "noSW2");
    h_sim_2_res_up_ch1_filtered->Draw("same");
    c_derivative->cd(5)->Modified();
    c_derivative->cd(5)->Update();
  
    c_derivative->cd(6);    
    h_cal_2_up_ch2_filtered->SetLineColor(kRed);
    h_cal_2_up_ch2_filtered->Draw();
    h_sim_2_res_up_ch2_filtered->Scale(scale_sim_2, "noSW2");
    h_sim_2_res_up_ch2_filtered->Draw("same");
    c_derivative->cd(6)->Modified();
    c_derivative->cd(6)->Update();
  
    c_derivative->cd(7);
    h_cal_2_filtered->Draw();
    h_sim_2_res_up_sig1_filtered->Scale(scale_sim_2, "noSW2");
    h_sim_2_res_up_sig1_filtered->Draw("same");
    c_derivative->cd(7)->Modified();
    c_derivative->cd(7)->Update();
  
    c_derivative->cd(8);
    h_cal_2_filtered->Draw();
    h_sim_2_res_up_sig2_filtered->Scale(scale_sim_2, "noSW2");
    h_sim_2_res_up_sig2_filtered->Draw("same");
    c_derivative->cd(8)->Modified();
    c_derivative->cd(8)->Update();
  
    if(iteration > 1){
      if (chi2_1 < thresh_delta_chi2_each){
      learning_rate_channel_1 = 0.;
      learning_rate_res_1 = 0.;
      }
  
      if (chi2_2 < thresh_delta_chi2_each){
      learning_rate_channel_2 = 0.;
      learning_rate_res_2 = 0.;
      }
    }
  
    if(delta_chi2 > thresh_delta_chi2){
      delete h_sim_1_res;
      delete h_cal_1;
      delete h_sim_1_res_up_ch1;
      delete h_cal_1_up_ch1;
      delete h_sim_1_res_up_ch2;
      delete h_cal_1_up_ch2;
      delete h_sim_1_res_up_sig1;
      delete h_sim_1_res_up_sig2;
  
      delete h_sim_1_res_filtered;
      delete h_cal_1_filtered;
      delete h_sim_1_res_up_ch1_filtered;
      delete h_cal_1_up_ch1_filtered;
      delete h_sim_1_res_up_ch2_filtered;
      delete h_cal_1_up_ch2_filtered;
      delete h_sim_1_res_up_sig1_filtered;
      delete h_sim_1_res_up_sig2_filtered;
  
      delete h_sim_2_res;
      delete h_cal_2;
      delete h_sim_2_res_up_ch1;
      delete h_cal_2_up_ch1;
      delete h_sim_2_res_up_ch2;
      delete h_cal_2_up_ch2;
      delete h_sim_2_res_up_sig1;
      delete h_sim_2_res_up_sig2;
  
      delete h_sim_2_res_filtered;
      delete h_cal_2_filtered;
      delete h_sim_2_res_up_ch1_filtered;
      delete h_cal_2_up_ch1_filtered;
      delete h_sim_2_res_up_ch2_filtered;
      delete h_cal_2_up_ch2_filtered;
      delete h_sim_2_res_up_sig1_filtered;
      delete h_sim_2_res_up_sig2_filtered;
  
    }
    else{}
  
    iteration++;
    }

}

int main(int argc, char **argv)
{
    if (argc != 2) {
        std::cerr << "Error: Expected exactly one argument (macro file) or --help.\n";
        return 1;
    }

    // If the user asks for help, display usage info and exit
    if (std::string(argv[1]) == "--help") {
        std::cout << "Usage: ./demo macro.mac or ./demo --help\n"
                  << "  macro.mac: A script file containing arguments\n"
                  << "  --help   : Show this help message\n";
        return 0;
    }

    // Otherwise, assume the argument is a macro file and process it
    std::vector<std::string> args = ReadMacroFile(argv[1]);

    // Convert vector of strings to array of C-style strings
    std::vector<char*> c_args;
    c_args.push_back(argv[0]); // Program name
    for (std::string& arg : args) {
        c_args.push_back(&arg[0]);
    }
    c_args.push_back(nullptr); // Null termination

    int new_argc = c_args.size() - 1;
    char** new_argv = c_args.data();

    if (ReadArgs(new_argc, new_argv) != 0) return 1;

    // Initialize application
    TApplication app("app", &argc, argv);

	std::cerr << "\nMain is running...\n";

	// Get Files and Trees
	GetTrees();
	
    // Create canvas and Handle window close event
    //TCanvas* c = new TCanvas("c", "Plot", 0, 0, 800, 600);

    //TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
    //rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	
    // Run the application
    std::cerr << "\nRun success.\n";
    app.Run();

    return 0;
}

