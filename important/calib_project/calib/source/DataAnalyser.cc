#include "DataAnalyser.hh"

double final_energy_calibration_coef_a;
double final_energy_calibration_coef_b;

double final_energy_calibration_coef_a_error;
double final_energy_calibration_coef_b_error;

double final_energy_resolution_coef_a;
double final_energy_resolution_coef_b;
double final_energy_resolution_coef_c;

double final_energy_resolution_coef_a_error;
double final_energy_resolution_coef_b_error;
double final_energy_resolution_coef_c_error;

// Function that performs fast fourier transformation (FFT) on a histogram h_channel.
TH1D* DataAnalyser::fft(TH1D* h_channel, double para_k, double para_c, std::string name_h_channel, std::string title_h_channel)
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
TH1D* DataAnalyser::Diff(TH1D* h_channel, std::string name_h_channel, std::string title_h_channel)
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
	double chi2_result = h_channel_1->Chi2Test(h_channel_2, "CHI2/NDF");
	return chi2_result;
}

// Function that performs linear calibration of energy and returns coefficients.
double fit_energy(int channel_1, int channel_2, TCanvas* c_fit, int i, const Config& config){
	TH1D* h_calibrate = new TH1D("h_calibrate", "Calibration fit", config.bin_mea, config.x_min_mea, config.x_max_mea);
	h_calibrate->SetBinContent(channel_1 - config.x_min_mea, config.energy_1);
	h_calibrate->SetBinContent(channel_2 - config.x_min_mea, config.energy_2);

	TF1* f_linear = new TF1("f_linear", "[0]*x + [1]", config.x_min_mea, config.x_max_mea);
	c_fit->cd(1);
	h_calibrate->Fit("f_linear", "Q");

	double coef_a = f_linear->GetParameter(0);
	double coef_b = f_linear->GetParameter(1);
	double coef_a_error = f_linear->GetParError(0);
	double coef_b_error = f_linear->GetParError(1);

	c_fit->Modified();
	c_fit->Update();

	delete h_calibrate;
	delete f_linear;

	if(i == 0){return coef_a;}
  else if(i == 1){return coef_b;}
  else if(i == 2){return coef_a_error;}
	else{return coef_b_error;}
}

// Function that performs calibration of energy resolution and returns coefficients.
double fit_energy_res(double res_sig_1, double res_sig_2, TH1D* h_1, TH1D* h_2, double x_min, double x_max, TCanvas* c_fit, int i, Config& config){
	TH1D* h_res = new TH1D("h_res", "Resolution fit", config.bin_mea, x_min, x_max);
	h_res->SetBinContent(h_1->FindBin(config.energy_1), res_sig_1/100.);
	h_res->SetBinContent(h_2->FindBin(config.energy_2), res_sig_2/100.);

	//! Choose fit function
	TF1* f_res = new TF1("f_res", "sqrt([0]*[0] + [1]*[1]/(x*x))");
	c_fit->cd(2);
	h_res->Fit("f_res", "Q", "", config.energy_1-0.1, config.energy_2+0.1);

	//! Correspond to chosen function
	double coef_a = std::abs(f_res->GetParameter(0));
	double coef_b = 0.;
	double coef_c = std::abs(f_res->GetParameter(1));
	double coef_a_error = std::abs(f_res->GetParError(0));
	double coef_b_error = 0.;
	double coef_c_error = std::abs(f_res->GetParError(1));

	c_fit->Modified();
	c_fit->Update();

	delete h_res;
	delete f_res;

	if(i == 0){return coef_a;}
	else if(i == 1){return coef_b;}
	else if(i == 2){return coef_c;}
	else if(i == 3){return coef_a_error;}
	else if(i == 4){return coef_b_error;}
	else{return coef_c_error;}
}

void DataAnalyser::Analyze(Config& config)
{
	std::cerr << "Sim1: " << config.name_f_sim_1 << "\n";
	std::cerr << "Sim2: " << config.name_f_sim_2 << "\n";
	std::cerr << "Mea1: " << config.name_f_mea_1 << "\n";
	std::cerr << "Mea2: " << config.name_f_mea_2 << "\n";
	
	//! Read the simulation and measurement files
	TFile* f_sim_1 = new TFile(config.name_f_sim_1.c_str(), "read");
	TFile* f_sim_2 = new TFile(config.name_f_sim_2.c_str(), "read");

  TCanvas* c_derivative = new TCanvas("c_derivative", "Derivative check", 1000, 500);
  c_derivative->Divide(2,2);

  int ch1;
  int ch2;

  //! Fill histograms for FFT to estimate Compton edge
  TH1D* h_mea_1_estimate = new TH1D();
  TH1D* h_mea_2_estimate = new TH1D();

  TFile* f_mea_1 = new TFile();
  TFile* f_mea_2 = new TFile();

  TTree* t_mea_1 = new TTree();
  TTree* t_mea_2 = new TTree();
  
  double event_mea_1;
  double event_mea_2;

  //! Variables used in case of ascii files
  int x1,x2;
  double y1,y2;
  std::vector<int> bin_input1, bin_input2;
  std::vector<double> content_input1, content_input2;

  if (config.ascii == 0){
    f_mea_1 = new TFile(config.name_f_mea_1.c_str(), "read");
    f_mea_2 = new TFile(config.name_f_mea_2.c_str(), "read");

    t_mea_1 = (TTree*) f_mea_1->Get("Events");
    t_mea_2 =  (TTree*) f_mea_2->Get("Events");

    
    t_mea_2->SetBranchAddress("Amplitude", &event_mea_2);
    t_mea_1->SetBranchAddress("Amplitude", &event_mea_1);

    h_mea_1_estimate = new TH1D("h_mea_1_estimate", "First measurement histogram", config.bin_mea, config.x_min_mea, config.x_max_mea);
    h_mea_2_estimate = new TH1D("h_mea_2_estimate", "Second measurement histogram", config.bin_mea, config.x_min_mea, config.x_max_mea);

    for(int i = 0; i < config.entries_mea_fft; i++){
      t_mea_1->GetEntry(i);
      h_mea_1_estimate->Fill(event_mea_1);

      t_mea_2->GetEntry(i);
      h_mea_2_estimate->Fill(event_mea_2);
    }
  }
  else{
    std::ifstream f_mea_1(config.name_f_mea_1.c_str());
    std::ifstream f_mea_2(config.name_f_mea_2.c_str());
    
    std::string line;
    while (std::getline(f_mea_1, line)) {
      if (!line.empty() && line[0] != '#') { // Ignore empty lines and comments
          std::istringstream iss(line);
          if (iss >> x1 >> y1) {  // Ensure valid data extraction
              bin_input1.push_back(x1);
              content_input1.push_back(y1);
          }
      }
    }
    
    while (std::getline(f_mea_2, line)) {
      if (!line.empty() && line[0] != '#') { // Ignore empty lines and comments
          std::istringstream iss(line);
          if (iss >> x2 >> y2) {  // Ensure valid data extraction
              bin_input2.push_back(x2);
              content_input2.push_back(y2);
          }
      }
    }

  //   //! change to read line by line
  //   while(f_mea_1 >> x1 >> y1){
  //     // std::cout << x1 << " " << y1 << "\n";
  //     //! if not #
  //     bin_input1.push_back(x1);
  //     content_input1.push_back(y1);
  //   }

    // while(f_mea_2 >> x2 >> y2){
    //   // std::cout << x2 << " " << y2 << "\n";
    //   bin_input2.push_back(x2);
    //   content_input2.push_back(y2);
    // }
  
    // std::cout << "bin_input.size() = " << bin_input.size() << "\n";
    // std::cout << "bin_input[0] = " << bin_input[0] << "\n";
    // std::cout << "bin_input[bin_input.size()-1] = " << bin_input[bin_input.size()-1] << "\n";

    config.bin_mea = bin_input1.size();
    config.bin_sim = config.bin_mea;
    config.x_min_mea = bin_input1[0];
    config.x_max_mea = bin_input1[bin_input1.size()-1];

    std::cout << "config.bin_mea = " << config.bin_mea << "\n";
    std::cout << "config.x_min_mea = " << config.x_min_mea << "\n";
    std::cout << "config.x_max_mea = " << config.x_max_mea << "\n";

    h_mea_1_estimate = new TH1D("h_mea_1_estimate", "First measurement histogram", config.bin_mea, config.x_min_mea, config.x_max_mea);
    h_mea_2_estimate = new TH1D("h_mea_2_estimate", "Second measurement histogram", config.bin_mea, config.x_min_mea, config.x_max_mea);
    for (int i=1; i<=config.bin_mea; i++){
      h_mea_1_estimate->SetBinContent(i, content_input1[i-1]);
      h_mea_2_estimate->SetBinContent(i, content_input2[i-1]);
      // std::cout << "i = " << i << "; content_input1[i-1] = " << content_input1[i-1] << "; content_input2[i-1] = " << content_input2[i-1] << "\n";
    }
  }

  TH1D* h_mea_1_estimate_filtered = fft(h_mea_1_estimate, config.rate_fft_ini_1, config.thresh_fft_ini_1,
  "h_mea_1_estimate_filtered", "First measurement histogram filtered");
  
  h_mea_1_estimate_filtered->GetXaxis()->SetTitle("Channel");
  h_mea_1_estimate_filtered->GetYaxis()->SetTitle("Count");
  TH1D* h_mea_1_estimate_diff = Diff(h_mea_1_estimate_filtered, "h_mea_1_estimate_diff", "Derivative of First measurement histogram");

  TH1D* h_mea_2_estimate_filtered = fft(h_mea_2_estimate, config.rate_fft_ini_2, config.thresh_fft_ini_2, "h_mea_2_estimate_filtered", "Second measurement histogram filtered");
  h_mea_2_estimate_filtered->GetXaxis()->SetTitle("Channel");
  h_mea_2_estimate_filtered->GetYaxis()->SetTitle("Count");
  TH1D* h_mea_2_estimate_diff = Diff(h_mea_2_estimate_filtered, "h_mea_2_estimate_diff", "Derivative of Second measurement histogram");

  //! Show derivatives on a canvas
  c_derivative->cd(1);
  h_mea_1_estimate_filtered->Draw();
  c_derivative->cd(3);
  h_mea_1_estimate_diff->Draw();
  h_mea_1_estimate_diff->GetXaxis()->SetRangeUser(config.x_min_fft_1, config.x_max_fft_1);
  ch1 = h_mea_1_estimate_diff->GetMinimumBin() + config.x_min_mea;
  std::cout << "ch1 = " << ch1 << "\n";

  c_derivative->cd(2);
  h_mea_2_estimate_filtered->Draw();
  c_derivative->cd(4);
  h_mea_2_estimate_diff->Draw();
  h_mea_2_estimate_diff->GetXaxis()->SetRangeUser(config.x_min_fft_2, config.x_max_fft_2);
  ch2 = h_mea_2_estimate_diff->GetMinimumBin() + config.x_min_mea;
  std::cout << "ch2 = " << ch2 << "\n";

  c_derivative->Modified();
  c_derivative->Update();

	TCanvas* c_fit = new TCanvas("c_fit", "c_fit", 1000, 500);
	c_fit->Divide(2,1);

	//! Select trees and branches from the files
	TTree* t_sim_1 = (TTree*) f_sim_1->Get("Events");
	TTree* t_sim_2 =  (TTree*) f_sim_2->Get("Events");

	double event_sim_1;
	double event_sim_2;

	t_sim_1->SetBranchAddress("Energy", &event_sim_1);
	t_sim_2->SetBranchAddress("Energy", &event_sim_2);

  //! Canvas to compare measurement and simulation spectrum while calibrating
	TCanvas* c_compare = new TCanvas("c_compare", "Compare Spectrum", 1200, 1000);
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
	double delta_chi2_1 = 0.;
	double delta_chi2_2 = 0.;
	double delta_chi2_3 = 0.;

  //! threshold for acceptance for changing of delta_chi2 for all spectrum 
	double thresh_delta_chi2 = 0.3;
  //! threshold of acceptance for changing of delta_chi2 for each spectrum
	double thresh_delta_chi2_each = 0.15;
	int iteration = 1;

  //! Option to use all the events in the file for the iterations
  int entries_mea_descent_used_1 = config.entries_mea_descent;
  int entries_mea_descent_used_2 = config.entries_mea_descent;  
  if (config.entries_mea_descent = 0){
    entries_mea_descent_used_1 = t_mea_1->GetEntries();
    entries_mea_descent_used_2 = t_mea_2->GetEntries();
    std::cerr << "All events are used for measurement file 1: " << entries_mea_descent_used_1 << "; and measurement file 2: " << entries_mea_descent_used_2 << "\n";
  }

  int entries_sim_descent_used_1 = config.entries_sim_descent;
  int entries_sim_descent_used_2 = config.entries_sim_descent;  
  if (config.entries_sim_descent = 0){
    entries_sim_descent_used_1 = t_sim_1->GetEntries();
    entries_sim_descent_used_2 = t_sim_2->GetEntries();
    std::cerr << "All events are used for simulation file 1: " << entries_sim_descent_used_1 << "; and simulation file 2: " << entries_sim_descent_used_2 << "\n";
  }

	while (delta_chi2>thresh_delta_chi2){
    std::cout << "\nIteration " << iteration << ":\n";
  
    //! Energy calibration coefficients
    double a = fit_energy(ch1, ch2, c_fit, 0, config);
    double b = fit_energy(ch1, ch2, c_fit, 1, config);
  
    double a_up_ch1 = fit_energy(ch1+config.step_channel_1, ch2, c_fit, 0, config);
    double b_up_ch1 = fit_energy(ch1+config.step_channel_1, ch2, c_fit, 1, config);
  
    double a_up_ch2 = fit_energy(ch1, ch2+config.step_channel_2, c_fit, 0, config);
    double b_up_ch2 = fit_energy(ch1, ch2+config.step_channel_2, c_fit, 1, config);
  
    //! Boundaries
    double x_min_cal = config.x_min_mea*a+b;
    double x_max_cal = config.x_max_mea*a+b;
  
    double x_min_cal_up_ch1 = config.x_min_mea*a_up_ch1+b_up_ch1;
    double x_max_cal_up_ch1 = config.x_max_mea*a_up_ch1+b_up_ch1;
  
    double x_min_cal_up_ch2 = config.x_min_mea*a_up_ch2+b_up_ch2;
    double x_max_cal_up_ch2 = config.x_max_mea*a_up_ch2+b_up_ch2;
  
    //! Histograms (1st file)
    TH1D* h_sim_1_res = new TH1D("h_sim_1_res", "Simulation 1 with resolution", config.bin_mea, x_min_cal, x_max_cal);
    TH1D* h_cal_1 = new TH1D("h_cal_1", "Measurement 1 Calibrated", config.bin_mea, x_min_cal, x_max_cal);
  
    TH1D* h_sim_1_res_up_ch1 = new TH1D("h_sim_1_res_up_ch1", "Simulation 1 with resolution", config.bin_mea, x_min_cal_up_ch1, x_max_cal_up_ch1);
    TH1D* h_cal_1_up_ch1 = new TH1D("h_cal_1_up_ch1", "Measurement 1 Calibrated", config.bin_mea, x_min_cal_up_ch1, x_max_cal_up_ch1);
  
    TH1D* h_sim_1_res_up_ch2 = new TH1D("h_sim_1_res_up_ch2", "Simulation 1 with resolution", config.bin_mea, x_min_cal_up_ch2, x_max_cal_up_ch2);
    TH1D* h_cal_1_up_ch2 = new TH1D("h_cal_1_up_ch2", "Measurement 1 Calibrated", config.bin_mea, x_min_cal_up_ch2, x_max_cal_up_ch2);
  
    TH1D* h_sim_1_res_up_sig1 = new TH1D("h_sim_1_res_up_sig1", "Simulation 1 with resolution", config.bin_mea, x_min_cal, x_max_cal);
  
    TH1D* h_sim_1_res_up_sig2 = new TH1D("h_sim_1_res_up_sig2", "Simulation 1 with resolution", config.bin_mea, x_min_cal, x_max_cal);
  
    //! Histograms (2nd file)
    TH1D* h_sim_2_res = new TH1D("h_sim_2_res", "Na22 Simulation with resolution", config.bin_mea, x_min_cal, x_max_cal);
    TH1D* h_cal_2 = new TH1D("h_cal_2", "Na22 Experiment Calibrated", config.bin_mea, x_min_cal, x_max_cal);
  
    TH1D* h_sim_2_res_up_ch1 = new TH1D("h_sim_2_res_up_ch1", "Na22 Simulation with resolution", config.bin_mea, x_min_cal_up_ch1, x_max_cal_up_ch1);
    TH1D* h_cal_2_up_ch1 = new TH1D("h_cal_2_up_ch1", "Na22 Experiment Calibrated", config.bin_mea, x_min_cal_up_ch1, x_max_cal_up_ch1);
  
    TH1D* h_sim_2_res_up_ch2 = new TH1D("h_sim_2_res_up_ch2", "Na22 Simulation with resolution", config.bin_mea, x_min_cal_up_ch2, x_max_cal_up_ch2);
    TH1D* h_cal_2_up_ch2 = new TH1D("h_cal_2_up_ch2", "Na22 Experiment Calibrated", config.bin_mea, x_min_cal_up_ch2, x_max_cal_up_ch2);
  
    TH1D* h_sim_2_res_up_sig1 = new TH1D("h_sim_2_res_up_sig1", "Na22 Simulation with resolution", config.bin_mea, x_min_cal, x_max_cal);
  
    TH1D* h_sim_2_res_up_sig2 = new TH1D("h_sim_2_res_up_sig2", "Na22 Simulation with resolution", config.bin_mea, x_min_cal, x_max_cal);

    //! Fill measurement histograms (1st and 2nd file)
    if (config.ascii == 0){
      for(int i = 0; i < entries_mea_descent_used_1; i++){
        t_mea_1->GetEntry(i);
        h_cal_1->Fill(a*(event_mea_1+0.5) + b);
        h_cal_1_up_ch1->Fill(a_up_ch1*(event_mea_1+0.5) + b_up_ch1);
        h_cal_1_up_ch2->Fill(a_up_ch2*(event_mea_1+0.5) + b_up_ch2);
      }
      for(int i = 0; i < entries_mea_descent_used_2; i++){
        t_mea_2->GetEntry(i);
        h_cal_2->Fill(a*(event_mea_2+0.5) + b);
        h_cal_2_up_ch1->Fill(a_up_ch1*(event_mea_2+0.5) + b_up_ch1);
        h_cal_2_up_ch2->Fill(a_up_ch2*(event_mea_2+0.5) + b_up_ch2);
      }
    }
    else {
      for (int i=1; i<=config.bin_mea; i++){
        h_cal_1->SetBinContent(i, content_input1[i-1]);
        h_cal_1_up_ch1->SetBinContent(i, content_input1[i-1]);
        h_cal_1_up_ch2->SetBinContent(i, content_input1[i-1]);

        h_cal_2->SetBinContent(i, content_input2[i-1]);
        h_cal_2_up_ch1->SetBinContent(i, content_input2[i-1]);
        h_cal_2_up_ch2->SetBinContent(i, content_input2[i-1]);
      }
    }
    //! Perform FFT to reduce noise (if needed)
    TH1D* h_cal_1_filtered = fft(h_cal_1, 0.1, config.thresh_fft_descent_mea, "h_cal_1_filtered", "h_cal_1_filtered");
    TH1D* h_cal_1_up_ch1_filtered = fft(h_cal_1_up_ch1, 0.1, config.thresh_fft_descent_mea, "h_cal_1_up_ch1_filtered", "h_cal_1_up_ch1_filtered");
    TH1D* h_cal_1_up_ch2_filtered = fft(h_cal_1_up_ch2, 0.1, config.thresh_fft_descent_mea, "h_cal_1_up_ch2_filtered", "h_cal_1_up_ch2_filtered");

    TH1D* h_cal_2_filtered = fft(h_cal_2, 0.1, config.thresh_fft_descent_mea, "h_cal_2_filtered", "h_cal_2_filtered");
    TH1D* h_cal_2_up_ch1_filtered = fft(h_cal_2_up_ch1, 0.1, config.thresh_fft_descent_mea, "h_cal_2_up_ch1_filtered", "h_cal_2_up_ch1_filtered");
    TH1D* h_cal_2_up_ch2_filtered = fft(h_cal_2_up_ch2, 0.1, config.thresh_fft_descent_mea, "h_cal_2_up_ch2_filtered", "h_cal_2_up_ch2_filtered");
  
    //! Prepare Energy Resolution coefficients
    double a_res = fit_energy_res(config.res_1, config.res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 0, config);
    double b_res = fit_energy_res(config.res_1, config.res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 1, config);
    double c_res = fit_energy_res(config.res_1, config.res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 2, config);
  
    double a_res_up_sig1 = fit_energy_res(config.res_1+config.step_res_1, config.res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 0, config);
    double b_res_up_sig1 = fit_energy_res(config.res_1+config.step_res_1, config.res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 1, config);
    double c_res_up_sig1 = fit_energy_res(config.res_1+config.step_res_1, config.res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 2, config);
  
    double a_res_up_sig2 = fit_energy_res(config.res_1, config.res_2+config.step_res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 0, config);
    double b_res_up_sig2 = fit_energy_res(config.res_1, config.res_2+config.step_res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 1, config);
    double c_res_up_sig2 = fit_energy_res(config.res_1, config.res_2+config.step_res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 2, config);
  
    TRandom3* ranGen = new TRandom3();
    //! Fill Simulation histograms (1st file)
    for (int i = 0; i < entries_sim_descent_used_1; i++)
    {
      t_sim_1->GetEntry(i);
      double sigma = event_sim_1*sqrt( pow(a_res,2) + pow(b_res/sqrt(event_sim_1),2) + pow(c_res/event_sim_1,2) );
      double energy = ranGen->Gaus(event_sim_1,sigma);
      h_sim_1_res->Fill(energy);
      h_sim_1_res_up_ch1->Fill(energy);
      h_sim_1_res_up_ch2->Fill(energy);
    }
    TH1D* h_sim_1_res_filtered = fft(h_sim_1_res, 0.1, config.thresh_fft_descent_sim, "h_sim_1_res_filtered", "h_sim_1_res_filtered");
    TH1D* h_sim_1_res_up_ch1_filtered = fft(h_sim_1_res_up_ch1, 0.1, config.thresh_fft_descent_sim, "h_sim_1_res_up_ch1_filtered", "h_sim_1_res_up_ch1_filtered");
    TH1D* h_sim_1_res_up_ch2_filtered = fft(h_sim_1_res_up_ch2, 0.1, config.thresh_fft_descent_sim, "h_sim_1_res_up_ch2_filtered", "h_sim_1_res_up_ch2_filtered");    
    for (int i = 0; i < entries_sim_descent_used_1; i++)
    {
      t_sim_1->GetEntry(i);
      double sigma = event_sim_1*sqrt( pow(a_res_up_sig1,2) + pow(b_res_up_sig1/sqrt(event_sim_1),2) + pow(c_res_up_sig1/event_sim_1,2) );
      h_sim_1_res_up_sig1->Fill(ranGen->Gaus(event_sim_1,sigma));
    }
    TH1D* h_sim_1_res_up_sig1_filtered = fft(h_sim_1_res_up_sig1, 0.1, config.thresh_fft_descent_sim, "h_sim_1_res_up_sig1_filtered", "h_sim_1_res_up_sig1_filtered");
  
    for (int i = 0; i < entries_sim_descent_used_1; i++)
    {
      t_sim_1->GetEntry(i);
      double sigma = event_sim_1*sqrt( pow(a_res_up_sig2,2) + pow(b_res_up_sig2/sqrt(event_sim_1),2) + pow((c_res_up_sig2)/event_sim_1,2) );
      h_sim_1_res_up_sig2->Fill(ranGen->Gaus(event_sim_1,sigma));
    }
    TH1D* h_sim_1_res_up_sig2_filtered = fft(h_sim_1_res_up_sig2, 0.1, config.thresh_fft_descent_sim, "h_sim_1_res_up_sig2_filtered", "h_sim_1_res_up_sig2_filtered");
  
    h_cal_1_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_1, config.x_max_descent_1);
    h_sim_1_res_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_1, config.x_max_descent_1);
  
    //! Fill Simulation histograms (2nd file)
    for (int i = 0; i < entries_sim_descent_used_2; i++)
    {
      t_sim_2->GetEntry(i);
      double sigma = event_sim_2*sqrt( pow(a_res,2) + pow(b_res/sqrt(event_sim_2),2) + pow(c_res/event_sim_2,2) );
      double energy = ranGen->Gaus(event_sim_2,sigma);
      h_sim_2_res->Fill(energy);
      h_sim_2_res_up_ch1->Fill(energy);
      h_sim_2_res_up_ch2->Fill(energy);
    }
    TH1D* h_sim_2_res_filtered = fft(h_sim_2_res, 0.1, config.thresh_fft_descent_sim, "h_sim_2_res_filtered", "h_sim_2_res_filtered");
    TH1D* h_sim_2_res_up_ch1_filtered = fft(h_sim_2_res_up_ch1, 0.1, config.thresh_fft_descent_sim, "h_sim_2_res_up_ch1_filtered", "h_sim_2_res_up_ch1_filtered");
    TH1D* h_sim_2_res_up_ch2_filtered = fft(h_sim_2_res_up_ch2, 0.1, config.thresh_fft_descent_sim, "h_sim_2_res_up_ch2_filtered", "h_sim_2_res_up_ch2_filtered");
  
    for (int i = 0; i < entries_sim_descent_used_2; i++)
    {
      t_sim_2->GetEntry(i);
      double sigma = event_sim_2*sqrt( pow(a_res_up_sig1,2) + pow(b_res_up_sig1/sqrt(event_sim_2),2) + pow(c_res_up_sig1/event_sim_2,2) );
      h_sim_2_res_up_sig1->Fill(ranGen->Gaus(event_sim_2,sigma));
    }
    TH1D* h_sim_2_res_up_sig1_filtered = fft(h_sim_2_res_up_sig1, 0.1, config.thresh_fft_descent_sim, "h_sim_2_res_up_sig1_filtered", "h_sim_2_res_up_sig1_filtered");
  
    for (int i = 0; i < entries_sim_descent_used_2; i++)
    {
      t_sim_2->GetEntry(i);
      double sigma = event_sim_2*sqrt( pow(a_res_up_sig2,2) + pow(b_res_up_sig2/sqrt(event_sim_2),2) + pow((c_res_up_sig2)/event_sim_2,2) );
      h_sim_2_res_up_sig2->Fill(ranGen->Gaus(event_sim_2,sigma));
    }
    TH1D* h_sim_2_res_up_sig2_filtered = fft(h_sim_2_res_up_sig2, 0.1, config.thresh_fft_descent_sim, "h_sim_2_res_up_sig2_filtered", "h_sim_2_res_up_sig2_filtered");
  
    h_cal_2_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_2, config.x_max_descent_2);
    h_sim_2_res_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_2, config.x_max_descent_2);
  
    //! Delta Chi2
    if(iteration == 2){
    	delta_chi2_1 = std::abs(chi2(h_cal_1_filtered, h_sim_1_res_filtered)) 
    	+ std::abs(chi2(h_cal_2_filtered, h_sim_2_res_filtered));
    }
    else if(iteration == 3){
    	delta_chi2_2 = std::abs(chi2(h_cal_1_filtered, h_sim_1_res_filtered) - chi2_1) 
    	+ std::abs(chi2(h_cal_2_filtered, h_sim_2_res_filtered) - chi2_2);
    }
    else if(iteration == 4){
    	delta_chi2_3 = std::abs(chi2(h_cal_1_filtered, h_sim_1_res_filtered) - chi2_1) 
    	+ std::abs(chi2(h_cal_2_filtered, h_sim_2_res_filtered) - chi2_2);
    	delta_chi2 = delta_chi2_1 + delta_chi2_2 + delta_chi2_3;
    }
    else if (iteration > 4){
    	delta_chi2_1 = delta_chi2_2;
    	delta_chi2_2 = delta_chi2_3;
    	delta_chi2_3 = std::abs(chi2(h_cal_1_filtered, h_sim_1_res_filtered) - chi2_1) 
    	+ std::abs(chi2(h_cal_2_filtered, h_sim_2_res_filtered) - chi2_2);
    	delta_chi2 = delta_chi2_1 + delta_chi2_2 + delta_chi2_3;
    }
  
    chi2_1 = chi2(h_cal_1_filtered, h_sim_1_res_filtered);
    chi2_2 = chi2(h_cal_2_filtered, h_sim_2_res_filtered);
    
    //std::cout << "chi2_1 = " << chi2_1 << "\n";
    //std::cout << "chi2_2 = " << chi2_2 << "\n";
  
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
    h_cal_1_up_ch1_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_1, config.x_max_descent_1);
    h_sim_1_res_up_ch1_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_1, config.x_max_descent_1);
    double chi2_1_up_ch1 = chi2(h_cal_1_up_ch1_filtered, h_sim_1_res_up_ch1_filtered);
  
    h_cal_1_up_ch2_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_1, config.x_max_descent_1);
    h_sim_1_res_up_ch2_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_1, config.x_max_descent_1);
    double chi2_1_up_ch2 = chi2(h_cal_1_up_ch2_filtered, h_sim_1_res_up_ch2_filtered);
  
    h_sim_1_res_up_sig1_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_1, config.x_max_descent_1);
    double chi2_1_up_sig1 = chi2(h_cal_1_filtered, h_sim_1_res_up_sig1_filtered);
  
    h_sim_1_res_up_sig2_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_1, config.x_max_descent_1);
    double chi2_1_up_sig2 = chi2(h_cal_1_filtered, h_sim_1_res_up_sig2_filtered);
  
    //! Chi2(2nd file)
    h_cal_2_up_ch1_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_2, config.x_max_descent_2);
    h_sim_2_res_up_ch1_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_2, config.x_max_descent_2);
    double chi2_2_up_ch1 = chi2(h_cal_2_up_ch1_filtered, h_sim_2_res_up_ch1_filtered);
  
    h_cal_2_up_ch2_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_2, config.x_max_descent_2);
    h_sim_2_res_up_ch2_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_2, config.x_max_descent_2);
    double chi2_2_up_ch2 = chi2(h_cal_2_up_ch2_filtered, h_sim_2_res_up_ch2_filtered);
  
    h_sim_2_res_up_sig1_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_2, config.x_max_descent_2);
    double chi2_2_up_sig1 = chi2(h_cal_2_filtered, h_sim_2_res_up_sig1_filtered);
  
    h_sim_2_res_up_sig2_filtered->GetXaxis()->SetRangeUser(config.x_min_descent_2, config.x_max_descent_2);
    double chi2_2_up_sig2 = chi2(h_cal_2_filtered, h_sim_2_res_up_sig2_filtered);
  
    //! Derivative
    double deri_chi2_up_ch1 = ((chi2_1_up_ch1+chi2_2_up_ch1)-(chi2_1+chi2_2))/config.step_channel_1;
    double deri_chi2_up_ch2 = ((chi2_1_up_ch2+chi2_2_up_ch2)-(chi2_1+chi2_2))/config.step_channel_2;
    double deri_chi2_up_sig1 = ((chi2_1_up_sig1+chi2_2_up_sig1)-(chi2_1+chi2_2))/config.step_res_1;
    double deri_chi2_up_sig2 = ((chi2_1_up_sig2+chi2_2_up_sig2)-(chi2_1+chi2_2))/config.step_res_2;
  
    // std::cout << "~~~~~~~~~~~CHECK~~~~~~~~~~~~\n";
    // std::cout << "deri_chi2_up_ch1 = " << deri_chi2_up_ch1 << "; deri_chi2_up_ch2 = " << deri_chi2_up_ch2 
    // << "\nderi_chi2_up_sig1 = " << deri_chi2_up_sig1 << "; deri_chi2_up_sig2 = " << deri_chi2_up_sig2 << "\n";
    // std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

    //! Moving along the gradient
    ch1 = ch1 - config.learning_rate_channel_1*deri_chi2_up_ch1;
    ch2 = ch2 - config.learning_rate_channel_2*deri_chi2_up_ch2;
  
    config.res_1 = config.res_1 - config.learning_rate_res_1*deri_chi2_up_sig1;
    config.res_2 = config.res_2 - config.learning_rate_res_2*deri_chi2_up_sig2;

    std::cout << "ch1 = " << ch1 << "; ch2 = " << ch2 
    << "\nSig1 = " << config.res_1 << "; Sig2 = " << config.res_2 << "\n";
    
     //std::cout << "\ndelta_Chi2_1 = " << delta_chi2_1 
     //<< "\ndelta_Chi2_2 = " << delta_chi2_2
     //<< "\ndelta_Chi2_3 = " << delta_chi2_3
     //<< "\nTotal delta_chi2 = " << delta_chi2 << "\n";
  
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

    double integral_mea_1 = h_cal_1_filtered->Integral(h_cal_1_filtered->FindBin(config.x_min_descent_1), h_cal_1_filtered->FindBin(config.x_max_descent_1));
    double integral_sim_1 = h_sim_1_res_filtered->Integral(h_cal_1_filtered->FindBin(config.x_min_descent_1), h_cal_1_filtered->FindBin(config.x_max_descent_1));
    if(integral_mea_1/integral_sim_1 < 1.){
      h_cal_1_filtered->Scale(integral_sim_1/integral_mea_1, "noSW2");
    }
    else{
      h_sim_1_res_filtered->Scale(integral_mea_1/integral_sim_1, "noSW2");
    }
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

    double integral_mea_2 = h_cal_2_filtered->Integral(h_cal_2_filtered->FindBin(config.x_min_descent_2), h_cal_2_filtered->FindBin(config.x_max_descent_2));
    double integral_sim_2 = h_sim_2_res_filtered->Integral(h_cal_2_filtered->FindBin(config.x_min_descent_2), h_cal_2_filtered->FindBin(config.x_max_descent_2));

    if(integral_mea_1/integral_sim_1 < 1.){
      h_cal_2_filtered->Scale(integral_sim_2/integral_mea_2, "noSW2");
    }
    else{
      h_sim_2_res_filtered->Scale(integral_mea_2/integral_sim_2, "noSW2");
    }
    
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
    h_sim_1_res_up_ch1_filtered->Scale(integral_mea_1/integral_sim_1, "noSW2");
    h_sim_1_res_up_ch1_filtered->Draw("same");
    c_derivative->cd(1)->Modified();
    c_derivative->cd(1)->Update();
  
    c_derivative->cd(2);    
    h_cal_2_up_ch1_filtered->SetLineColor(kRed);
    h_cal_2_up_ch1_filtered->Draw();
    h_sim_2_res_up_ch1_filtered->Scale(integral_mea_2/integral_sim_2, "noSW2");
    h_sim_2_res_up_ch1_filtered->Draw("same");
    c_derivative->cd(2)->Modified();
    c_derivative->cd(2)->Update();
  
    c_derivative->cd(3);
    h_cal_1_filtered->Draw();
    h_sim_1_res_up_ch1_filtered->Draw("same");
    c_derivative->cd(3)->Modified();
    c_derivative->cd(3)->Update();
  
    c_derivative->cd(4);
    h_cal_2_filtered->Draw();
    h_sim_2_res_up_ch1_filtered->Draw("same");
    c_derivative->cd(4)->Modified();
    c_derivative->cd(4)->Update();
  /*
    c_derivative->cd(5);    
    h_cal_2_up_ch1_filtered->SetLineColor(kRed);
    h_cal_2_up_ch1_filtered->Draw();
    h_sim_2_res_up_ch1_filtered->Scale(integral_mea_2/integral_sim_2, "noSW2");
    h_sim_2_res_up_ch1_filtered->Draw("same");
    c_derivative->cd(5)->Modified();
    c_derivative->cd(5)->Update();
  
    c_derivative->cd(6);    
    h_cal_2_up_ch2_filtered->SetLineColor(kRed);
    h_cal_2_up_ch2_filtered->Draw();
    h_sim_2_res_up_ch2_filtered->Scale(integral_mea_2/integral_sim_2, "noSW2");
    h_sim_2_res_up_ch2_filtered->Draw("same");
    c_derivative->cd(6)->Modified();
    c_derivative->cd(6)->Update();
  
    c_derivative->cd(7);
    h_cal_2_filtered->Draw();
    h_sim_2_res_up_sig1_filtered->Scale(integral_mea_2/integral_sim_2, "noSW2");
    h_sim_2_res_up_sig1_filtered->Draw("same");
    c_derivative->cd(7)->Modified();
    c_derivative->cd(7)->Update();
  
    c_derivative->cd(8);
    h_cal_2_filtered->Draw();
    h_sim_2_res_up_sig2_filtered->Scale(integral_mea_2/integral_sim_2, "noSW2");
    h_sim_2_res_up_sig2_filtered->Draw("same");
    c_derivative->cd(8)->Modified();
    c_derivative->cd(8)->Update();
  */
    if(iteration > 1){
      if (chi2_1 < thresh_delta_chi2_each){
      config.learning_rate_channel_1 = 0.3*config.learning_rate_channel_1;
      config.learning_rate_res_1 = 0.3*config.learning_rate_res_1;
      }
  
      if (chi2_2 < thresh_delta_chi2_each){
      config.learning_rate_channel_2 = 0.3*config.learning_rate_channel_2;
      config.learning_rate_res_2 = 0.3*config.learning_rate_res_2;
      }
    }
  
    if(delta_chi2 > thresh_delta_chi2){
      
      final_energy_calibration_coef_a = a;
      final_energy_calibration_coef_b = b;

      final_energy_calibration_coef_a_error = fit_energy(ch1, ch2, c_fit, 2, config);
      final_energy_calibration_coef_b_error = fit_energy(ch1, ch2, c_fit, 3, config);

      final_energy_resolution_coef_a = a_res;
      final_energy_resolution_coef_b = b_res;
      final_energy_resolution_coef_c = c_res;

      final_energy_resolution_coef_a_error = fit_energy_res(config.res_1, config.res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 3, config);
      final_energy_resolution_coef_b_error = fit_energy_res(config.res_1, config.res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 4, config);
      final_energy_resolution_coef_c_error = fit_energy_res(config.res_1, config.res_2, h_cal_1_filtered, h_cal_2_filtered, x_min_cal, x_max_cal, c_fit, 5, config);

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
  
  std::cout << "ch1 = " << ch1 << "; ch2 = " << ch2 
  << "\nSig1 = " << config.res_1 << "; Sig2 = " << config.res_2 << "\n";
  
  std::cout << "\nEnergy calibration function: E = a*Ch + b\n";
  std::cout << "Coefficients:\na = " << final_energy_calibration_coef_a << " +/- " << final_energy_calibration_coef_a_error << 
  ";\nb = " << final_energy_calibration_coef_b << " +/- " << final_energy_calibration_coef_b_error << "\n";
  std::cout << "\nEnergy resolution curve: sigma_E = E*sqrt(a^2 + b^2/x + c^2/x^2)\n";
  std::cout << "Coefficients:\na = " << final_energy_resolution_coef_a << " +/- " << final_energy_resolution_coef_a_error <<
  ";\nb = " << final_energy_resolution_coef_b << " +/- " << final_energy_resolution_coef_b_error <<
  ";\nc = " << final_energy_resolution_coef_c << " +/- " << final_energy_resolution_coef_c_error << "\n";
  std::cout << "\nChi-square 1 = " << chi2_1 << "\n";
  std::cout << "Chi-square 2 = " << chi2_2 << "\n";

  std::string file_to_be_saved = "../coefficients.txt";
  std::ofstream fileSave(file_to_be_saved);
  fileSave << "The function used to perform energy calibration is a linear equation E = a*Ch + b\n";
  fileSave << "a = " << final_energy_calibration_coef_a << "\n";
  fileSave << "b = " << final_energy_calibration_coef_b << "\n";
  fileSave << "a_error = " << final_energy_calibration_coef_a_error << "\n";
  fileSave << "b_error = " << final_energy_calibration_coef_b_error << "\n";
  fileSave << "\nThe function used for the energy resolution curve is sigma_E = E*sqrt(a^2 + b^2/x + c^2/x^2)\n";
  fileSave << "a = " << final_energy_resolution_coef_a << "\n";
  fileSave << "b = " << final_energy_resolution_coef_b << "\n";
  fileSave << "c = " << final_energy_resolution_coef_c << "\n";
  fileSave << "a_error = " << final_energy_resolution_coef_a_error << "\n";
  fileSave << "b_error = " << final_energy_resolution_coef_b_error << "\n";
  fileSave << "c_error = " << final_energy_resolution_coef_c_error << "\n";
  //! Include chi-square distance output (only within the descent range)
  fileSave << "\nChi-square distance for 1st histogram from " << config.x_min_descent_1 << " to " << config.x_max_descent_1 << " (MeV) = " << chi2_1 << "\n";
  fileSave << "Chi-square distance for 2nd histogram from " << config.x_min_descent_2 << " to " << config.x_max_descent_2 << " (MeV) = " << chi2_2 << "\n";

  fileSave.close();
  std::cout << "\nOutput file is saved to" << file_to_be_saved;
}
