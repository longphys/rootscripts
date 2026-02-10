void calib(){
    TCanvas* c_fit = new TCanvas("c_fit", "Calibration Fit", 1200, 600);
	TH1D* h_calibrate = new TH1D("h_calibrate", "Calibration fit", 4096, 0, 4096);
	// h_calibrate->SetBinContent(1012, 4.6);
	// h_calibrate->SetBinContent(1066, 4.784);
	// h_calibrate->SetBinContent(1234, 5.304);
	// h_calibrate->SetBinContent(1294, 5.49);
	// h_calibrate->SetBinContent(1455, 6.002);
	// h_calibrate->SetBinContent(1973, 7.687);
	
	// h_calibrate->SetBinContent(1012, 4.6);
	h_calibrate->SetBinContent(939, 4.784);
	// h_calibrate->SetBinContent(1234, 5.304);
	h_calibrate->SetBinContent(1200, 5.49);
	h_calibrate->SetBinContent(1371, 6.002);
	h_calibrate->SetBinContent(1909, 7.687);

	TF1* f_linear = new TF1("f_linear", "[0]*x + [1]", 0, 2000);
	c_fit->cd();
	h_calibrate->Fit("f_linear", "Q");


	double coef_a = f_linear->GetParameter(0);
	double coef_b = f_linear->GetParameter(1);
	double coef_a_error = f_linear->GetParError(0);
	double coef_b_error = f_linear->GetParError(1);

    std::cout << "Calibration Coefficient a: " << coef_a << " ± " << coef_a_error << "\n";
    std::cout << "Calibration Coefficient b: " << coef_b << " ± " << coef_b_error << "\n";

	c_fit->Modified();
	c_fit->Update();
}