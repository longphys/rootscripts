#include "TRandom3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraph.h"

void circle()
{
	TRandom3* rnd = new TRandom3(0);
	double radius1 = 1.0;
	double radius2 = 4.0;
	double phi = 2.0*TMath::Pi();
	std::vector<double> v_x, v_y;
	TCanvas* c1 = new TCanvas("c1", "c1", 700, 700);
	int N = (int)5.0e6;
	for (int i = 0; i < N; i++)
	{
		double r1 = rnd->Uniform(0.0*radius1, radius1);
		double r2 = rnd->Uniform(0.0*radius2, radius2);
		double angle = rnd->Uniform(0.0, phi);
		double x = r1*TMath::Cos(angle);
		double y = r2*TMath::Sin(angle);
		v_x.push_back(x);
		v_y.push_back(y);
	}
	TGraph* gr = new TGraph(N, &v_x[0], &v_y[0]);
	gr->SetMarkerStyle(7);
	gr->SetMarkerColor(kBlue + 2);
	c1->cd();
	gr->Draw("AP");
}
