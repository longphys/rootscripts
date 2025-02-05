#include "TF1.h"
#include "TH1.h"
#include "TVirtualFFT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include <cstdlib>  // For `atof`

int main(int argc, char **argv)
{
    // Initialize application
    TApplication app("app", &argc, argv);

    // Default values
    double xmin = -5, xmax = 5;
    std::string func_expr = "sin(x)";
    double func_np = 100;

    // Handling arguments (basic example)
    if (argc > 2) {
        xmin = atof(argv[1]);  // Convert first argument to double
        xmax = atof(argv[2]);  // Convert second argument to double
    }
    if (argc > 3) {
        func_expr = argv[3];   // Use third argument as the function expression
    }
    if (argc > 4) {
        func_np = atof(argv[4]);   // Use fourth argument as the number of points used to draw the function
    }

    // Create canvas
    TCanvas* c = new TCanvas("c", "Dynamic Function Plot", 0, 0, 800, 600);

    // Define the function with user-defined range and expression
    TF1 *f1 = new TF1("f1", func_expr.c_str(), xmin, xmax);
    f1->SetNpx(func_np);
    f1->SetLineColor(kBlue+1);
    f1->SetTitle("User-defined Function;x; f(x)");
    f1->Draw();

    // Update canvas
    c->Modified(); 
    c->Update();

    // Handle window close event
    TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
    rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

    // Run the application
    app.Run();

    return 0;
}

