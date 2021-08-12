#include <TCanvas.h>
#include <TH1.h>
#include <TFile.h>
#include "FitFunctions.h"
#include "TStopwatch.h"


void FitSpectra();

int main()
{
    TStopwatch timer;
    
    FitSpectra();
    
    timer.Print();

    return 0;
}

void FitSpectra()
{
    TFile *fin = TFile::Open("b_spect.root");
    
    TH1F *hEx_singles = (TH1F*)fin->Get("Ex_sinlges");
    
    TCanvas *c1 = new TCanvas();
    
    hEx_singles->SetStats(0);
    hEx_singles->Draw("E");
    
    LoadPeakData();
    
    fLineshape = new TF1("fLineshape",Lineshape,Sp+0.05,7,2+3*vEx.size());
    fLineshape->SetNpx(1.e6);
    fLineshape->SetParameter(0,0.);   // Parametro background
    fLineshape->SetParameter(1,1.25);  // Parametro r0
    
    printf("Wigner check: %f\n",Wigner_width(4,1,8,1,1,0. - Sp, r0 * (pow(8.,1./3.)+pow(1.,1./3.))));
    printf("vEx.size: %d\n",(int)vEx.size());
    

    for(unsigned int i=0;i<vEx.size();i++)
    {
        fLineshape->SetParameter(2+3*i,vEx.at(i));
        fLineshape->SetParameter(3+3*i,vgamma2_proton.at(i));
        fLineshape->SetParameter(4+3*i,1000.);
    }
    
    std::cout << "Initial Parameters :" << std::endl;
    std::cout << "p0 :" << fLineshape->GetParameter(0) << std::endl;
    std::cout << "p1 :" << fLineshape->GetParameter(1) << std::endl;
    std::cout << "p2 :" << fLineshape->GetParameter(2) << std::endl;
    std::cout << "p3 :" << fLineshape->GetParameter(3) << std::endl;
    std::cout << "p4 :" << fLineshape->GetParameter(4) << std::endl;
    std::cout << "p5 :" << fLineshape->GetParameter(5) << std::endl;
    std::cout << "p6 :" << fLineshape->GetParameter(6) << std::endl;
    std::cout << "p7 :" << fLineshape->GetParameter(7) << std::endl;
    std::cout << "p8 :" << fLineshape->GetParameter(8) << std::endl;
    std::cout << "p9 :" << fLineshape->GetParameter(9) << std::endl;
    std::cout << "p10 :" << fLineshape->GetParameter(10) << std::endl;
    std::cout << "p11 :" << fLineshape->GetParameter(11) << std::endl;
    std::cout << "p12 :" << fLineshape->GetParameter(12) << std::endl;
    std::cout << "p13 :" << fLineshape->GetParameter(13) << std::endl;
    std::cout << "p14 :" << fLineshape->GetParameter(14) << std::endl;
    std::cout << "p15 :" << fLineshape->GetParameter(15) << std::endl;
    std::cout << "p16 :" << fLineshape->GetParameter(16) << std::endl;
    std::cout << "p17 :" << fLineshape->GetParameter(17) << std::endl;
    std::cout << "p18 :" << fLineshape->GetParameter(18) << std::endl;
    std::cout << "p19 :" << fLineshape->GetParameter(19) << std::endl;
    std::cout << "END reading initial parameters.." << std::endl;
    std::cout << " " << std::endl;

    printf("fLineshape->Eval(0): %f\n",fLineshape->Eval(0));
    std::cout << " " << std::endl;
    
    fLineshape->Draw("same");
    
    c1->SaveAs("figures/initial/1.lineshape_raw.png");
 //   c1->SaveAs("figures/initial/lineshape_raw.pdf");
 //   c1->SaveAs("figures/initial/lineshape_raw.C");
    
    
    fResponse = new TF1("fResponse",Response,-5,5,1);
    fResponse->SetNpx(1.e6);
    fResponse->SetParameter(0,2.48794e-02);
    fResponse->Draw();
    
    c1->SaveAs("figures/initial/2.response_raw.png");
 //   c1->SaveAs("figures/initial/response_raw.pdf");
 //   c1->SaveAs("figures/initial/response_raw.C");
    
    hEx_singles->Draw("E");
    
    double *params = new double[fLineshape->GetNpar()+fResponse->GetNpar()];
    
    fConvolution = new TF1Convolution(fLineshape,fResponse,Sp+0.05,7,true);
    fConvolution->SetNofPointsFFT(1000);
    params[0] = 0.;
    
    for(unsigned int i=0;i<vEx.size();i++)
    {
        params[2+3*i] = vEx.at(i);
        params[3+3*i] = vgamma2_proton.at(i);
        params[4+3*i] = InitialAmps.at(i);
    }
    
    params[fLineshape->GetNpar()+fResponse->GetNpar()-1] = fResponse->GetParameter(0);
    
    TF1 *fConvOutput = new TF1("fConvOutput",*fConvolution,Sp+0.075,6.5,fConvolution->GetNpar()); 
    fConvOutput->SetNpx(1.e6);
    fConvOutput->SetParameters(params);
    
    for(int i=0;i<vEx.size();i++)
        fConvOutput->SetParLimits(4+3*i,0.,2000.);
    
    fConvOutput->FixParameter(0,0.);
    fConvOutput->FixParameter(1,1.25);
//     fConvOutput->SetParLimits(1,1.0,2.5);
    fConvOutput->FixParameter(fConvOutput->GetNpar()-1,fConvOutput->GetParameter(fConvOutput->GetNpar()-1));
    
    for(int i=0;i<fConvOutput->GetNpar();i++)
        printf("Parameter %d: %f\n",i,fConvOutput->GetParameter(i));
    
    fConvOutput->SetLineColor(3);
    fConvOutput->Draw("same");

// El problema con el fit pueden ser los parametros iniciales  
//    for(int i=0;i<100;i++){ r = hEx_singles->Fit(fConvOutput,"SR");}   // multiples fit
  
    TFitResultPtr r = hEx_singles->Fit(fConvOutput,"BRLME");   // agregue V y S,  -> default: BRLME

    int status = int(r);
    printf("Fit Status %d\n",status);

    std::cout << " " << std::endl;
    std::cout << " Chi-Square of Fit: " << fConvOutput->GetChisquare() << std::endl;
    std::cout << " Degrees of Freedom: " << fConvOutput->GetNDF() << std::endl;
    
//     printf("Conv Output %f\n",fConvOutput->Eval(0.0));
    
    c1->SaveAs("figures/initial/3.conv_raw.png");
 //   c1->SaveAs("figures/initial/conv_raw.pdf");
 //   c1->SaveAs("figures/initial/conv_raw.C");
}
