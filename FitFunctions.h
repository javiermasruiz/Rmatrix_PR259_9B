#include "function_library.h"
#include <vector>
#include <TF1.h>
#include <TF1Convolution.h>

double Sp = -0.1859; //MeV - unbound to proton emission
double r0 = 1.25; //fm 

TF1 *fLineshape;
TF1 *fResponse;
TF1Convolution *fConvolution;

std::vector<double> vEx;
std::vector<double> vgamma2_proton; //reduced widths
std::vector<double> ell_proton;
std::vector<double> InitialAmps;
//Wigner_width(int Z1, int Z2, double A1, double A2, int L, double E, double r)

void LoadPeakData()
{
    vEx.push_back(0);
    vgamma2_proton.push_back(0.54e-3/Wigner_width(4,1,8,1,1,0.0-Sp, r0 * (pow(8.,1./3.)+pow(1.,1./3.)))); // See Wheldon 2015  
    ell_proton.push_back(1);
    InitialAmps.push_back(1000);
    
    vEx.push_back(2.345);
    vgamma2_proton.push_back(81.e-3/Wigner_width(4,1,8,1,3,2.345-Sp, r0 * (pow(8.,1./3.)+pow(1.,1./3.)))); // See Wheldon 2015  
    ell_proton.push_back(3);
    InitialAmps.push_back(300);
   
    vEx.push_back(2.751);
    vgamma2_proton.push_back(0.614/Wigner_width(4,1,8,1,2,2.751-Sp, r0 * (pow(8.,1./3.)+pow(1.,1./3.)))); 
    ell_proton.push_back(2);
    InitialAmps.push_back(100);
   /* 
    vEx.push_back(2.780);
    vgamma2_proton.push_back(3.13/Wigner_width(4,1,8,1,1,2.780-Sp, r0 * (pow(8.,1./3.)+pow(1.,1./3.)))); // See Wheldon 2015  
    ell_proton.push_back(1);
    InitialAmps.push_back(100);
    
    vEx.push_back(4.8);
    vgamma2_proton.push_back(1.2/Wigner_width(4,1,8,1,1,4.8-Sp, r0 * (pow(8.,1./3.)+pow(1.,1./3.)))); // See Wheldon 2015  
    ell_proton.push_back(1);
    InitialAmps.push_back(10);
*/
    printf("Loaded peak data\n");
}

double Lineshape(double *x, double *pars)
{
    double result = 0;
    
    //pars[0] - background
    //pars[1] - r0
    //pars[2+3*i] - Ex
    //pars[3+3*i] - gamma2
    //pars[4+3*i] - Amplitude
    result = pars[0];
    
    for(unsigned int i=0;i<vEx.size();i++)
    {
        double Gamma_p = pars[3+3*i] * Wigner_width(4,1,8,1,ell_proton.at(i),x[0] - Sp, pars[1] * (pow(8.,1./3.)+pow(1.,1./3.)));
        result += 4*pars[4+3*i] * pow(Gamma_p,2.) / (pow(x[0] - pars[2+3*i],2.) + 0.25*pow(Gamma_p,2.));  
        // Convolucion por suma. vEx.size son los pasos de la convolucion ??
    }
    
//     printf("x = %f\tresult = %f\n",x[0],result);
    return result;
}

double Response(double *x, double *pars)
{
    //response function for the spectrometer - currently a Gaussian but may become more complicated - talk to Kevin about Landau    
    //pars[0] - resolution - everything else is just "normal" - 1 to use normalised Gaussian (hopefully!)
    
    return TMath::Gaus(x[0],0,pars[0],1);
}

double FitFunctionSingles(double *x, double *pars)
{
  //  return fConvolution->EvalNumConv(x[0]);
}
