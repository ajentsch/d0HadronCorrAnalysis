/////ROOT LIBRARIES///////////////
#include "TROOT.h"
#include "TSystem.h" 
#include "TApplication.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <Riostream.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <iomanip>
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TPad.h>
#include <TMath.h>
#include <TF1.h>
#include <THnSparse.h>
#include <TLeaf.h>
#include <TLatex.h>
#include "TError.h"
//gErrorIgnoreLevel = 4000;  //default is kInfo

///////C++ Libraries/////////////////
#include <iostream>
#include <iomanip>
#include <math.h>
#include <map>
#include <vector>
#include <utility>
#include <climits>
#include <sstream>
#include <string>
#include <fstream>
#include <TAttMarker.h>
#include <algorithm>
#include <memory>
#include <iterator>
#include <ctype.h>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>



using namespace std;


int dipoleLimit(){

    double phiBinShift = (TMath::Pi()/12.0);

    //myfit[2][0]->SetParameters(-.01, -.01, .03, .11, .3, .3, .03, .8);
    
    double offset    = .005;
    double dipoleAmp = -0.005;
    double nonPeriodicGaussAmp    = .02;
    double nonPeriodicGuassWidth  = .6;
    double periodicGaussAmp       = .01;  
    double periodicGaussWidth     = 1.25;
    
    
    TH1D * displayHistogram = new TH1D("","", 100, -2*TMath::TwoPi(), 2*TMath::TwoPi());
    
    TF1* dipole = new TF1("dipole","[0]*cos(x) + [1]",  -TMath::TwoPi(), TMath::TwoPi());
    
    
    TF1* nonPeriodicGauss = new TF1("etaGauss","[0]*exp(-0.5*((x-3.14159)*(x-3.14159))/([1]*[1]))", -TMath::TwoPi(), TMath::TwoPi());
    TF1* periodicGauss = new TF1("ASGauss","[0]*(exp(-0.5*((x-3.14159)*(x-3.14159))/([1]*[1]))+exp(-0.5*((x+3.14159)*(x+3.14159))/([1]*[1])))",-TMath::TwoPi(), TMath::TwoPi());
    

    dipole->SetParameter(0,dipoleAmp); 
    dipole->SetParameter(1,offset); 
    periodicGauss->SetParameter(0, periodicGaussAmp);
    periodicGauss->SetParameter(1, periodicGaussWidth);

    TCanvas * can = new TCanvas("", "", 1100, 1000);
    
    dipole->SetLineColor(2);
    periodicGauss->SetLineColor(1);
    
    displayHistogram->SetMaximum(.03);
    displayHistogram->Draw();
    dipole->Draw("SAME");
    periodicGauss->Draw("SAME");
    
    can->SaveAs("dipole_Gaussian_comparison.png");
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   

}