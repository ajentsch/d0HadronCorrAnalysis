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


int fitFunctionTester(){

    double phiBinShift = (TMath::Pi()/12.0);

    //myfit[2][0]->SetParameters(-.01, -.01, .03, .11, .3, .3, .03, .8);
    
    double A0       = -0.002;
    double AD       = -.006;
    double AQ       = .018;
    double A1       = .075;  
    double phiWidth = .15;
    double etaWidth = .18;   
    
    //stuff for D*
    double A2       = .09;  
    //double A2PhiWidth = .2;
    double A2EtaWidth = .0002;
    
    double etaGaussPeak  = .06;
    double etaGaussWidth = .02;
    
    double ASGaussPeak  = .02 ;
    double ASGaussWidth = 1.0;

    TH2D * offsetPlot = new TH2D("offset", "offset", 7, -1.556, 1.556, 12, -TMath::PiOver2()+phiBinShift , 3*TMath::PiOver2()+phiBinShift);
    TH2D * dipolePlot = new TH2D("dipole", "dipole", 7, -1.556, 1.556, 12, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D * quadrupolePlot = new TH2D("quadrupole", "quadrupole", 7, -1.556, 1.556, 12, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D * jetPeakPlot = new TH2D("Jet Peak", "Jet Peak", 7, -1.556, 1.556, 12, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D * dStarPeakPlot = new TH2D("D* Peak", "D* Peak", 7, -1.556, 1.556, 12, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D * etaGaussPlot = new TH2D("eta gauss", "eta gauss", 7, -1.556, 1.556, 12, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D * ASGaussPlot = new TH2D("Away side gauss", "Away side gauss", 7, -1.556, 1.556, 12, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D * fullFitPlot = new TH2D("Full fit", "Full fit", 7, -1.556, 1.556, 12, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D * fullFitASGaussPlot = new TH2D("Full fit with AS Gauss", "Full fit with AS Gauss", 7, -1.556, 1.556, 12, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D * dipolePlusQuadrupolePlot = new TH2D("Dipole plus Quadrupole", "Dipole plus Quadrupole", 7, -1.556, 1.556, 12, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);

    TF2* offsetFit = new TF2("offset","[0]", -1.556, 1.556, -TMath::PiOver2(), 3*TMath::PiOver2());
    TF2* dipoleFit = new TF2("dipole","[0]*cos(y)", -1.556, 1.556,-TMath::PiOver2(), 3*TMath::PiOver2());
    TF2* quadrupoleFit = new TF2("quadrupole","[0]*2*cos(2*y)", -1.556, 1.556,-TMath::PiOver2(), 3*TMath::PiOver2());
    TF2* jetPeakFit = new TF2("jetPeak","[0]*exp(-0.5*(x*x)/([1]*[1]))*(exp(-0.5*(y*y)/([2]*[2]))+exp(-0.5*((y-6.28)*(y-6.28))/([2]*[2])))", -1.556, 1.556, -TMath::PiOver2(), 3*TMath::PiOver2());
    TF2* dStarPeakFit = new TF2("dStarPeak","[0]*exp(-0.5*(x*x)/([1]*[1]))*(exp(-0.5*(y*y)/([2]*[2]))+exp(-0.5*((y-6.28)*(y-6.28))/([2]*[2])))", -1.556, 1.556, -TMath::PiOver2(), 3*TMath::PiOver2());
    
    TF2* etaGaussFit = new TF2("etaGauss","[0]*exp(-0.5*(x*x)/([1]*[1]))", -1.556, 1.556,-TMath::PiOver2(), 3*TMath::PiOver2());
    TF2* ASGaussFit = new TF2("ASGauss","[0]*exp(-0.5*(x*x)/([1]*[1]))*(exp(-0.5*((y-3.14)*(y-3.14))/([2]*[2]))+exp(-0.5*((y+3.14)*(y+3.14))/([2]*[2])))", -1.556, 1.556,-TMath::PiOver2(), 3*TMath::PiOver2());
    
    TF2* ASEtaGaussAndDipole = new TF2("AS_Eta_Gauss_and_Dipole", " [0]*cos(y)*exp(-0.5*(x*x)/([1]*[1]))", -1.556, 1.556,-TMath::PiOver2(), 3*TMath::PiOver2());
	
    TF2* fullFitDStarPeak = new TF2("fullFit","[0] + [1]*cos(y) + [2]*2*cos(2*y) + ([3]*exp(-0.5*(x*x)/([4]*[4]))+([5]*exp(-0.5*(x*x)/([6]*[6]))))*(exp(-0.5*(y*y)/([7]*[7]))+exp(-0.5*((y-6.28)*(y-6.28))/([7]*[7]))) ",-1.556, 1.556,-TMath::PiOver2(), 3*TMath::PiOver2());
    TF2* fullFit6Parameter = new TF2("fullFit6","[0] + [1]*cos(y) + [2]*2*cos(2*y) + [3]*exp(-0.5*(x*x)/([4]*[4]))*(exp(-0.5*(y*y)/([5]*[5]))+exp(-0.5*((y-6.28)*(y-6.28))/([5]*[5]))) ",-1.556, 1.556,-TMath::PiOver2(), 3*TMath::PiOver2());
    TF2* fullFitASGauss = new TF2("fullFitASGauss","[0] + [1]*2*cos(2*y) + ([2]*exp(-0.5*(x*x)/([3]*[3]))+([4]*exp(-0.5*(x*x)/([5]*[5]))))*(exp(-0.5*(y*y)/([6]*[6]))+exp(-0.5*((y-6.28)*(y-6.28))/([6]*[6]))) + [7]*exp(-0.5*(y*y)/([8]*[8])) + [9]*(exp(-0.5*((y-3.14)*(y-3.14))/([10]*[10]))+exp(-0.5*((y+3.14)*(y+3.14))/([10]*[10])))",-1.556, 1.556,-TMath::PiOver2(), 3*TMath::PiOver2());
    //    TF2* fullFitASGauss = new TF2("fullFitASGauss","[0] + [1]*2*cos(2*y) + [2]*((exp(-0.5*(y*y)/([3]*[3]))+exp(-0.5*((y-6.28)*(y-6.28))/([3]*[3])))*exp(-0.5*(x*x)/([4]*[4]))) + [5]*exp(-0.5*(y*y)/([6]*[6])) + [7]*(exp(-0.5*((y-3.14)*(y-3.14))/([8]*[8]))+exp(-0.5*((y+3.14)*(y+3.14))/([8]*[8])))",-1.556, 1.556,-TMath::PiOver2(), 3*TMath::PiOver2());

    TF2* dipolePlusQuadrupole = new TF2("dipolePlusQuadrupole", "[0]*cos(y) + [1]*2*cos(2*y)", -1.556, 1.556, -TMath::PiOver2(), 3*TMath::PiOver2());

    offsetFit->FixParameter(0, A0);
    dipoleFit->FixParameter(0, AD); 
    quadrupoleFit->FixParameter(0, AQ);
    dipolePlusQuadrupole->FixParameter(0, AD);
    dipolePlusQuadrupole->FixParameter(1, AQ);
    jetPeakFit->FixParameter(0, A1);
    jetPeakFit->FixParameter(1, etaWidth);
    jetPeakFit->FixParameter(2, phiWidth);
    dStarPeakFit->FixParameter(0, A2);
    dStarPeakFit->FixParameter(1, A2EtaWidth);
    dStarPeakFit->FixParameter(2, phiWidth);
    
   
    etaGaussFit->FixParameter(0, etaGaussPeak);
    etaGaussFit->FixParameter(1, etaGaussWidth); 
    ASGaussFit->FixParameter(0, ASGaussPeak);
    ASGaussFit->FixParameter(1, ASGaussWidth); 
    
    cout << "Test: " << jetPeakFit->Eval(0, 0, 0, 0) << endl;
    
    fullFit6Parameter->FixParameter(0, A0);    //offset A0
    fullFit6Parameter->FixParameter(1, AD);    //v1     AD
    fullFit6Parameter->FixParameter(2, AQ);    //v2     AQ
    fullFit6Parameter->FixParameter(3, A1);    //A1 jet (tall) amp
    //fullFit->FixParameter(4, A1SigEta);     //eta tall width
    //fullFit->FixParameter(5, A2);    //A2 jet (tall) amp
    //fullFit->FixParameter(6, A2SigEta);     //eta short width
    fullFit6Parameter->FixParameter(5, phiWidth);    //jet peak phi width
    //fullFit->FixParameter(8, etaPeak);     //eta_Gauss_peak
    fullFit6Parameter->FixParameter(4, etaWidth);    //eta_Gauss width
    
    
    fullFitDStarPeak->FixParameter(0, A0);    //offset A0
    fullFitDStarPeak->FixParameter(1, AD);    //v1     AD
    fullFitDStarPeak->FixParameter(2, AQ);    //v2     AQ
    fullFitDStarPeak->FixParameter(3, A1);    //A1 jet (tall) amp
    fullFitDStarPeak->FixParameter(4, etaWidth);     //eta tall width
    fullFitDStarPeak->FixParameter(5, A2);    //A2 D* amp
    fullFitDStarPeak->FixParameter(6, A2EtaWidth);     //eta short width
    fullFitDStarPeak->FixParameter(5, phiWidth);    //jet peak phi width
    //fullFit->FixParameter(8, etaPeak);     //eta_Gauss_peak
    fullFitDStarPeak->FixParameter(7, phiWidth);    //eta_Gauss width
    
    //fullFitASGauss->FixParameter(0, A0);    //offset A0
    //fullFitASGauss->FixParameter(1, AQ);    //v2     AQ
    //fullFitASGauss->FixParameter(2, A1);    //A1 jet (tall) amp
    //fullFitASGauss->FixParameter(3, A1SigEta);     //eta tall width
    //fullFitASGauss->FixParameter(4, A2);    //A2 jet (tall) amp
    //fullFitASGauss->FixParameter(5, A2SigEta); //eta short width
    //fullFitASGauss->FixParameter(6, sigPhi); //jet peak phi width
    //fullFitASGauss->FixParameter(7, etaPeak);     //eta_Gauss_peak
    //fullFitASGauss->FixParameter(8, etaWidth);    //eta_Gauss width
    //fullFitASGauss->FixParameter(9, ASGaussPeak);
    //fullFitASGauss->FixParameter(10, ASGaussPeak);
    
    TCanvas * canvas = new TCanvas("", "", 4800, 2400);
    canvas->Divide(4,2);
    
    //canvas->cd(1);
    //offsetPlot->Eval(offsetFit);
    //offsetPlot->SetMaximum(A1);
    //offsetPlot->Draw("SURF1");
    canvas->cd(1);
    dipolePlot->Eval(dipoleFit);
    dipolePlot->SetMaximum(A1);
    dipolePlot->Draw("SURF1");
    canvas->cd(2);
    quadrupolePlot->Eval(quadrupoleFit);
    quadrupolePlot->SetMaximum(A1);
    quadrupolePlot->Draw("SURF1");
    canvas->cd(3);
    dipolePlusQuadrupolePlot->Eval(dipolePlusQuadrupole);
    dipolePlusQuadrupolePlot->SetMaximum(A1);
    dipolePlusQuadrupolePlot->Draw("SURF1");
    canvas->cd(4);
    jetPeakPlot->Eval(jetPeakFit);
    jetPeakPlot->SetMaximum(A1);
    jetPeakPlot->Draw("SURF1");
    gPad->SetTheta(0); // default is 30
    gPad->SetPhi(90); // default is 30
    gPad->Update();
    
    canvas->cd(5);
    //jetPeakPlot->Eval(jetPeakFit);
    //jetPeakPlot->SetMaximum(A1);
    jetPeakPlot->Draw("SURF1");
    gPad->SetTheta(20); // default is 30
    gPad->SetPhi(0); // default is 30
    gPad->Update();
    
    canvas->cd(6);
    dStarPeakPlot->Eval(dStarPeakFit);
    dStarPeakPlot->SetMaximum(A1);
    dStarPeakPlot->Draw("SURF1");
    gPad->SetTheta(20); // default is 30
    gPad->SetPhi(0); // default is 30
    gPad->Update();
    
    //canvas->cd(5);
    //etaGaussPlot->Eval(etaGaussFit);
    //etaGaussPlot->SetMaximum(A1);
    //etaGaussPlot->Draw("SURF1");
    canvas->cd(7);
    fullFitPlot->Eval(fullFitDStarPeak);
    fullFitPlot->Draw("SURF1");
    gPad->SetTheta(20); // default is 30
    gPad->SetPhi(35); // default is 30
    gPad->Update();
    //ASGaussPlot->Eval(ASGaussFit);
    //ASGaussPlot->SetMaximum(A1);
    //ASGaussPlot->Draw("SURF1");
    canvas->cd(8);
    fullFitPlot->Draw("SURF1");
    gPad->SetTheta(0); // default is 30
    gPad->SetPhi(90); // default is 30
    gPad->Update();
    
    //canvas->cd(8);
    //fullFitPlot->Draw("SURF1");
    //fullFitASGaussPlot->Eval(fullFitASGauss);
    //fullFitASGaussPlot->Draw("SURF1");
    //gPad->SetTheta(10); // default is 30
    //gPad->SetPhi(5); // default is 30
    //gPad->Update();
    
    canvas->SaveAs("fits.png");
    
    TCanvas * single = new TCanvas("", "", 1300, 1000);
    offsetFit->SetParameter(0, 0);
    offsetFit->SetNameTitle("","");
    offsetFit->Draw("SURF1");
    single->SaveAs("offset.png");
    dipoleFit->SetParameter(0, -.5);
    dipoleFit->SetMaximum(2.0);
    dipoleFit->SetNameTitle("","");
    dipoleFit->Draw("SURF1");
    single->SaveAs("dipole.png");
    quadrupoleFit->SetParameter(0, .5);
    quadrupoleFit->SetMaximum(2.0);
    quadrupoleFit->SetNameTitle("","");
    quadrupoleFit->Draw("SURF1");
    single->SaveAs("Quadrupole.png");
    jetPeakFit->SetParameter(0, 1.0);
    jetPeakFit->SetParameter(1, .8);
    jetPeakFit->SetParameter(2, .8);
    jetPeakFit->SetNameTitle("","");
    jetPeakFit->Draw("SURF1");
    single->SaveAs("SSPeak.png");
    ASGaussFit->SetParameter(0, 1.0);
    ASGaussFit->SetParameter(1, 1.2);
    ASGaussFit->SetParameter(2, 1.0);
    ASGaussFit->SetMaximum(2.0);
    ASGaussFit->SetNameTitle("","");
    ASGaussFit->Draw("SURF1");
    single->SaveAs("ASPeak.png");
	ASEtaGaussAndDipole->SetParameter(0, -.01);
	ASEtaGaussAndDipole->SetParameter(1, .3);
	ASEtaGaussAndDipole->SetMaximum(.015);
    ASEtaGaussAndDipole->SetNameTitle("","");
    ASEtaGaussAndDipole->Draw("SURF1");
	single->SaveAs("AS_Eta_Gauss_And_dipole.png");
    
    TFile * output = new TFile("output.root", "RECREATE");
    
    fullFitPlot->Write();
    
    output->Close();

}