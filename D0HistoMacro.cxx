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


int D0HistoMacro(){

    
    
    
    TFile* input = new TFile("C:/Users/ajentsch/Desktop/output files from analysis/best one to date.root", "READ");
    
    
    TH1D* invMass            = new TH1D("invariant mass", "D0 Invariant Mass", 45, 1.6, 2.1);
    TH2D* AngCorr2D          = new TH2D("2D Ang Dist.", "2D Ang Dist.", 500, -2, 2, 500, -TMath::PiOver2(), 3*TMath::PiOver2());
    //TH1D* signal             = new TH1D("signal", "signal", 45, 1.6, 2.1);
    //TH1D* backgroundFitHist  = new TH1D("bkgd", "bkgd", 45, 1.6, 2.1);
    
    invMass = (TH1D*)input->Get("D0 minus LS BG");
    invMass->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    invMass->Sumw2();
    
    AngCorr2D = (TH2D*)input->Get("2D Ang Corr");
    AngCorr2D->GetXaxis()->SetTitle("#Delta#eta");
    AngCorr2D->GetYaxis()->SetTitle("#Delta#phi");
    AngCorr2D->SetTitle("Sibling Dist.");
    
    //TF1* invMassFit = new TF1("invMassFit", "gaus", 1.82, 1.91);
    //TF1* background = new TF1("background", "pol3", 1.6, 2.1);
    //invMassFit->SetParLimits(0,580,630);
    //invMassFit->SetParLimits(1,1.82,1.89);
    
    //TCanvas *c = new TCanvas("c2", "Histograms", 850, 850);
    
    //histogram->Draw();
    
    //c->SaveAs("C:/Users/ajentsch/Desktop/CeHistogram.png"); 

    //invMass->Fit("invMassFit","R");
    //CreateHistogram();
    //invMass->Fit("background", "R");
    //backgroundFitHist = (TH1D*)background->GetHistogram();
    TFile* output = new TFile("C:/Users/ajentsch/Desktop/outputTest.root", "RECREATE");

    //signal->Add(invMass, backgroundFitHist, 1, -1);
    
    invMass->Write();
    AngCorr2D->Write();
    //backgroundFitHist->Write();
    //signal->Write();
    
    input->Close();
    output->Close();
    
}