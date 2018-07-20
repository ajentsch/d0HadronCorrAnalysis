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
gErrorIgnoreLevel = 4000;  //default is kInfo

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
#include <istream>

/****
    
    Code to extract error histograms from delRho over Rho histograms
    Author: Alex Jentsch -- alex.jentsch@utexas.edu
    
    Inputs: delRho over Rho correlation hists from d0HadronCorrMaker
    
    Outputs: Histograms of errors
    
    Last update: 3/29/2018
    
    
*****/


using namespace std;


int D0ErrorExtractionFromDelRhoHistograms(){

//constants///////////////////////////////////////
    int NUM_ETA_CORR_BINS = 9;
    int NUM_PHI_CORR_BINS = 12;
    int NUM_ETA_PRIM_BINS = 9;
    int NUM_PHI_PRIM_BINS = 24;
    
    double phiBinShift = TMath::Pi()/12.0;
    double phiBinShiftPrimary = TMath::Pi()/24.0;
    
    TH2D* delRhoOverRhoHists[3][5][10][16];   //final errors stored here
    TH2D* errorHists[3][5][10][16];   //final errors stored here
    
    ////////////////////////////////////////////////
    
    TFile * inputFile = new TFile("corrHistogramMaker_best_nominal_new_eff_correction.root");
    TFile * outputFile = new TFile("extracted_error_output_nominal_root_errors.root", "RECREATE");
    
    TString vzLabel          = "_VzBin_";
    TString multLabel        = "_CentBin_";
   
    TString bandLabel[3]     = {"SideBandLeft_", "US_", "SideBandRight_"};
    TString ptBinLabel[5]    = {"0", "1", "2", "3", "4"};
    TString vzBinLabel[10]   = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
    TString multBinLabel[16] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"};
    
    
    
    TString mainLabel = "delRho_over_RhoRef_";
    TString corrLabel = "correlation_PtBin_";
    TString errorHistLabel = "ErrorHistogram";
  
    TString str1;
    
    for(int band = 0; band < 3; band++){
        for(int ptBin = 0; ptBin < 5; ptBin++){
            for(int vzBin = 0; vzBin < 10; vzBin++){
                for(int multBin = 0; multBin < 16; multBin++){
                
                    str1 = mainLabel + bandLabel[band] + corrLabel + ptBinLabel[ptBin] + vzLabel + vzBinLabel[vzBin] + multLabel + multBinLabel[multBin];
                    delRhoOverRhoHists[band][ptBin][vzBin][multBin] = new TH2D(str1, str1, NUM_ETA_CORR_BINS, -2, 2, NUM_PHI_CORR_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                    delRhoOverRhoHists[band][ptBin][vzBin][multBin]->SetDirectory(0);
                    
                    str1 = errorHistLabel + ptBinLabel[ptBin] + vzLabel + vzBinLabel[vzBin] + multLabel + multBinLabel[multBin];
                    errorHists[band][ptBin][vzBin][multBin] = new TH2D(str1, str1, NUM_ETA_CORR_BINS, -2, 2, NUM_PHI_CORR_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                    errorHists[band][ptBin][vzBin][multBin]->SetDirectory(0);
                
                }
            }
        }
    }