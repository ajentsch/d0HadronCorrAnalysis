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



using namespace std;



int dataPointHistToTextFile(){

    //TString mainPath       = "C:/Users/ajentsch/Desktop/";    
    //TString lannyFolder     = "Data_Files_for_Lanny/";

    double phiBinShift = (TMath::Pi()/12.0);
    
    double ETA_BIN = 2.0;
    int NUM_ETA_BINS = 13;
    int NUM_PHI_BINS = 12;

    TFile * file = new TFile("d0HadronCorrMaker_OUTPUT_pt_2_10_GeV_nominal_cuts.root");
    
    TH2D * peripheral = new TH2D("","", NUM_ETA_BINS, -ETA_BIN,ETA_BIN, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D * midCentral = new TH2D("","", NUM_ETA_BINS, -ETA_BIN, ETA_BIN, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D * central = new TH2D("","", NUM_ETA_BINS, -ETA_BIN, ETA_BIN, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);

    peripheral = (TH2D*) file->Get("delRho_over_RhoRef_SideBandLeft_correlation_PtBin_2_VzBin_5_CentBin_10");
    midCentral = (TH2D*) file->Get("delRho_over_RhoRef_US_correlation_PtBin_2_VzBin_5_CentBin_10");
    central    = (TH2D*) file->Get("delRho_over_RhoRef_SideBandRight_correlation_PtBin_2_VzBin_5_CentBin_10");
    
    ofstream dataPointsForLanny;
    TString dataPointsForLannyFile = "Data_Points_12Phi_13Eta_Pt_2_10_LeftSideBand_Middle_TPC_MultBin_10.txt";
    TString dataPointsForLannyOutput = dataPointsForLannyFile;
    dataPointsForLanny.open (dataPointsForLannyOutput);
    
    dataPointsForLanny << "Eta Bin" << "\t" << "Phi Bin" << "\t" << "Content" << "\t" << "Error" << endl;
    
    for(int etaBin = 1; etaBin < NUM_ETA_BINS+1; etaBin++){
        for(int phiBin = 1; phiBin < NUM_PHI_BINS+1; phiBin++){
    
            if(TMath::Abs(peripheral->GetXaxis()->GetBinCenter(etaBin)) < .0001) { dataPointsForLanny << "0.0" << "\t" <<  peripheral->GetYaxis()->GetBinCenter(phiBin) << "\t" << peripheral->GetBinContent(etaBin, phiBin) << "\t" <<
            peripheral->GetBinError(etaBin, phiBin) << endl; }
            
            else if(TMath::Abs(peripheral->GetYaxis()->GetBinCenter(phiBin)) < .0001) { dataPointsForLanny <<  peripheral->GetXaxis()->GetBinCenter(etaBin) << "\t" << "0.0" << "\t" << peripheral->GetBinContent(etaBin, phiBin) << "\t" <<
            peripheral->GetBinError(etaBin, phiBin) << endl; }
            
            
            else dataPointsForLanny << peripheral->GetXaxis()->GetBinCenter(etaBin) << "\t" <<  peripheral->GetYaxis()->GetBinCenter(phiBin) << "\t" << peripheral->GetBinContent(etaBin, phiBin) << "\t" <<
            peripheral->GetBinError(etaBin, phiBin) << endl;
        }
    }        
    
    
    dataPointsForLanny.close();
 

    ////////////////////////////////////////////////////////////////////////////////////
   

    ofstream dataPointsForLanny;
    TString dataPointsForLannyFile = "Data_Points_12Phi_13Eta_Pt_2_10_SignalRegion_Middle_TPC_MultBin_10.txt";
    TString dataPointsForLannyOutput = dataPointsForLannyFile;
    dataPointsForLanny.open (dataPointsForLannyOutput);
    
    dataPointsForLanny << "Eta Bin" << "\t" << "Phi Bin" << "\t" << "Content" << "\t" << "Error" << endl;
    
    for(int etaBin = 1; etaBin < NUM_ETA_BINS+1; etaBin++){
        for(int phiBin = 1; phiBin < NUM_PHI_BINS+1; phiBin++){
    
            if(TMath::Abs(midCentral->GetXaxis()->GetBinCenter(etaBin)) < .0001) { dataPointsForLanny << "0.0" << "\t" <<  midCentral->GetYaxis()->GetBinCenter(phiBin) << "\t" << midCentral->GetBinContent(etaBin, phiBin) << "\t" <<
            midCentral->GetBinError(etaBin, phiBin) << endl; }
            
            else if(TMath::Abs(midCentral->GetYaxis()->GetBinCenter(phiBin)) < .0001) { dataPointsForLanny <<  midCentral->GetXaxis()->GetBinCenter(etaBin) << "\t" << "0.0" << "\t" << midCentral->GetBinContent(etaBin, phiBin) << "\t" <<
            midCentral->GetBinError(etaBin, phiBin) << endl; }
            
            
            else dataPointsForLanny << midCentral->GetXaxis()->GetBinCenter(etaBin) << "\t" <<  midCentral->GetYaxis()->GetBinCenter(phiBin) << "\t" << midCentral->GetBinContent(etaBin, phiBin) << "\t" <<
            midCentral->GetBinError(etaBin, phiBin) << endl;
        }
    }        
    
    
    dataPointsForLanny.close();
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    ofstream dataPointsForLanny;
    TString dataPointsForLannyFile = "Data_Points_12Phi_13Eta_Pt_2_10_RightSideBand_Middle_TPC_MultBin_10.txt";
    TString dataPointsForLannyOutput = dataPointsForLannyFile;
    dataPointsForLanny.open (dataPointsForLannyOutput);
    
    dataPointsForLanny << "Eta Bin" << "\t" << "Phi Bin" << "\t" << "Content" << "\t" << "Error" << endl;
    
    for(int etaBin = 1; etaBin < NUM_ETA_BINS+1; etaBin++){
        for(int phiBin = 1; phiBin < NUM_PHI_BINS+1; phiBin++){
    
            if(TMath::Abs(central->GetXaxis()->GetBinCenter(etaBin)) < .0001) { dataPointsForLanny << "0.0" << "\t" <<  central->GetYaxis()->GetBinCenter(phiBin) << "\t" << central->GetBinContent(etaBin, phiBin) << "\t" <<
            central->GetBinError(etaBin, phiBin) << endl; }
            
            else if(TMath::Abs(central->GetYaxis()->GetBinCenter(phiBin)) < .0001) { dataPointsForLanny <<  central->GetXaxis()->GetBinCenter(etaBin) << "\t" << "0.0" << "\t" << central->GetBinContent(etaBin, phiBin) << "\t" <<
            central->GetBinError(etaBin, phiBin) << endl; }
            
            
            else dataPointsForLanny << central->GetXaxis()->GetBinCenter(etaBin) << "\t" <<  central->GetYaxis()->GetBinCenter(phiBin) << "\t" << central->GetBinContent(etaBin, phiBin) << "\t" <<
            central->GetBinError(etaBin, phiBin) << endl;
        }
    }        
    
    
    dataPointsForLanny.close();
    
    
    file->Close();
}    