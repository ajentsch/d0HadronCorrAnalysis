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

plotFinalFitResultsUpdated(){


    
    bool USE_ETA_GAP = false;
    
    

    double NDF = 1;
    
    double projectionScaleFactor;
    
   
    int NUM_PT_BINS = 3;
    int FIRST_PT_BIN = 2;
    int NUM_PHI_BINS = 12;
    int NUM_ETA_BINS = 13;
    double ETA_RANGE = 1.05;
    double phiBinShift = (TMath::Pi()/12.0);
   
    double nBarPrime[3] = {74.4, 351, 1027};
    
    //double nBarPrime[3] = {100, 351, 1027};
   
    TString inputFileName  = "d0HadronCorrMaker_OUTPUT_pt_2_10_GeV_nominal_cuts.root";
	TString outputFileName = "fittingOutput_OUTPUT_12_phi_bins_pt_2_10_symmeterized_PROJECTIONS.root";
    TString outputFolder   = "deltaPhiProjectionFits/";
    TString path = "C:/Users/ajentsch/desktop/D0_Analysis_Code_Local_Copies/D0CorrelationFitter/";
	
	TString fullySubtractedLabelCustomCentrality = "Full_DStar_Correction_FullSubtractedCorr_SideBand__PtBin_";
	TString ptBinLabel[5] = {"0", "1", "2", "3", "4"};
    TString centralityBin[3] = {"Peripheral", "MidCentral", "Central"};
    TString fileType = ".png";
    
    TString lightFlavorPtLabel[10] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
    
    TString ptLabel     = "PtBin_";
    TString dipoleLabel = "Dipole_AD_";
    
    TString offsetLabel = "Offset_A0_";
    TString quadLabel   = "Quadrupole_AQ_";
    
    TString jetAmpLabel = "Same_Side_Amplitude_A_{SS}_";
    TString sigEtaLabel = "Same_Side_Eta_Width_";
    TString sigPhiLabel = "Same_Side_Phi_Width";
    
    TString lightFlavorLabel = "Light_Flavor_";
    
    TString awaySideAmpLabel = "AwaySideAmp_";
    TString awaySidePhiWidthLabel = "AwaySide_Phi_Width_";
    
    TString pythiaLabel = "Pythia_";
    
    TString leftLabel = "_left_";
    TString rightLabel = "_right_";
	
	TString phiProj = "Phi_projection_";
	TString str1;
   
	TH2D* fullySubtractedCorrCustomCentralitySideBand[5][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandPhiProj[5][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandPhiProjLeft[5][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandPhiProjRight[5][3];
	TH1D* fullySubtractedCorrCustomCentralitySideBandFitPhiProj[5][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandResPhiProj[5][3];
   
    TH1D* fitParameterOffsetPlotVsCent[5];
    TH1D* fitParameterDipolePlotVsCent[5];
    TH1D* fitParameterQuadrupolePlotVsCent[5];
    TH1D* fitParameterJetAmpPlotVsCent[5];
    TH1D* fitParameterSigmaPhiPlotVsCent[5];
    TH1D* fitParameterSigmaEtaPlotVsCent[5];
    TH1D* fitParameterASAmpPlotVsCent[5];
    TH1D* fitParameterASSigmaPhiPlotVsCent[5];
    
    TH1D* lightFlavorResultsJetAmpPt1[5];
    TH1D* lightFlavorResultsEtaWidthPt1[5];
    TH1D* lightFlavorResultsPhiWidthPt1[5];
    TH1D* lightFlavorResultsJetAmpPt2[5];
    TH1D* lightFlavorResultsEtaWidthPt2[5];
    TH1D* lightFlavorResultsPhiWidthPt2[5];
    TH1D* lightFlavorResultsJetAmpPt3[5];
    TH1D* lightFlavorResultsEtaWidthPt3[5];
    TH1D* lightFlavorResultsPhiWidthPt3[5];
    TH1D* lightFlavorResultsJetAmpPt4[5];
    TH1D* lightFlavorResultsEtaWidthPt4[5];
    TH1D* lightFlavorResultsPhiWidthPt4[5];
    TH1D* lightFlavorResultsJetAmpPt5[5];
    TH1D* lightFlavorResultsEtaWidthPt5[5];
    TH1D* lightFlavorResultsPhiWidthPt5[5];
    TH1D* lightFlavorResultsJetAmpPt6[5];
    TH1D* lightFlavorResultsEtaWidthPt6[5];
    TH1D* lightFlavorResultsPhiWidthPt6[5];
    
    TH1D* pythiaResultsJetAmp[5];
    TH1D* pythiaResultsEtaWidth[5];
    TH1D* pythiaResultsPhiWidth[5];
    
    TH1D* fitParameterPublishedQuadVsCent[5];
    
    TH1D * jetVolumePlot[5];
    
   // TH2D* errorOnBin[3];
    TH1D* sigmaPerBin[5][3];
    
   
    
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
    
        
        str1 = quadLabel + ptLabel + ptBinLabel[k];
        fitParameterQuadrupolePlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        fitParameterPublishedQuadVsCent[k] = new TH1D("","",3, 0, 3);
        str1 = jetAmpLabel + ptLabel + ptBinLabel[k];
        fitParameterJetAmpPlotVsCent[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = sigPhiLabel + ptLabel + ptBinLabel[k];
        fitParameterSigmaPhiPlotVsCent[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = sigEtaLabel + ptLabel + ptBinLabel[k];
        fitParameterSigmaEtaPlotVsCent[k] = new TH1D(str1, str1, 4, 0, 4);
       
        str1 = lightFlavorLabel + lightFlavorPtLabel[0] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt1[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[0] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt1[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[0] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt1[k] = new TH1D(str1, str1, 4, 0, 4);
        
        str1 = lightFlavorLabel + lightFlavorPtLabel[1] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt2[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[1] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt2[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[1] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt2[k] = new TH1D(str1, str1, 4, 0, 4);
        
        str1 = lightFlavorLabel + lightFlavorPtLabel[2] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt3[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[2] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt3[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[2] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt3[k] = new TH1D(str1, str1, 4, 0, 4);
        
        str1 = lightFlavorLabel + lightFlavorPtLabel[3] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt4[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[3] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt4[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[3] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt4[k] = new TH1D(str1, str1, 4, 0, 4);
        
        str1 = lightFlavorLabel + lightFlavorPtLabel[4] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt5[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[4] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt5[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[4] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt5[k] = new TH1D(str1, str1, 4, 0, 4);
        
        str1 = lightFlavorLabel + lightFlavorPtLabel[5] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt6[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[5] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt6[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[5] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt6[k] = new TH1D(str1, str1, 4, 0, 4);
        
        str1 = pythiaLabel + jetAmpLabel + ptLabel + ptBinLabel[k];
        pythiaResultsJetAmp[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = pythiaLabel + sigEtaLabel + ptLabel + ptBinLabel[k];
        pythiaResultsEtaWidth[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = pythiaLabel + sigPhiLabel + ptLabel + ptBinLabel[k];
        pythiaResultsPhiWidth[k] = new TH1D(str1, str1, 4, 0, 4);
        
        pythiaResultsJetAmp[k]->GetXaxis()->SetLabelSize(.06);
        pythiaResultsJetAmp[k]->GetXaxis()->SetBinLabel(1, "pythia");
        pythiaResultsJetAmp[k]->GetXaxis()->SetBinLabel(2, "50-100%");
        pythiaResultsJetAmp[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        pythiaResultsJetAmp[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        pythiaResultsEtaWidth[k]->GetXaxis()->SetLabelSize(.06);
        pythiaResultsEtaWidth[k]->GetXaxis()->SetBinLabel(1, "pythia");
        pythiaResultsEtaWidth[k]->GetXaxis()->SetBinLabel(2, "50-100%");
        pythiaResultsEtaWidth[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        pythiaResultsEtaWidth[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        pythiaResultsPhiWidth[k]->GetXaxis()->SetLabelSize(.06);
        pythiaResultsPhiWidth[k]->GetXaxis()->SetBinLabel(1, "pythia");
        pythiaResultsPhiWidth[k]->GetXaxis()->SetBinLabel(2, "50-100%");
        pythiaResultsPhiWidth[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        pythiaResultsPhiWidth[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        
        
        fitParameterQuadrupolePlotVsCent[k]->SetLabelSize(.06);
        fitParameterQuadrupolePlotVsCent[k]->GetYaxis()->SetTitle("Quadrupole Amplitude, A_{Q}");
        fitParameterQuadrupolePlotVsCent[k]->GetYaxis()->SetTitleSize(.045);
        fitParameterQuadrupolePlotVsCent[k]->GetYaxis()->SetTitleOffset(1.5);
        fitParameterQuadrupolePlotVsCent[k]->SetTitle("");
        
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetTitle("Centrality (%)");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetTitleOffset(1.3);
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetLabelSize(.06);
        fitParameterJetAmpPlotVsCent[k]->GetYaxis()->SetLabelSize(.05);
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "pythia");
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "50-100%");
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetLabelSize(.06);
        fitParameterSigmaEtaPlotVsCent[k]->GetYaxis()->SetLabelSize(.04);
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "pythia");
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "50-100%");
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetTitle("Centrality (%)");
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetTitleOffset(1.3);
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetLabelSize(.06);
        fitParameterSigmaPhiPlotVsCent[k]->GetYaxis()->SetLabelSize(.04);
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "pythia");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "50-100%");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetTitle("Centrality (%)");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetTitleOffset(1.3);
        
        lightFlavorResultsJetAmpPt1[k]->GetXaxis()->SetLabelSize(.06);
        lightFlavorResultsJetAmpPt1[k]->GetYaxis()->SetLabelSize(.05);
        lightFlavorResultsJetAmpPt1[k]->GetXaxis()->SetBinLabel(1, "pythia");
        lightFlavorResultsJetAmpPt1[k]->GetXaxis()->SetBinLabel(2, "50-100%");
        lightFlavorResultsJetAmpPt1[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        lightFlavorResultsJetAmpPt1[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        lightFlavorResultsEtaWidthPt1[k]->GetXaxis()->SetLabelSize(.07);
        lightFlavorResultsEtaWidthPt1[k]->GetYaxis()->SetLabelSize(.05);
        lightFlavorResultsEtaWidthPt1[k]->GetXaxis()->SetBinLabel(1, "pythia");
        lightFlavorResultsEtaWidthPt1[k]->GetXaxis()->SetBinLabel(2, "50-100%");
        lightFlavorResultsEtaWidthPt1[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        lightFlavorResultsEtaWidthPt1[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        lightFlavorResultsPhiWidthPt1[k]->GetXaxis()->SetLabelSize(.06);
        lightFlavorResultsPhiWidthPt1[k]->GetYaxis()->SetLabelSize(.05);
        lightFlavorResultsPhiWidthPt1[k]->GetXaxis()->SetBinLabel(1, "pythia");
        lightFlavorResultsPhiWidthPt1[k]->GetXaxis()->SetBinLabel(2, "50-100%");
        lightFlavorResultsPhiWidthPt1[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        lightFlavorResultsPhiWidthPt1[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        
        
        
        fitParameterJetAmpPlotVsCent[k]->SetStats(0);
        fitParameterSigmaEtaPlotVsCent[k]->SetStats(0);
        fitParameterSigmaPhiPlotVsCent[k]->SetStats(0);
        lightFlavorResultsJetAmpPt1[k]->SetStats(0);
        lightFlavorResultsEtaWidthPt1[k]->SetStats(0);
        lightFlavorResultsPhiWidthPt1[k]->SetStats(0);
        fitParameterQuadrupolePlotVsCent[k]->SetStats(0);
        
        pythiaResultsJetAmp[k]->SetStats(0);
        pythiaResultsEtaWidth[k]->SetStats(0);
        pythiaResultsPhiWidth[k]->SetStats(0);
    }

    TString binLabels[11] = {"PYTHIA", "95%", "85%", "75%", "65%", "55%", "45%", "35%", "25%", "15%", "5%"};
    
    TH1D* fullPhiWidthPlot = new TH1D("", "", 11, 1, 12);
    
    for(int i = 0; i < 11; i++){
    
        fullPhiWidthPlot->GetXaxis()->SetBinLabel(i+1, binLabels[i]);
    }
    
    
    
                             // PYTHIA, peripheral, mid-central, central
    double centrality[4]       = {1.5, 5.5, 8.5, 11.5};  
    double phiWidthAlexData[4] = {.310772, .349718, .669247, .754439};    
    
    //double phiErrors[4]        = {.0343975, 
     
    double eCentLow[4] = {.5, 2.5, 1.5, 1.5};
    double eCentHigh[4] = {.5, 1.5, 1.5, .5}; 
   
    double ephiWidthAlexLow[4] = { 0.05, 0.05, 0.05, 0.05};
    double ephiWidthAlexHigh[4] = { 0.05, 0.05, 0.05, 0.05};
    
    
  
    TGraphAsymmErrors *AlexData = new TGraphAsymmErrors(4, centrality, phiWidthAlexData, eCentLow, eCentHigh, ephiWidthAlexLow, ephiWidthAlexHigh);
  
    AlexData->SetMarkerColor(2);
    AlexData->SetMarkerStyle(21);
    AlexData->SetMarkerSize(2);
    AlexData->SetLineColor(2);
    AlexData->SetLineWidth(3);
    
    TCanvas * testing = new TCanvas("yourmom", "", 1200, 1100);
    
    fullPhiWidthPlot->SetStats(0);
    fullPhiWidthPlot->Draw();
    AlexData->Draw("SAME P");
    
    testing->SaveAs("testing.png");

   
    
}


double getJetYieldError(double A, double sigEta, double sigPhi, double eA, double eSigEta, double eSigPhi){

    cout << endl <<"error on A, sigPhi, sigEta - " << eA << " , " << eSigPhi << " , " << eSigEta << endl;

    double arg = ((eA/A)*(eA/A)) + ((eSigEta/sigEta)*(eSigEta/sigEta)) + ((eSigPhi/sigPhi)*(eSigPhi/sigPhi));
    double product = TMath::Abs((A*sigEta*sigPhi));
    
    //cout << "prduct value: " << product << endl;
    //cout << "arg value: " << arg << endl;
    //cout << "sqrt of arg: " << TMath::Sqrt(arg) << endl;
    
    return TMath::Abs((A*sigEta*sigPhi))*TMath::Sqrt(arg);
    
}








