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

plotFinalFitResults(){


    bool USE_ETA_GAP = false;
    
    

    double NDF = 1;
    
    double projectionScaleFactor;
    
   
    int NUM_PT_BINS = 3;
    int FIRST_PT_BIN = 2;
    int NUM_PHI_BINS = 12;
    int NUM_ETA_BINS = 13;
    double ETA_RANGE = 1.05;
    double phiBinShift = (TMath::Pi()/12.0);
   
    //double nBarPrime[3] = {74.4, 351, 1027};
    
    double nBarPrime[3] = {80, 351, 1027};
   
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
    
		double xbinsCent[4] = {0.0, 1.0, 2.0, 2.67};
		double xbinsCentWithPythia[5] = {0.0, 1.0, 2.0, 3.0, 3.67};
	
        str1 = offsetLabel + ptLabel + ptBinLabel[k];
        fitParameterOffsetPlotVsCent[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = dipoleLabel + ptLabel + ptBinLabel[k];
        fitParameterDipolePlotVsCent[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = quadLabel + ptLabel + ptBinLabel[k];
		
		
		
        fitParameterQuadrupolePlotVsCent[k] = new TH1D(str1, str1, 3, xbinsCent);
        fitParameterPublishedQuadVsCent[k] = new TH1D("","",3, 0, 3);
        str1 = jetAmpLabel + ptLabel + ptBinLabel[k];
        fitParameterJetAmpPlotVsCent[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = sigPhiLabel + ptLabel + ptBinLabel[k];
        fitParameterSigmaPhiPlotVsCent[k] = new TH1D(str1, str1, 4, xbinsCentWithPythia);
        str1 = sigEtaLabel + ptLabel + ptBinLabel[k];
        fitParameterSigmaEtaPlotVsCent[k] = new TH1D(str1, str1, 4, xbinsCentWithPythia);
        
        str1 = awaySideAmpLabel + ptLabel + ptBinLabel[k];
        fitParameterASAmpPlotVsCent[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = awaySidePhiWidthLabel + ptLabel + ptBinLabel[k];
        fitParameterASSigmaPhiPlotVsCent[k] = new TH1D(str1, str1, 4, 0, 4);
        
        str1 = lightFlavorLabel + lightFlavorPtLabel[0] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt1[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[0] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt1[k] = new TH1D(str1, str1, 4, xbinsCentWithPythia);
        str1 = lightFlavorLabel + lightFlavorPtLabel[0] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt1[k] = new TH1D(str1, str1, 4, xbinsCentWithPythia);
        
        str1 = lightFlavorLabel + lightFlavorPtLabel[1] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt2[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[1] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt2[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[1] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt2[k] = new TH1D(str1, str1, 4, 0, 4);
        
        str1 = lightFlavorLabel + lightFlavorPtLabel[2] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt3[k] = new TH1D(str1, str1, 4, 0, 4);
        str1 = lightFlavorLabel + lightFlavorPtLabel[2] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt3[k] = new TH1D(str1, str1, 4, xbinsCentWithPythia);
        str1 = lightFlavorLabel + lightFlavorPtLabel[2] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt3[k] = new TH1D(str1, str1, 4, xbinsCentWithPythia);
        
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
        pythiaResultsEtaWidth[k] = new TH1D(str1, str1, 4, xbinsCentWithPythia);
        str1 = pythiaLabel + sigPhiLabel + ptLabel + ptBinLabel[k];
        pythiaResultsPhiWidth[k] = new TH1D(str1, str1, 4, xbinsCentWithPythia);
        
        pythiaResultsJetAmp[k]->GetXaxis()->SetLabelSize(.06);
        pythiaResultsJetAmp[k]->GetXaxis()->SetBinLabel(1, "pythia");
        pythiaResultsJetAmp[k]->GetXaxis()->SetBinLabel(2, "50-80%");
        pythiaResultsJetAmp[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        pythiaResultsJetAmp[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        pythiaResultsEtaWidth[k]->GetXaxis()->SetLabelSize(.06);
        pythiaResultsEtaWidth[k]->GetXaxis()->SetBinLabel(1, "pythia");
        pythiaResultsEtaWidth[k]->GetXaxis()->SetBinLabel(2, "50-80%");
        pythiaResultsEtaWidth[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        pythiaResultsEtaWidth[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        pythiaResultsPhiWidth[k]->GetXaxis()->SetLabelSize(.06);
        pythiaResultsPhiWidth[k]->GetXaxis()->SetBinLabel(1, "pythia");
        pythiaResultsPhiWidth[k]->GetXaxis()->SetBinLabel(2, "50-80%");
        pythiaResultsPhiWidth[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        pythiaResultsPhiWidth[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetLabelSize(.06);
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "pythia");
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "50-80%");
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetLabelSize(.07);
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetTitle("Centrality (%)");
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetTitleOffset(1.3);
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetBinLabel(1, "pythia");
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetBinLabel(2, "50-80%");
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        fitParameterQuadrupolePlotVsCent[k]->SetLabelSize(.06);
        fitParameterQuadrupolePlotVsCent[k]->GetYaxis()->SetTitle("Quadrupole Amplitude, A_{Q}");
        fitParameterQuadrupolePlotVsCent[k]->GetYaxis()->SetTitleSize(.045);
        fitParameterQuadrupolePlotVsCent[k]->GetYaxis()->SetTitleOffset(1.5);
        fitParameterQuadrupolePlotVsCent[k]->SetTitle("");
        
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-80%");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetTitle("Centrality (%)");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetTitleOffset(1.3);
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetLabelSize(.06);
        fitParameterJetAmpPlotVsCent[k]->GetYaxis()->SetLabelSize(.05);
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "pythia");
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "50-80%");
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetLabelSize(.06);
        fitParameterSigmaEtaPlotVsCent[k]->GetYaxis()->SetLabelSize(.04);
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "pythia");
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "50-80%");
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetTitle("Centrality (%)");
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetTitleOffset(1.3);
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetLabelSize(.06);
        fitParameterSigmaPhiPlotVsCent[k]->GetYaxis()->SetLabelSize(.04);
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "pythia");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "50-80%");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetTitle("Centrality (%)");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetTitleOffset(1.3);
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetLabelSize(.06);
        fitParameterASAmpPlotVsCent[k]->GetYaxis()->SetLabelSize(.05);
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-80%");
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        
        fitParameterASSigmaPhiPlotVsCent[k]->GetXaxis()->SetLabelSize(.07);
        fitParameterASSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-80%");
        fitParameterASSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterASSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        
        lightFlavorResultsJetAmpPt1[k]->GetXaxis()->SetLabelSize(.06);
        lightFlavorResultsJetAmpPt1[k]->GetYaxis()->SetLabelSize(.05);
        lightFlavorResultsJetAmpPt1[k]->GetXaxis()->SetBinLabel(1, "pythia");
        lightFlavorResultsJetAmpPt1[k]->GetXaxis()->SetBinLabel(2, "50-80%");
        lightFlavorResultsJetAmpPt1[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        lightFlavorResultsJetAmpPt1[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        lightFlavorResultsEtaWidthPt1[k]->GetXaxis()->SetLabelSize(.07);
        lightFlavorResultsEtaWidthPt1[k]->GetYaxis()->SetLabelSize(.05);
        lightFlavorResultsEtaWidthPt1[k]->GetXaxis()->SetBinLabel(1, "pythia");
        lightFlavorResultsEtaWidthPt1[k]->GetXaxis()->SetBinLabel(2, "50-80%");
        lightFlavorResultsEtaWidthPt1[k]->GetXaxis()->SetBinLabel(3, "20-50%");
        lightFlavorResultsEtaWidthPt1[k]->GetXaxis()->SetBinLabel(4, "0-20%");
        lightFlavorResultsPhiWidthPt1[k]->GetXaxis()->SetLabelSize(.06);
        lightFlavorResultsPhiWidthPt1[k]->GetYaxis()->SetLabelSize(.05);
        lightFlavorResultsPhiWidthPt1[k]->GetXaxis()->SetBinLabel(1, "pythia");
        lightFlavorResultsPhiWidthPt1[k]->GetXaxis()->SetBinLabel(2, "50-80%");
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



    //Hard-Code numbers here
    
    //fitParameterPublishedQuadVsCent[2]->SetBinContent(1, 0);
    //fitParameterPublishedQuadVsCent[2]->SetBinError(1, .000860232);
    fitParameterPublishedQuadVsCent[2]->SetBinContent(2, .00816);
    fitParameterPublishedQuadVsCent[2]->SetBinContent(1, -1000);
    fitParameterPublishedQuadVsCent[2]->SetBinContent(3, -1000);
    //fitParameterPublishedQuadVsCent[2]->SetBinError(2, .001);
    //fitParameterPublishedQuadVsCent[2]->SetBinContent(3, 0);
    
    /* ///TEMPORARY NUMBERS--------------------
    
    
    //temporary numbers from different fits
    lightFlavorResultsJetAmpPt1[2]->SetBinContent(1, .119999); 
    lightFlavorResultsJetAmpPt1[2]->SetBinError(1,  .0597798);
    lightFlavorResultsJetAmpPt1[2]->SetBinContent(2, .043);
    lightFlavorResultsJetAmpPt1[2]->SetBinContent(3, .033);
    
    lightFlavorResultsEtaWidthPt1[2]->SetBinContent(1, .262319);
    lightFlavorResultsEtaWidthPt1[2]->SetBinError(1, .0348749);
    lightFlavorResultsEtaWidthPt1[2]->SetBinContent(2, 1.41);
    lightFlavorResultsEtaWidthPt1[2]->SetBinContent(3, 2.24);
    
    lightFlavorResultsPhiWidthPt1[2]->SetBinContent(1, .216455);
    lightFlavorResultsPhiWidthPt1[2]->SetBinError(1, .0190255);
    lightFlavorResultsPhiWidthPt1[2]->SetBinContent(2, .53);
    lightFlavorResultsPhiWidthPt1[2]->SetBinContent(3, .57);
     */
    //---------------------------------------------
    
    //----------------------pythia values------------------------
    
    pythiaResultsJetAmp[2]->SetBinContent(1, 2.21967);
    pythiaResultsJetAmp[2]->SetBinError(1, .764791);
       
    pythiaResultsEtaWidth[2]->SetBinContent(1, .330681);
    pythiaResultsEtaWidth[2]->SetBinError(1,    .0310280);
        
    pythiaResultsPhiWidth[2]->SetBinContent(1, .310772);
    pythiaResultsPhiWidth[2]->SetBinError(1,   .0343975);
    
    
    
    //------------------------------------------------------------------
    //------------------from symmetric space with correctly scaled errors------------------------
    //------------------------------------------------------------------
    
    
    //nominal values
    
    /*fitParameterQuadrupolePlotVsCent[2]->SetBinContent(1, .00402219);  //-- from AS 1D Gaussian fit
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(1, .00249265);   //-- from AS 1D Gaussian fit
	//fitParameterQuadrupolePlotVsCent[2]->SetBinContent(1, .00867822);   //from 6-parameter
    //fitParameterQuadrupolePlotVsCent[2]->SetBinError(1, .000925368);    //from 6-parameter
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(2, .0905741);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(2, .0188205);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(2, .349718);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(2, .0644555);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(2, .312136);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(2, .113712);*/
    
    fitParameterQuadrupolePlotVsCent[2]->SetBinContent(1, .00463036);  //-- from AS 1D Gaussian fit
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(1, .00268507);   //-- from AS 1D Gaussian fit
	//fitParameterQuadrupolePlotVsCent[2]->SetBinContent(1, .00867822);   //from 6-parameter
    //fitParameterQuadrupolePlotVsCent[2]->SetBinError(1, .000925368);    //from 6-parameter
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(2, .0851668);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(2, .0189719);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(2, .355024);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(2, .0675623);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(2, .306846);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(2, .112368);
    
   
   
    fitParameterQuadrupolePlotVsCent[2]->SetBinContent(2, .00647402);
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(2, .00320154);
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(3, .0389386);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(3, .00344479);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(3, .669247);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(3, .0618335);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(3, 1.40393);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(3, .353493);
    
   
    
    /*fitParameterQuadrupolePlotVsCent[2]->SetBinContent(3, .00264);
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(3, .000272316);
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(3, .00751125);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(3, .00174301);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(3, .433461);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(3, .0522924);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(3, 1.8);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(3, .258073);*/ //non-zero quadrupole
    
    fitParameterQuadrupolePlotVsCent[2]->SetBinContent(3, 0.0);
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(3,   .000866073);
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(4,     .0439447);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(4,        .00666497);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(4,    .754439);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(4,       .0807210);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(4,     1.24592);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(4,       .221779);
    
    
    
    
    //////////pt = 2.56
    //lightFlavorResultsJetAmpPt1[2]->SetBinContent(2, .653);
    lightFlavorResultsJetAmpPt1[2]->SetBinContent(3, 4.76);
    lightFlavorResultsJetAmpPt1[2]->SetBinContent(4, 12.17);
    
    //lightFlavorResultsEtaWidthPt1[2]->SetBinContent(2, .57);
    lightFlavorResultsEtaWidthPt1[2]->SetBinContent(3, 1.41);
    lightFlavorResultsEtaWidthPt1[2]->SetBinContent(4, 2.24);
	//lightFlavorResultsEtaWidthPt1[2]->SetBinError(2, .57*.05);
    lightFlavorResultsEtaWidthPt1[2]->SetBinError(3, 1.41*.09);
    lightFlavorResultsEtaWidthPt1[2]->SetBinError(4, 2.24*.09);
    
    //lightFlavorResultsPhiWidthPt1[2]->SetBinContent(2, .43);
    lightFlavorResultsPhiWidthPt1[2]->SetBinContent(3, .53);
    lightFlavorResultsPhiWidthPt1[2]->SetBinContent(4, .57);
	//lightFlavorResultsPhiWidthPt1[2]->SetBinError(2, .43*.05);
    lightFlavorResultsPhiWidthPt1[2]->SetBinError(3, .53*.05);
    lightFlavorResultsPhiWidthPt1[2]->SetBinError(4, .57*.05);
   
    ///////////////////pt = 5.7
    //lightFlavorResultsJetAmpPt3[2]->SetBinContent(2, 1.103);
    lightFlavorResultsJetAmpPt3[2]->SetBinContent(3, 2.08);
    lightFlavorResultsJetAmpPt3[2]->SetBinContent(4, 6.34);
    
   // lightFlavorResultsEtaWidthPt3[2]->SetBinContent(2, .49);
    lightFlavorResultsEtaWidthPt3[2]->SetBinContent(3, .63);
    lightFlavorResultsEtaWidthPt3[2]->SetBinContent(4, 1.08);
	//lightFlavorResultsEtaWidthPt3[2]->SetBinError(2, .49*.1);
    lightFlavorResultsEtaWidthPt3[2]->SetBinError(3, .63*.15);
    lightFlavorResultsEtaWidthPt3[2]->SetBinError(4, 1.08*.17);
    
    //lightFlavorResultsPhiWidthPt3[2]->SetBinContent(2, .48);
    lightFlavorResultsPhiWidthPt3[2]->SetBinContent(3, .37);
    lightFlavorResultsPhiWidthPt3[2]->SetBinContent(4, .44);
	//lightFlavorResultsPhiWidthPt3[2]->SetBinError(2, .48*.12);
    lightFlavorResultsPhiWidthPt3[2]->SetBinError(3, .37*.11);
    lightFlavorResultsPhiWidthPt3[2]->SetBinError(4, .44*.1);
    
    double centrality[1]  = {1.67};  
    double eCentLow[1]    = {.33};
    double eCentHigh[1]   = {.33}; 
    
    
    double JetYieldKettlerPt1[1] = {.653};    
    double eJetYieldKettlerPt1Low[4] = {.653*.15};
    double eJetYieldKettlerPt1High[4] = {.653*.15};
    
    double phiWidthKettlerPt1[1] = {.43};    
    double ePhiWidthKettlerPt1Low[1] = {.43*.05};
    double ePhiWidthKettlerPt1High[1] = {.43*.05};
    
    double etaWidthKettlerPt1[1] = {.57};    
    double eEtaWidthKettlerPt1Low[1] = {.57*.05};
    double eEtaWidthKettlerPt1High[1] = {.57*.05};
    
    double JetYieldKettlerPt3[1] = {1.103};    
    double eJetYieldKettlerPt3Low[1] = {1.103*.15};
    double eJetYieldKettlerPt3High[1] = {1.103*.15};
    
    double phiWidthKettlerPt3[1] = {.48};    
    double ePhiWidthKettlerPt3Low[1] = {.48*.12};
    double ePhiWidthKettlerPt3High[1] = {.48*.12};
    
    double etaWidthKettlerPt3[1] = {.49};    
    double eEtaWidthKettlerPt3Low[1] = {.49*.1};
    double eEtaWidthKettlerPt3High[1] = {.49*.1};
    
    //double centrality[1]  = {1.67};  
    //double eCentLow[1]    = {.33};
    //double eCentHigh[1]   = {.33};
    
    //double kettlerQuadPeripheralPt1[1] = {.0145};    
    //double ekettlerQuadPeripheralPt1Low[1] = {.0005};
    //double ekettlerQuadPeripheralPt1High[1] = {.0005};
    
    //double centralityMidCentral[1]  = {2.165, 2.495, 2.825};  
    //double eCentMidCentLow[1]    = {.165,.165,.165};
    //double eCentMidCentHigh[1]   = {.165,.165,.165};
    
    //double kettlerQuadMidCentralPt1[1] = {.0135, .012, .008};    
    //double ekettlerQuadMidCentralPt1Low[1] = {.0005, .0005, .0005};
    //double ekettlerQuadMidCentralPt1High[1] = {.0005, .0005, .0005};
    
	//Quadrupole data
	
	fitParameterQuadrupolePlotVsCent[2]->SetBinContent(1, .00463036);  //-- from AS 1D Gaussian fit
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(1, .00268507);   //-- from AS 1D Gaussian fit
	fitParameterQuadrupolePlotVsCent[2]->SetBinContent(2, .00647402);
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(2, .00320154);
	fitParameterQuadrupolePlotVsCent[2]->SetBinContent(3, 0.0);
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(3,   .000866073);
	
	double alexDataQuadrupole[3] = {.00463036, .00647402, 0.0};
	double ealexDataQuadrupole[3] = {.00268507, .00320154 , .000866073};
	//double alexDataQuadrupole[3] = {.00463036, .00647402, 0.0};
	
	double alexQuadCentrality[3] = {4, 10, 14.5};
	double alexQuadCentralityLow[3] = {3, 3, 2.5};
	double alexQuadCentralityHigh[3] = {3, 3, 2.5};
    //           2	2   2   2   2   2   2   2   2   2   
	//PYTHIA, 100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0
	
	TH1D* quadTemplateHistogram = new TH1D("","", 16, 1, 17);
	quadTemplateHistogram->GetXaxis()->SetBinLabel(1, "80");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(2, "");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(3, "70");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(4, "");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(5, "60");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(6, "");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(7, "50");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(8, "");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(9, "40");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(10, "");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(11, "30");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(12, "");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(13, "20");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(14, "");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(15, "10");
	quadTemplateHistogram->GetXaxis()->SetBinLabel(16, "");
	
	
	
    double kettlerCentrality[7] = {.67, 1.165, 1.495, 1.825, 2.1675, 2.41875, 2.58625};
    double kettlerCentLow[7]    = {.33, .165, .165, .165, .1675, .08375, .08375};
    double kettlerCentHigh[7]   = {.33, .165, .165, .165, .1675, .08375, .08375}; 
    
    double kettlerQuadValuesPt1[7]      = {.0145, .0135, .012, .008, .0038, .001, 0.0};    
    double ekettlerQuadValuesPt1Low[7]  = {.0005, .0005,.0005, .0005, .0005, .0005,.0005};
    double ekettlerQuadValuesPt1High[7] = {.0005, .0005,.0005, .0005, .0005, .0005,.0005};
    
    double kettlerQuadValuesPt3[7]      = {.0125, .0135, .0116, .0095, .0057, 0.0, 0.0};    
    double ekettlerQuadValuesPt3Low[7]  = {.0015, .0015,.0015, .0015, .0015,.0015, .0015};
    double ekettlerQuadValuesPt3High[7] = {.0015, .0015,.0015, .0015, .0015,.0015, .0015};
    
    TGraphAsymmErrors *kettlerPeripheralPt1JetYield = new TGraphAsymmErrors(1, centrality, JetYieldKettlerPt1, eCentLow, eCentHigh, eJetYieldKettlerPt1Low, eJetYieldKettlerPt1High);
    TGraphAsymmErrors *kettlerPeripheralPt1Phi = new TGraphAsymmErrors(1, centrality, phiWidthKettlerPt1, eCentLow, eCentHigh, ePhiWidthKettlerPt1Low, ePhiWidthKettlerPt1High);
    TGraphAsymmErrors *kettlerPeripheralPt1Eta = new TGraphAsymmErrors(1, centrality, etaWidthKettlerPt1, eCentLow, eCentHigh, eEtaWidthKettlerPt1Low, eEtaWidthKettlerPt1High);
    
    TGraphAsymmErrors *kettlerPt1Quad = new TGraphAsymmErrors(7, kettlerCentrality, kettlerQuadValuesPt1, kettlerCentLow, kettlerCentHigh, ekettlerQuadValuesPt1Low, ekettlerQuadValuesPt1High);
    
   
    TGraphAsymmErrors *kettlerPeripheralPt3JetYield = new TGraphAsymmErrors(1, centrality, JetYieldKettlerPt3, eCentLow, eCentHigh, eJetYieldKettlerPt3Low, eJetYieldKettlerPt3High);
    TGraphAsymmErrors *kettlerPeripheralPt3Phi = new TGraphAsymmErrors(1, centrality, phiWidthKettlerPt3, eCentLow, eCentHigh, ePhiWidthKettlerPt3Low, ePhiWidthKettlerPt3High);
    TGraphAsymmErrors *kettlerPeripheralPt3Eta = new TGraphAsymmErrors(1, centrality, etaWidthKettlerPt3, eCentLow, eCentHigh, eEtaWidthKettlerPt3Low, eEtaWidthKettlerPt3High);
  
    TGraphAsymmErrors *kettlerPt3Quad = new TGraphAsymmErrors(7, kettlerCentrality, kettlerQuadValuesPt3, kettlerCentLow, kettlerCentHigh, ekettlerQuadValuesPt3Low, ekettlerQuadValuesPt3High);
  
    //TGraphAsymmErrors *kettlerPeripheralPt3Quad = new TGraphAsymmErrors(1, centrality, etaWidthKettlerPt1, eCentLow, eCentHigh, eEtaWidthKettlerPt1Low, eEtaWidthKettlerPt1High);
    //TGraphAsymmErrors *kettlerMidCentralPt3Quad = new TGraphAsymmErrors(1, centrality, etaWidthKettlerPt1, eCentLow, eCentHigh, eEtaWidthKettlerPt1Low, eEtaWidthKettlerPt1High);
    //TGraphAsymmErrors *kettlerCentralPt3Quad = new TGraphAsymmErrors(1, centrality, etaWidthKettlerPt1, eCentLow, eCentHigh, eEtaWidthKettlerPt1Low, eEtaWidthKettlerPt1High);
    
	double kettlerEtaWidthValuesPt1[8] = {.5, .7, 1.0, 1.5, 1.9, 2.2, 2.5, 2.2};
	double ekettlerEtaWidthValuesPt1Low[8] = {.025, .025, .03, .05, .05, .1, .2, .2}; 
	double ekettlerEtaWidthValuesPt1High[8] = {.025, .025, .03, .05, .05, .1, .2, .2}; 
    
    TH1D* kettlerJetYieldResults[3];
	kettlerJetYieldResults[0] = new TH1D("","", 4, xbinsCentWithPythia); //1.15
	kettlerJetYieldResults[1] = new TH1D("","", 4, xbinsCentWithPythia); //2.56
	kettlerJetYieldResults[2] = new TH1D("","", 4, xbinsCentWithPythia); //5.7
	
	//kettlerJetYieldResults[0]->SetBinContent(2, .630);
	kettlerJetYieldResults[0]->SetBinContent(3, 3.45);
	kettlerJetYieldResults[0]->SetBinContent(4, 8.32);
	//kettlerJetYieldResults[1]->SetBinContent(2, .653);
	kettlerJetYieldResults[1]->SetBinContent(3, 4.76);
	kettlerJetYieldResults[1]->SetBinContent(4, 12.17);
	//kettlerJetYieldResults[2]->SetBinContent(2, 1.103);
	kettlerJetYieldResults[2]->SetBinContent(3, 2.08);
	kettlerJetYieldResults[2]->SetBinContent(4, 6.34);
    
    //kettlerJetYieldResults[0]->SetBinError(2, .630*.15);
	kettlerJetYieldResults[0]->SetBinError(3, 3.45*.15);
	kettlerJetYieldResults[0]->SetBinError(4, 8.32*.15);
	//kettlerJetYieldResults[1]->SetBinError(2, .653*.15);
	kettlerJetYieldResults[1]->SetBinError(3, 4.76*.15);
	kettlerJetYieldResults[1]->SetBinError(4, 12.17*.15);
	//kettlerJetYieldResults[2]->SetBinError(2, 1.103*.15);
	kettlerJetYieldResults[2]->SetBinError(3, 2.08*.15);
	kettlerJetYieldResults[2]->SetBinError(4, 6.34*.15);
    
    kettlerPeripheralPt1JetYield->SetMarkerStyle(9);
    kettlerPeripheralPt1JetYield->SetMarkerSize(2);
    kettlerPeripheralPt1JetYield->SetMarkerColor(38);
    kettlerPeripheralPt1JetYield->SetLineColor(38);
    kettlerPeripheralPt1JetYield->SetLineWidth(3);
    
    kettlerPt1Quad->SetMarkerStyle(9);
    kettlerPt1Quad->SetMarkerSize(2);
    kettlerPt1Quad->SetMarkerColor(38);
    kettlerPt1Quad->SetLineColor(38);
    kettlerPt1Quad->SetLineWidth(3);
    
   
    kettlerPeripheralPt3JetYield->SetMarkerStyle(9);
    kettlerPeripheralPt3JetYield->SetMarkerSize(2);
    kettlerPeripheralPt3JetYield->SetMarkerColor(4);
    kettlerPeripheralPt3JetYield->SetLineColor(4);
    kettlerPeripheralPt3JetYield->SetLineWidth(3);
    
    kettlerPt3Quad->SetMarkerStyle(9);
    kettlerPt3Quad->SetMarkerSize(2);
    kettlerPt3Quad->SetMarkerColor(4);
    kettlerPt3Quad->SetLineColor(4);
    kettlerPt3Quad->SetLineWidth(3);
    
    
    kettlerJetYieldResults[0]->SetMarkerStyle(9);
	kettlerJetYieldResults[0]->SetMarkerSize(2);
	kettlerJetYieldResults[0]->SetMarkerColor(8);
	kettlerJetYieldResults[1]->SetMarkerStyle(9);
	kettlerJetYieldResults[1]->SetMarkerSize(2);
	kettlerJetYieldResults[1]->SetMarkerColor(38);
	kettlerJetYieldResults[2]->SetMarkerStyle(9);
	kettlerJetYieldResults[2]->SetMarkerSize(2);
	kettlerJetYieldResults[2]->SetMarkerColor(4);
    
    kettlerJetYieldResults[0]->SetLineColor(8);
    kettlerJetYieldResults[0]->SetLineWidth(3);
    kettlerJetYieldResults[1]->SetLineColor(38);
    kettlerJetYieldResults[1]->SetLineWidth(3);
	kettlerJetYieldResults[2]->SetLineColor(4);
    kettlerJetYieldResults[2]->SetLineWidth(3);
	
	
    double jetVolumeYield;
    double jetVolumeYieldError;
    TString jetVolumeYieldLabel = "Jet per-trigger yield ";
   
    double integralEtaDep[3] = { fitParameterSigmaEtaPlotVsCent[2]->GetBinContent(2), 2.71437, 2.57664};
   
    for(int ptBin = 2; ptBin < NUM_PT_BINS; ptBin++){
        
            str1 = jetVolumeYieldLabel + ptBinLabel[ptBin];
            jetVolumePlot[ptBin] = new TH1D(str1, str1, 4, xbinsCentWithPythia);
            jetVolumePlot[ptBin]->GetXaxis()->SetBinLabel(1, "pythia");
            jetVolumePlot[ptBin]->GetXaxis()->SetBinLabel(2, "50-80%");
            jetVolumePlot[ptBin]->GetXaxis()->SetBinLabel(3, "20-50%");
            jetVolumePlot[ptBin]->GetXaxis()->SetBinLabel(4, "0-20%");
			
			jetVolumePlot[ptBin]->GetXaxis()->SetLabelSize(.06);
            jetVolumePlot[ptBin]->GetXaxis()->SetTitle("Centrality (%)");
            jetVolumePlot[ptBin]->GetXaxis()->SetTitleOffset(1.3);
			jetVolumePlot[ptBin]->GetYaxis()->SetTitle("NS Associated Yield");
			jetVolumePlot[ptBin]->GetYaxis()->SetTitleOffset(1.0);
			jetVolumePlot[ptBin]->GetYaxis()->SetTitleSize(.05);
			jetVolumePlot[ptBin]->GetYaxis()->SetLabelSize(.04);
            jetVolumePlot[ptBin]->SetMarkerStyle(9);
            
        
        for(int centBin = 0; centBin < 3; centBin++){
         
            if(centBin == 0) { jetVolumeYield = ((nBarPrime[centBin]*fitParameterJetAmpPlotVsCent[2]->GetBinContent(centBin+2)*fitParameterSigmaPhiPlotVsCent[2]->GetBinContent(centBin+2)*integralEtaDep[centBin])/2.0); }
            else jetVolumeYield = ((nBarPrime[centBin]*fitParameterJetAmpPlotVsCent[2]->GetBinContent(centBin+2)*fitParameterSigmaPhiPlotVsCent[2]->GetBinContent(centBin+2)*integralEtaDep[centBin])/(2.0*TMath::Sqrt(TMath::TwoPi())));
            jetVolumePlot[ptBin]->SetBinContent(centBin+2, jetVolumeYield);
            
            jetVolumeYieldError = (nBarPrime[centBin]/(TMath::Sqrt(TMath::TwoPi())))*getJetYieldError(fitParameterJetAmpPlotVsCent[2]->GetBinContent(centBin+2), integralEtaDep[centBin],fitParameterSigmaPhiPlotVsCent[2]->GetBinContent(centBin+2),
                                                   fitParameterJetAmpPlotVsCent[2]->GetBinError(centBin+2), fitParameterSigmaEtaPlotVsCent[2]->GetBinError(centBin+2),fitParameterSigmaPhiPlotVsCent[2]->GetBinError(centBin+2));
            
            //cout << jetVolumeYieldError << endl;
            
            jetVolumePlot[ptBin]->SetBinError(centBin+2, jetVolumeYieldError);
            
            cout << "# jet tracks per D0 in Centrality Bin: " << centBin <<"  = " << jetVolumeYield << " +/- " << jetVolumeYieldError << endl;
        
        }  
    }   
    
    
    //pythia -- 
    
    
    double pythiaNBar = 6.10505;
    //double pythiaNBar = 5;
    double pythiaIntEta = .520653;
    
    jetVolumeYield = (pythiaNBar*pythiaResultsJetAmp[2]->GetBinContent(1)*pythiaResultsPhiWidth[2]->GetBinContent(1)*pythiaResultsEtaWidth[2]->GetBinContent(1))/2.0);
    //jetVolumePlot[2]->SetBinContent(1, jetVolumeYield);
            
    jetVolumeYieldError = (pythiaNBar/2.0)*getJetYieldError(pythiaResultsJetAmp[2]->GetBinContent(1), pythiaResultsEtaWidth[2]->GetBinContent(1),pythiaResultsPhiWidth[2]->GetBinContent(1),
                                           pythiaResultsJetAmp[2]->GetBinError(1), pythiaResultsEtaWidth[2]->GetBinError(1),pythiaResultsPhiWidth[2]->GetBinError(1));
            
    //cout << jetVolumeYieldError << endl;

    TH1D* jetVolumePythiaPlot = new TH1D("","", 4, xbinsCentWithPythia);    
       
    jetVolumePythiaPlot->SetMarkerStyle(22);
    jetVolumePythiaPlot->SetMarkerSize(3);
    jetVolumePythiaPlot->SetMarkerColor(8);
    jetVolumePythiaPlot->SetLineWidth(5);
    jetVolumePythiaPlot->SetLineColor(8);
       
    jetVolumePythiaPlot->SetBinContent(1, jetVolumeYield);      
    jetVolumePythiaPlot->SetBinError(1, jetVolumeYieldError);
    //jetVolumePlot[2]->SetBinError(1, jetVolumeYieldError);
            
    cout << "# per-trigg yield PYTHIA  = " << jetVolumeYield << " +/- " << jetVolumeYieldError <<endl;
    
    TLegend * leg = new TLegend(.17, .88, .89, .6);
    leg->AddEntry(lightFlavorResultsJetAmpPt3[2], "di-Hadron, Mean Trigger p_{T} =  5.7 GeV/c [3]");
    leg->AddEntry(lightFlavorResultsJetAmpPt1[2], "di-Hadron, Mean Trigger p_{T} =  2.56 GeV/c [3]");
    //leg->AddEntry(lightFlavorResultsJetAmpPt2[2], "LF 3.82 GeV/c");
    
    //leg->AddEntry(lightFlavorResultsJetAmpPt4[2], "LF 1.71 GeV/c");
    //leg->AddEntry(lightFlavorResultsJetAmpPt5[2], "diHadron, 1.15 GeV/c Trigger [2]");
    //leg->AddEntry(lightFlavorResultsJetAmpPt6[2], "LF .77 GeV/c");
    leg->AddEntry(pythiaResultsJetAmp[2],         "Pythia D^{0}-Hadron, Mean D^{0} p_{T} = 3 GeV/c");
    leg->AddEntry(fitParameterJetAmpPlotVsCent[2], "D^{0}-Hadron AuAu 200 GeV, Mean D^{0} p_{T} = 3 GeV/c");
    
	leg->SetTextAlign(12);
	leg->SetFillColor(0);
    leg->SetMargin(.15);
    leg->SetBorderSize(0);
    
    
    
    TLegend * legQuad = new TLegend(.12, .88, .80, .80);
    legQuad->AddEntry(fitParameterPublishedQuadVsCent[2], "Average v_{2} in D^{0} P_{t} 2-10 GeV/c [3]");
    
    
    //finalFitCorrGridCanvas->cd(1);
    
    fitParameterQuadrupolePlotVsCent[2]->SetMarkerStyle(29);
    fitParameterQuadrupolePlotVsCent[2]->SetMarkerSize(5);
    fitParameterQuadrupolePlotVsCent[2]->SetLineWidth(5);
	fitParameterQuadrupolePlotVsCent[2]->SetMarkerColor(2);
    fitParameterQuadrupolePlotVsCent[2]->SetLineColor(2);
    fitParameterQuadrupolePlotVsCent[2]->SetMinimum(-.001);
    fitParameterQuadrupolePlotVsCent[2]->SetMaximum(.025);
    //fitParameterQuadrupolePlotVsCent[2]->Draw("E1");
    fitParameterPublishedQuadVsCent[2]->SetMarkerStyle(34);
    fitParameterPublishedQuadVsCent[2]->SetMarkerSize(3);
    fitParameterPublishedQuadVsCent[2]->SetMarkerColor(2);
    
    //fitParameterPublishedQuadVsCent[2]->Draw("SAME P");
    //legQuad->Draw("SAME");
    
    TCanvas * finalJetYieldCorrGridCanvas = new TCanvas("", "", 1200, 1100);
    //finalJetYieldCorrGridCanvas->Divide(2,1, 0.000001,0.000001);
    
    finalJetYieldCorrGridCanvas->cd(1);
    
    //jetVolumePlot[2]->SetMarkerStyle(20);
    //jetVolumePlot[2]->SetMarkerSize(2);
    //jetVolumePlot[2]->SetLineWidth(2);
    //jetVolumePlot[2]->SetLineColor(1);
    //jetVolumePlot[2]->SetMarkerColor(1);
    //jetVolumePlot[2]->SetStats(0);
    //jetVolumePlot[2]->Draw("P E1");
    
    lightFlavorResultsJetAmpPt1[2]->SetMinimum(0);
    lightFlavorResultsJetAmpPt1[2]->SetMaximum(1.8);
    lightFlavorResultsJetAmpPt1[2]->SetMarkerStyle(9);
    lightFlavorResultsJetAmpPt1[2]->SetMarkerSize(2);
    lightFlavorResultsJetAmpPt1[2]->SetMarkerColor(38);
    //lightFlavorResultsJetAmpPt1[2]->Draw("P");
    pythiaResultsJetAmp[2]->SetMarkerStyle(22);
    pythiaResultsJetAmp[2]->SetMarkerSize(2);
    pythiaResultsJetAmp[2]->SetMarkerColor(8);
    pythiaResultsJetAmp[2]->SetLineWidth(2);
    pythiaResultsJetAmp[2]->SetLineColor(8);
   // pythiaResultsJetAmp[2]->Draw("SAME P E1");
    lightFlavorResultsJetAmpPt2[2]->SetMarkerStyle(9);
    lightFlavorResultsJetAmpPt2[2]->SetMarkerSize(2);
    lightFlavorResultsJetAmpPt2[2]->SetMarkerColor(3);
    //lightFlavorResultsJetAmpPt2[2]->Draw("SAME P");
    lightFlavorResultsJetAmpPt3[2]->SetMarkerStyle(9);
    lightFlavorResultsJetAmpPt3[2]->SetMarkerSize(2);
    lightFlavorResultsJetAmpPt3[2]->SetMarkerColor(4);
    //lightFlavorResultsJetAmpPt3[2]->Draw("SAME P");
    lightFlavorResultsJetAmpPt4[2]->SetMarkerStyle(9);
    lightFlavorResultsJetAmpPt4[2]->SetMarkerSize(2);
    lightFlavorResultsJetAmpPt4[2]->SetMarkerColor(6);
    //lightFlavorResultsJetAmpPt4[2]->Draw("SAME P");
    lightFlavorResultsJetAmpPt5[2]->SetMarkerStyle(9);
    lightFlavorResultsJetAmpPt5[2]->SetMarkerSize(2);
    lightFlavorResultsJetAmpPt5[2]->SetMarkerColor(8);
    //lightFlavorResultsJetAmpPt5[2]->Draw("SAME P");
    lightFlavorResultsJetAmpPt6[2]->SetMarkerStyle(9);
    lightFlavorResultsJetAmpPt6[2]->SetMarkerSize(2);
    lightFlavorResultsJetAmpPt6[2]->SetMarkerColor(12);
    //lightFlavorResultsJetAmpPt6[2]->Draw("SAME P");
    fitParameterJetAmpPlotVsCent[2]->SetMarkerStyle(29);
    fitParameterJetAmpPlotVsCent[2]->SetMarkerSize(4);
    fitParameterJetAmpPlotVsCent[2]->SetLineWidth(2);
	fitParameterJetAmpPlotVsCent[2]->SetMarkerColor(2);
	fitParameterJetAmpPlotVsCent[2]->SetLineColor(2);
    //fitParameterJetAmpPlotVsCent[2]->Draw("SAME E1"); 
   
    TLegend * legJetYield = new TLegend(.12, .89, .77, .73);
    legJetYield->AddEntry(lightFlavorResultsJetAmpPt3[2], "di-Hadron, Mean Trigger p_{T} =  5.7 GeV/c [3]");
    legJetYield->AddEntry(lightFlavorResultsJetAmpPt1[2], "di-Hadron, Mean Trigger p_{T} =  2.56 GeV/c [3]");
    //leg->AddEntry(lightFlavorResultsJetAmpPt2[2], "LF 3.82 GeV/c");
    
    //leg->AddEntry(lightFlavorResultsJetAmpPt4[2], "LF 1.71 GeV/c");
    //leg->AddEntry(lightFlavorResultsJetAmpPt5[2], "diHadron, 1.15 GeV/c Trigger [2]");
    //leg->AddEntry(lightFlavorResultsJetAmpPt6[2], "LF .77 GeV/c");
    legJetYield->AddEntry(pythiaResultsJetAmp[2],         "Pythia D^{0}-Hadron, Mean D^{0} p_{T} = 3 GeV/c");
    legJetYield->AddEntry(fitParameterJetAmpPlotVsCent[2], "D^{0}-Hadron AuAu 200 GeV, Mean D^{0} p_{T} = 3 GeV/c");
    
	legJetYield->SetTextAlign(12);
	legJetYield->SetFillColor(0);
    legJetYield->SetMargin(.15);
	legJetYield->SetBorderSize(0);
   
   
    jetVolumePlot[2]->SetMaximum(99);
    jetVolumePlot[2]->SetMarkerStyle(29);
    jetVolumePlot[2]->SetMarkerSize(5);
    jetVolumePlot[2]->SetLineWidth(5);
    jetVolumePlot[2]->SetLineColor(2);
    jetVolumePlot[2]->SetMarkerColor(2);
    jetVolumePlot[2]->SetStats(0);
	jetVolumePlot[2]->SetTitle("");
    //jetVolumePlot[2]
    finalJetYieldCorrGridCanvas->cd(1)->SetLogy(1);
    
    jetVolumePlot[2]->Draw("P E1");
    
    jetVolumePythiaPlot->Draw("SAME P E1");
    
    double perTrigYieldValues[3];
    perTrigYieldValues[0] = jetVolumePlot[2]->GetBinContent(2);
    perTrigYieldValues[1] = jetVolumePlot[2]->GetBinContent(3);
    perTrigYieldValues[2] = jetVolumePlot[2]->GetBinContent(4);
    
    double perTrigYieldSystematics[3] = { 0.03429352, 0.58216638 , 1.688138436};
    
    //perTrigYieldSystematics[0] = perTrigYieldValues[0]*perTrigYieldSystematics[0];
    //perTrigYieldSystematics[1] = perTrigYieldValues[1]*perTrigYieldSystematics[1];
    //perTrigYieldSystematics[2] = perTrigYieldValues[2]*perTrigYieldSystematics[2];
    
    
    //SYSTEMATICS ERRORS HERE
	
	double specialXValuesWithPythia[3] = {1.5, 2.5, 3.335};
    
    for(int i = 0; i < 3; i++) {
        double x1 = specialXValuesWithPythia[i]-.06;
        double x2 = specialXValuesWithPythia[i]+0.06;
        double y1 = perTrigYieldValues[i]-perTrigYieldSystematics[i];
        double y2 = perTrigYieldValues[i]+perTrigYieldSystematics[i];

   
    
        TLine *la = new TLine(x1, y1, x1, y1+0.01);
        la->SetLineColor(2);
        la->SetLineWidth(5);
        la->Draw("same");
        TLine *lb = new TLine(x2, y1, x2, y1+0.01);
        lb->SetLineColor(2);
        lb->SetLineWidth(5);
        lb->Draw("same");
        TLine *lc = new TLine(x1, y2, x1, y2-0.01);
        lc->SetLineColor(2);
        lc->SetLineWidth(5);
        lc->Draw("same");
        TLine *ld = new TLine(x2, y2, x2, y2-0.01);
        ld->SetLineColor(2);
        ld->SetLineWidth(5);
        ld->Draw("same");
        TLine *le = new TLine(x1, y1, x2, y1);
        le->SetLineWidth(4);
        le->SetLineColor(2);
        le->Draw("same");
        TLine *lf = new TLine(x1, y2, x2, y2);
        lf->SetLineWidth(4);
        lf->SetLineColor(2);
        lf->Draw("same");
    }
    
	//kettlerJetYieldResults[0]->Draw("SAME P");
	kettlerJetYieldResults[1]->Draw("SAME P");
	kettlerJetYieldResults[2]->Draw("SAME P");
	kettlerPeripheralPt1JetYield->Draw("SAME P");
    kettlerPeripheralPt3JetYield->Draw("SAME P");
    
    legJetYield->Draw("SAME");
	
	finalJetYieldCorrGridCanvas->SaveAs("Jet_Yield_Plot.png");
	
	TCanvas * finalFitCorrGridCanvas = new TCanvas("", "", 2400, 1100);
    finalFitCorrGridCanvas->Divide(2,1, 0.000001,0.000001);
	
    finalFitCorrGridCanvas->cd(1);
	finalFitCorrGridCanvas->cd(1)->SetLeftMargin(.15);
    fitParameterSigmaPhiPlotVsCent[2]->SetMarkerStyle(29);
    fitParameterSigmaPhiPlotVsCent[2]->SetMaximum(2.0);
    fitParameterSigmaPhiPlotVsCent[2]->SetMarkerSize(5);
    fitParameterSigmaPhiPlotVsCent[2]->SetLineWidth(4);
	fitParameterSigmaPhiPlotVsCent[2]->SetMarkerColor(2);
    fitParameterSigmaPhiPlotVsCent[2]->SetLineColor(2);
	fitParameterSigmaPhiPlotVsCent[2]->SetTitle("");
	fitParameterSigmaPhiPlotVsCent[2]->GetYaxis()->SetTitle("#sigma_{NS, #Delta#phi}");
	fitParameterSigmaPhiPlotVsCent[2]->GetYaxis()->SetTitleSize(.06);
	fitParameterSigmaPhiPlotVsCent[2]->GetYaxis()->SetTitleOffset(1.0);
    fitParameterSigmaPhiPlotVsCent[2]->Draw("E1");
    
    double phiWidthValues[3] = {.349718, .669247, .754439}; //these are the nominal values
    
    
    double phiWidthSystematics[3] = { .037051958, .075790561 , .077030402};
    
   // phiWidthSystematics[0] = phiWidthValues[0]*phiWidthSystematics[0];
    //phiWidthSystematics[1] = phiWidthValues[1]*phiWidthSystematics[1];
    //phiWidthSystematics[2] = phiWidthValues[2]*phiWidthSystematics[2];
    
    //SYSTEMATICS ERRORS HERE
    
	
	
    for(int i = 0; i < 3; i++) {
        double x1 = specialXValuesWithPythia[i]-.06;
        double x2 = specialXValuesWithPythia[i]+0.06;
        double y1 = phiWidthValues[i]-phiWidthSystematics[i];
        double y2 = phiWidthValues[i]+phiWidthSystematics[i];

   
    
        TLine *la = new TLine(x1, y1, x1, y1+0.01);
        la->SetLineColor(2);
        la->SetLineWidth(5);
        la->Draw("same");
        TLine *lb = new TLine(x2, y1, x2, y1+0.01);
        lb->SetLineColor(2);
        lb->SetLineWidth(5);
        lb->Draw("same");
        TLine *lc = new TLine(x1, y2, x1, y2-0.01);
        lc->SetLineColor(2);
        lc->SetLineWidth(5);
        lc->Draw("same");
        TLine *ld = new TLine(x2, y2, x2, y2-0.01);
        ld->SetLineColor(2);
        ld->SetLineWidth(5);
        ld->Draw("same");
        TLine *le = new TLine(x1, y1, x2, y1);
        le->SetLineWidth(4);
        le->SetLineColor(2);
        le->Draw("same");
        TLine *lf = new TLine(x1, y2, x2, y2);
        lf->SetLineWidth(4);
        lf->SetLineColor(2);
        lf->Draw("same");
    }
    
    ////////////////////////////////////////////////
    
    lightFlavorResultsPhiWidthPt1[2]->SetMinimum(0);
    lightFlavorResultsPhiWidthPt1[2]->SetMaximum(2.0);
    lightFlavorResultsPhiWidthPt1[2]->SetMarkerStyle(9);
    lightFlavorResultsPhiWidthPt1[2]->SetMarkerSize(2);
    lightFlavorResultsPhiWidthPt1[2]->SetMarkerColor(38);
    lightFlavorResultsPhiWidthPt1[2]->SetLineWidth(4);
	lightFlavorResultsPhiWidthPt1[2]->SetLineColor(38);
    lightFlavorResultsPhiWidthPt1[2]->Draw("SAME P");
    kettlerPeripheralPt1Phi->SetMarkerStyle(9);
    kettlerPeripheralPt1Phi->SetMarkerSize(2);
    kettlerPeripheralPt1Phi->SetMarkerColor(38);
    kettlerPeripheralPt1Phi->SetLineWidth(4);
    kettlerPeripheralPt1Phi->SetLineColor(38);
    kettlerPeripheralPt1Phi->Draw("SAME P");

    pythiaResultsPhiWidth[2]->SetMarkerStyle(22);
    pythiaResultsPhiWidth[2]->SetMarkerSize(3);
    pythiaResultsPhiWidth[2]->SetMarkerColor(8);
    pythiaResultsPhiWidth[2]->SetLineWidth(4);
    pythiaResultsPhiWidth[2]->SetLineColor(8);
    pythiaResultsPhiWidth[2]->Draw("SAME P E1");
    lightFlavorResultsPhiWidthPt2[2]->SetMarkerStyle(9);
    lightFlavorResultsPhiWidthPt2[2]->SetMarkerSize(2);
    lightFlavorResultsPhiWidthPt2[2]->SetMarkerColor(3);
    //lightFlavorResultsPhiWidthPt2[2]->Draw("SAME P");
    lightFlavorResultsPhiWidthPt3[2]->SetMarkerStyle(9);
    lightFlavorResultsPhiWidthPt3[2]->SetMarkerSize(2);
    lightFlavorResultsPhiWidthPt3[2]->SetMarkerColor(4);
    lightFlavorResultsPhiWidthPt3[2]->SetLineWidth(4);
	lightFlavorResultsPhiWidthPt3[2]->SetLineColor(4);
    lightFlavorResultsPhiWidthPt3[2]->Draw("SAME P");
    kettlerPeripheralPt3Phi->SetMarkerStyle(9);
    kettlerPeripheralPt3Phi->SetMarkerSize(2);
    kettlerPeripheralPt3Phi->SetMarkerColor(4);
    kettlerPeripheralPt3Phi->SetLineWidth(4);
    kettlerPeripheralPt3Phi->SetLineColor(4);
    kettlerPeripheralPt3Phi->Draw("SAME P");
    //lightFlavorResultsPhiWidthPt4[2]->SetMarkerStyle(9);
    //lightFlavorResultsPhiWidthPt4[2]->SetMarkerSize(2);
    //lightFlavorResultsPhiWidthPt4[2]->SetMarkerColor(6);
    //lightFlavorResultsPhiWidthPt4[2]->Draw("SAME P");
    //lightFlavorResultsPhiWidthPt5[2]->SetMarkerStyle(9);
    //lightFlavorResultsPhiWidthPt5[2]->SetMarkerSize(2);
    //lightFlavorResultsPhiWidthPt5[2]->SetMarkerColor(8);
	//lightFlavorResultsPhiWidthPt5[2]->SetLineColor(8);
    //lightFlavorResultsPhiWidthPt5[2]->Draw("SAME P");
    //lightFlavorResultsPhiWidthPt6[2]->SetMarkerStyle(9);
    //lightFlavorResultsPhiWidthPt6[2]->SetMarkerSize(2);
    //lightFlavorResultsPhiWidthPt6[2]->SetMarkerColor(12);
    //lightFlavorResultsPhiWidthPt6[2]->Draw("SAME P");
    
    leg->Draw("SAME");
    
    
    finalFitCorrGridCanvas->cd(2);
    fitParameterSigmaEtaPlotVsCent[2]->SetMarkerStyle(29);
    fitParameterSigmaEtaPlotVsCent[2]->SetMarkerSize(5);
    fitParameterSigmaEtaPlotVsCent[2]->SetLineWidth(5);
	fitParameterSigmaEtaPlotVsCent[2]->SetMarkerColor(2);
	fitParameterSigmaEtaPlotVsCent[2]->SetLineColor(2);
    fitParameterSigmaEtaPlotVsCent[2]->SetMaximum(2.5);
	fitParameterSigmaEtaPlotVsCent[2]->SetTitle("");
	fitParameterSigmaEtaPlotVsCent[2]->GetYaxis()->SetTitle("#sigma_{NS, #Delta#eta}");
	fitParameterSigmaEtaPlotVsCent[2]->GetYaxis()->SetTitleSize(.06);
	fitParameterSigmaEtaPlotVsCent[2]->GetYaxis()->SetTitleOffset(1.0);
    fitParameterSigmaEtaPlotVsCent[2]->Draw("E1");
    
    double etaWidthValues[3] = {.312136, 1.40393, 1.24592};
    double etaWidthSystematics[3] = { .032325296, .200152368 , .129029439};
    
    
    //etaWidthSystematics[0] = etaWidthValues[0]*etaWidthSystematics[0];
    //etaWidthSystematics[1] = etaWidthValues[1]*etaWidthSystematics[1];
    //etaWidthSystematics[2] = etaWidthValues[2]*etaWidthSystematics[2];
    
    //SYSTEMATICS ERRORS HERE
    
    for(int i = 0; i < 3; i++) {
        double x1 = specialXValuesWithPythia[i]-.06;
        double x2 = specialXValuesWithPythia[i]+0.06;
        double y1 = etaWidthValues[i]-etaWidthSystematics[i];
        double y2 = etaWidthValues[i]+etaWidthSystematics[i];

   
    
        TLine *la = new TLine(x1, y1, x1, y1+0.01);
        la->SetLineColor(2);
        la->SetLineWidth(5);
        la->Draw("same");
        TLine *lb = new TLine(x2, y1, x2, y1+0.01);
        lb->SetLineColor(2);
        lb->SetLineWidth(5);
        lb->Draw("same");
        TLine *lc = new TLine(x1, y2, x1, y2-0.01);
        lc->SetLineColor(2);
        lc->SetLineWidth(5);
        lc->Draw("same");
        TLine *ld = new TLine(x2, y2, x2, y2-0.01);
        ld->SetLineColor(2);
        ld->SetLineWidth(5);
        ld->Draw("same");
        TLine *le = new TLine(x1, y1, x2, y1);
        le->SetLineWidth(4);
        le->SetLineColor(2);
        le->Draw("same");
        TLine *lf = new TLine(x1, y2, x2, y2);
        lf->SetLineWidth(4);
        lf->SetLineColor(2);
        lf->Draw("same");
    }
    
    
    lightFlavorResultsEtaWidthPt1[2]->SetMinimum(0);
    lightFlavorResultsEtaWidthPt1[2]->SetMaximum(3.0);
    lightFlavorResultsEtaWidthPt1[2]->SetMarkerStyle(9);
    lightFlavorResultsEtaWidthPt1[2]->SetMarkerSize(2);
    lightFlavorResultsEtaWidthPt1[2]->SetMarkerColor(38);
    lightFlavorResultsEtaWidthPt1[2]->SetLineWidth(4);
	lightFlavorResultsEtaWidthPt1[2]->SetLineColor(38);
    lightFlavorResultsEtaWidthPt1[2]->Draw("SAME P");
    pythiaResultsEtaWidth[2]->SetMarkerStyle(22);
    pythiaResultsEtaWidth[2]->SetMarkerSize(3);
    pythiaResultsEtaWidth[2]->SetMarkerColor(8);
    pythiaResultsEtaWidth[2]->SetLineWidth(4);
    pythiaResultsEtaWidth[2]->SetLineColor(8);
    pythiaResultsEtaWidth[2]->Draw("SAME P E1");
    kettlerPeripheralPt1Eta->SetMarkerStyle(9);
    kettlerPeripheralPt1Eta->SetMarkerSize(2);
    kettlerPeripheralPt1Eta->SetMarkerColor(38);
    kettlerPeripheralPt1Eta->SetLineWidth(4);
    kettlerPeripheralPt1Eta->SetLineColor(38);
    kettlerPeripheralPt1Eta->Draw("SAME P");
    lightFlavorResultsEtaWidthPt2[2]->SetMarkerStyle(9);
    lightFlavorResultsEtaWidthPt2[2]->SetMarkerSize(2);
    lightFlavorResultsEtaWidthPt2[2]->SetMarkerColor(3);
    //lightFlavorResultsEtaWidthPt2[2]->Draw("SAME P");
    lightFlavorResultsEtaWidthPt3[2]->SetMarkerStyle(9);
    lightFlavorResultsEtaWidthPt3[2]->SetMarkerSize(2);
    lightFlavorResultsEtaWidthPt3[2]->SetMarkerColor(4);
	lightFlavorResultsEtaWidthPt3[2]->SetLineColor(4);
    lightFlavorResultsEtaWidthPt3[2]->SetLineWidth(4);
    lightFlavorResultsEtaWidthPt3[2]->Draw("SAME P");
    kettlerPeripheralPt3Eta->SetMarkerStyle(9);
    kettlerPeripheralPt3Eta->SetMarkerSize(2);
    kettlerPeripheralPt3Eta->SetMarkerColor(4);
    kettlerPeripheralPt3Eta->SetLineWidth(4);
    kettlerPeripheralPt3Eta->SetLineColor(4);
    kettlerPeripheralPt3Eta->Draw("SAME P");
    lightFlavorResultsEtaWidthPt4[2]->SetMarkerStyle(9);
    lightFlavorResultsEtaWidthPt4[2]->SetMarkerSize(2);
    lightFlavorResultsEtaWidthPt4[2]->SetMarkerColor(6);
    //lightFlavorResultsEtaWidthPt4[2]->Draw("SAME P");
    lightFlavorResultsEtaWidthPt5[2]->SetMarkerStyle(9);
    lightFlavorResultsEtaWidthPt5[2]->SetMarkerSize(2);
    lightFlavorResultsEtaWidthPt5[2]->SetMarkerColor(8);
    //lightFlavorResultsEtaWidthPt5[2]->Draw("SAME P");
    lightFlavorResultsEtaWidthPt6[2]->SetMarkerStyle(9);
    lightFlavorResultsEtaWidthPt6[2]->SetMarkerSize(2);
    lightFlavorResultsEtaWidthPt6[2]->SetMarkerColor(12);
    //lightFlavorResultsEtaWidthPt6[2]->Draw("SAME P");
    
    
    //finalFitCorrGridCanvas->cd(5);
    //jetVolumePlot[2]->SetMarkerStyle(20);
    //jetVolumePlot[2]->Draw("P");
    
    finalFitCorrGridCanvas->SaveAs("NS_Widths_plots.png");
    
    TCanvas * finalFitCorrGridCanvasQuadrupole = new TCanvas("", "", 1200, 1100);
    
  
    fitParameterQuadrupolePlotVsCent[2]->Draw("E1");
    finalFitCorrGridCanvasQuadrupole->cd(1);
    gPad->SetRightMargin(.01);
    gPad->SetLeftMargin(.15);
     double quadValues[3];
    quadValues[0] = fitParameterQuadrupolePlotVsCent[2]->GetBinContent(1);
    quadValues[1] = fitParameterQuadrupolePlotVsCent[2]->GetBinContent(2);
    quadValues[2] = fitParameterQuadrupolePlotVsCent[2]->GetBinContent(3);
    
    double quadSystematics[3] = { .000373544, 0.000601248 , .00092870878};
    
    
    //SYSTEMATICS ERRORS HERE
    
	double specialXValuesNoPythia[3] = {.5, 1.5, 2.335};
	
    for(int i = 0; i < 3; i++) {
        double x1 = specialXValuesNoPythia[i]-.06;
        double x2 = specialXValuesNoPythia[i]+0.06;
        double y1 = quadValues[i]-quadSystematics[i];
        double y2 = quadValues[i]+quadSystematics[i];

   
    
        TLine *la = new TLine(x1, y1, x1, y1+0.0001);
        la->SetLineColor(2);
        la->SetLineWidth(5);
        la->Draw("same");
        TLine *lb = new TLine(x2, y1, x2, y1+0.0001);
        lb->SetLineColor(2);
        lb->SetLineWidth(5);
        lb->Draw("same");
        TLine *lc = new TLine(x1, y2, x1, y2-0.0001);
        lc->SetLineColor(2);
        lc->SetLineWidth(5);
        lc->Draw("same");
        TLine *ld = new TLine(x2, y2, x2, y2-0.0001);
        ld->SetLineColor(2);
        ld->SetLineWidth(5);
        ld->Draw("same");
        TLine *le = new TLine(x1, y1, x2, y1);
        le->SetLineWidth(4);
        le->SetLineColor(2);
        le->Draw("same");
        TLine *lf = new TLine(x1, y2, x2, y2);
        lf->SetLineWidth(4);
        lf->SetLineColor(2);
        lf->Draw("same");
    }
    
    kettlerPt1Quad->Draw("SAME P");
    kettlerPt3Quad->Draw("SAME P");
    
	//gStyle->SetLegendFillColor(0);
	
    TLegend * leg2 = new TLegend(.18, .88, .98, .7);
    leg2->AddEntry(lightFlavorResultsJetAmpPt3[2], "di-Hadron, Mean Trigger p_{T} =  5.7 GeV/c [3]");
    leg2->AddEntry(lightFlavorResultsJetAmpPt1[2], "di-Hadron, Mean Trigger p_{T} =  2.56 GeV/c [3]");
    leg2->AddEntry(fitParameterJetAmpPlotVsCent[2], "D^{0}-Hadron Au+Au 200 GeV, Mean D^{0} p_{T} = 3 GeV/c");
    leg2->SetTextAlign(12);
	leg2->SetFillColor(0);
    leg2->SetMargin(.095);
	leg2->SetBorderSize(0);
    
    leg2->Draw("SAME");
    
    finalFitCorrGridCanvasQuadrupole->SaveAs("final_output_fit_parameters_QUADRUPOLE.png");
    
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


/*///////////////////pt = 1.71
    lightFlavorResultsJetAmpPt4[2]->SetBinContent(2, .044);
    lightFlavorResultsJetAmpPt4[2]->SetBinContent(3, .035);
    lightFlavorResultsJetAmpPt4[2]->SetBinContent(4, .026);
    
    lightFlavorResultsEtaWidthPt4[2]->SetBinContent(2, .66);
    lightFlavorResultsEtaWidthPt4[2]->SetBinContent(3, 1.55);
    lightFlavorResultsEtaWidthPt4[2]->SetBinContent(4, 2.39);
    
    lightFlavorResultsPhiWidthPt4[2]->SetBinContent(2, .5);
    lightFlavorResultsPhiWidthPt4[2]->SetBinContent(3, .59);
    lightFlavorResultsPhiWidthPt4[2]->SetBinContent(4, .62);
    
    ///////////////////pt = 1.15
    lightFlavorResultsJetAmpPt5[2]->SetBinContent(2, .630);
    lightFlavorResultsJetAmpPt5[2]->SetBinContent(3, 3.45);
    lightFlavorResultsJetAmpPt5[2]->SetBinContent(4, 8.32);
    
    lightFlavorResultsEtaWidthPt5[2]->SetBinContent(2, .73);
    lightFlavorResultsEtaWidthPt5[2]->SetBinContent(3, 1.56);
    lightFlavorResultsEtaWidthPt5[2]->SetBinContent(4, 2.51);
    
    lightFlavorResultsPhiWidthPt5[2]->SetBinContent(2, .59);
    lightFlavorResultsPhiWidthPt5[2]->SetBinContent(3, .63);
    lightFlavorResultsPhiWidthPt5[2]->SetBinContent(4, .66);
    
    ///////////////////pt = .77
    lightFlavorResultsJetAmpPt6[2]->SetBinContent(2, .0243);
    lightFlavorResultsJetAmpPt6[2]->SetBinContent(3, .020);
    lightFlavorResultsJetAmpPt6[2]->SetBinContent(4, .016);
    
    lightFlavorResultsEtaWidthPt6[2]->SetBinContent(2, .8);
    lightFlavorResultsEtaWidthPt6[2]->SetBinContent(3, 1.60);
    lightFlavorResultsEtaWidthPt6[2]->SetBinContent(4, 2.74);
    
    lightFlavorResultsPhiWidthPt6[2]->SetBinContent(2, .69);
    lightFlavorResultsPhiWidthPt6[2]->SetBinContent(3, .70);
    lightFlavorResultsPhiWidthPt6[2]->SetBinContent(4, .72);
    
     
    
    
    ///////////pt = 3.82
    lightFlavorResultsJetAmpPt2[2]->SetBinContent(2, .076);
    lightFlavorResultsJetAmpPt2[2]->SetBinContent(3, .046);
    lightFlavorResultsJetAmpPt2[2]->SetBinContent(4, .036);
    
    lightFlavorResultsEtaWidthPt2[2]->SetBinContent(2, .49);
    lightFlavorResultsEtaWidthPt2[2]->SetBinContent(3, .93);
    lightFlavorResultsEtaWidthPt2[2]->SetBinContent(4, 1.86);
    
    lightFlavorResultsPhiWidthPt2[2]->SetBinContent(2, .40);
    lightFlavorResultsPhiWidthPt2[2]->SetBinContent(3, .44);
    lightFlavorResultsPhiWidthPt2[2]->SetBinContent(4, .54);
    
    
    
    
    
    
    */
    
     //SPECIAL TEST WITH ONLY 50-70% centrality
    
    /*fitParameterQuadrupolePlotVsCent[2]->SetBinContent(1, .00465721);  //-- from AS 1D Gaussian fit
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(1, .00267145);   //-- from AS 1D Gaussian fit
	//fitParameterQuadrupolePlotVsCent[2]->SetBinContent(1, .00867822);   //from 6-parameter
    //fitParameterQuadrupolePlotVsCent[2]->SetBinError(1, .000925368);    //from 6-parameter
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(2, .0806960);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(2, .0198731);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(2, .315697);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(2, .0733015);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(2, .299356);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(2, .109332);*/
    
    
    //--------------------------------------
    
    //-------------from symmetric fits----------
    //----------------------------
    //---------------------------
    /*fitParameterQuadrupolePlotVsCent[2]->SetBinContent(1, .00501008);  //-- from AS 1D Gaussian fit
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(1, .00134967);   //-- from AS 1D Gaussian fit
	//fitParameterQuadrupolePlotVsCent[2]->SetBinContent(1, .00867822);   //from 6-parameter
    //fitParameterQuadrupolePlotVsCent[2]->SetBinError(1, .000925368);    //from 6-parameter
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(1, .0886923);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(1, .00970244);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(1, .332103);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(1, .0661020);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(1, .289871);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(1, .0523706);
    
   
    fitParameterQuadrupolePlotVsCent[2]->SetBinContent(2, .00605651);
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(2, .00162442);
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(2, .0411934);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(2, .00174167);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(2, .673615);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(2, .0282700);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(2, 1.95394);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(2, .474411);
    
   
    
    /*fitParameterQuadrupolePlotVsCent[2]->SetBinContent(3, .00264);
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(3, .000272316);
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(3, .00751125);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(3, .00174301);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(3, .433461);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(3, .0522924);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(3, 1.8);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(3, .258073); //non-zero quadrupole
    
    fitParameterQuadrupolePlotVsCent[2]->SetBinContent(3, 0.0);
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(3,   .000312240);
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(3,     .0403940);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(3,        .00265756);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(3,    .71155);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(3,       .0363435);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(3,     1.06623);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(3,       .0950461);*/
    
    
    //------------------------------------------------------------------
    //------------------from absolute value fits------------------------
    //------------------------------------------------------------------
    
   /*  fitParameterQuadrupolePlotVsCent[2]->SetBinContent(1, .00345693);  //-- from AS 1D Gaussian fit
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(1, .00388678);   //-- from AS 1D Gaussian fit
	//fitParameterQuadrupolePlotVsCent[2]->SetBinContent(1, .00867822);   //from 6-parameter
    //fitParameterQuadrupolePlotVsCent[2]->SetBinError(1, .000925368);    //from 6-parameter
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(1, .0954529);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(1, .0115554);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(1, .331899);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(1, .116135);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(1, .292285);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(1, .0561364);
    
   
    fitParameterQuadrupolePlotVsCent[2]->SetBinContent(2, .00439827);
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(2, .00116807);
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(2, .0678038);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(2, .00537518);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(2, .771187);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(2, .0559492);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(2, 1.86604);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(2, .402917);
    
  
    fitParameterQuadrupolePlotVsCent[2]->SetBinContent(3, 0.0);
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(3,   .00282599);
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(3,     .0402207);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(3,        .00435868);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(3,    .705670);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(3,       .0580549);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(3,     1.01926);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(3,       .144411); */
    