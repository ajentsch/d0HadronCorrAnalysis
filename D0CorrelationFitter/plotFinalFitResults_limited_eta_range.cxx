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
   
    double nBarPrime[3] = {74.4, 351, 1027};
   
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
    
        str1 = offsetLabel + ptLabel + ptBinLabel[k];
        fitParameterOffsetPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = dipoleLabel + ptLabel + ptBinLabel[k];
        fitParameterDipolePlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = quadLabel + ptLabel + ptBinLabel[k];
        fitParameterQuadrupolePlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        fitParameterPublishedQuadVsCent[k] = new TH1D("","",3,0,3);
        str1 = jetAmpLabel + ptLabel + ptBinLabel[k];
        fitParameterJetAmpPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = sigPhiLabel + ptLabel + ptBinLabel[k];
        fitParameterSigmaPhiPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = sigEtaLabel + ptLabel + ptBinLabel[k];
        fitParameterSigmaEtaPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        
        str1 = awaySideAmpLabel + ptLabel + ptBinLabel[k];
        fitParameterASAmpPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = awaySidePhiWidthLabel + ptLabel + ptBinLabel[k];
        fitParameterASSigmaPhiPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        
        str1 = lightFlavorLabel + lightFlavorPtLabel[0] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt1[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = lightFlavorLabel + lightFlavorPtLabel[0] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt1[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = lightFlavorLabel + lightFlavorPtLabel[0] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt1[k] = new TH1D(str1, str1, 3, 0, 3);
        
        str1 = lightFlavorLabel + lightFlavorPtLabel[1] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt2[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = lightFlavorLabel + lightFlavorPtLabel[1] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt2[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = lightFlavorLabel + lightFlavorPtLabel[1] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt2[k] = new TH1D(str1, str1, 3, 0, 3);
        
        str1 = lightFlavorLabel + lightFlavorPtLabel[2] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt3[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = lightFlavorLabel + lightFlavorPtLabel[2] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt3[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = lightFlavorLabel + lightFlavorPtLabel[2] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt3[k] = new TH1D(str1, str1, 3, 0, 3);
        
        str1 = lightFlavorLabel + lightFlavorPtLabel[3] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt4[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = lightFlavorLabel + lightFlavorPtLabel[3] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt4[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = lightFlavorLabel + lightFlavorPtLabel[3] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt4[k] = new TH1D(str1, str1, 3, 0, 3);
        
        str1 = lightFlavorLabel + lightFlavorPtLabel[4] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt5[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = lightFlavorLabel + lightFlavorPtLabel[4] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt5[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = lightFlavorLabel + lightFlavorPtLabel[4] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt5[k] = new TH1D(str1, str1, 3, 0, 3);
        
        str1 = lightFlavorLabel + lightFlavorPtLabel[5] + jetAmpLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsJetAmpPt6[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = lightFlavorLabel + lightFlavorPtLabel[5] + sigEtaLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsEtaWidthPt6[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = lightFlavorLabel + lightFlavorPtLabel[5] + sigPhiLabel + ptLabel + ptBinLabel[k];
        lightFlavorResultsPhiWidthPt6[k] = new TH1D(str1, str1, 3, 0, 3);
        
        str1 = pythiaLabel + jetAmpLabel + ptLabel + ptBinLabel[k];
        pythiaResultsJetAmp[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = pythiaLabel + sigEtaLabel + ptLabel + ptBinLabel[k];
        pythiaResultsEtaWidth[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = pythiaLabel + sigPhiLabel + ptLabel + ptBinLabel[k];
        pythiaResultsPhiWidth[k] = new TH1D(str1, str1, 3, 0, 3);
        
        pythiaResultsJetAmp[k]->GetXaxis()->SetLabelSize(.07);
        pythiaResultsJetAmp[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        pythiaResultsJetAmp[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        pythiaResultsJetAmp[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        pythiaResultsEtaWidth[k]->GetXaxis()->SetLabelSize(.07);
        pythiaResultsEtaWidth[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        pythiaResultsEtaWidth[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        pythiaResultsEtaWidth[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        pythiaResultsPhiWidth[k]->GetXaxis()->SetLabelSize(.07);
        pythiaResultsPhiWidth[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        pythiaResultsPhiWidth[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        pythiaResultsPhiWidth[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetLabelSize(.07);
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetLabelSize(.07);
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterQuadrupolePlotVsCent[k]->SetLabelSize(.07);
        
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetLabelSize(.07);
        fitParameterJetAmpPlotVsCent[k]->GetYaxis()->SetLabelSize(.05);
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetLabelSize(.07);
        fitParameterSigmaEtaPlotVsCent[k]->GetYaxis()->SetLabelSize(.05);
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetLabelSize(.07);
        fitParameterSigmaPhiPlotVsCent[k]->GetYaxis()->SetLabelSize(.05);
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetLabelSize(.07);
        fitParameterASAmpPlotVsCent[k]->GetYaxis()->SetLabelSize(.05);
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        
        fitParameterASSigmaPhiPlotVsCent[k]->GetXaxis()->SetLabelSize(.07);
        fitParameterASSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterASSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterASSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        
        lightFlavorResultsJetAmpPt1[k]->GetXaxis()->SetLabelSize(.07);
        lightFlavorResultsJetAmpPt1[k]->GetYaxis()->SetLabelSize(.05);
        lightFlavorResultsJetAmpPt1[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        lightFlavorResultsJetAmpPt1[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        lightFlavorResultsJetAmpPt1[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        lightFlavorResultsEtaWidthPt1[k]->GetXaxis()->SetLabelSize(.07);
        lightFlavorResultsEtaWidthPt1[k]->GetYaxis()->SetLabelSize(.05);
        lightFlavorResultsEtaWidthPt1[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        lightFlavorResultsEtaWidthPt1[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        lightFlavorResultsEtaWidthPt1[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        lightFlavorResultsPhiWidthPt1[k]->GetXaxis()->SetLabelSize(.07);
        lightFlavorResultsPhiWidthPt1[k]->GetYaxis()->SetLabelSize(.05);
        lightFlavorResultsPhiWidthPt1[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        lightFlavorResultsPhiWidthPt1[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        lightFlavorResultsPhiWidthPt1[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        
        
        
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
    
    
    /* ///TEMPORARY NUMBERS--------------------
    fitParameterPublishedQuadVsCent[2]->SetBinContent(1, .00973648);
    fitParameterPublishedQuadVsCent[2]->SetBinError(1, .000860232);
    fitParameterPublishedQuadVsCent[2]->SetBinContent(2, .0072);
    fitParameterPublishedQuadVsCent[2]->SetBinContent(3, .0036);
    
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
    
    pythiaResultsJetAmp[2]->SetBinContent(1, 1.6);
    pythiaResultsJetAmp[2]->SetBinError(1, .280516);
       
    pythiaResultsEtaWidth[2]->SetBinContent(1, .415422);
    pythiaResultsEtaWidth[2]->SetBinError(1,    .0331115);
        
    pythiaResultsPhiWidth[2]->SetBinContent(1, .406196);
    pythiaResultsPhiWidth[2]->SetBinError(1,   .0247533);
    
    
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
    
    //------------------------------------------------------------------
    //------------------from symmetric space with correctly scaled errors------------------------
    //------------------------------------------------------------------
    
    fitParameterQuadrupolePlotVsCent[2]->SetBinContent(1, .00500994);  //-- from AS 1D Gaussian fit
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(1, .00248577);   //-- from AS 1D Gaussian fit
	//fitParameterQuadrupolePlotVsCent[2]->SetBinContent(1, .00867822);   //from 6-parameter
    //fitParameterQuadrupolePlotVsCent[2]->SetBinError(1, .000925368);    //from 6-parameter
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(1, .0886915);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(1, .0191773);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(1, .332108);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(1, .135612);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(1, .289873);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(1, .0466348);
    
   
    fitParameterQuadrupolePlotVsCent[2]->SetBinContent(2, .00605651);
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(2, .00324884);
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(2, .0411947);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(2, .00347766);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(2, .673604);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(2, .0565171);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(2, 1.95405);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(2, .922342);
    
   
    
    /*fitParameterQuadrupolePlotVsCent[2]->SetBinContent(3, .00264);
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(3, .000272316);
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(3, .00751125);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(3, .00174301);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(3, .433461);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(3, .0522924);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(3, 1.8);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(3, .258073);*/ //non-zero quadrupole
    
    fitParameterQuadrupolePlotVsCent[2]->SetBinContent(3, 0.0);
    fitParameterQuadrupolePlotVsCent[2]->SetBinError(3,   .00257986);
    fitParameterJetAmpPlotVsCent[2]->SetBinContent(3,     .043952);
    fitParameterJetAmpPlotVsCent[2]->SetBinError(3,        .005208);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(3,    .711567);
    fitParameterSigmaPhiPlotVsCent[2]->SetBinError(3,       .0699830);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(3,     1.06627);
    fitParameterSigmaEtaPlotVsCent[2]->SetBinError(3,       .184991);
    
    //////////pt = 2.56
    lightFlavorResultsJetAmpPt1[2]->SetBinContent(1, .058);
    lightFlavorResultsJetAmpPt1[2]->SetBinContent(2, .043);
    lightFlavorResultsJetAmpPt1[2]->SetBinContent(3, .033);
    
    lightFlavorResultsEtaWidthPt1[2]->SetBinContent(1, .57);
    lightFlavorResultsEtaWidthPt1[2]->SetBinContent(2, 1.41);
    lightFlavorResultsEtaWidthPt1[2]->SetBinContent(3, 2.24);
    
    lightFlavorResultsPhiWidthPt1[2]->SetBinContent(1, .43);
    lightFlavorResultsPhiWidthPt1[2]->SetBinContent(2, .53);
    lightFlavorResultsPhiWidthPt1[2]->SetBinContent(3, .57);
    
    
    
    ///////////pt = 3.82
    lightFlavorResultsJetAmpPt2[2]->SetBinContent(1, .076);
    lightFlavorResultsJetAmpPt2[2]->SetBinContent(2, .046);
    lightFlavorResultsJetAmpPt2[2]->SetBinContent(3, .036);
    
    lightFlavorResultsEtaWidthPt2[2]->SetBinContent(1, .49);
    lightFlavorResultsEtaWidthPt2[2]->SetBinContent(2, .93);
    lightFlavorResultsEtaWidthPt2[2]->SetBinContent(3, 1.86);
    
    lightFlavorResultsPhiWidthPt2[2]->SetBinContent(1, .40);
    lightFlavorResultsPhiWidthPt2[2]->SetBinContent(2, .44);
    lightFlavorResultsPhiWidthPt2[2]->SetBinContent(3, .54);
    
    ///////////////////pt = 5.7
    lightFlavorResultsJetAmpPt3[2]->SetBinContent(1, .102);
    lightFlavorResultsJetAmpPt3[2]->SetBinContent(2, .051);
    lightFlavorResultsJetAmpPt3[2]->SetBinContent(3, .031);
    
    lightFlavorResultsEtaWidthPt3[2]->SetBinContent(1, .49);
    lightFlavorResultsEtaWidthPt3[2]->SetBinContent(2, .63);
    lightFlavorResultsEtaWidthPt3[2]->SetBinContent(3, 1.08);
    
    lightFlavorResultsPhiWidthPt3[2]->SetBinContent(1, .48);
    lightFlavorResultsPhiWidthPt3[2]->SetBinContent(2, .37);
    lightFlavorResultsPhiWidthPt3[2]->SetBinContent(3, .44);
    
    ///////////////////pt = 1.71
    lightFlavorResultsJetAmpPt4[2]->SetBinContent(1, .044);
    lightFlavorResultsJetAmpPt4[2]->SetBinContent(2, .035);
    lightFlavorResultsJetAmpPt4[2]->SetBinContent(3, .026);
    
    lightFlavorResultsEtaWidthPt4[2]->SetBinContent(1, .66);
    lightFlavorResultsEtaWidthPt4[2]->SetBinContent(2, 1.55);
    lightFlavorResultsEtaWidthPt4[2]->SetBinContent(3, 2.39);
    
    lightFlavorResultsPhiWidthPt4[2]->SetBinContent(1, .5);
    lightFlavorResultsPhiWidthPt4[2]->SetBinContent(2, .59);
    lightFlavorResultsPhiWidthPt4[2]->SetBinContent(3, .62);
    
    ///////////////////pt = 1.15
    lightFlavorResultsJetAmpPt5[2]->SetBinContent(1, .032);
    lightFlavorResultsJetAmpPt5[2]->SetBinContent(2, .025);
    lightFlavorResultsJetAmpPt5[2]->SetBinContent(3, .019);
    
    lightFlavorResultsEtaWidthPt5[2]->SetBinContent(1, .73);
    lightFlavorResultsEtaWidthPt5[2]->SetBinContent(2, 1.56);
    lightFlavorResultsEtaWidthPt5[2]->SetBinContent(3, 2.51);
    
    lightFlavorResultsPhiWidthPt5[2]->SetBinContent(1, .59);
    lightFlavorResultsPhiWidthPt5[2]->SetBinContent(2, .63);
    lightFlavorResultsPhiWidthPt5[2]->SetBinContent(3, .66);
    
    ///////////////////pt = .77
    lightFlavorResultsJetAmpPt6[2]->SetBinContent(1, .0243);
    lightFlavorResultsJetAmpPt6[2]->SetBinContent(2, .020);
    lightFlavorResultsJetAmpPt6[2]->SetBinContent(3, .016);
    
    lightFlavorResultsEtaWidthPt6[2]->SetBinContent(1, .8);
    lightFlavorResultsEtaWidthPt6[2]->SetBinContent(2, 1.60);
    lightFlavorResultsEtaWidthPt6[2]->SetBinContent(3, 2.74);
    
    lightFlavorResultsPhiWidthPt6[2]->SetBinContent(1, .69);
    lightFlavorResultsPhiWidthPt6[2]->SetBinContent(2, .70);
    lightFlavorResultsPhiWidthPt6[2]->SetBinContent(3, .72);
    
    
    
    double jetVolumeYield;
    double jetVolumeYieldError;
    TString jetVolumeYieldLabel = "Jet per-trigger yield ";
   
    double integralEtaDep[3] = { .362012, .957845, .870826};
   
    for(int ptBin = 2; ptBin < NUM_PT_BINS; ptBin++){
        
            str1 = jetVolumeYieldLabel + ptBinLabel[ptBin];
            jetVolumePlot[ptBin] = new TH1D(str1, str1, 4, 0, 4);
            jetVolumePlot[ptBin]->GetXaxis()->SetBinLabel(1, "pythia");
            jetVolumePlot[ptBin]->GetXaxis()->SetBinLabel(2, "50-100%");
            jetVolumePlot[ptBin]->GetXaxis()->SetBinLabel(3, "20-50%");
            jetVolumePlot[ptBin]->GetXaxis()->SetBinLabel(4, "0-20%");
        
        for(int centBin = 0; centBin < 3; centBin++){
         
            jetVolumeYield = ((nBarPrime[centBin]*fitParameterJetAmpPlotVsCent[2]->GetBinContent(centBin+1)*fitParameterSigmaPhiPlotVsCent[2]->GetBinContent(centBin+1)*integralEtaDep[centBin])/2.0)/(TMath::Sqrt(TMath::TwoPi()));
            jetVolumePlot[ptBin]->SetBinContent(centBin+2, jetVolumeYield);
            
            jetVolumeYieldError = (nBarPrime[centBin]/(2.0*TMath::Sqrt(TMath::TwoPi())))*getJetYieldError(fitParameterJetAmpPlotVsCent[2]->GetBinContent(centBin+1), integralEtaDep[centBin],fitParameterSigmaPhiPlotVsCent[2]->GetBinContent(centBin+1),
                                                   fitParameterJetAmpPlotVsCent[2]->GetBinError(centBin+1), fitParameterSigmaEtaPlotVsCent[2]->GetBinError(centBin+1),fitParameterSigmaPhiPlotVsCent[2]->GetBinError(centBin+1));
            
            //cout << jetVolumeYieldError << endl;
            
            jetVolumePlot[ptBin]->SetBinError(centBin+2, jetVolumeYieldError);
            
            cout << "# jet tracks per D0 in Centrality Bin: " << centBin <<"  = " << jetVolumeYield << " +/- " << jetVolumeYieldError << endl;
        
        }  
    }   
    
    
    //pythia -- 
    
    
    double pythiaNBar = 8.726;
    
    jetVolumeYield = (pythiaNBar*pythiaResultsJetAmp[2]->GetBinContent(1)*pythiaResultsPhiWidth[2]->GetBinContent(1)*pythiaResultsEtaWidth[2]->GetBinContent(1))/2.0;
    jetVolumePlot[2]->SetBinContent(1, jetVolumeYield);
            
    jetVolumeYieldError = (pythiaNBar/2.0)*getJetYieldError(pythiaResultsJetAmp[2]->GetBinContent(1), pythiaResultsEtaWidth[2]->GetBinContent(1),pythiaResultsPhiWidth[2]->GetBinContent(1),
                                           pythiaResultsJetAmp[2]->GetBinError(1), pythiaResultsEtaWidth[2]->GetBinError(1),pythiaResultsPhiWidth[2]->GetBinError(1));
            
    //cout << jetVolumeYieldError << endl;
            
    jetVolumePlot[2]->SetBinError(1, jetVolumeYieldError);
            
    cout << "# per-trigg yield PYTHIA  = " << jetVolumeYield << " +/- " << jetVolumeYieldError <<endl;
    
    TLegend * leg = new TLegend(.12, .88, .7, .6);
    leg->AddEntry(lightFlavorResultsJetAmpPt3[2], "LF 5.7 GeV/c");
    leg->AddEntry(lightFlavorResultsJetAmpPt1[2], "LF 2.56 GeV/c");
    //leg->AddEntry(lightFlavorResultsJetAmpPt2[2], "LF 3.82 GeV/c");
    
    //leg->AddEntry(lightFlavorResultsJetAmpPt4[2], "LF 1.71 GeV/c");
    leg->AddEntry(lightFlavorResultsJetAmpPt5[2], "LF 1.15 GeV/c");
    //leg->AddEntry(lightFlavorResultsJetAmpPt6[2], "LF .77 GeV/c");
    leg->AddEntry(pythiaResultsJetAmp[2],         "pythia D0-h 2-10 GeV/c");
    leg->AddEntry(fitParameterJetAmpPlotVsCent[2], "STAR AuAu 200 GeV");
    
    TCanvas * finalFitCorrGridCanvas = new TCanvas("", "", 4800, 1100);
    finalFitCorrGridCanvas->Divide(4,1);
    
    
    finalFitCorrGridCanvas->cd(1);
    
    fitParameterQuadrupolePlotVsCent[2]->SetMarkerStyle(20);
    fitParameterQuadrupolePlotVsCent[2]->SetMarkerSize(2);
    fitParameterQuadrupolePlotVsCent[2]->SetLineWidth(2);
    fitParameterQuadrupolePlotVsCent[2]->SetLineColor(1);
    fitParameterQuadrupolePlotVsCent[2]->Draw("E1");
    fitParameterPublishedQuadVsCent[2]->SetMarkerStyle(9);
    fitParameterPublishedQuadVsCent[2]->SetMarkerSize(2);
    fitParameterPublishedQuadVsCent[2]->SetMarkerColor(2);
    fitParameterPublishedQuadVsCent[2]->Draw("SAME P");
    finalFitCorrGridCanvas->cd(2);
    
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
    lightFlavorResultsJetAmpPt1[2]->SetMarkerColor(2);
    //lightFlavorResultsJetAmpPt1[2]->Draw("P");
    pythiaResultsJetAmp[2]->SetMarkerStyle(29);
    pythiaResultsJetAmp[2]->SetMarkerSize(5);
    pythiaResultsJetAmp[2]->SetMarkerColor(2);
    pythiaResultsJetAmp[2]->SetLineWidth(2);
    pythiaResultsJetAmp[2]->SetLineColor(2);
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
    lightFlavorResultsJetAmpPt5[2]->SetMarkerColor(11);
    //lightFlavorResultsJetAmpPt5[2]->Draw("SAME P");
    lightFlavorResultsJetAmpPt6[2]->SetMarkerStyle(9);
    lightFlavorResultsJetAmpPt6[2]->SetMarkerSize(2);
    lightFlavorResultsJetAmpPt6[2]->SetMarkerColor(12);
    //lightFlavorResultsJetAmpPt6[2]->Draw("SAME P");
    fitParameterJetAmpPlotVsCent[2]->SetMarkerStyle(22);
    fitParameterJetAmpPlotVsCent[2]->SetMarkerSize(3);
    fitParameterJetAmpPlotVsCent[2]->SetLineWidth(2);
    fitParameterJetAmpPlotVsCent[2]->SetLineColor(1);
    //fitParameterJetAmpPlotVsCent[2]->Draw("SAME E1"); 
   
    jetVolumePlot[2]->SetMarkerStyle(20);
    jetVolumePlot[2]->SetMarkerSize(2);
    jetVolumePlot[2]->SetLineWidth(2);
    jetVolumePlot[2]->SetLineColor(1);
    jetVolumePlot[2]->SetMarkerColor(1);
    jetVolumePlot[2]->SetStats(0);
    jetVolumePlot[2]->Draw("P E1");
    
    
    finalFitCorrGridCanvas->cd(3);
    fitParameterSigmaPhiPlotVsCent[2]->SetMarkerStyle(22);
    fitParameterSigmaPhiPlotVsCent[2]->SetMaximum(2.0);
    fitParameterSigmaPhiPlotVsCent[2]->SetMarkerSize(3);
    fitParameterSigmaPhiPlotVsCent[2]->SetLineWidth(2);
    fitParameterSigmaPhiPlotVsCent[2]->SetLineColor(1);
    fitParameterSigmaPhiPlotVsCent[2]->Draw("E1");
    lightFlavorResultsPhiWidthPt1[2]->SetMinimum(0);
    lightFlavorResultsPhiWidthPt1[2]->SetMaximum(2.0);
    lightFlavorResultsPhiWidthPt1[2]->SetMarkerStyle(9);
    lightFlavorResultsPhiWidthPt1[2]->SetMarkerSize(2);
    lightFlavorResultsPhiWidthPt1[2]->SetMarkerColor(2);
    lightFlavorResultsPhiWidthPt1[2]->Draw("SAME P");
    pythiaResultsPhiWidth[2]->SetMarkerStyle(29);
    pythiaResultsPhiWidth[2]->SetMarkerSize(5);
    pythiaResultsPhiWidth[2]->SetMarkerColor(2);
    pythiaResultsPhiWidth[2]->SetLineWidth(2);
    pythiaResultsPhiWidth[2]->SetLineColor(2);
    pythiaResultsPhiWidth[2]->Draw("SAME P E1");
    lightFlavorResultsPhiWidthPt2[2]->SetMarkerStyle(9);
    lightFlavorResultsPhiWidthPt2[2]->SetMarkerSize(2);
    lightFlavorResultsPhiWidthPt2[2]->SetMarkerColor(3);
    //lightFlavorResultsPhiWidthPt2[2]->Draw("SAME P");
    lightFlavorResultsPhiWidthPt3[2]->SetMarkerStyle(9);
    lightFlavorResultsPhiWidthPt3[2]->SetMarkerSize(2);
    lightFlavorResultsPhiWidthPt3[2]->SetMarkerColor(4);
    lightFlavorResultsPhiWidthPt3[2]->Draw("SAME P");
    lightFlavorResultsPhiWidthPt4[2]->SetMarkerStyle(9);
    lightFlavorResultsPhiWidthPt4[2]->SetMarkerSize(2);
    lightFlavorResultsPhiWidthPt4[2]->SetMarkerColor(6);
    //lightFlavorResultsPhiWidthPt4[2]->Draw("SAME P");
    lightFlavorResultsPhiWidthPt5[2]->SetMarkerStyle(9);
    lightFlavorResultsPhiWidthPt5[2]->SetMarkerSize(2);
    lightFlavorResultsPhiWidthPt5[2]->SetMarkerColor(11);
    lightFlavorResultsPhiWidthPt5[2]->Draw("SAME P");
    lightFlavorResultsPhiWidthPt6[2]->SetMarkerStyle(9);
    lightFlavorResultsPhiWidthPt6[2]->SetMarkerSize(2);
    lightFlavorResultsPhiWidthPt6[2]->SetMarkerColor(12);
    //lightFlavorResultsPhiWidthPt6[2]->Draw("SAME P");
    
    leg->Draw("SAME");
    
    
    finalFitCorrGridCanvas->cd(4);
    fitParameterSigmaEtaPlotVsCent[2]->SetMarkerStyle(22);
    fitParameterSigmaEtaPlotVsCent[2]->SetMarkerSize(3);
    fitParameterSigmaEtaPlotVsCent[2]->SetLineWidth(2);
    fitParameterSigmaEtaPlotVsCent[2]->SetLineColor(1);
    fitParameterSigmaEtaPlotVsCent[2]->SetMaximum(3.0);
    fitParameterSigmaEtaPlotVsCent[2]->Draw("E1");
    lightFlavorResultsEtaWidthPt1[2]->SetMinimum(0);
    lightFlavorResultsEtaWidthPt1[2]->SetMaximum(3.0);
    lightFlavorResultsEtaWidthPt1[2]->SetMarkerStyle(9);
    lightFlavorResultsEtaWidthPt1[2]->SetMarkerSize(2);
    lightFlavorResultsEtaWidthPt1[2]->SetMarkerColor(2);
    lightFlavorResultsEtaWidthPt1[2]->Draw("SAME P");
    pythiaResultsEtaWidth[2]->SetMarkerStyle(29);
    pythiaResultsEtaWidth[2]->SetMarkerSize(5);
    pythiaResultsEtaWidth[2]->SetMarkerColor(2);
    pythiaResultsEtaWidth[2]->SetLineWidth(2);
    pythiaResultsEtaWidth[2]->SetLineColor(2);
    pythiaResultsEtaWidth[2]->Draw("SAME P E1");
    lightFlavorResultsEtaWidthPt2[2]->SetMarkerStyle(9);
    lightFlavorResultsEtaWidthPt2[2]->SetMarkerSize(2);
    lightFlavorResultsEtaWidthPt2[2]->SetMarkerColor(3);
    //lightFlavorResultsEtaWidthPt2[2]->Draw("SAME P");
    lightFlavorResultsEtaWidthPt3[2]->SetMarkerStyle(9);
    lightFlavorResultsEtaWidthPt3[2]->SetMarkerSize(2);
    lightFlavorResultsEtaWidthPt3[2]->SetMarkerColor(4);
    lightFlavorResultsEtaWidthPt3[2]->Draw("SAME P");
    lightFlavorResultsEtaWidthPt4[2]->SetMarkerStyle(9);
    lightFlavorResultsEtaWidthPt4[2]->SetMarkerSize(2);
    lightFlavorResultsEtaWidthPt4[2]->SetMarkerColor(6);
    //lightFlavorResultsEtaWidthPt4[2]->Draw("SAME P");
    lightFlavorResultsEtaWidthPt5[2]->SetMarkerStyle(9);
    lightFlavorResultsEtaWidthPt5[2]->SetMarkerSize(2);
    lightFlavorResultsEtaWidthPt5[2]->SetMarkerColor(11);
    lightFlavorResultsEtaWidthPt5[2]->Draw("SAME P");
    lightFlavorResultsEtaWidthPt6[2]->SetMarkerStyle(9);
    lightFlavorResultsEtaWidthPt6[2]->SetMarkerSize(2);
    lightFlavorResultsEtaWidthPt6[2]->SetMarkerColor(12);
    //lightFlavorResultsEtaWidthPt6[2]->Draw("SAME P");
    
    
    finalFitCorrGridCanvas->cd(5);
    //jetVolumePlot[2]->SetMarkerStyle(20);
    //jetVolumePlot[2]->Draw("P");
    
    finalFitCorrGridCanvas->SaveAs("final_output_fit_parameters_from_different_functions.png");

}


double getJetYieldError(double A, double sigEta, double sigPhi, double eA, double eSigEta, double eSigPhi){

    double arg = ((eA/A)*(eA/A)) + ((eSigEta/sigEta)*(eSigEta/sigEta)) + ((eSigPhi/sigPhi)*(eSigPhi/sigPhi));
    double product = TMath::Abs((A*sigEta*sigPhi));
    
    //cout << "prduct value: " << product << endl;
    //cout << "arg value: " << arg << endl;
    //cout << "sqrt of arg: " << TMath::Sqrt(arg) << endl;
    
    return TMath::Abs((A*sigEta*sigPhi))*TMath::Sqrt(arg);
    
}