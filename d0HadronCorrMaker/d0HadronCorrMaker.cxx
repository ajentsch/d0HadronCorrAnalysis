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

#include "d0HadronCorrMaker.h"
#include "histogramExtractAndSave.cxx"


/******************************

Author: Alex Jentsch

Most Recent Update: 2/22/2018

The original code is named corrHistogramMaker. It can still be used to produce the standard results, without the D* correction.

This code will replace that older named version, as it contains a cleaner interface, less overall functionality (for speed and for ease of use), and contains the D* correction.

To be added: easier functionality for error input.


******************************/


Bool_t reject;
Double_t fline(Double_t *x, Double_t *par)
{
    if ( (reject && (x[0] > 1.81 && x[0] < 1.91)) || (reject && (x[0] > 1.6 && x[0] < 1.7)) ){     //|| (reject && (x[0] > 2.0 && x[0] < 2.1))
    //if ( (reject && (x[0] > 1.81 && x[0] < 1.91)) ){
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0];
}

Double_t fexpo(Double_t *x, Double_t *par)
{
    if ( (reject && (x[0] > 1.81 && x[0] < 1.91)) || (reject && (x[0] > 1.6 && x[0] < 1.7)) ){     //|| (reject && (x[0] > 2.0 && x[0] < 2.1))
      TF1::RejectPoint();
      return 0;
   }
   return exp(par[0] + par[1]*x[0]);
}

double calcScaleFactorforFit(double sigma)
{
    double preFactor = TMath::Sqrt(TMath::PiOver2());
    
    return preFactor*sigma*.5;
    
}
     
//Double_t SecH(     

int d0HadronCorrMaker(){
  

  TString inputFileName   = "systematic_test_13"; //MAIN FILE NAME FOR INPUT

//-------------------------------------------------------------
//BEGIN VARIABLES SECTION
//-------------------------------------------------------------  
    int NUM_PHI_BINS = 12;
    int NUM_ETA_BINS = 13;
    int NUM_CENT_BINS = 16;
    int NUM_VZ_BINS   = 10;
    int NUM_PT_BINS   = 3;                // the last pt-bin is currently the integrated bin. This needs to be reconfigured to have the integrated bin first so I am not pulling a bin out from the middle.
    double NUM_FILES_COMBINED = 1;
    double ETA_RANGE = 1.69;
    double phiBinShift = (TMath::Pi()/12.0);
    
    int FIRST_PT_BIN = 2;
    
    double SBLScaleFactor = .5;
    double SBRScaleFactor = .5;
    
    //-----------------------
    //file naming and input/output
    //-----------------------

    TString rootLabel          = ".root";
    TString inputFileNameFull  = inputFileName + rootLabel;  
    TString outputFileName     = "d0HadronCorrMaker_OUTPUT_";
    outputFileName     += inputFileNameFull;
    
    TFile *file = new TFile(inputFileNameFull); //Root file with raw histograms to use for calculations
    
    cout << "Root Data File open. " << endl << endl;
    
    //histogramExtractAndSave(file, path);    //--THIS PRINTS A BUNCH OF QA STUFF -- TURNED OFF IN STPICOD0ANAMAKER FOR NOW
    
    TFile *output = new TFile(outputFileName, "RECREATE"); //root file to store calculated and scaled histograms
    
    cout << "Root Output File open. " << endl << endl;
    
    bool PRINT_SUB_HISTOS_RAW = false;
    bool PRINT_SUB_HISTOS_SCALED = false;
    bool PRINT_SUB_HISTOS_DEL_RHO = false;
    bool PRINT_DELRHO_OVER_REF_HISTOS = true;
    
    TString fileType    = ".png"; //file type for output histogram pictures
    TString fileTypeEps = ".eps"; //file type for progess report
    TString fileTypePDF = ".pdf";
    
    TString binLabelVz[10] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
    TString binLabelCent[16] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"};
    
    TString binLabelPt[6]    = {"0", "1", "2", "3", "4", ""};
    TString PtBinLabel[6]   = {"_PtBin_", "_PtBin_", "_PtBin_", "_PtBin_", "_PtBin_", ""};
    TString PtSubFolders[6] = {"PtBin0/", "PtBin1/", "PtBin2/", "PtBin3/", "PtBin4/", "PtIntegrated/"};
    
    TH1D*   oneDUSHistos[6][4];     
    TString oneDUSStrings[6] = {"D0_US_invMass_Pt_Bin_0", "D0_US_invMass_Pt_Bin_1", "D0_US_invMass_Pt_Bin_2", "D0_US_invMass_Pt_Bin_3", "D0_US_invMass_Pt_Bin_4", "unlikeSign"};
                                              
    TH1D*   oneDLSHistos[6][4];     
    TString oneDLSStrings[6] = {"LS_invMass_Pt_Bin_0", "LS_invMass_Pt_Bin_1", "LS_invMass_Pt_Bin_2", "LS_invMass_Pt_Bin_3", "LS_invMass_Pt_Bin_4", "LikeSignBG"};
                                              
    TH1D*   oneDScaledLSHistos[6][4];     
    TString oneDScaledLSStrings[6] = {"ScaledLS_PtBin_0", "ScaledLS_PtBin_1", "ScaledLS_PtBin_2", "ScaledLS_PtBin_3", "ScaledLS_PtBin_4", "ScaledLS"};
                                                  
    TH1D*   oneDSubtractedInvMassHistos[6][4];                                              
    TString oneDSubtractedInvMassStrings[6] = {"D0_Minus_Scaled_LS_BG_PtBin_0", "D0_Minus_Scaled_LS_BG_PtBin_1", "D0_Minus_Scaled_LS_BG_PtBin_2", "D0_Minus_Scaled_LS_BG_PtBin_3", "D0_Minus_Scaled_LS_BG_PtBin_4","D0_Minus_Scaled_LS_BG"}; 
    
    TH1D*   oneDSubtractedInvMassHistosFunction[6][4];                                              
    TString oneDSubtractedInvMassStringsFunction[6] = {"US_Minus_Expo_Fit_BG_PtBin_0", "US_Minus_Expo_Fit_BG_PtBin_1", "US_Minus_Expo_Fit_BG_PtBin_2", "US_Minus_Expo_Fit_BG_PtBin_3", "US_Minus_Expo_Fit_BG_PtBin_4", "US_Minus_Expo_Fit_BG"}; 
                                                                                                                          
    TString binLabelCentralityClass[3] = {"_Peripheral", "_MidCentral", "_Central"};
    
    
    double integralSideBandCounts[2];
    double integralCounts[2];
   
    //-----------------------------------
    //Correlation variables and arrays
    //-----------------------------------
    
    // US histograms
    TH2D* sibCorrBin[4][6][10][16];    //Raw US sibling histograms binned by centrality
    TH2D* mixCorrBin[4][6][10][16];    //Raw US mixed histograms binned by centrality
    TH2D* scaledMixCorrBin[4][6][10][16];  //storage for the eventual scaled mixed US histograms  
    TH2D* sibMinusScaledMix[4][6][10][16]; // Raw US Sibling minus scaled mixed                   
    TH2D* delRhoOverRhoRefBin[4][6][10][16];
    
    TH2D* rootErrorsOnDelRhoOverRhoRef[4][6][10][16];
    
    //ERRORS STORED HERE--ALL ERRORS STORED AS VARIANCES -- SQUARE ROOT MUST BE TAKEN FOR ACTUAL ERROR
    
    double sibCorrBinStatErrors[4][6][10][16][9][12];    
    double mixCorrBinStatErrors[4][6][10][16][9][12];
    double scaledMixCorrBinStatErrors[4][6][10][16][9][12];
    double sibMinusScaledMixStatErrors[4][6][10][16][9][12];
    double delRhoOverRhoRefBinStatErrors[4][6][10][16][9][12];
    //------------------------------------------------------------------
    
    TH2D* sibCorrBinVzInt[4][6][16];
    TH2D* mixCorrBinVzInt[4][6][16];
    TH2D* scaledMixCorrBinVzInt[4][6][16];  
    TH2D* sibMinusScaledMixVzInt[4][6][16]; 
    TH2D* delRhoOverRhoRefBinVzInt[4][6][16];
    
    TH1D* sibCorrBinPhiProj[4][6][10][16];    //Raw US sibling histograms binned by centrality
    TH1D* mixCorrBinPhiProj[4][6][10][16];    //Raw US mixed histograms binned by centrality
    TH1D* scaledMixCorrBinPhiProj[4][6][10][16];  //storage for the eventual scaled mixed US histograms
    TH1D* sibMinusScaledMixPhiProj[4][6][10][16]; // Raw US Sibling minus scaled mixed
    TH1D* delRhoOverRhoRefBinPhiProj[4][6][10][16];
    
    TH1D* sibCorrBinPhiProjVzInt[4][6][16];
    TH1D* mixCorrBinPhiProjVzInt[4][6][16];
    TH1D* scaledMixCorrBinPhiProjVzInt[4][6][16];  
    TH1D* sibMinusScaledMixPhiProjVzInt[4][6][16]; 
    TH1D* delRhoOverRhoRefBinPhiProjVzInt[4][6][16];
  
    
    TH1D* sibCorrUSPhiProjMinusLSPhiProjVzInt[6][16];
   
   
    TH2D* sideBandAverage[6][16];
    TH1D* sideBandAveragePhiProj[6][16];    
    TH2D* fullySubtractedCorrSideBandCent[6][16];
    TH1D* fullySubtractedCorrSideBandCentPhiProj[6][16];
    
    TH2D* delRhoOverRhoRefBinVzIntCentInt[4][6][3];
    TH1D* delRhoOverRhoRefBinVzIntCentIntPhiProj[4][6][3];
    
    TH2D* fullySubtractedCorrLSCent[6][16];
    TH1D* fullySubtractedCorrLSCentPhiProj[6][16];
    
    //dStarStuff
    TH2D* mixCorrCentrality[5][3];
	TH2D* sibCorrCentrality[5][3];
    
    //other testing
    TH2D* SBLeftCorrCentrality[5][3];
    TH2D* SBRightCorrCentrality[5][3];
    
    //TH2D* mixCorrUSVzInt[3][16];
    //TH2D* mixCorrLSVzInt[3][16];
    
    //Correlations with integrated centralities
    TH2D* fullySubtractedCorrCustomCentralityUS[6][3];
    TH2D* fullySubtractedCorrCustomCentralityLS[6][3];
    TH1D* fullySubtractedCorrCustomCentralityLSPhiProj[6][3];
    TH2D* fullySubtractedCorrCustomCentralitySideBand[6][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandPhiProj[6][3];
    double VzWeightFactorCustom[3][10];
    double VzWeightFactorCustomTotal;
    
    TString fullySubtractedLabelCustomCentrality = "FullSubtractedCorr_";
    TString centralityBin[3] = {"Peripheral", "MidCentral", "Central"};
    
    TString SibCorrLabels[4]          = {"Sibling_SideBandLeft_correlation", "Sibling_US_correlation", "Sibling_SideBandRight_correlation", "Sibling_LS_correlation"};
    TString MixCorrLabels[4]          = {"Mixed_SideBandLeft_correlation", "Mixed_US_correlation", "Mixed_SideBandRight_correlation", "Mixed_LS_correlation"};
    TString ScaledMixCorrLabels[4]     = {"Scaled_Mixed_SideBandLeft_correlation", "Scaled_Mixed_US_correlation", "Scaled_Mixed_SideBandRight_correlation", "Scaled_Mixed_LS_correlation"};
    TString SibMinusScaledMixLabels[4] = {"Sib_Minus_Scaled_Mix_SideBandLeft_correlation", "Sib_Minus_Scaled_Mix_US_correlation", "Sib_Minus_Scaled_Mix_SideBandRight_correlation", "Sib_Minus_Scaled_Mix_LS_correlation"};
    TString delRhoOverRhoRefLabels[4]  = {"delRho_over_RhoRef_SideBandLeft_correlation", "delRho_over_RhoRef_US_correlation", "delRho_over_RhoRef_SideBandRight_correlation", "delRho_over_RhoRef_LS_correlation"};
    TString errorsFromRootDelRhoOverRhoRef[4]  = {"Errors_delRho_over_RhoRef_SideBandLeft_correlation", "Errors_delRho_over_RhoRef_US_correlation", "Errors_delRho_over_RhoRef_SideBandRight_correlation", "Errors_delRho_over_RhoRef_LS_correlation"};
   
    TString sibLabel = "Sibling_";
    TString mixLabel = "Mixed_";
    TString diffLabel = "Diff_";
    
    TString subtractedCorrLabels[2] = {"US_Minus_SidebandAverage_Corr", "US_Minus_LS_Corr"};
    
    TString sideBandAverageLabel = "sideBandAverage";
    
    TString phiProj = "Phi_projection_";
    TString SSEtaProj = "Same_Side_Eta_Projection_";
    TString ASEtaProj = "Away_Side_Eta_Projection_";
    
    
    TString USMinusLSSubtractedLabel1 = "US_Minus_LS_Subtracted_CentBin_";
    
    TString withFitLabel = "_With_Fit";
    
    //TString delRhoOverRhoRefUSMinusLSLabel1 = "delRho_over_RhoRef_US_Minus_Scaled_LS_";
    
    TString centBinLabel = "_CentBin_";
    TString VzBinLabel = "_VzBin_";
    TString VzIntLabel = "_Vz_Integrated_";
    TString bandLabel[4] = {"SideBandLeft", "US", "SideBandRight", "LS"};

    double numSib = 0;
    double numMix = 0;
    double numMixUS = 0;
    double numMixLS = 0;
    
    double integralRawMixHistos[4][6][10][16];
    double integralRawMixHistosVzInt[4][6][16];
    //double integralRawMixHistos[4][4][10][16];
    double integralRawMixHistosVzIntCentInt[4][6][3];
   
    double scaleFactorVzCent = 0;
    
    double scaleFactorVz = 0;
    double scaleFactorLS = 0;
    
    
    double LSScaleFactor[6][6] = {{1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, {1.0,1.0,1.0,1.0,1.0,1.0}};
    double ScalingRatio[3][16];
    
    double BOverSPlusBLS[6][4];
    double BOverSPlusBFit[6][4];
    double SPlusBOverSLS[6][4];
    double SPlusBOverSFit[6][4];
    
    double totalPairsSibling = 0;
    double integralError[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
    //fully subtraced US minus LS
    
    TString delRhoOverRhoRefUSMinusLSLabel1 = "US_delRhoOverRho_Minus_LS_delRhoOverRho_";
    
    double v2Array[5][3];
    double v2ArrayStatError[5][3];
    
    //----------------------------
    //General use variables
    //----------------------------
   
    TString str1;
    TString str2;
    TString str3;
    TString str4;
    

//-------------------------------------------------------------
//END VARIABLES SECTION
//-------------------------------------------------------------    
    
    

    double countsPeakRegionUS = 0;
    double countsPeakRegionLS = 0;
    double countsSideBandUS   = 0;
    double countsSideBandLS   = 0;
    double sideBandLow      = 2.0;
    double sideBandHigh     = 2.1;
    double massLow          = 1.825;   //1
    double massHigh         = 1.90;
    
    double USSideBandLeftLow   = 1.72;  //1.62
    double USSideBandLeftHigh  = 1.80;  //1.7
    double USSideBandRightLow  = 1.92;   //2.0
    double USSideBandRightHigh = 2.0;   //2.1 
    
    double integralLeftSB  = 0;
    double integralRightSB = 0;
    double SBAverage       = 0;
	double sqrtSPlusB       = 1;
    double integralFit      = 0;
	double DMesonSig        = 0;
    double integral         = 0;
    double signalOverBG     = 0;
    double optimizationFactor = 0;
    //double LSScaleFactor    = 0;
    
    Int_t binMassLowUS;
    Int_t binMassHighUS;
    Int_t binMassLowLS;
    Int_t binMassHighLS;
        
    TString title;
    TString tmp;
        
    Double_t par[5];
    Double_t parBG[2];
        
    TF1 * g1 = new TF1("m1", "gaus", 1.8, 1.9);
    TF1 * e1 = new TF1("m2", "expo", 1.6, 1.78);
    
    //TF1 * eL = new TF1("m3", "expo", 1.69, 1.81);
    //TF1 * eR = new TF1("m4", "expo", 1.92, 2.0);
    
    TF1 * BGFunc = new TF1("BGfunc", "expo", 1.6, 2.1);
    
    TF1 * fun = new TF1("fun","gaus(0)+expo(3)",1.6,2.1);
    
    TF1 * finalGaussian = new TF1("final Guassian", "gaus", 1.82, 1.9);
    
    double par0;
    double par1;
    double par2;
    double par3;

    TAxis* xAxisUS;
    TAxis* xAxisLS;
    
    int LabelVersion = 0;
    
//------------------------------------------------------------------
//Code section to write information about this dataset to file
//------------------------------------------------------------------    
    /* //Fill Cuts histogram here for producing output text file
    histOfCuts->SetBinContent(1, D0InvMassLow); //D0InvMassLow
    histOfCuts->SetBinContent(2, D0InvMassHigh); //D0InvMassHigh
    histOfCuts->SetBinContent(3, USSideBandLeftLow); //USSideBandLeftLow
    histOfCuts->SetBinContent(4, USSideBandLeftHigh); //USSideBandLeftHigh
    histOfCuts->SetBinContent(5, USSideBandRightLow); //USSideBandRightLow
    histOfCuts->SetBinContent(6, USSideBandRightHigh); //USSideBandRightHigh
    histOfCuts->SetBinContent(7, d0PtLow); //d0PtLow 
    histOfCuts->SetBinContent(8, d0PtHigh); //d0PtHigh
    histOfCuts->SetBinContent(9, d0DaughterPionPtMin); //d0DaughterPionPtMin
    histOfCuts->SetBinContent(10, d0DaughterKaonPtMin); //d0DaughterKaonPtMin
    histOfCuts->SetBinContent(11, d0DecayLengthMax); //d0DecayLengthMax
    
    histOfCuts->SetBinContent(12, d0TopologicalCutArray[0][0]); //d0DecayLengthMin
    histOfCuts->SetBinContent(13, d0TopologicalCutArray[0][1]); //daughterDCA
    histOfCuts->SetBinContent(14, d0TopologicalCutArray[0][2]); //d0DCAtoPV
    histOfCuts->SetBinContent(15, d0TopologicalCutArray[0][3]); //pionDCA
    histOfCuts->SetBinContent(16, d0TopologicalCutArray[0][4]); //kaonDCA
    
    histOfCuts->SetBinContent(17, d0TopologicalCutArray[1][0]); //d0DecayLengthMin
    histOfCuts->SetBinContent(18, d0TopologicalCutArray[1][1]); //daughterDCA
    histOfCuts->SetBinContent(19, d0TopologicalCutArray[1][2]); //d0DCAtoPV
    histOfCuts->SetBinContent(20, d0TopologicalCutArray[1][3]); //pionDCA
    histOfCuts->SetBinContent(21, d0TopologicalCutArray[1][4]); //kaonDCA
    
    histOfCuts->SetBinContent(22, d0TopologicalCutArray[2][0]); //d0DecayLengthMin
    histOfCuts->SetBinContent(23, d0TopologicalCutArray[2][1]); //daughterDCA
    histOfCuts->SetBinContent(24, d0TopologicalCutArray[2][2]); //d0DCAtoPV
    histOfCuts->SetBinContent(25, d0TopologicalCutArray[2][3]); //pionDCA
    histOfCuts->SetBinContent(26, d0TopologicalCutArray[2][4]); //kaonDCA
    
    histOfCuts->SetBinContent(27, d0TopologicalCutArray[3][0]); //d0DecayLengthMin
    histOfCuts->SetBinContent(28, d0TopologicalCutArray[3][1]); //daughterDCA
    histOfCuts->SetBinContent(29, d0TopologicalCutArray[3][2]); //d0DCAtoPV
    histOfCuts->SetBinContent(30, d0TopologicalCutArray[3][3]); //pionDCA
    histOfCuts->SetBinContent(31, d0TopologicalCutArray[3][4]); //kaonDCA
    
    histOfCuts->SetBinContent(32, d0TopologicalCutArray[4][0]); //d0DecayLengthMin
    histOfCuts->SetBinContent(33, d0TopologicalCutArray[4][1]); //daughterDCA
    histOfCuts->SetBinContent(34, d0TopologicalCutArray[4][2]); //d0DCAtoPV
    histOfCuts->SetBinContent(35, d0TopologicalCutArray[4][3]); //pionDCA
    histOfCuts->SetBinContent(36, d0TopologicalCutArray[4][4]); //kaonDCA
    
    
    histOfCuts->SetBinContent(37, hadronPtMin); //hadronPtMin
    histOfCuts->SetBinContent(38, hadronPtMax); //hadronPtMax
    histOfCuts->SetBinContent(39, BUFFER_SIZE); //BUFFER_SIZE
    histOfCuts->SetBinContent(40, NUM_PHI_BINS); //NUM_PHI_BINS
    histOfCuts->SetBinContent(41, NUM_ETA_BINS); //NUM_ETA_BINS
    histOfCuts->SetBinContent(42, 0); // Require HFT
    histOfCuts->SetBinContent(43, trackChi2max); // chi2 cut
    histOfCuts->SetBinContent(44, trackDCAtoPvtx); // DCA cut
    histOfCuts->SetBinContent(49, 1); //Scale factor counter to produce text file
    histOfCuts->SetBinContent(50, 0); //USE TOF STRICTLY
    histOfCuts->SetBinContent(51, 0); //USE TOF HYBRID
    histOfCuts->SetBinContent(52, 0); //USE DOUBLE MIS-PID PROTECTION
    
    
    */
    
        double numEvents = 0;
        
        numEvents = (((TH1I*) file->Get("number_of_events_used"))->GetBinContent(2))/(1000000);
        
        ofstream cutFile;
        TString cutFileName = "Cuts_and_data_information.txt";
        //TString cutOutputFile = path + cutFileName;
        cutFile.open (cutFileName);
        
        cutFile << "Important information about this data run" << endl << endl;
        cutFile << "Number of Events: " << numEvents << "M Events" << endl << endl;
        
        NUM_FILES_COMBINED = ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(49);
        
        ((TH1D*) file->Get("HistOfCuts"))->Scale(1/NUM_FILES_COMBINED);
        
        cutFile << "Trigger D0 Cuts" << endl << endl;
        cutFile << "D0 InvMass Signal Band -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(1) << "\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(2) << endl;
        cutFile << "D0 SideBandLeft   Band -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(3) << "\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(4) << endl;
        cutFile << "D0 SideBandRight  Band -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(5) << "\t\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(6) << endl;
        cutFile << "D0 Pt Cuts             -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(7) << "\t\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(8) << endl;
        cutFile << "D0 DecayLengthMax(no upper bound)      -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(11) << endl;
        
        cutFile << "D0 DaughterPtCut       -\t" << "Pi : " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(9) << "\t" << "K  : " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(10) << endl;
        
        cutFile << "Use TOF STRICT: -\t" <<  ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(50) << endl; 
        cutFile << "Use TOF HYBRID: -\t" <<  ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(51) << endl; 
        cutFile << "Use DOUBLE MIS_PID_PROTECTION: -\t" <<  ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(52) << endl; 
        
        cutFile << endl << "_____________________________________________________________" << endl;
        
        cutFile << endl << "D0 pt 0-1 GeV/c cuts" << endl;
        
        cutFile << "D0 DecayLengthMinimum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(12) << endl;
        cutFile << "D0 DaughterDCAMaximum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(13) << endl;
        cutFile << "D0 DCA to PV Maximum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(14) << endl;
        cutFile << "Daughter Pion DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(15) << endl;
        cutFile << "Daughter Kaon DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(16) << endl;
        
        cutFile << endl << "_____________________________________________________________" << endl;
        
        cutFile << endl << "D0 pt 1.2-2 GeV/c cuts" << endl;
        
        cutFile << "D0 DecayLengthMinimum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(17) << endl;
        cutFile << "D0 DaughterDCAMaximum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(18) << endl;
        cutFile << "D0 DCA to PV Maximum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(19) << endl;
        cutFile << "Daughter Pion DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(20) << endl;
        cutFile << "Daughter Kaon DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(21) << endl;
        
        cutFile << endl <<"_____________________________________________________________"<< endl;
        
        cutFile << endl << "D0 pt 2-3 GeV/c cuts" << endl;
        
        cutFile << "D0 DecayLengthMinimum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(22) << endl;
        cutFile << "D0 DaughterDCAMaximum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(23) << endl;
        cutFile << "D0 DCA to PV Maximum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(24) << endl;
        cutFile << "Daughter Pion DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(25) << endl;
        cutFile << "Daughter Kaon DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(26) << endl;
        
        cutFile << endl <<"_____________________________________________________________"<< endl;
        
        cutFile << endl << "D0 pt 3-4 GeV/c cuts" << endl;
        
        cutFile << "D0 DecayLengthMinimum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(27) << endl;
        cutFile << "D0 DaughterDCAMaximum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(28) << endl;
        cutFile << "D0 DCA to PV Maximum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(29) << endl;
        cutFile << "Daughter Pion DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(30) << endl;
        cutFile << "Daughter Kaon DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(31) << endl;
        
        cutFile << endl << "_____________________________________________________________" << endl;
        
        cutFile << endl << "D0 pt 4-10 GeV/c cuts" << endl;
        
        cutFile << "D0 DecayLengthMinimum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(32) << endl;
        cutFile << "D0 DaughterDCAMaximum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(33) << endl;
        cutFile << "D0 DCA to PV Maximum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(34) << endl;
        cutFile << "Daughter Pion DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(35) << endl;
        cutFile << "Daughter Kaon DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(36) << endl;
        
        
        
        cutFile << endl << "------------------------------------------" << endl << "Associated hadron Cuts" << endl << endl;
        cutFile << "associated hadron Pt   -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(37) << "\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(38) << endl;
        cutFile << "Num of events used to mix  -\t" << ": " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(39) << endl;
        cutFile << "HFT Tracks only??  -\t" << "YES(1), NO(0): " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(42) << endl;
        cutFile << "Track Chi2 Max:  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(43) << endl;
        cutFile << "Track DCA Max:  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(44) << endl;
        
        cutFile.close();
        
    
//---------------------------------------------    
//BEGIN INVARIANT MASS HISTOGRAM INFORMATION   
//---------------------------------------------
    
    ofstream invariantMassInfo;
    TString invariantMassInfoName = "Signal_and_background_information.txt";
    invariantMassInfo.open (invariantMassInfoName);
    
    double sigPerBin[5][3];
    double dMesonSigArray[3];
    
    
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){  //Initialize the histograms 
    
        for(int i = 0; i < 3; i++){
        
            //if(k == 3) { LabelVersion = 1; }
       
            tmp = oneDUSStrings[k] + binLabelCentralityClass[i];
            oneDUSHistos[k][i] = (TH1D*)file->Get(tmp);     //US spectrum
            tmp = oneDLSStrings[k] + binLabelCentralityClass[i];
            oneDLSHistos[k][i] = (TH1D*)file->Get(tmp);     //LS spectrum
        
            oneDScaledLSHistos[k][i] = (TH1D*)oneDLSHistos[k][i]->Clone();  //clone of LS spectrum to be scaled   
            tmp = oneDScaledLSStrings[k] + binLabelCentralityClass[i];
            oneDScaledLSHistos[k][i]->SetTitle(oneDScaledLSStrings[k]);        
            tmp = oneDSubtractedInvMassStrings[k] + binLabelCentralityClass[i];
            oneDSubtractedInvMassHistos[k][i] = new TH1D(tmp, tmp, 50, 1.6, 2.1); //final spectrum -- US minus LS (scaled)
       
            oneDUSHistos[k][i]->GetXaxis()->SetTitle("Invariant Mass [GeV/c^{2}]");
            oneDUSHistos[k][i]->GetYaxis()->SetTitle("counts");
            oneDUSHistos[k][i]->GetYaxis()->SetTitleOffset(1.4);
            oneDLSHistos[k][i]->GetXaxis()->SetTitle("Invariant Mass [GeV/c^{2}]");
            oneDLSHistos[k][i]->GetYaxis()->SetTitle("counts");
            oneDLSHistos[k][i]->GetYaxis()->SetTitleOffset(1.4);
            oneDScaledLSHistos[k][i]->GetXaxis()->SetTitle("Invariant Mass [GeV/c^{2}]");
            oneDScaledLSHistos[k][i]->GetYaxis()->SetTitle("counts");
            oneDScaledLSHistos[k][i]->GetYaxis()->SetTitleOffset(1.4);
            oneDSubtractedInvMassHistos[k][i]->GetXaxis()->SetTitle("Invariant Mass [GeV/c^{2}]");
            oneDSubtractedInvMassHistos[k][i]->GetYaxis()->SetTitle("counts");
            oneDSubtractedInvMassHistos[k][i]->GetYaxis()->SetTitleOffset(1.4);
        
            xAxisUS = oneDUSHistos[k][i]->GetXaxis();
            xAxisLS = oneDLSHistos[k][i]->GetXaxis();
        
            binMassLowUS  = oneDUSHistos[k][i]->GetXaxis()->FindBin(sideBandLow);            //get normalization factors to scale the LS distribution
            binMassHighUS = oneDUSHistos[k][i]->GetXaxis()->FindBin(sideBandHigh);
            binMassLowLS = oneDLSHistos[k][i]->GetXaxis()->FindBin(sideBandLow);            //get normalization factors to scale the LS distribution
            binMassHighLS = oneDLSHistos[k][i]->GetXaxis()->FindBin(sideBandHigh);
           
            countsSideBandUS = oneDUSHistos[k][i]->Integral(binMassLowUS, binMassHighUS);
            
            binMassLowUS  = oneDUSHistos[k][i]->GetXaxis()->FindBin(USSideBandLeftLow);        //calculate the integral number in the LEFT SB   
            binMassHighUS = oneDUSHistos[k][i]->GetXaxis()->FindBin(USSideBandLeftHigh);
            integralLeftSB   = oneDUSHistos[k][i]->Integral(binMassLowUS, binMassHighUS);
            
            binMassLowUS  = oneDUSHistos[k][i]->GetXaxis()->FindBin(USSideBandRightLow);       //calculate the integral number in the RIGHT SB     
            binMassHighUS = oneDUSHistos[k][i]->GetXaxis()->FindBin(USSideBandRightHigh);
            integralRightSB  = oneDUSHistos[k][i]->Integral(binMassLowUS, binMassHighUS);
            
            SBAverage = .5*(integralLeftSB + integralRightSB);
            
            countsSideBandLS = oneDLSHistos[k][i]->Integral(binMassLowLS, binMassHighLS);
        
            //Calculate peak region values
        
            binMassLowUS  = xAxisUS->FindBin(massLow);            
            binMassHighUS = xAxisUS->FindBin(massHigh);
            
            cout << "bin min: " << binMassLowUS << "    bin max: " << binMassHighUS << endl;
            binMassLowLS  = xAxisLS->FindBin(massLow);            
            binMassHighLS = xAxisLS->FindBin(massHigh);
        
            oneDUSHistos[k][i]->Sumw2();
            oneDUSHistos[k][i]->SetMarkerStyle(20);
            oneDLSHistos[k][i]->Sumw2();
            oneDLSHistos[k][i]->SetMarkerStyle(20);
            oneDSubtractedInvMassHistos[k][i]->Sumw2();
            oneDSubtractedInvMassHistos[k][i]->SetMarkerStyle(20);
        
           /*  oneDUSHistos[k][i]->Draw();
                
            tmp = oneDUSStrings[k] + binLabelCentralityClass[i] + fileType;
            c->SaveAs(tmp); */
        
            //oneDUSHistos[i]->SetTitle("");
            //oneDUSHistos[i]->Draw();
            //tmp = path + oneDUSStrings[i] + fileTypeEps;
            //c->SaveAs(tmp);
        
            LSScaleFactor[k][i] = countsSideBandUS/countsSideBandLS;
        
            oneDScaledLSHistos[k][i]->Scale(LSScaleFactor[k][i]);            //normalize LS spectrum here
        
            countsPeakRegionUS = oneDUSHistos[k][i]->Integral(binMassLowUS, binMassHighUS);
            countsPeakRegionLS = oneDScaledLSHistos[k][i]->Integral(binMassLowLS, binMassHighLS);  //this gets an estimate for B from the LS scaled to the BG in the US, using a single sideband
        
            invariantMassInfo << endl << "____________________________________________________________________________________________________" << endl;
            invariantMassInfo << "PtBin (0: 0-1, 1: 1-2, 2: 2-3, 3: 3-5, 4: 5-10, 5: all): " << k << endl;
            invariantMassInfo << "CentBin (0: per, 1: mid, 2: cent, 3: all): " << i << endl;
            //cout << "LS Scale Factor: " << LSScaleFactor[i] << endl;
            invariantMassInfo << endl << "S+B (from integral in mass window): " << countsPeakRegionUS << endl;
        
            oneDSubtractedInvMassHistos[k][i]->Add(oneDUSHistos[k][i], oneDScaledLSHistos[k][i], 1, -1); // form subtracted spectrum here
        
            oneDSubtractedInvMassHistosFunction[k][i] = (TH1D*)oneDSubtractedInvMassHistos[k][i]->Clone();
            tmp = oneDSubtractedInvMassStringsFunction[k] + binLabelCentralityClass[i];
            oneDSubtractedInvMassHistosFunction[k][i]->SetTitle(tmp);
        
            oneDSubtractedInvMassHistos[k][i]->Fit(g1, "qR0");
            oneDSubtractedInvMassHistos[k][i]->Fit(e1, "qR0+");
       
            g1->GetParameters(&par[0]);
            e1->GetParameters(&par[3]);
        
            fun->SetParameters(par);
            oneDSubtractedInvMassHistos[k][i]->Fit(fun, "qR+");
            oneDSubtractedInvMassHistos[k][i]->SetMarkerColor(2);
        
            BGFunc->SetParameter(0, fun->GetParameter(3));
            BGFunc->SetParameter(1, fun->GetParameter(4));
        
            //------------------------This is where things are controversial -- how do I fit the residual background????-------------------
        
        
            TF1 *fl = new TF1("fl",fline,1.6,2.1,2);
            fl->SetParameters(2,-1);
            //fit only the linear background excluding the signal area
            reject = kTRUE;
            oneDSubtractedInvMassHistos[k][i]->Fit(fl,"qR");
            reject = kFALSE;
        
            oneDSubtractedInvMassHistosFunction[k][i]->Add(fl, -1);
        
        
            //-----------------------------------------------------------------------------------------------------------------------------
        
            oneDSubtractedInvMassHistosFunction[k][i]->Fit(finalGaussian, "qR");
        
            integralFit = oneDSubtractedInvMassHistosFunction[k][i]->Integral(binMassLowUS, binMassHighUS);
        
            //integralFit is the estimate for "S"
            //countsPeakRegionUS is S+B
            //B = countsPeakRegionUS - integralFit
        
            BOverSPlusBFit[k][i] = (countsPeakRegionUS - integralFit)/countsPeakRegionUS;    //  B/S+B scale factor from background subtraction from fit
            SPlusBOverSFit[k][i] = countsPeakRegionUS/integralFit;                           //  S+B/S used as a final scale factor on the final correlation quantity -- from fit
            signalOverBG = integralFit/(countsPeakRegionUS-integralFit);
            optimizationFactor = ((signalOverBG*signalOverBG)*((countsPeakRegionUS-integralFit)*(countsPeakRegionUS-integralFit)))/((countsPeakRegionUS-integralFit)*(signalOverBG+2));
        
            invariantMassInfo << "S (from fit estimate): " << integralFit << endl;
            invariantMassInfo << "B (from fit estimate): " << (countsPeakRegionUS-integralFit) << endl;
            invariantMassInfo << "S+B/S from LS subtraction + expo Fit:  " << SPlusBOverSFit[k][i] << endl;
            invariantMassInfo << "B/S+B from LS subtraction + expo Fit: " << BOverSPlusBFit[k][i] << endl;
        
            invariantMassInfo << endl << "Signal over background (S/B): " << signalOverBG << endl;
            invariantMassInfo << "Optimization Factor: " << optimizationFactor << endl << endl;
        
            //invariantMassInfo << endl << "Left SB: " << integralLeftSB << endl;
            //invariantMassInfo <<  "Right SB: " << integralRightSB << endl;
            //invariantMassInfo << endl << "SB Average (how close is it to B?): " << SBAverage << endl << endl;
        
            /*  oneDUSHistos[k][i]->Draw();
            oneDScaledLSHistos[k][i]->Draw("SAME");
            tmp = oneDUSStrings[k] + binLabelCentralityClass[i] + withFitLabel + fileType;
            c->SaveAs(tmp);
            oneDLSHistos[k][i]->Draw();
            tmp = oneDLSStrings[k] + binLabelCentralityClass[i] + fileType;
            c->SaveAs(tmp);
            oneDSubtractedInvMassHistos[k][i]->Draw();
            tmp = oneDSubtractedInvMassStrings[k] + binLabelCentralityClass[i] + fileType;
            c->SaveAs(tmp); */
			
			sqrtSPlusB = TMath::Sqrt(countsPeakRegionUS);
			DMesonSig = integralFit/sqrtSPlusB;
            sigPerBin[k][i] = DMesonSig;
			tmp = "#frac{S}{#sqrt{S+B}} = ";
			str1 = Form("%.2f",DMesonSig);
			str2 = tmp + str1;
			
			invariantMassInfo << "D0 Significance:   " << DMesonSig << endl << endl;
			
			dMesonSigArray[i] = DMesonSig;
            
           /*  pad_number = (3*i)+1;
            invMassCanvas->cd(pad_number);
            oneDUSHistos[k][i]->Draw();
            
            pad_number = (3*i)+2;
            invMassCanvas->cd(pad_number);
            oneDSubtractedInvMassHistos[k][i]->Draw(); */
            
            //pad_number = (3*i)+3;
            //invMassCanvas->cd(pad_number);
            //oneDSubtractedInvMassHistosFunction[k][i]->Draw();
            //finalGaussian->Draw("SAME");
			//textBoxSignifigance->Draw("SAME");
            /*tmp = oneDSubtractedInvMassStringsFunction[k] + binLabelCentralityClass[i] + fileType;
            c->SaveAs(tmp);*/
        }
        
       
    }
    
    
    TCanvas * invMassCanvas = new TCanvas("invMassCanvas", "invMassCanvas", 4800, 3600);
    invMassCanvas->Divide(3,3);

    double pad_number;
    
    for(int i = 0; i < 3; i++){
    
        tmp = "#frac{S}{#sqrt{S+B}} = ";
        str1 = Form("%.2f",dMesonSigArray[i]);
		str2 = tmp + str1;
    
        TPaveText *textBoxSignifigance = new TPaveText(0.15, 0.7, .35, .85, "NB NDC");
        textBoxSignifigance->SetFillColor(0);
        textBoxSignifigance->AddText(str2);
        textBoxSignifigance->GetLine(0)->SetTextSize(.041);
    
        pad_number = (3*i)+1;
        invMassCanvas->cd(pad_number);
        oneDUSHistos[2][i]->Draw();
        
        pad_number = (3*i)+2;
        invMassCanvas->cd(pad_number);
        oneDSubtractedInvMassHistos[2][i]->Draw();
        
        pad_number = (3*i)+3;
        invMassCanvas->cd(pad_number);
        oneDSubtractedInvMassHistosFunction[2][i]->Draw();
        textBoxSignifigance->Draw("SAME");
    }
    
    invMassCanvas->SaveAs("invariant_mass_histograms.png");
    
    invariantMassInfo.close();
    
    cout << endl << "Invariant Mass histograms processed and S & B calculated for each Pt/Centrality bin" << endl << endl;
    
    //Print EPS file for Lanny DOE review
    
    TCanvas *c = new TCanvas("c2", "Histograms", 1100, 850);
   
//------------------------------------------
//BEGIN CORRELATION HISTOGRAM INFORMATION
//------------------------------------------

    /**********************************************************************************
        EXTRACT SIBLING AND MIXED HISTOGRAMS FROM RAW DATA FILE HERE
    **********************************************************************************/
    for(int band = 0; band < 4; band++){ //begin loop to get raw sibling and mixed histograms in both US and LS
        for(int i = 0; i < NUM_CENT_BINS; i++){
            for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
    
                //Sibling histogram information
                str1 = SibCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                sibCorrBinVzInt[band][k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            
                str1 = SibCorrLabels[band] + phiProj + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                sibCorrBinPhiProjVzInt[band][k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                formatCorrHist(sibCorrBinVzInt[band][k][i]);
            
                //mixed histogram information
                str1 = MixCorrLabels[band]  + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                mixCorrBinVzInt[band][k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            
                str1 = MixCorrLabels[band]  + phiProj + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                mixCorrBinPhiProjVzInt[band][k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                formatCorrHist(mixCorrBinPhiProjVzInt[band][k][i]);
            
                for(int j = 0; j < 10; j++){
                
                    str1 = SibCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
                    str2 = MixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
                    
                    sibCorrBin[band][k][j][i] = (TH2D*)file->Get(str1);
                    mixCorrBin[band][k][j][i] = (TH2D*)file->Get(str2);
                    
                    //sibCorrBin[band][k][j][i]->Sumw2();
                    //mixCorrBin[band][k][j][i]->Sumw2();
                    
                    //This gets the stat errors for the raw sibling and mixed histograms
                    for(int etaBin = 0; etaBin < NUM_ETA_BINS; etaBin++){
                        for(int phiBin = 0; phiBin < NUM_PHI_BINS; phiBin++){
                            
                            sibCorrBinStatErrors[band][k][j][i][etaBin][phiBin] = sibCorrBin[band][k][j][i]->GetBinContent(etaBin, phiBin);    
                            mixCorrBinStatErrors[band][k][j][i][etaBin][phiBin] = mixCorrBin[band][k][j][i]->GetBinContent(etaBin, phiBin);
                        }
                    }//Stat errors
                    
                    sibCorrBinVzInt[band][k][i]->Add(sibCorrBin[band][k][j][i]);
                    formatCorrHist(sibCorrBin[band][k][j][i]);  //histogram formatting for axes
                    mixCorrBinVzInt[band][k][i]->Add(mixCorrBin[band][k][j][i]);
                    formatCorrHist(mixCorrBin[band][k][j][i]);  //histogram formatting for axes
                    
                    
                            
                    
                    //sib
                    
                    str1 = phiProj + SibCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];         //sibling phi projection US
                    sibCorrBinPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);  
                    sibCorrBinPhiProj[band][k][j][i] = (TH1D*)sibCorrBin[band][k][j][i]->ProjectionY();  
                    sibCorrBinPhiProjVzInt[band][k][i]->Add(sibCorrBinPhiProj[band][k][j][i]);
                             
                    //mix
                    
                    str1 = phiProj + MixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];     //mixling phi projection US
                    mixCorrBinPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);  
                    mixCorrBinPhiProj[band][k][j][i] = (TH1D*)mixCorrBin[band][k][j][i]->ProjectionY();  
                    mixCorrBinPhiProjVzInt[band][k][i]->Add(mixCorrBinPhiProj[band][k][j][i]);
                             
                }
            }
        }
    }//end loop to get raw sibling and mixed histograms in both US and LS
    
  
    
  

    cout << endl << "Sibling and mixed histograms extracted" << endl;
   
    /**********************************************************************************
        NORMALIZE THE MIXED HISTOGRAMS HERE
    **********************************************************************************/
   
    for(int band = 0; band < 4; band++){ // begin loop to make scaled mixed histograms
        for(int i = 0; i < NUM_CENT_BINS; i++){ 
            for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
            
                str1 = ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                scaledMixCorrBinVzInt[band][k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                
                str1 = phiProj + ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                scaledMixCorrBinPhiProjVzInt[band][k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                      
                for(int j = 0; j < 10; j++){
            
                    scaledMixCorrBin[band][k][j][i] = (TH2D*) mixCorrBin[band][k][j][i]->Clone();
                    
                    //US histograms
                    numSib = sibCorrBin[band][k][j][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
                    numMixUS = mixCorrBin[band][k][j][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
                    integralRawMixHistos[band][k][j][i] = numMixUS;
                    scaledMixCorrBin[band][k][j][i]->Scale(numSib/numMixUS);
                    //VzScaleFactorUS[j] = numMixUS;
                    
                    for(int etaBin = 0; etaBin < NUM_ETA_BINS; etaBin++){
                        for(int phiBin = 0; phiBin < NUM_PHI_BINS; phiBin++){
                            
                            scaledMixCorrBinStatErrors[band][k][j][i][etaBin][phiBin] = mixCorrBinStatErrors[band][k][j][i][etaBin][phiBin]/((numSib/numMixUS)*(numSib/numMixUS));    
                            
                        }
                    }//Stat errors
             
             
                    //scaledMixCorrBinVzInt[band][k][i]->Add(
             
                    //US scaled histograms
                    
                    str1 = ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
            
                    formatCorrHist(scaledMixCorrBin[band][k][j][i], str1);  //formatting
            
                    str1 = phiProj + ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];       //scaled mix phi projection US 
                    scaledMixCorrBinPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);  
                    scaledMixCorrBinPhiProj[band][k][j][i] = (TH1D*)scaledMixCorrBin[band][k][j][i]->ProjectionY();  
                }
            }
        }  
    }//end loop to make scaled mixed histograms  */
   //--------------------------------------------------------calculate some useful quantities here--------------------------------------
   
   cout << endl << "Mixed Histograms Normalized." << endl << endl;
   
    for(int band = 0; band < 4; band++){                               //build arrays to store information for adding the various alpha's (Vz, centrality, etc.)
       for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
           
           integralRawMixHistosVzIntCentInt[band][k][0] = 0;
           integralRawMixHistosVzIntCentInt[band][k][1] = 0;
           integralRawMixHistosVzIntCentInt[band][k][2] = 0;
   
            for(int i = 0; i < NUM_CENT_BINS; i++){
            
                integralRawMixHistosVzInt[band][k][i] = 0;
                //cout << band << "   " << k << "   "<< i << endl;
               
                for(int j = 0; j < 10; j++){
            
                    integralRawMixHistosVzInt[band][k][i] = integralRawMixHistosVzInt[band][k][i] + integralRawMixHistos[band][k][j][i];
                    
                }   
            }
        }
    }   
 //------------------------------------------------------------------------------------------------------------------------------------
 
   cout << "Weight Factors Calculated." << endl << endl;
    /**********************************************************************************
        CALCULATE DELTA RHO HERE
    **********************************************************************************/
 
    for(int band = 0; band < 4; band++){ // begin loop to make sib minus scaledMixed histos (delRho)
        for(int i = 0; i < NUM_CENT_BINS; i++){ 
            for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
            
                str1 = SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                sibMinusScaledMixVzInt[band][k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                str1 = phiProj + SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                sibMinusScaledMixPhiProjVzInt[band][k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
           
                for(int j = 0; j < 10; j++){
        
                    str1 = SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
                    sibMinusScaledMix[band][k][j][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                    sibMinusScaledMix[band][k][j][i]->Add(sibCorrBin[band][k][j][i], scaledMixCorrBin[band][k][j][i], 1, -1);
                    
                    formatCorrHist(sibMinusScaledMix[band][k][j][i]); //formatting
        
                    str1 = phiProj + SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
                    sibMinusScaledMixPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);  
                    sibMinusScaledMixPhiProj[band][k][j][i] = (TH1D*)sibMinusScaledMix[band][k][j][i]->ProjectionY();  
                
                    sibMinusScaledMixVzInt[band][k][i]->Add(sibMinusScaledMix[band][k][j][i]);
                    sibMinusScaledMixPhiProjVzInt[band][k][i]->Add(sibMinusScaledMixPhiProj[band][k][j][i]);
                }
            }
        }// end loop to make sibling minus scaled mixed histograms (delRho)
    }
  
   cout << "deltaRho for each band (US, LS, SBL, SBR) calculated." << endl << endl;
  
   //-----------------------------------------------------------------------------------------------------------------------------------------------------------
   //-----------------------------------------------------------------------------------------------------------------------------------------------------------
   //-----------------------------------------BEGIN SECTION TO CALCULATE ACTUAL CORRELATIONS--------------------------------------------------------------------
   //-----------------------------------------------------------------------------------------------------------------------------------------------------------
   //-----------------------------------------------------------------------------------------------------------------------------------------------------------
   
    TFile * statErrorFile = new TFile("output_errors_correct_2_10GeV.root");
    TH2D * statErrorHistTemp; 
    output->cd();
    bool errorFile = false;
    if(!statErrorFile) { cout << "No error file found!!!!!!" << endl; errorFile = false;}
    //Example: ErrorHistogram_SBR_PtBin_0_VzBin_0_CentBin_0
    //Example: ErrorHistogram_US_Signal_Region_PtBin_0_VzBin_0_CentBin_0
    TString errorHistLabel = "ErrorHistogram_";
    TString errorBand[3] = {"SBL_", "US_Signal_Region_", "SBR_"};
    TString errorPtBin[5] = {"PtBin_0", "PtBin_1", "PtBin_2", "PtBin_3", "PtBin_4"};
    TString errorVzBin = "_VzBin_";
    //TString centBinLabel = "_CentBin_";
    //TString binLabelVz[10] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
    //TString binLabelCent[16] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"};
    
   
    for(int band = 0; band < 4; band++){     // begin loop to make delRho/RhoRef US histos -- this will calculate the delRho/Rho histograms for all 4 bands -- LS, US, SBR, SBL
        for(int i = 0; i < NUM_CENT_BINS; i++){ 
            for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
            
                str1 = delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                delRhoOverRhoRefBinVzInt[band][k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                str1 = phiProj + delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                delRhoOverRhoRefBinPhiProjVzInt[band][k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            
                for(int j = 0; j < 10; j++){

                    str1 = delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
                    delRhoOverRhoRefBin[band][k][j][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            
                    str1 = errorsFromRootDelRhoOverRhoRef[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
                    rootErrorsOnDelRhoOverRhoRef[band][k][j][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);;
            
                    if(scaledMixCorrBin[band][k][j][i]->GetEntries() > 0){
                    
                        delRhoOverRhoRefBin[band][k][j][i]->Divide(sibMinusScaledMix[band][k][j][i], scaledMixCorrBin[band][k][j][i], 1, 1);
                        delRhoOverRhoRefBin[band][k][j][i]->Write();
                        for(int eta = 1; eta < NUM_ETA_BINS+1; eta++){
                            for(int phi = 1; phi < NUM_PHI_BINS+1; phi++){
                
                                rootErrorsOnDelRhoOverRhoRef[band][k][j][i]->SetBinContent(eta, phi, delRhoOverRhoRefBin[band][k][j][i]->GetBinError(eta, phi));
                            }
                        }
                        
                        rootErrorsOnDelRhoOverRhoRef[band][k][j][i]->Write();
                    }
                    
                    //STORE REAL STAT ERROR VALUES HERE!!!!!!!!
                    if(band != 3 && k > 1 && errorFile){
                     
                        str2 = errorHistLabel + errorBand[band] + errorPtBin[k] + errorVzBin + binLabelVz[j] + centBinLabel + binLabelCent[i];
                        cout << str2 << endl;
                        statErrorHistTemp = (TH2D*) statErrorFile->Get(str2);
                        
                        for(int etabin = 1; etabin < NUM_ETA_BINS+1; etabin++){
                            for(int phibin = 1; phibin < NUM_PHI_BINS+1; phibin++){
                            
                                //LEFT AND RIGHT SIDEBAND ARE FLIPPED IN THE ACTUAL ERROR FILE!!!!!
                                //FLIP THEM BACK RIGHT HERE!!!!
                                delRhoOverRhoRefBin[band][k][j][i]->SetBinError(etabin, phibin, statErrorHistTemp->GetBinContent(etabin,phibin));
                            }
                        }    
                    }
                    
                    //ERRORS INPUT AND DONE!!!!
                    
                    formatCorrHist(delRhoOverRhoRefBin[band][k][j][i]);
                    //VzScaleFactorUS[j] = delRhoOverRhoRefUSBin[k][j][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
                    delRhoOverRhoRefBin[band][k][j][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                
                    str1 = phiProj + delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];        //sib - scaled mix phi projection US 
                    delRhoOverRhoRefBinPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);  
                    delRhoOverRhoRefBinPhiProj[band][k][j][i] = (TH1D*)delRhoOverRhoRefBin[band][k][j][i]->ProjectionY();  
                
                    scaleFactorVz = integralRawMixHistos[band][k][j][i]/integralRawMixHistosVzInt[band][k][i];
                    delRhoOverRhoRefBinVzInt[band][k][i]->Add(delRhoOverRhoRefBin[band][k][j][i], scaleFactorVz);
                
                }
 
                
                delRhoOverRhoRefBinVzInt[band][k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                
                delRhoOverRhoRefBinVzInt[band][k][i]->Write();
                
                delRhoOverRhoRefBinPhiProjVzInt[band][k][i] = (TH1D*)delRhoOverRhoRefBinVzInt[band][k][i]->ProjectionY();
            
            }
        } 
    }//end loop to make delRho/Rho_ref histograms*/
    
    statErrorFile->Close();
    
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------THIS SECTION DOES THE FINAL SUBTRACTION IN THE 16 MULTIPLICITY BINS------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------
    
        for(int i = 0; i < NUM_CENT_BINS; i++){// begin side band subtraction here -- this will take the US delRho/Rho and subtract the SB average delRho/Rho from it
            for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
           
                str1 = sideBandAverageLabel + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                sideBandAverage[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                str1 = phiProj + sideBandAverageLabel + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                sideBandAveragePhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
               
                
                sideBandAverage[k][i]->Add(delRhoOverRhoRefBinVzInt[0][k][i], delRhoOverRhoRefBinVzInt[2][k][i], SBLScaleFactor, SBRScaleFactor);
                sideBandAveragePhiProj[k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[0][k][i], delRhoOverRhoRefBinPhiProjVzInt[2][k][i], SBLScaleFactor, SBRScaleFactor);
                
                str1 = subtractedCorrLabels[0] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                fullySubtractedCorrSideBandCent[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                str1 = phiProj + subtractedCorrLabels[0] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                fullySubtractedCorrSideBandCentPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
           
                fullySubtractedCorrSideBandCent[k][i]->Add(delRhoOverRhoRefBinVzInt[1][k][i], sideBandAverage[k][i], 1, -BOverSPlusBFit[k][3]);
               
                sideBandAverage[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                
                fullySubtractedCorrSideBandCent[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                fullySubtractedCorrSideBandCent[k][i]->Write();
                
                fullySubtractedCorrSideBandCentPhiProj[k][i] = (TH1D*)fullySubtractedCorrSideBandCent[k][i]->ProjectionY();
            }    
        }// end side band subtraction here
        
        
        /* for(int i = 0; i < NUM_CENT_BINS; i++){// begin LS subtraction here  -- this will take the US delRho/Rho and subtract the LS delRho/Rho from it
            for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
           
                str1 = subtractedCorrLabels[1] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                fullySubtractedCorrLSCent[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                str1 = phiProj + subtractedCorrLabels[1] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                fullySubtractedCorrLSCentPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
           
                fullySubtractedCorrLSCent[k][i]->Add(delRhoOverRhoRefBinVzInt[1][k][i], delRhoOverRhoRefBinVzInt[3][k][i], 1, -BOverSPlusBFit[k][3]);
                
                //fullySubtractedCorrLSCentPhiProj[k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[1][k][i], delRhoOverRhoRefBinPhiProjVzInt[3][k][i], 1, -BOverSPlusB[k]);
            
                str1 = path + outputFolders[5] + outputFolders[6] + PtSubFolders[k] + subtractedCorrLabels[1] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i] + fileType; 
                fullySubtractedCorrLSCent[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                fullySubtractedCorrLSCent[k][i]->Draw("SURF1");
                c->SaveAs(str1);
                
                str1 = path + outputFolders[5] + outputFolders[6] + PtSubFolders[k] + phiProj + subtractedCorrLabels[1] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i] + fileType; 
                fullySubtractedCorrLSCentPhiProj[k][i] = (TH1D*) fullySubtractedCorrLSCent[k][i]->ProjectionY();
                fullySubtractedCorrLSCentPhiProj[k][i]->Draw();
                c->SaveAs(str1);
           
           
            }    
        }// end LS subtraction here */
          

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------      
    //---------------------------THIS SECTION DOES THE FINAL SUBTRACTION IN THE 3 centrality BINS----------------------------------------------------------------      
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------

    /* if(nTracks >= 1   && nTracks < 42)  { return 0;  }
    if(nTracks >= 42  && nTracks < 86)  { return 1;  }
    if(nTracks >= 86  && nTracks < 131)  { return 2;  }
    if(nTracks >= 131  && nTracks < 183) { return 3;  }
    if(nTracks >= 183 && nTracks < 235) { return 4;  }
    if(nTracks >= 235 && nTracks < 288) { return 5;  }
    if(nTracks >= 288 && nTracks < 340) { return 6;  }
    if(nTracks >= 340 && nTracks < 392) { return 7;  }
    if(nTracks >= 392 && nTracks < 440) { return 8;  }
    if(nTracks >= 440 && nTracks < 491) { return 9;  }
    if(nTracks >= 491 && nTracks < 542) { return 10;  }
    if(nTracks >= 542 && nTracks < 593) { return 11;  }
    if(nTracks >= 593 && nTracks < 644) { return 12;  }
    if(nTracks >= 644 && nTracks < 695) { return 13;  }
    if(nTracks >= 695 && nTracks < 746) { return 14;  }
    if(nTracks >= 746)                  { return 15; } */
    
    
    TString bandSubtract[3]          = {"LS_","SideBand_", "US"};
    int temp = 0;
    int centralityCombinations[6] = {0, 3, 4, 7, 8, 15}; //This was calculated for combining multiplicity bins into the proper centrality bins
    //int centralityCombinations[6] = {0, 2, 3, 15}; //This was calculated for combining multiplicity bins into the proper centrality bins
    int numBinsInCentrality = 0;
   
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
        for(int i = 0; i < 3; i++){
           
            str1 =  fullySubtractedLabelCustomCentrality + bandSubtract[2] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            fullySubtractedCorrCustomCentralityUS[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
           
            //str1 =  fullySubtractedLabelCustomCentrality + bandSubtract[0] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            //fullySubtractedCorrCustomCentralityLS[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            str1 =  fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            fullySubtractedCorrCustomCentralitySideBand[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            //str1 =  phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[0] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            // fullySubtractedCorrCustomCentralityLSPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            str1 =  phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
           
		    str1 =  SibCorrLabels[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            sibCorrCentrality[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            str1 =  MixCorrLabels[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            mixCorrCentrality[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            str1 =  SibCorrLabels[1] + bandLabel[0] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            SBLeftCorrCentrality[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            str1 =  SibCorrLabels[1] + bandLabel[2] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            SBRightCorrCentrality[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
           
            for(int band = 0; band < 4; band++){
           
                delRhoOverRhoRefBinVzIntCentInt[band][k][i] = new TH2D("","",  NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                delRhoOverRhoRefBinVzIntCentIntPhiProj[band][k][i] = new TH1D("","", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            }
           
            numBinsInCentrality = centralityCombinations[(2*i)+1] - centralityCombinations[(2*i)];
           
            for(int band = 0; band < 4; band++){
                for(int bin = 0; bin < numBinsInCentrality+1; bin++){
               
                    temp = centralityCombinations[(2*i)]+bin;
                    integralRawMixHistosVzIntCentInt[band][k][i] = integralRawMixHistosVzIntCentInt[band][k][i] + integralRawMixHistosVzInt[band][k][temp];
           
                }
            }
        }    
       
       cout << endl << endl;
        
        for(int i = 0; i < 3; i++){        
           
            numBinsInCentrality = centralityCombinations[(2*i)+1] - centralityCombinations[(2*i)];
           
            for(int bin = 0; bin < numBinsInCentrality+1; bin++){
           
                temp = centralityCombinations[(2*i)]+bin;
                
                scaleFactorVzCent = integralRawMixHistosVzInt[1][k][temp]/integralRawMixHistosVzIntCentInt[1][k][i];                 //US
                delRhoOverRhoRefBinVzIntCentInt[1][k][i]->Add(delRhoOverRhoRefBinVzInt[1][k][temp], scaleFactorVzCent);
                delRhoOverRhoRefBinVzIntCentIntPhiProj[1][k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[1][k][temp], scaleFactorVzCent);
                
                cout << "US cent Scale Factor (cent bin: " << i << ") :" << scaleFactorVzCent << endl;
               
                
                scaleFactorVzCent = integralRawMixHistosVzInt[2][k][temp]/integralRawMixHistosVzIntCentInt[2][k][i];                 //Right SB
                delRhoOverRhoRefBinVzIntCentInt[2][k][i]->Add(delRhoOverRhoRefBinVzInt[2][k][temp], scaleFactorVzCent);
                delRhoOverRhoRefBinVzIntCentIntPhiProj[2][k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[2][k][temp], scaleFactorVzCent);
                
                scaleFactorVzCent = integralRawMixHistosVzInt[0][k][temp]/integralRawMixHistosVzIntCentInt[0][k][i];                 //Left SB
                delRhoOverRhoRefBinVzIntCentInt[0][k][i]->Add(delRhoOverRhoRefBinVzInt[0][k][temp], scaleFactorVzCent);
                delRhoOverRhoRefBinVzIntCentIntPhiProj[0][k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[0][k][temp], scaleFactorVzCent);
            
                //summing the sib and mixed histograms
                sibCorrCentrality[k][i]->Add(sibCorrBinVzInt[1][k][temp]);
                mixCorrCentrality[k][i]->Add(mixCorrBinVzInt[1][k][temp]);
                
                SBLeftCorrCentrality[k][i]->Add(sibCorrBinVzInt[0][k][temp]);
                SBRightCorrCentrality[k][i]->Add(sibCorrBinVzInt[2][k][temp]);
            
            }
        }    
    }   
      
      
    //TESTING
    
    //BOverSPlusBFit[k][i] = (countsPeakRegionUS - integralFit)/countsPeakRegionUS;    //  B/S+B scale factor from background subtraction from fit
    //SPlusBOverSFit[k][i] = countsPeakRegionUS/integralFit;  

    TH1D* testProjUS = new TH1D("DEEPA_SIBLING_TEST_US", "DEEPA_SIBLING_TEST_US", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH1D* testProjSBL = new TH1D("DEEPA_SIBLING_TEST_SBL", "DEEPA_SIBLING_TEST_SBL", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH1D* testProjSBR = new TH1D("DEEPA_SIBLING_TEST_SBR", "DEEPA_SIBLING_TEST_SBR", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH1D* testProjSBAvg = new TH1D("DEEPA_SIBLING_TEST_SB_AVG", "DEEPA_SIBLING_TEST_SB_AVG", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH1D* testProjFinal = new TH1D("DEEPA_SIBLING_TEST_FINAL", "DEEPA_SIBLING_TEST_FINAL", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    testProjUS = (TH1D*) sibCorrCentrality[2][0]->ProjectionY("DEEPA_SIBLING_TEST_US", 2, 12, "");
    testProjSBL = (TH1D*) SBLeftCorrCentrality[2][0]->ProjectionY("DEEPA_SIBLING_TEST_SBL", 2, 12, "");
    testProjSBR = (TH1D*) SBRightCorrCentrality[2][0]->ProjectionY("DEEPA_SIBLING_TEST_SBR", 2, 12, "");
    //testProj->SetNameTitle("SUPERTEST_DEEPA_COMPARE_PYTHIA", "SUPERTEST_DEEPA_COMPARE_PYTHIA");
    double testScaleNumber = 1.0/2768;
    testProjSBAvg->Add(testProjSBL, testProjSBR, .5, .5);
    testProjFinal->Add(testProjUS, testProjSBAvg, 1, -BOverSPlusBFit[2][0]);
    testProjFinal->Scale(testScaleNumber);
    
    testProjUS->Write();
    testProjSBL->Write();
    testProjSBR->Write();
    testProjFinal->Write();
    
    //testProj->Scale(testScaleNumber); 
    //testProj->Write(); 
      
    cout << "WORKS GOOD" << endl;  
      
       double phiProjectionScaleFactor = 1.0/9.0;
       
        for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){    
            for(int i = 0; i < 3; i++){
      
                cout << "k, i: " << k << " , " << i << endl;
               
                //fullySubtractedCorrCustomCentralityLS[k][i]->Add(delRhoOverRhoRefBinVzIntCentInt[1][k][i], delRhoOverRhoRefBinVzIntCentInt[3][k][i], 1, -BOverSPlusBFit[k][i]); //LS
                //fullySubtractedCorrCustomCentralityLS[k][i]->Scale(SPlusBOverSFit[k][i]); //LS
                
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Add(delRhoOverRhoRefBinVzIntCentInt[1][k][i]);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Add(delRhoOverRhoRefBinVzIntCentInt[0][k][i], -SBLScaleFactor*BOverSPlusBFit[k][i]); //SBLScaleFactor*
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Add(delRhoOverRhoRefBinVzIntCentInt[2][k][i], -SBRScaleFactor*BOverSPlusBFit[k][i]);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Scale(SPlusBOverSFit[k][i]);
                
                //fullySubtractedCorrCustomCentralityLSPhiProj[k][i] = (TH1D*) fullySubtractedCorrCustomCentralityLS[k][i]->ProjectionY(); //LS
                //fullySubtractedCorrCustomCentralityLSPhiProj[k][i]->Scale(phiProjectionScaleFactor);                                     //LS
                
                fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i] = (TH1D*) fullySubtractedCorrCustomCentralitySideBand[k][i]->ProjectionY(fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetTitle(),3,11, "");
                fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Scale(phiProjectionScaleFactor);
       
                //fullySubtractedCorrCustomCentralityLS[k][i]->Write(); //LS
           
                fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Write();
                fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Write();
           
                str1 = fullySubtractedCorrCustomCentralitySideBand[k][i]->GetTitle() + fileType;
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Draw("SURF1");
                c->SaveAs(str1);
                //fullySubtractedCorrCustomCentralityLSPhiProj[k][i]->Draw(); //LS
           
             }    
        }
 
 
    TCanvas * finalCorrCanvasNoDStar = new TCanvas("", "", 5100, 1200);
    finalCorrCanvasNoDStar->Divide(3,1);
    int pad = 1;
    
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){    //plot projections -- phi, SS eta, AS eta
        for(int i = 0; i < 3; i++){
        
            finalCorrCanvasNoDStar->cd(pad);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->SetStats(0);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetTitle("#Delta#eta");
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->CenterTitle();
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetTitle("#Delta#phi");
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->CenterTitle();
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetTitle("#Delta#rho/#rho_{mix}");
            
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetTitleSize(.1);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetTitleSize(.1); 
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetTitleSize(.07);
    
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetTitleOffset(1.15);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetTitleOffset(1.15); 
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetTitleOffset(1.0);
			fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetLabelSize(.06);
			fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetLabelSize(.06);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetLabelSize(.05);
			fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetNdivisions(5,3,0, 1);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->Draw("SURF1");
            pad++;
        }
    }
    
    finalCorrCanvasNoDStar->SaveAs("final_correlations_on_one_plot_No_DStar.png");
    
    
    
 
 
 
    cout << endl <<"Starting D* correction now...." << endl;
    
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
 
    //D* Correction
    //Algorithm: Read-in the D* (delEta, delPhi) distributions
    //           Take same-event minus mix-event to get yield (normalization factors calculated separately, hardcoded below)
    //           Take result, divide it by the (kpi) mixed event, normalized by #pairs, and multiplied by S/S+B
 
    TString DStarLabel = "KPiPi_minus_KPi_";
    TString siblingDStarLabel      = "Sibling_";
    TString mixDStarLabel          = "Mixed_";
    TString diffDStarLabel         = "Diff_";
    TString correctionLabel        = "Full_DStar_Correction_";
    TString angularLabel = "_Angular_Distribution";
    TString DStarLabelAng = "DStarHistograms_";
   
    TH2D* dStarAngularDistSib[5][3];
    TH2D* dStarAngularDistMix[5][3];
    TH2D* dStarAngularDistDiff[5][3];
    TH2D* dStarCorrection[5][3];
    
    TH1D* finalCorrDStarCorrectedPhiProj[5][3];
    TH1D* finalCorrDStarCorrectedSSEtaProj[5][3];
    TH1D* finalCorrDStarCorrectedASEtaProj[5][3];
    
    double scaleFactorOnMix;
    
    for(int k = 2; k < 5; k++){
        for(int j = 0; j < 3; j++){
    
            str1 = siblingDStarLabel + angularLabel + DStarLabelAng + PtBinLabel[0] + binLabelPt[k] + binLabelCentralityClass[j];
            dStarAngularDistSib[k][j] = (TH2D*) file->Get(str1);
            //dStarAngularDistSib[k][j]->Sumw2();
            str1 = mixDStarLabel + angularLabel + DStarLabelAng + PtBinLabel[0] + binLabelPt[k] + binLabelCentralityClass[j];
            dStarAngularDistMix[k][j] = (TH2D*) file->Get(str1);
            //dStarAngularDistMix[k][j]->Sumw2();
            str1 = diffDStarLabel + angularLabel + DStarLabelAng + PtBinLabel[0] + binLabelPt[k] + binLabelCentralityClass[j];
            dStarAngularDistDiff[k][j] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, (3*TMath::PiOver2())+phiBinShift);
            
            str1 = correctionLabel + angularLabel + DStarLabelAng + PtBinLabel[0] + binLabelPt[k] + binLabelCentralityClass[j];
            dStarCorrection[k][j] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, (3*TMath::PiOver2())+phiBinShift);
    
        }
    }
 
    double dStarNormFactor[5][3]; // = {{.2, .2, .2}, {.2, .2, .2}, {0.205836, 0.204568, 0.200935}, {0.216633, 0.206026, 0.20196}, {0.241584, 0.208599, 0.201621}};
    dStarNormFactor[2][0] = .212336; //0.205836;
    dStarNormFactor[2][1] = .203867; //0.204568;
    dStarNormFactor[2][2] = .20176; //0.200935;
    dStarNormFactor[3][0] = 0.216633;
    dStarNormFactor[3][1] = 0.206026;
    dStarNormFactor[3][2] = 0.20196;
    dStarNormFactor[4][0] = 0.241584;
    dStarNormFactor[4][1] = 0.208599;
    dStarNormFactor[4][2] = 0.201621;
    double nEventsInBin[5][3]; // = {{1, 1, 1}, {1, 1, 1}, {2317, 22238, 48808}, {977, 7459, 11881}, {452, 3043, 3715}};
    //nEventsInBin[2][0] = 2317;
    //nEventsInBin[2][1] = 22238;
    //nEventsInBin[2][2] = 48808;
    nEventsInBin[3][0] = 977;
    nEventsInBin[3][1] = 7459;
    nEventsInBin[3][2] = 11881;
    nEventsInBin[4][0] = 452;
    nEventsInBin[4][1] = 3043;
    nEventsInBin[4][2] = 3715;
    
    nEventsInBin[2][0] = 3538;
    nEventsInBin[2][1] = 30993;
    nEventsInBin[2][2] = 61168;
    //nEventsInBin[3][0] = 977;
    //nEventsInBin[3][1] = 7459;
    //nEventsInBin[3][2] = 11881;
    //nEventsInBin[4][0] = 452;
    //nEventsInBin[4][1] = 3043;
    //nEventsInBin[4][2] = 3715;
	
	//double dStarInvMassNorm;
    
    double binValueSib;
    double binValueMix, binValueDStarMix;
    double binValueCorrection;
    
    double OneOverNEvent;
    double OneOverNEventTimesNMix;
    double integralDStarSib;
	double integralSib;
    double integralMix, integralDStarMix;
    double correctionScaleFactor, dStarPairNorm;
    
    double errorEstimateCenter;
    double errorEstimate1;
    double errorEstimate2;
    double errorOnMixZeroBin, errorOnDStarZeroBin;
    double binValueMix, binValueDStar, ratioDStarToMix;
    
    double ratioOfMixed;
    double intOfSibDStar, intOfMixDStar, ratioDStarSibToMix;
    
    cout << endl;
 
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){    
        for(int i = 0; i < 3; i++){
      
            integralDStarSib = 0;
            integralMix = 0;
			integralSib = 0;
            integralDStarMix = 0;
            
            integralMix = mixCorrCentrality[k][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
            integralSib = sibCorrCentrality[k][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
            
            OneOverNEvent = 1/nEventsInBin[k][i];
            OneOverNEventTimesNMix = OneOverNEvent*(1/5.0);
            //calculate the D0-softPion sibling dist here
            //dStarNormFactor[k][i] = .2*dStarNormFactor[k][i];
			dStarAngularDistSib[k][i]->Write();
            dStarAngularDistMix[k][i]->Write();
            
            str1 = sibLabel + angularLabel + DStarLabelAng + PtBinLabel[0] + binLabelPt[k] + binLabelCentralityClass[i] + fileType;
            dStarAngularDistSib[k][i]->Draw("SURF1");
            c->SaveAs(str1);
            str1 = mixLabel + angularLabel + DStarLabelAng + PtBinLabel[0] + binLabelPt[k] + binLabelCentralityClass[i] + fileType;
            dStarAngularDistMix[k][i]->Draw("SURF1");
            c->SaveAs(str1);
            //Get the total number of sibling D0-SoftPion pairs
            dStarAngularDistDiff[k][i]->Add(dStarAngularDistSib[k][i], dStarAngularDistMix[k][i], 1, -dStarNormFactor[k][i]);
            dStarAngularDistDiff[k][i]->Write();
            //
            
            str1 = diffLabel + angularLabel + DStarLabelAng + PtBinLabel[0] + binLabelPt[k] + binLabelCentralityClass[i] + fileType;
            dStarAngularDistDiff[k][i]->Draw("SURF1");
            c->SaveAs(str1);
            
            scaleFactorOnMix = (1/SPlusBOverSFit[k][i])*(integralSib/integralMix);    //this factor multiplies the TOTAL signal-region mized event to estimate the pair-count
			cout << "S/S+B - " << 1/SPlusBOverSFit[k][i] << endl;
            cout << "Int of Sib / Int of Mix = " << (integralSib/integralMix) << endl;
            cout << "mix Scale Factor: " << scaleFactorOnMix << endl;                 // that comes from real D0-hadron mixed pairs
            
            //mixCorrCentrality[k][i]->Write();
            mixCorrCentrality[k][i]->Scale(scaleFactorOnMix);  
            mixCorrCentrality[k][i]->Write();
           
            
            
            //dStarAngularDistMix[k][i]->Scale(ratioDStarSibToMix);   //scales the D* mixed-event to normalize the counts
            //----------------------------
           
                        //scales the signal-region mixed events dist.
            
            dStarCorrection[k][i]->Divide(dStarAngularDistDiff[k][i], dStarAngularDistMix[k][i], 1, 1); 
            dStarAngularDistMix[k][i]->Divide(mixCorrCentrality[k][i]);
            dStarCorrection[k][i]->Multiply(dStarAngularDistMix[k][i]);
            //dStarCorrection[k][i]->Scale(ratioOfMixed);
            
            //errorOnMixZeroBin = mixCorrCentrality[k][i]->GetBinError(7, 3);
            //cout << "Error on Mix Zero Bin: " << errorOnMixZeroBin << endl;
            //errorOnDStarZeroBin = dStarAngularDistDiff[k][i]->GetBinError(7, 3);
            //cout << "Error on D* Zero Bin: " << errorOnDStarZeroBin << endl;
            //errorEstimate1 = TMath::Sqrt((errorOnMixZeroBin*errorOnMixZeroBin) + (errorOnDStarZeroBin*errorOnDStarZeroBin));
            
            //cout << "error on (0,0) bin: " << errorEstimate1 << endl;
            
            dStarCorrection[k][i]->Write();
            str1 = correctionLabel + angularLabel + DStarLabelAng + PtBinLabel[0] + binLabelPt[k] + binLabelCentralityClass[i] + fileType;
            dStarCorrection[k][i]->Draw("SURF1");
            c->SaveAs(str1);
            errorEstimateCenter = 1.41*fullySubtractedCorrCustomCentralitySideBand[k][i]->GetBinError(7,4);
            errorEstimate1 = 1.41*fullySubtractedCorrCustomCentralitySideBand[k][i]->GetBinError(8,4);
            errorEstimate2 = 1.41*fullySubtractedCorrCustomCentralitySideBand[k][i]->GetBinError(7,5);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->Add(dStarCorrection[k][i], -1);
            
            
            str1 = correctionLabel + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            fullySubtractedCorrCustomCentralitySideBand[k][i]->SetNameTitle(str1, str1);
            //fullySubtractedCorrCustomCentralitySideBand[k][i]->SetBinError(7, 4, errorEstimateCenter);
            //fullySubtractedCorrCustomCentralitySideBand[k][i]->SetBinError(8, 4, errorEstimate1);
            //fullySubtractedCorrCustomCentralitySideBand[k][i]->SetBinError(6, 4, errorEstimate1);
            //fullySubtractedCorrCustomCentralitySideBand[k][i]->SetBinError(7, 5, errorEstimate2);
            //fullySubtractedCorrCustomCentralitySideBand[k][i]->SetBinError(7, 3, errorEstimate2);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->Write();
            
            str1 = correctionLabel + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fileType;
            fullySubtractedCorrCustomCentralitySideBand[k][i]->Draw("SURF1");
            c->SaveAs(str1);
            

            //TString phiProj = "Phi_projection_";
            //TString SSEtaProj = "Same_Side_Eta_Projection_";
            //TString ASEtaProj = "Away_Side_Eta_Projection_";
            
        }    
    }
 
    TCanvas * finalCorrCanvas = new TCanvas("", "", 5100, 1200);
    finalCorrCanvas->Divide(3,1);
    int pad = 1;
    
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){    //plot projections -- phi, SS eta, AS eta
        for(int i = 0; i < 3; i++){
        
            finalCorrCanvas->cd(pad);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->SetStats(0);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetTitle("#Delta#eta");
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->CenterTitle();
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetTitle("#Delta#phi");
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->CenterTitle();
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetTitle("#Delta#rho/#rho_{mix}");
            
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetTitleSize(.1);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetTitleSize(.1); 
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetTitleSize(.07);
    
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetTitleOffset(1.15);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetTitleOffset(1.15); 
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetTitleOffset(1.0);
			fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetLabelSize(.06);
			fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetLabelSize(.06);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetLabelSize(.05);
			fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetNdivisions(5,3,0, 1);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->Draw("SURF1");
            pad++;
        }
    }
    
    finalCorrCanvas->SaveAs("final_correlations_on_one_plot.png");
        
 
 
    TCanvas * projectionGridCanvas = new TCanvas("", "", 3600, 1100);
    projectionGridCanvas->Divide(3,1);
    int pad = 1;
    double projScale;
    TString corrProjLabel = "Fully_DStar_Corrected_Projections_";
 
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){    //plot projections -- phi, SS eta, AS eta
        for(int i = 0; i < 3; i++){
  
            finalCorrDStarCorrectedPhiProj[k][i] = (TH1D*)fullySubtractedCorrCustomCentralitySideBand[k][i]->ProjectionY(phiProj, 4, 10, "");
            str1 = phiProj + correctionLabel + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            finalCorrDStarCorrectedPhiProj[k][i]->SetTitle(str1);
            projScale = 1.0/7.0;
            finalCorrDStarCorrectedPhiProj[k][i]->Scale(projScale);
            
            finalCorrDStarCorrectedSSEtaProj[k][i] = (TH1D*)fullySubtractedCorrCustomCentralitySideBand[k][i]->ProjectionX(SSEtaProj, 1, 6, ""); 
            str1 = SSEtaProj + correctionLabel + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            finalCorrDStarCorrectedSSEtaProj[k][i]->SetTitle(str1);
            projScale = 1.0/6.0;
            finalCorrDStarCorrectedSSEtaProj[k][i]->Scale(projScale);
            finalCorrDStarCorrectedSSEtaProj[k][i]->SetMaximum(finalCorrDStarCorrectedPhiProj[k][i]->GetMaximum());
            
            finalCorrDStarCorrectedASEtaProj[k][i] = (TH1D*)fullySubtractedCorrCustomCentralitySideBand[k][i]->ProjectionX(ASEtaProj, 7, 12, ""); 
            str1 = ASEtaProj + correctionLabel + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            finalCorrDStarCorrectedASEtaProj[k][i]->SetTitle(str1);
            projScale = 1.0/6.0;
            finalCorrDStarCorrectedASEtaProj[k][i]->Scale(projScale);
            finalCorrDStarCorrectedASEtaProj[k][i]->SetMaximum(finalCorrDStarCorrectedPhiProj[k][i]->GetMaximum());
            
            
            projectionGridCanvas->cd(pad);
            finalCorrDStarCorrectedPhiProj[k][i]->Draw();
            pad++;
            projectionGridCanvas->cd(pad);
            finalCorrDStarCorrectedSSEtaProj[k][i]->Draw();
            pad++;
            projectionGridCanvas->cd(pad);
            finalCorrDStarCorrectedASEtaProj[k][i]->Draw();
            
            str1 = corrProjLabel + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fileType;
            
            projectionGridCanvas->SaveAs(str1);
            pad = 1;
        }
    }
  
  
  
    sibCorrBinVzInt[0][2][7]->SetNameTitle("SE Left Sideband", "SE Left Sideband");
    sibCorrBinVzInt[0][2][7]->SetStats(0);
    sibCorrBinVzInt[1][2][7]->SetNameTitle("SE Signal Region", "SE Signal Region");
    sibCorrBinVzInt[1][2][7]->SetStats(0);
    sibCorrBinVzInt[2][2][7]->SetNameTitle("SE Right Sideband", "SE Right Sideband");
    sibCorrBinVzInt[2][2][7]->SetStats(0);
    dStarAngularDistSib[2][1]->SetNameTitle("SE D^{0}+\pi_{soft} Correlation", "SE D^{0} + \pi_{soft} Correlation");
    dStarAngularDistSib[2][1]->SetStats(0);
  
    TCanvas * SECanvas = new TCanvas("", "", 4800, 1100);
    SECanvas->Divide(4,1);
   
    SECanvas->cd(1);
    sibCorrBinVzInt[1][2][7]->GetZaxis()->SetTitle("");
    sibCorrBinVzInt[1][2][7]->Draw("SURF1");
    SECanvas->cd(2);
    sibCorrBinVzInt[0][2][7]->GetZaxis()->SetTitle("");
    sibCorrBinVzInt[0][2][7]->Draw("SURF1");
    SECanvas->cd(3);
    sibCorrBinVzInt[2][2][7]->GetZaxis()->SetTitle("");
    sibCorrBinVzInt[2][2][7]->Draw("SURF1");
    SECanvas->cd(4);
    dStarAngularDistSib[2][1]->GetXaxis()->SetTitle("#Delta#eta");
    dStarAngularDistSib[2][1]->GetYaxis()->SetTitle("#Delta#phi");
    dStarAngularDistSib[2][1]->GetXaxis()->SetTitleSize(.065);
    dStarAngularDistSib[2][1]->GetYaxis()->SetTitleSize(.065); 
    dStarAngularDistSib[2][1]->GetXaxis()->SetTitleOffset(1.5);
    dStarAngularDistSib[2][1]->GetYaxis()->SetTitleOffset(1.5); 
    
    dStarAngularDistSib[2][1]->GetXaxis()->CenterTitle();
    dStarAngularDistSib[2][1]->GetYaxis()->CenterTitle();
    dStarAngularDistSib[2][1]->Draw("SURF1");
    
    SECanvas->SaveAs("SE_Distributions_All_Bands.png");
  
    file->Close();
    output->Close();
   
}


           /*  integralMix = 0;
            integralDStarMix = 0;
            
            for(int binDelEta = 1; binDelEta < NUM_ETA_BINS+1; binDelEta++){                   //this calculates the correct scale factor for the ratio of mixed-events histos
                for(int binDelPhi = 1; binDelPhi < NUM_PHI_BINS+1; binDelPhi++){
                
                    binValueMix = mixCorrCentrality[k][i]->GetBinContent(binDelEta, binDelPhi);
                    binValueDStarMix = dStarAngularDistMix[k][i]->GetBinContent(binDelEta, binDelPhi);
                    
                    if(binValueDStarMix == 0 || binValueMix == 0) { continue; }
                    
                    integralMix += binValueMix;
                    integralDStarMix += binValueDStarMix;
                }
            }
            
            
            
            ratioOfMixed = integralDStarMix/integralMix;    //This gives us the scale factor that comes from the ratio of mixed-event distributions
            cout << "ratio mix to mix: " << ratioOfMixed << endl;
            
            for(int binDelEta = 1; binDelEta < NUM_ETA_BINS+1; binDelEta++){
                for(int binDelPhi = 1; binDelPhi < NUM_PHI_BINS+1; binDelPhi++){
            
                    binValueSib = dStarAngularDistDiff[k][i]->GetBinContent(binDelEta, binDelPhi);
                    //binValueMix = mixCorrCentrality[k][i]->GetBinContent(binDelEta, binDelPhi);
                    binValueDStarMix = dStarAngularDistMix[k][i]->GetBinContent(binDelEta, binDelPhi);
                    if(binValueSib == 0 || binValueDStarMix == 0) {
                    
                        binValueCorrection = 0.0;
                        dStarCorrection[k][i]->SetBinContent(binDelEta, binDelPhi, binValueCorrection);
                        dStarCorrection[k][i]->SetBinError(binDelEta, binDelPhi, 0);
                        continue;
                    }
                    binValueCorrection = (binValueSib/binValueDStarMix);
                    
                    cout << "eta bin: " << binDelEta << "   phi bin: " << binDelPhi << "     binValueCorr: " << binValueCorrection << endl;
                    
                    errorOnMixZeroBin = dStarAngularDistMix[k][i]->GetBinError(binDelEta, binDelPhi);
                    errorOnDStarZeroBin = dStarAngularDistDiff[k][i]->GetBinError(binDelEta, binDelPhi);
                    binValueDStarMix = dStarAngularDistMix[k][i]->GetBinContent(binDelEta, binDelPhi);
                    binValueDStar = dStarAngularDistDiff[k][i]->GetBinContent(binDelEta, binDelPhi);
                    
                    ratioDStarToMix = (binValueDStar/binValueDStarMix);
                    
                    errorEstimate1 = ratioDStarToMix*TMath::Sqrt(((errorOnMixZeroBin*errorOnMixZeroBin)/(binValueDStarMix*binValueDStarMix)) 
                                                               + ((errorOnDStarZeroBin*errorOnDStarZeroBin)/(binValueDStar*binValueDStar)));
                    
                    
                    
                    //dStarCorrection[k][i]->SetBinContent(binDelEta, binDelPhi, binValueCorrection);
                    //dStarCorrection[k][i]->SetBinError(binDelEta, binDelPhi, errorEstimate1);
                    //dStarCorrection[k][i]->Scale(ratioOfMixed);
                }
            } */

            
            /*  integralSib = sibCorrCentrality[k][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
			integralMix = mixCorrCentrality[k][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
            integralDStarSib = dStarAngularDistSib[k][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
            integralDStarMix = dStarAngularDistMix[k][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
            
            correctionScaleFactor = integralDStarMix/integralMix;
            //cout << "ratio mix to mix: " << correctionScaleFactor << endl;
            dStarPairNorm = integralDStarSib/integralDStarMix;
            
            //cout << integralSib/integralMix << endl;
            //cout << "S/S+B: " << 1/SPlusBOverSFit[k][i] << endl;
            scaleFactorOnMix = (1/SPlusBOverSFit[k][i])*(integralSib/integralMix);    //this factor multiplies the TOTAL signal-region mized event to estimate the pair-count
			cout << "mix Scale Factor: " << scaleFactorOnMix << endl;                 // that comes from real D0-hadron mixed pairs
            
            //mixCorrCentrality[k][i]->Scale(scaleFactorOnMix);             */
            
            /*intOfSibDStar = dStarAngularDistSib[k][i]->Integral(1,NUM_ETA_BINS,1,NUM_PHI_BINS);
            intOfMixDStar = dStarAngularDistMix[k][i]->Integral(1,NUM_ETA_BINS,1,NUM_PHI_BINS);
            cout << "int D* sib: " << intOfSibDStar << endl;
            cout << "int D* Mix: " << intOfMixDStar << endl;
            ratioDStarSibToMix = intOfSibDStar/(intOfMixDStar*dStarNormFactor[k][i]);
            cout << "ratio thing mah jig: " << ratioDStarSibToMix << endl;*/