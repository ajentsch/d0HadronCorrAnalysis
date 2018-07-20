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

#include "C:/Users/ajentsch/Desktop/corrHistogramMaker.h"



Bool_t reject;
Double_t fline(Double_t *x, Double_t *par)
{
    if ( (reject && (x[0] > 1.81 && x[0] < 1.91)) || (reject && (x[0] > 1.6 && x[0] < 1.7)) ){     //|| (reject && (x[0] > 2.0 && x[0] < 2.1))
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0];
}

double calcScaleFactorforFit(double sigma)
{
    double preFactor = TMath::Sqrt(TMath::PiOver2());
    
    return preFactor*sigma*.5;
    
}
     
//Double_t SecH(     

int corrHistogramMaker(){
  

    

//-------------------------------------------------------------
//BEGIN VARIABLES SECTION
//-------------------------------------------------------------  
    int NUM_PHI_BINS = 12;
    int NUM_ETA_BINS = 9;
    
    int NUM_CENT_BINS = 16;
    int NUM_VZ_BINS   = 10;
    int NUM_PT_BINS   = 4;
    
    double NUM_FILES_COMBINED = 1;
    
    double ETA_RANGE = 1.55;
	
	double phiBinShift = (TMath::Pi()/12.0);
    
    bool oldFileWithOldCuts = false;
    
    //-----------------------
    //file naming and input/output
    //-----------------------
    
    TString numPhi = "12";
    TString numEta = "9";
    
    TString inputRootFolder = "D0_Hadron_Correlation_Root_Files/";    
    //TString subFolder1 = "corrHistograms_";
    //TString subFolder2 = "_Eta_";
    //TString subFolder3 = "_Phi_bins";
    TString rootFile   = ".root";
    TString outputRoot = "corrHistogramMaker";
    TString slash      = "/";
    
    TString mainPath       = "C:/Users/ajentsch/Desktop/";    
    TString mainFolder     = "D0 Correlation Output Histograms/";
    //TString binSubFolder   = subFolder1 + numEta + subFolder2 + numPhi + subFolder3;
    TString binSubFolder   = "currentDataFile";
    
    //TString inputFileName  = mainPath + inputRootFolder + numEta + subFolder2 + numPhi + subFolder3 + rootFile;
    TString inputFileName  = mainPath + inputRootFolder + binSubFolder + rootFile;
    
    TString path           = mainPath + mainFolder + binSubFolder + slash;
    TString outputFileName = path + outputRoot + rootFile;
    
    TFile *file = new TFile(inputFileName); //Root file with raw histograms to use for calculations
    
    TFile *output = new TFile(outputFileName, "RECREATE"); //root file to store calculated and scaled histograms
    
    
    bool PRINT_SUB_HISTOS_RAW = false;
    bool PRINT_SUB_HISTOS_SCALED = false;
    bool PRINT_SUB_HISTOS_DEL_RHO = false;
    bool PRINT_DELRHO_OVER_REF_HISTOS = true;
    
    
    
    
    TString fileType    = ".png"; //file type for output histogram pictures
    TString fileTypeEps = ".eps"; //file type for progess report
    TString fileTypePDF = ".pdf";
    
    TString binLabelVz[10] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
    TString binLabelCent[16] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"};
    TString binLabelPt[4]    = {"0", "1", "2", ""};
    
    TString PtBinLabel[4]   = {"_PtBin_", "_PtBin_", "_PtBin_", ""};
    
    TString outputFolders[9] = {"SibCorr/", "MixedCorr/", "ScaledMixedCorr/", "SibMinusScaledMixCorr/", 
                                 "delRho_over_rhoRefCorr/", "FullySubtractedCorr/", "USMinusLS/", 
                                 "SideBandSubtraction/", "Centralities/"};
                                
    TString invariantMassOutputFolders[4] = {"PtBin0/", "PtBin1/", "PtBin2/", "PtIntegrated/"};
    TString invariantMassHistogramsFolder = "invariantMassHistograms/";
        
                                
    TString VzSubFolders[11] = {"VzBin0/", "VzBin1/", "VzBin2/", "VzBin3/", "VzBin4/", "VzBin5/", "VzBin6/", "VzBin7/", "VzBin8/", "VzBin9/", "VzIntegrated/"};
    TString PtSubFolders[4] = {"PtBin0/", "PtBin1/", "PtBin2/", "PtIntegrated/"};
    TString bandFolders[5]   = {"LeftSideBand/","UnlikeSign/", "RightSideBand/", "LikeSign/", "SideBandAverage/"};    
    
    
    
    
   
    TH1D*   oneDUSHistos[4][4];     
    TString oneDUSStrings[4] = {"D0_US_invMass_Pt_Bin_0", "D0_US_invMass_Pt_Bin_1", "D0_US_invMass_Pt_Bin_2", "unlikeSign"};
                                              
    TH1D*   oneDLSHistos[4][4];     
    TString oneDLSStrings[4] = {"LS_invMass_Pt_Bin_0", "LS_invMass_Pt_Bin_1", "LS_invMass_Pt_Bin_2", "LikeSignBG"};
                                              
    TH1D*   oneDScaledLSHistos[4][4];     
    TString oneDScaledLSStrings[4] = {"ScaledLS_PtBin_0", "ScaledLS_PtBin_1", "ScaledLS_PtBin_2", "ScaledLS"};
                                                  
    TH1D*   oneDSubtractedInvMassHistos[4][4];                                              
    TString oneDSubtractedInvMassStrings[4] = {"D0_Minus_Scaled_LS_BG_PtBin_0", "D0_Minus_Scaled_LS_BG_PtBin_1", "D0_Minus_Scaled_LS_BG_PtBin_2", "D0_Minus_Scaled_LS_BG"}; 
    
    TH1D*   oneDSubtractedInvMassHistosFunction[4][4];                                              
    TString oneDSubtractedInvMassStringsFunction[4] = {"US_Minus_Expo_Fit_BG_PtBin_0", "US_Minus_Expo_Fit_BG_PtBin_1", "US_Minus_Expo_Fit_BG_PtBin_2", "US_Minus_Expo_Fit_BG"}; 
                                                                                                                          
    TString binLabelCentralityClass[3] = {"_Peripheral", "_MidCentral", "_Central"};
    
    
    double integralSideBandCounts[2];
    double integralCounts[2];
   
    //-----------------------------------
    //Correlation variables and arrays
    //-----------------------------------
    
    // US histograms
    TH2D* sibCorrBin[4][4][10][16];    //Raw US sibling histograms binned by centrality
    TH2D* mixCorrBin[4][4][10][16];    //Raw US mixed histograms binned by centrality
    TH2D* scaledMixCorrBin[4][4][10][16];  //storage for the eventual scaled mixed US histograms  
    TH2D* sibMinusScaledMix[4][4][10][16]; // Raw US Sibling minus scaled mixed                   
    TH2D* delRhoOverRhoRefBin[4][4][10][16];
    
    //ERRORS STORED HERE--ALL ERRORS STORED AS VARIANCES -- SQUARE ROOT MUST BE TAKEN FOR ACTUAL ERROR
    
    double sibCorrBinStatErrors[4][4][10][16][9][12];    
    double mixCorrBinStatErrors[4][4][10][16][9][12];
    double scaledMixCorrBinStatErrors[4][4][10][16][9][12];
    double sibMinusScaledMixStatErrors[4][4][10][16][9][12];
    double delRhoOverRhoRefBinStatErrors[4][4][10][16][9][12];
    //------------------------------------------------------------------
    
    TH2D* sibCorrBinVzInt[4][4][16];
    TH2D* mixCorrBinVzInt[4][4][16];
    TH2D* scaledMixCorrBinVzInt[4][4][16];  
    TH2D* sibMinusScaledMixVzInt[4][4][16]; 
    TH2D* delRhoOverRhoRefBinVzInt[4][4][16];
    
    TH1D* sibCorrBinPhiProj[4][4][10][16];    //Raw US sibling histograms binned by centrality
    TH1D* mixCorrBinPhiProj[4][4][10][16];    //Raw US mixed histograms binned by centrality
    TH1D* scaledMixCorrBinPhiProj[4][4][10][16];  //storage for the eventual scaled mixed US histograms
    TH1D* sibMinusScaledMixPhiProj[4][4][10][16]; // Raw US Sibling minus scaled mixed
    TH1D* delRhoOverRhoRefBinPhiProj[4][4][10][16];
    
    TH1D* sibCorrBinPhiProjVzInt[4][4][16];
    TH1D* mixCorrBinPhiProjVzInt[4][4][16];
    TH1D* scaledMixCorrBinPhiProjVzInt[4][4][16];  
    TH1D* sibMinusScaledMixPhiProjVzInt[4][4][16]; 
    TH1D* delRhoOverRhoRefBinPhiProjVzInt[4][4][16];
  
    
    TH1D* sibCorrUSPhiProjMinusLSPhiProjVzInt[3][16];
   
   
    TH2D* sideBandAverage[4][16];
    TH1D* sideBandAveragePhiProj[4][16];    
    TH2D* fullySubtractedCorrSideBandCent[4][16];
    TH1D* fullySubtractedCorrSideBandCentPhiProj[4][16];
    
    TH2D* delRhoOverRhoRefBinVzIntCentInt[4][4][3];
    TH1D* delRhoOverRhoRefBinVzIntCentIntPhiProj[4][4][3];
    
    TH2D* fullySubtractedCorrLSCent[4][16];
    TH1D* fullySubtractedCorrLSCentPhiProj[4][16];
    
    //TH2D* mixCorrUSVzInt[3][16];
    //TH2D* mixCorrLSVzInt[3][16];
    
    //Correlations with integrated centralities
    TH2D* fullySubtractedCorrCustomCentralityUS[4][3];
    TH2D* fullySubtractedCorrCustomCentralityLS[4][3];
    TH1D* fullySubtractedCorrCustomCentralityLSPhiProj[4][3];
    TH2D* fullySubtractedCorrCustomCentralitySideBand[4][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandPhiProj[4][3];
    double VzWeightFactorCustom[3][10];
    double VzWeightFactorCustomTotal;
    
    TString fullySubtractedLabelCustomCentrality = "FullSubtractedCorr_";
    TString centralityBin[3] = {"Peripheral", "MidCentral", "Central"};
    
    TString SibCorrLabels[4]          = {"Sibling_SideBandLeft_correlation", "Sibling_US_correlation", "Sibling_SideBandRight_correlation", "Sibling_LS_correlation"};
    TString MixCorrLabels[4]          = {"Mixed_SideBandLeft_correlation", "Mixed_US_correlation", "Mixed_SideBandRight_correlation", "Mixed_LS_correlation"};
    TString ScaledMixCorrLabels[4]     = {"Scaled_Mixed_SideBandLeft_correlation", "Scaled_Mixed_US_correlation", "Scaled_Mixed_SideBandRight_correlation", "Scaled_Mixed_LS_correlation"};
    TString SibMinusScaledMixLabels[4] = {"Sib_Minus_Scaled_Mix_SideBandLeft_correlation", "Sib_Minus_Scaled_Mix_US_correlation", "Sib_Minus_Scaled_Mix_SideBandRight_correlation", "Sib_Minus_Scaled_Mix_LS_correlation"};
    TString delRhoOverRhoRefLabels[4]  = {"delRho_over_RhoRef_SideBandLeft_correlation", "delRho_over_RhoRef_US_correlation", "delRho_over_RhoRef_SideBandRight_correlation", "delRho_over_RhoRef_LS_correlation"};
   
    TString subtractedCorrLabels[2] = {"US_Minus_SidebandAverage_Corr", "US_Minus_LS_Corr"};
    
    TString sideBandAverageLabel = "sideBandAverage";
    
    TString phiProj = "Phi_projection_";
    
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
    
    double integralRawMixHistos[4][4][10][16];
    double integralRawMixHistosVzInt[4][4][16];
    //double integralRawMixHistos[4][4][10][16];
    double integralRawMixHistosVzIntCentInt[4][4][3];
   
    double scaleFactorVzCent = 0;
    
    double scaleFactorVz = 0;
    double scaleFactorLS = 0;
    
    
    double LSScaleFactor[4][4] = {{1.0, 1.0, 1.0, 1.0}, {1.0,1.0,1.0,1.0}};
    double ScalingRatio[3][16];
    
    double BOverSPlusBLS[4][4];
    double BOverSPlusBFit[4][4];
    double SPlusBOverSLS[4][4];
    double SPlusBOverSFit[4][4];
    
    double totalPairsSibling = 0;
    double integralError[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
    //fully subtraced US minus LS
    
    TString delRhoOverRhoRefUSMinusLSLabel1 = "US_delRhoOverRho_Minus_LS_delRhoOverRho_";
    
    double v2Array[4][3];
    double v2ArrayStatError[4][3];
    
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
    
    TCanvas *c = new TCanvas("c2", "Histograms", 1100, 850);

    double countsPeakRegionUS = 0;
    double countsPeakRegionLS = 0;
    double countsSideBandUS   = 0;
    double countsSideBandLS   = 0;
    double sideBandLow      = 2.0;
    double sideBandHigh     = 2.1;
    double massLow          = 1.82;
    double massHigh         = 1.90;
	
	double sqrtSPlusB       = 1;
    
    double integralFit      = 0;
	
	double DMesonSig        = 0;
    
    double integral         = 0;
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
    histOfCuts->SetBinContent(49, 1); //Scale factor counter to produce text file*/
    
    if(!oldFileWithOldCuts){ 
    
        double numEvents = 0;
        
        numEvents = (((TH1I*) file->Get("number_of_events_used"))->GetBinContent(2))/(1000000);
        
        ofstream cutFile;
        TString cutFileName = "Cuts_and_data_information.txt";
        TString cutOutputFile = path + cutFileName;
        cutFile.open (cutOutputFile);
        
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
        
        cutFile << endl << "D0 pt 0-1 GeV/c cuts" << endl;
        
        cutFile << "D0 DecayLengthMinimum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(12) << endl;
        cutFile << "D0 DaughterDCAMaximum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(13) << endl;
        cutFile << "D0 DCA to PV Maximum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(14) << endl;
        cutFile << "Daughter Pion DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(15) << endl;
        cutFile << "Daughter Kaon DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(16) << endl;
        
        cutFile << endl << "D0 pt 1-2 GeV/c cuts" << endl;
        
        cutFile << "D0 DecayLengthMinimum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(17) << endl;
        cutFile << "D0 DaughterDCAMaximum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(18) << endl;
        cutFile << "D0 DCA to PV Maximum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(19) << endl;
        cutFile << "Daughter Pion DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(20) << endl;
        cutFile << "Daughter Kaon DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(21) << endl;
        
        cutFile << endl << "D0 pt 2-3 GeV/c cuts" << endl;
        
        cutFile << "D0 DecayLengthMinimum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(22) << endl;
        cutFile << "D0 DaughterDCAMaximum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(23) << endl;
        cutFile << "D0 DCA to PV Maximum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(24) << endl;
        cutFile << "Daughter Pion DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(25) << endl;
        cutFile << "Daughter Kaon DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(26) << endl;
        
        cutFile << endl << "D0 pt 3-5 GeV/c cuts" << endl;
        
        cutFile << "D0 DecayLengthMinimum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(27) << endl;
        cutFile << "D0 DaughterDCAMaximum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(28) << endl;
        cutFile << "D0 DCA to PV Maximum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(29) << endl;
        cutFile << "Daughter Pion DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(30) << endl;
        cutFile << "Daughter Kaon DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(31) << endl;
        
        cutFile << endl << "D0 pt 5-10 GeV/c cuts" << endl;
        
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
        
    }
    
    else if(oldFileWithOldCuts){ 
    
        double numEvents = 0;
    
        numEvents = (((TH1I*) file->Get("number of events used"))->GetBinContent(2))/(1000000);
    
        ofstream cutFile;
        TString cutFileName = "Cuts_and_data_information.txt";
        TString cutOutputFile = path + cutFileName;
        cutFile.open (cutOutputFile);
    
        cutFile << "Important information about this data run" << endl << endl;
        cutFile << "Number of Events: " << numEvents << "M Events" << endl << endl;
    
        NUM_FILES_COMBINED = ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(29);
    
        ((TH1D*) file->Get("HistOfCuts"))->Scale(1/NUM_FILES_COMBINED);
    
        cutFile << "Trigger D0 Cuts" << endl << endl;
        cutFile << "D0 InvMass Signal Band -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(1) << "\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(2) << endl;
        cutFile << "D0 SideBandLeft   Band -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(3) << "\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(4) << endl;
        cutFile << "D0 SideBandRight  Band -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(5) << "\t\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(6) << endl;
        cutFile << "D0 Pt Cuts             -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(7) << "\t\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(8) << endl;
        cutFile << "D0 DecayLengthCut      -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(9) << "\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(10) << endl;
        cutFile << "D0 DaughterPtCut       -\t" << "Pi : " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(12) << "\t" << "K  : " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(13) << endl;
        cutFile << "D0 K/Pi DCA to PV      -\t" << "Pi : " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(15) << "\t" << "K  : " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(14) << endl;
        cutFile << "DaughterDCA            -\t" << "Less Than : " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(11) << endl;
        cutFile << "D0 DCA to PV           -\t" << "Greater Than: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(16) << endl;
    
        cutFile << "Associated hadron Cuts" << endl << endl;
        cutFile << "associated hadron Pt   -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(17) << "\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(18) << endl;
        cutFile << "Num of events used to mix  -\t" << ": " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(19) << endl;
        cutFile << "HFT Tracks only??  -\t" << "YES(1), NO(0): " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(22) << endl;
        cutFile << "Track Chi2 Max:  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(23) << endl;
        cutFile << "Track DCA Max:  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(24) << endl;
    
        cutFile.close();
        
    }


//---------------------------------------------    
//BEGIN INVARIANT MASS HISTOGRAM INFORMATION   
//---------------------------------------------
    
    ofstream invariantMassInfo;
    TString invariantMassInfoName = "Signal_and_background_information.txt";
    TString invariantMassInfoFile = path + invariantMassInfoName;
    invariantMassInfo.open (invariantMassInfoFile);
    
    
    for(int k = 0; k < NUM_PT_BINS; k++){  //Initialize the histograms 
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
       
            oneDUSHistos[k][i]->GetXaxis()->SetTitle("Invariant Mass GeV/c^{2}");
            oneDUSHistos[k][i]->GetYaxis()->SetTitle("counts");
            oneDUSHistos[k][i]->GetYaxis()->SetTitleOffset(1.4);
            oneDLSHistos[k][i]->GetXaxis()->SetTitle("Invariant Mass GeV/c^{2}");
            oneDLSHistos[k][i]->GetYaxis()->SetTitle("counts");
            oneDLSHistos[k][i]->GetYaxis()->SetTitleOffset(1.4);
            oneDScaledLSHistos[k][i]->GetXaxis()->SetTitle("Invariant Mass GeV/c^{2}");
            oneDScaledLSHistos[k][i]->GetYaxis()->SetTitle("counts");
            oneDScaledLSHistos[k][i]->GetYaxis()->SetTitleOffset(1.4);
            oneDSubtractedInvMassHistos[k][i]->GetXaxis()->SetTitle("Invariant Mass GeV/c^{2}");
            oneDSubtractedInvMassHistos[k][i]->GetYaxis()->SetTitle("counts");
            oneDSubtractedInvMassHistos[k][i]->GetYaxis()->SetTitleOffset(1.4);
        
            xAxisUS = oneDUSHistos[k][i]->GetXaxis();
            xAxisLS = oneDLSHistos[k][i]->GetXaxis();
        
            binMassLowUS  = oneDUSHistos[k][i]->GetXaxis()->FindBin(sideBandLow);            //get normalization factors to scale the LS distribution
            binMassHighUS = oneDUSHistos[k][i]->GetXaxis()->FindBin(sideBandHigh);
            binMassLowLS = oneDLSHistos[k][i]->GetXaxis()->FindBin(sideBandLow);            //get normalization factors to scale the LS distribution
            binMassHighLS = oneDLSHistos[k][i]->GetXaxis()->FindBin(sideBandHigh);
        
            countsSideBandUS = oneDUSHistos[k][i]->Integral(binMassLowUS, binMassHighUS);
            countsSideBandLS = oneDLSHistos[k][i]->Integral(binMassLowLS, binMassHighLS);
        
            //Calculate peak region values
        
            binMassLowUS  = xAxisUS->FindBin(massLow);            
            binMassHighUS = xAxisUS->FindBin(massHigh);
            binMassLowLS  = xAxisLS->FindBin(massLow);            
            binMassHighLS = xAxisLS->FindBin(massHigh);
        
            oneDUSHistos[k][i]->Sumw2();
            oneDUSHistos[k][i]->SetMarkerStyle(20);
            oneDLSHistos[k][i]->Sumw2();
            oneDLSHistos[k][i]->SetMarkerStyle(20);
            oneDSubtractedInvMassHistos[k][i]->Sumw2();
            oneDSubtractedInvMassHistos[k][i]->SetMarkerStyle(20);
        
            oneDUSHistos[k][i]->Draw();
                
            tmp = path + invariantMassHistogramsFolder + invariantMassOutputFolders[k] + oneDUSStrings[k] + binLabelCentralityClass[i] + fileType;
            c->SaveAs(tmp);
        
            //oneDUSHistos[i]->SetTitle("");
            //oneDUSHistos[i]->Draw();
            //tmp = path + oneDUSStrings[i] + fileTypeEps;
            //c->SaveAs(tmp);
        
            LSScaleFactor[k][i] = countsSideBandUS/countsSideBandLS;
        
            oneDScaledLSHistos[k][i]->Scale(LSScaleFactor[k][i]);            //normalize LS spectrum here
        
            countsPeakRegionUS = oneDUSHistos[k][i]->Integral(binMassLowUS, binMassHighUS);
            countsPeakRegionLS = oneDScaledLSHistos[k][i]->Integral(binMassLowLS, binMassHighLS);  //this gets an estimate for B from the LS scaled to the BG in the US, using a single sideband
        
            invariantMassInfo << "____________________________________________________________________________________________________" << endl;
            invariantMassInfo << "PtBin (0: 0-1, 1: 1-4, 2: 4-20, 3: all): " << k << endl;
            invariantMassInfo << "CentBin (0: per, 1: mid, 2: cent, 3: all): " << i << endl;
            //cout << "LS Scale Factor: " << LSScaleFactor[i] << endl;
            invariantMassInfo << "S+B (from integral in mass window): " << countsPeakRegionUS << endl;
        
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
        
            invariantMassInfo << "S (from fit estimate): " << integralFit << endl;
            invariantMassInfo << "B (from fit estimate): " << (countsPeakRegionUS-integralFit) << endl;
            invariantMassInfo << "S+B/S from LS subtraction + expo Fit:  " << SPlusBOverSFit[k][i] << endl;
            invariantMassInfo << "B/S+B from LS subtraction + expo Fit: " << BOverSPlusBFit[k][i] << endl;
        
            oneDUSHistos[k][i]->Draw();
            oneDScaledLSHistos[k][i]->Draw("SAME");
            tmp = path + invariantMassHistogramsFolder + invariantMassOutputFolders[k] + oneDUSStrings[k] + binLabelCentralityClass[i] + withFitLabel + fileType;
            c->SaveAs(tmp);
            oneDLSHistos[k][i]->Draw();
            tmp = path + invariantMassHistogramsFolder + invariantMassOutputFolders[k] + oneDLSStrings[k] + binLabelCentralityClass[i] + fileType;
            c->SaveAs(tmp);
            oneDSubtractedInvMassHistos[k][i]->Draw();
            tmp = path + invariantMassHistogramsFolder + invariantMassOutputFolders[k] + oneDSubtractedInvMassStrings[k] + binLabelCentralityClass[i] + fileType;
            c->SaveAs(tmp);
			
			sqrtSPlusB = TMath::Sqrt(countsPeakRegionUS);
			DMesonSig = integralFit/sqrtSPlusB;
			tmp = "#frac{S}{#sqrt{S+B}} = ";
			str1 = Form("%.2f",DMesonSig);
			str2 = tmp + str1;
			
			invariantMassInfo << "D0 Significance:   " << DMesonSig << endl;
			
			TPaveText *textBoxSignifigance = new TPaveText(0.25, 0.75, .45, .95, "NB NDC");
            textBoxSignifigance->AddText(str2);
            textBoxSignifigance->GetLine(0)->SetTextSize(.041);
            
			oneDSubtractedInvMassHistosFunction[k][i]->Draw();
			textBoxSignifigance->Draw("SAME");
            tmp = path + invariantMassHistogramsFolder + invariantMassOutputFolders[k] + oneDSubtractedInvMassStringsFunction[k] + binLabelCentralityClass[i] + fileType;
            c->SaveAs(tmp);
        }
    }
    
    invariantMassInfo.close();

    cout << endl << "Invariant Mass histograms process and S & B calculated for each Pt/Centrality bin" << endl << endl;
    
//------------------------------------------
//BEGIN CORRELATION HISTOGRAM INFORMATION
//------------------------------------------

    /**********************************************************************************
        EXTRACT SIBLING AND MIXED HISTOGRAMS FROM RAW DATA FILE HERE
    **********************************************************************************/
    for(int band = 0; band < 4; band++){ //begin loop to get raw sibling and mixed histograms in both US and LS
        for(int i = 0; i < NUM_CENT_BINS; i++){
            for(int k = 0; k < NUM_PT_BINS; k++){
    
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
                
                
                if(k==3){totalPairsSibling = totalPairsSibling + sibCorrBin[band][k][j][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);}
                
                sibCorrBinVzInt[band][k][i]->Add(sibCorrBin[band][k][j][i]);
                formatCorrHist(sibCorrBin[band][k][j][i]);  //histogram formatting for axes
                mixCorrBinVzInt[band][k][i]->Add(mixCorrBin[band][k][j][i]);
                formatCorrHist(mixCorrBin[band][k][j][i]);  //histogram formatting for axes
                
            
                //Sibling US
                if(PRINT_SUB_HISTOS_RAW){
                    sibCorrBin[band][k][j][i]->Draw("SURF1");                              //sibling 2D US histogram
                    str1 = path + outputFolders[0] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + SibCorrLabels[band] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;
                    c->SaveAs(str1);
                }    
        
                str1 = phiProj + SibCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];         //sibling phi projection US
                sibCorrBinPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);  
                sibCorrBinPhiProj[band][k][j][i] = (TH1D*)sibCorrBin[band][k][j][i]->ProjectionY();  
                if(PRINT_SUB_HISTOS_RAW){
                    str1 = path + outputFolders[0] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + phiProj + SibCorrLabels[band] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;        
                    sibCorrBinPhiProj[band][k][j][i]->Draw();
                    c->SaveAs(str1);
                }
                
                sibCorrBinPhiProjVzInt[band][k][i]->Add(sibCorrBinPhiProj[band][k][j][i]);
                         
                //Mixed 
                if(PRINT_SUB_HISTOS_RAW){
                    mixCorrBin[band][k][j][i]->Draw("SURF1");                              //mixling 2D US histogram
                    str1 = path + outputFolders[1] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + MixCorrLabels[band] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;
                    c->SaveAs(str1);
                }
                
                str1 = phiProj + MixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];     //mixling phi projection US
                mixCorrBinPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);  
                mixCorrBinPhiProj[band][k][j][i] = (TH1D*)mixCorrBin[band][k][j][i]->ProjectionY();  
                
                if(PRINT_SUB_HISTOS_RAW){
                    str1 = path + outputFolders[1] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + phiProj + MixCorrLabels[band] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;        
                    mixCorrBinPhiProj[band][k][j][i]->Draw();
                    c->SaveAs(str1);
                }
        
                mixCorrBinPhiProjVzInt[band][k][i]->Add(mixCorrBinPhiProj[band][k][j][i]);
                         
            }
        
            if(PRINT_SUB_HISTOS_RAW){
                
                //sibling
                sibCorrBinVzInt[band][k][i]->Draw("SURF1");
                str1 = path + outputFolders[0] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + SibCorrLabels[band] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;
                c->SaveAs(str1);
        
                sibCorrBinPhiProjVzInt[band][k][i]->Draw();
                str1 = path + outputFolders[0] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + phiProj + SibCorrLabels[band] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;
                c->SaveAs(str1);
        
                //mixed
                mixCorrBinVzInt[band][k][i]->Draw("SURF1");
                str1 = path + outputFolders[1] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + MixCorrLabels[band] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;
                c->SaveAs(str1);
        
                mixCorrBinPhiProjVzInt[band][k][i]->Draw();
                str1 = path + outputFolders[1] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + phiProj + MixCorrLabels[band] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;
                c->SaveAs(str1);
        
           }
       }
    }//end loop to get raw sibling and mixed histograms in both US and LS
    
    cout << "Total sibling pairs for band " << band << " : " << totalPairsSibling << endl;
    
   }

   cout << endl << "Sibling and mixed histograms extracted" << endl;
   
    /**********************************************************************************
        NORMALIZE THE MIXED HISTOGRAMS HERE
    **********************************************************************************/
   
  for(int band = 0; band < 4; band++){ // begin loop to make scaled mixed histograms
    for(int i = 0; i < NUM_CENT_BINS; i++){ 
        for(int k = 0; k < NUM_PT_BINS; k++){
            
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
         
         
                //US scaled histograms
                
                str1 = ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
        
                formatCorrHist(scaledMixCorrBin[band][k][j][i], str1);  //formatting
        
                if(PRINT_SUB_HISTOS_SCALED){
                    scaledMixCorrBin[band][k][j][i]->Draw("SURF1");
                    str1 = path + outputFolders[2] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;
                    c->SaveAs(str1);
                }
                
                str1 = phiProj + ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];       //scaled mix phi projection US 
                scaledMixCorrBinPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);  
                scaledMixCorrBinPhiProj[band][k][j][i] = (TH1D*)scaledMixCorrBin[band][k][j][i]->ProjectionY();  
                
                if(PRINT_SUB_HISTOS_SCALED){
                    str1 = path + outputFolders[2] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + phiProj + ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;        
                    scaledMixCorrBinPhiProj[band][k][j][i]->Draw();
                    c->SaveAs(str1);
                }
            }
                
            if(PRINT_SUB_HISTOS_SCALED){
            
                str1 = path + outputFolders[2] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;        
                scaledMixCorrBinVzInt[band][k][i]->Draw("SURF1");
                c->SaveAs(str1);
        
                str1 = path + outputFolders[2] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + phiProj + ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;        
                scaledMixCorrBinPhiProjVzInt[band][k][i]->Draw();
                c->SaveAs(str1);
        
           }
        }
      }  //end loop to make scaled mixed histograms  */
   }
   //--------------------------------------------------------calculate some useful quantities here--------------------------------------
   
   cout << endl << "Mixed Histograms Normalized." << endl << endl;
   
   for(int band = 0; band < 4; band++){                               //build arrays to store information for adding the various alpha's (Vz, centrality, etc.)
       for(int k = 0; k < NUM_PT_BINS; k++){
           
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
            for(int k = 0; k < NUM_PT_BINS; k++){
            
                str1 = SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                sibMinusScaledMixVzInt[band][k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                str1 = phiProj + SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                sibMinusScaledMixPhiProjVzInt[band][k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
           
                for(int j = 0; j < 10; j++){
        
                
                    str1 = SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
                
                    sibMinusScaledMix[band][k][j][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                    sibMinusScaledMix[band][k][j][i]->Add(sibCorrBin[band][k][j][i], scaledMixCorrBin[band][k][j][i], 1, -1);
                    //integralError = 0.0;
                    for(int etaBin = 0; etaBin < NUM_ETA_BINS; etaBin++){
                        for(int phiBin = 0; phiBin < NUM_PHI_BINS; phiBin++){
                        
                            sibMinusScaledMixStatErrors[band][k][j][i][etaBin][phiBin] = sibCorrBinStatErrors[band][k][j][i][etaBin][phiBin] + scaledMixCorrBinStatErrors[band][k][j][i][etaBin][phiBin];    
                            integralError[phiBin] = integralError[phiBin] + sibMinusScaledMixStatErrors[band][k][j][i][etaBin][phiBin];
                        }
                    }//Stat errors
                
                      
                
                    formatCorrHist(sibMinusScaledMix[band][k][j][i]); //formatting
        
                    if(PRINT_SUB_HISTOS_DEL_RHO){
                        sibMinusScaledMix[band][k][j][i]->Draw("SURF1");
                        str1 = path + outputFolders[3] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + SibMinusScaledMixLabels[band] + binLabelCent[i] + fileType;
                        c->SaveAs(str1);
                   
                    }
                   
                    str1 = phiProj + SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
                    sibMinusScaledMixPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);  
                    sibMinusScaledMixPhiProj[band][k][j][i] = (TH1D*)sibMinusScaledMix[band][k][j][i]->ProjectionY();  
                
                    if(PRINT_SUB_HISTOS_DEL_RHO){
                        str1 = path + outputFolders[3] +  bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + phiProj + SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;        
                        
                        for(int phiBin = 0; phiBin < NUM_PHI_BINS; phiBin++){
                                                
                            sibMinusScaledMixPhiProj[band][k][j][i]->SetBinError(phiBin+1, integralError[phiBin]);
                        }
                        
                        sibMinusScaledMixPhiProj[band][k][j][i]->Draw("E0");
                        c->SaveAs(str1);
                    }
                
                    sibMinusScaledMixVzInt[band][k][i]->Add(sibMinusScaledMix[band][k][j][i]);
                    sibMinusScaledMixPhiProjVzInt[band][k][i]->Add(sibMinusScaledMixPhiProj[band][k][j][i]);
                }
     
                if(PRINT_SUB_HISTOS_DEL_RHO){ 
                    str1 = path+outputFolders[3]+bandFolders[band]+PtSubFolders[k]+VzSubFolders[10]+SibMinusScaledMixLabels[band]+PtBinLabel[k]+binLabelPt[k]+VzIntLabel + centBinLabel + binLabelCent[i] + fileType;
                    sibMinusScaledMixVzInt[band][k][i]->Draw("SURF1");
                    c->SaveAs(str1);
             
                    str1 = path + outputFolders[3] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + phiProj + SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;
                    sibMinusScaledMixPhiProjVzInt[band][k][i]->Draw();
                    c->SaveAs(str1);
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
   
   
    for(int band = 0; band < 4; band++){     // begin loop to make delRho/RhoRef US histos -- this will calculate the delRho/Rho histograms for all 4 bands -- LS, US, SBR, SBL
        for(int i = 0; i < NUM_CENT_BINS; i++){ 
            for(int k = 0; k < NUM_PT_BINS; k++){
            
                str1 = delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                delRhoOverRhoRefBinVzInt[band][k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                str1 = phiProj + delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                delRhoOverRhoRefBinPhiProjVzInt[band][k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            
                for(int j = 0; j < 10; j++){

                    str1 = delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
                    delRhoOverRhoRefBin[band][k][j][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                    
                    if(scaledMixCorrBin[band][k][j][i]->GetEntries() > 0){
                    
                        delRhoOverRhoRefBin[band][k][j][i]->Divide(sibMinusScaledMix[band][k][j][i], scaledMixCorrBin[band][k][j][i], 1, 1);
                    }
                    
                    formatCorrHist(delRhoOverRhoRefBin[band][k][j][i]);
                    //VzScaleFactorUS[j] = delRhoOverRhoRefUSBin[k][j][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
                    delRhoOverRhoRefBin[band][k][j][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                
                    if(PRINT_DELRHO_OVER_REF_HISTOS){ 
                        delRhoOverRhoRefBin[band][k][j][i]->Draw("SURF1");
                        str1 = path + outputFolders[4] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;
                        c->SaveAs(str1);
                    }    
        
                    str1 = phiProj + delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];        //sib - scaled mix phi projection US 
                    delRhoOverRhoRefBinPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);  
                    delRhoOverRhoRefBinPhiProj[band][k][j][i] = (TH1D*)delRhoOverRhoRefBin[band][k][j][i]->ProjectionY();  
                
                    if(PRINT_DELRHO_OVER_REF_HISTOS){               
                        str1 = path + outputFolders[4] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + phiProj + delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;      
                        delRhoOverRhoRefBinPhiProj[band][k][j][i]->Draw();
                        c->SaveAs(str1);
                    }
                
                    scaleFactorVz = integralRawMixHistos[band][k][j][i]/integralRawMixHistosVzInt[band][k][i];
                    delRhoOverRhoRefBinVzInt[band][k][i]->Add(delRhoOverRhoRefBin[band][k][j][i], scaleFactorVz);
                
                }
 
                str1 = path + outputFolders[4] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType; 
                delRhoOverRhoRefBinVzInt[band][k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                if(PRINT_DELRHO_OVER_REF_HISTOS){ 
                    delRhoOverRhoRefBinVzInt[band][k][i]->Draw("SURF1");
                    c->SaveAs(str1);
                    
                }
            
                delRhoOverRhoRefBinPhiProjVzInt[band][k][i] = (TH1D*)delRhoOverRhoRefBinVzInt[band][k][i]->ProjectionY();
            
                str1 = path + outputFolders[4] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + phiProj + delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;
                if(PRINT_DELRHO_OVER_REF_HISTOS){ 
                    delRhoOverRhoRefBinPhiProjVzInt[band][k][i]->Draw();
                    c->SaveAs(str1);
                }    
            }
        } //end loop to make delRho/Rho_ref histograms*/
    }
    
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------THIS SECTION DOES THE FINAL SUBTRACTION IN THE 16 MULTIPLICITY BINS------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------
    
        for(int i = 0; i < NUM_CENT_BINS; i++){// begin side band subtraction here -- this will take the US delRho/Rho and subtract the SB average delRho/Rho from it
            for(int k = 0; k < NUM_PT_BINS; k++){
           
                str1 = sideBandAverageLabel + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                sideBandAverage[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                str1 = phiProj + sideBandAverageLabel + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                sideBandAveragePhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                
                sideBandAverage[k][i]->Add(delRhoOverRhoRefBinVzInt[0][k][i], delRhoOverRhoRefBinVzInt[2][k][i], .5, .5);
                sideBandAveragePhiProj[k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[0][k][i], delRhoOverRhoRefBinPhiProjVzInt[2][k][i], .5, .5);
                
                str1 = subtractedCorrLabels[0] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                fullySubtractedCorrSideBandCent[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                str1 = phiProj + subtractedCorrLabels[0] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                fullySubtractedCorrSideBandCentPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
           
                fullySubtractedCorrSideBandCent[k][i]->Add(delRhoOverRhoRefBinVzInt[1][k][i], sideBandAverage[k][i], 1, -BOverSPlusBFit[k][3]);
               
                str1 = path + outputFolders[4] + bandFolders[4] + PtSubFolders[k] + VzSubFolders[10] + sideBandAverageLabel + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType; 
                sideBandAverage[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                sideBandAverage[k][i]->Draw("SURF1");
                c->SaveAs(str1);
                
                str1 = path + outputFolders[5] + outputFolders[7] + PtSubFolders[k] + subtractedCorrLabels[0] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i] + fileType; 
                fullySubtractedCorrSideBandCent[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                fullySubtractedCorrSideBandCent[k][i]->Draw("SURF1");
                c->SaveAs(str1);
                fullySubtractedCorrSideBandCent[k][i]->Write();
                
                str1 = path + outputFolders[5] + outputFolders[7] + PtSubFolders[k] + phiProj + subtractedCorrLabels[0] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i] + fileType; 
                fullySubtractedCorrSideBandCentPhiProj[k][i] = (TH1D*)fullySubtractedCorrSideBandCent[k][i]->ProjectionY();
                fullySubtractedCorrSideBandCentPhiProj[k][i]->Draw();
                c->SaveAs(str1);
           
           
            }    
        }// end side band subtraction here
        
        
        for(int i = 0; i < NUM_CENT_BINS; i++){// begin LS subtraction here  -- this will take the US delRho/Rho and subtract the LS delRho/Rho from it
            for(int k = 0; k < NUM_PT_BINS; k++){
           
                str1 = subtractedCorrLabels[1] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                fullySubtractedCorrLSCent[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                str1 = phiProj + subtractedCorrLabels[1] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                fullySubtractedCorrLSCentPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
           
                fullySubtractedCorrLSCent[k][i]->Add(delRhoOverRhoRefBinVzInt[1][k][i], delRhoOverRhoRefBinVzInt[3][k][i], 1, -BOverSPlusBFit[k][3]);
                
                //fullySubtractedCorrLSCentPhiProj[k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[1][k][i], delRhoOverRhoRefBinPhiProjVzInt[3][k][i], 1, -BOverSPlusB[k]);
            
                str1 = path + outputFolders[5] + outputFolders[6] + PtSubFolders[k] + subtractedCorrLabels[0] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i] + fileType; 
                fullySubtractedCorrLSCent[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                fullySubtractedCorrLSCent[k][i]->Draw("SURF1");
                c->SaveAs(str1);
                
                str1 = path + outputFolders[5] + outputFolders[6] + PtSubFolders[k] + phiProj + subtractedCorrLabels[0] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i] + fileType; 
                fullySubtractedCorrLSCentPhiProj[k][i] = (TH1D*) fullySubtractedCorrLSCent[k][i]->ProjectionY();
                fullySubtractedCorrLSCentPhiProj[k][i]->Draw();
                c->SaveAs(str1);
           
           
            }    
        }// end LS subtraction here
          

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------      
    //---------------------------THIS SECTION DOES THE FINAL SUBTRACTION IN THE 3 centrality BINS----------------------------------------------------------------      
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------

    
          
      
    
    TString bandSubtract[3]          = {"LS_","SideBand_", "US"};
    int temp = 0;
    int centralityCombinations[6] = {0, 2, 3, 8, 9, 15}; //This was calculated for combining multiplicity bins into the proper centrality bins
    int numBinsInCentrality = 0;
   
    for(int k = 0; k < NUM_PT_BINS; k++){
        for(int i = 0; i < 3; i++){
           
           
           str1 =  fullySubtractedLabelCustomCentrality + bandSubtract[2] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
           fullySubtractedCorrCustomCentralityUS[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
           
           str1 =  fullySubtractedLabelCustomCentrality + bandSubtract[0] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
           fullySubtractedCorrCustomCentralityLS[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
           str1 =  fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
           fullySubtractedCorrCustomCentralitySideBand[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
           str1 =  phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[0] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
           fullySubtractedCorrCustomCentralityLSPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
           str1 =  phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
           fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
           
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
                
                scaleFactorVzCent = integralRawMixHistosVzInt[3][k][temp]/integralRawMixHistosVzIntCentInt[3][k][i];                 //LS
                delRhoOverRhoRefBinVzIntCentInt[3][k][i]->Add(delRhoOverRhoRefBinVzInt[3][k][temp], scaleFactorVzCent);
                delRhoOverRhoRefBinVzIntCentIntPhiProj[3][k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[3][k][temp], scaleFactorVzCent);
                
                scaleFactorVzCent = integralRawMixHistosVzInt[2][k][temp]/integralRawMixHistosVzIntCentInt[2][k][i];                 //Right SB
                delRhoOverRhoRefBinVzIntCentInt[2][k][i]->Add(delRhoOverRhoRefBinVzInt[2][k][temp], scaleFactorVzCent);
                delRhoOverRhoRefBinVzIntCentIntPhiProj[2][k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[2][k][temp], scaleFactorVzCent);
                
                scaleFactorVzCent = integralRawMixHistosVzInt[0][k][temp]/integralRawMixHistosVzIntCentInt[0][k][i];                 //Left SB
                delRhoOverRhoRefBinVzIntCentInt[0][k][i]->Add(delRhoOverRhoRefBinVzInt[0][k][temp], scaleFactorVzCent);
                delRhoOverRhoRefBinVzIntCentIntPhiProj[0][k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[0][k][temp], scaleFactorVzCent);
            
           }
       }    
    }   
      
       double phiProjectionScaleFactor = 1.0/9.0;
       
        for(int k = 0; k < 3; k++){    //k = 3 will give the LS, which we don't need right now
            for(int i = 0; i < 3; i++){
      
                cout << "k, i: " << k << " , " << i << endl;
               
                fullySubtractedCorrCustomCentralityLS[k][i]->Add(delRhoOverRhoRefBinVzIntCentInt[1][k][i], delRhoOverRhoRefBinVzIntCentInt[3][k][i], 1, -BOverSPlusBFit[k][i]);
                fullySubtractedCorrCustomCentralityLS[k][i]->Scale(SPlusBOverSFit[k][i]);
                
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Add(delRhoOverRhoRefBinVzIntCentInt[1][k][i]);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Add(delRhoOverRhoRefBinVzIntCentInt[0][k][i], -.5*BOverSPlusBFit[k][i]);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Add(delRhoOverRhoRefBinVzIntCentInt[2][k][i], -.5*BOverSPlusBFit[k][i]);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Scale(SPlusBOverSFit[k][i]);
                
                fullySubtractedCorrCustomCentralityLSPhiProj[k][i] = (TH1D*) fullySubtractedCorrCustomCentralityLS[k][i]->ProjectionY();
                fullySubtractedCorrCustomCentralityLSPhiProj[k][i]->Scale(phiProjectionScaleFactor);
                
                fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i] = (TH1D*) fullySubtractedCorrCustomCentralitySideBand[k][i]->ProjectionY();
                fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Scale(phiProjectionScaleFactor);
       
                str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + fullySubtractedLabelCustomCentrality + bandSubtract[0] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fileType;
                fullySubtractedCorrCustomCentralityLS[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                fullySubtractedCorrCustomCentralityLS[k][i]->Draw("SURF1");
                c->SaveAs(str1);
                
                fullySubtractedCorrCustomCentralityLS[k][i]->Write();
           
                str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fileType;
                fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Draw("SURF1");
                c->SaveAs(str1);
                
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Write();
           
                str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[0] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fileType;
                fullySubtractedCorrCustomCentralityLSPhiProj[k][i]->Draw();
                c->SaveAs(str1);
           
                str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fileType;
                fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Draw();
                c->SaveAs(str1);
            }    
        }
 
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
 
    ofstream dataPointsForLanny;
    TString dataPointsForLannyFile = "Data_Points_PtIntegrated_Central.txt";
    TString dataPointsForLannyOutput = path + dataPointsForLannyFile;
    dataPointsForLanny.open (dataPointsForLannyOutput);
    
    dataPointsForLanny << "Eta Bin" << "\t" << "Phi Bin" << "\t" << "Content" << "\t" << "Error" << endl;
    
    for(int etaBin = 2; etaBin < 9; etaBin++){
        for(int phiBin = 1; phiBin < 13; phiBin++){
    
            dataPointsForLanny << (etaBin-1) << "\t" << phiBin << "\t" << fullySubtractedCorrCustomCentralitySideBand[3][2]->GetBinContent(etaBin, phiBin) << "\t" <<
            fullySubtractedCorrCustomCentralitySideBand[3][2]->GetBinError(etaBin, phiBin) << endl;
        }
    }        
    
    
    dataPointsForLanny.close();
    
    cout << endl << "Data file for Lanny's fitting code produced." << endl;
 
 
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------FITTING AND PARAMETER EXTRACTION DONE HERE---------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
 
    TCanvas *c = new TCanvas("c3", "Histograms", 1100, 850);
 
    TString fitFunctionLabel = "parameterFitFunction";
 
    //TF2 *myfit   = new TF2("parameterFitFunction","[0] + [1]*cos(y) + [2]*2*cos(2*y) + [3]*(exp(-0.5*(y*y)/([4]*[4]))+exp(-0.5*((y-6.28)*(y-6.28))/([4]*[4])))*exp(-0.5*(x*x)/([5]*[5]))", -ETA_RANGE, ETA_RANGE,-1.57,4.71);
    
    TF2 *myfit[4][3];  

    TF1 *myfit1D = new TF1("parameterFitFunction","[0] + [1]*cos(x) + [2]*2*cos(2*x)" , -1.57, 4.71);
    TF1 *fitProjection = new TF1("parameterFitFunction","[0] + [1]*cos(x) + [2]*2*cos(2*x) + [3]*[4]*(exp(-0.5*(x*x)/([5]*[5]))+exp(-0.5*((x-6.28)*(x-6.28))/([5]*[5])))" , -1.57, 4.71);
    TF1 *quadrupole1    = new TF1("Quadrupole", "[0]*2*cos(2*x)", -1.57, 4.71);
    TF1 *quadrupole2    = new TF1("Quadrupole", "[0]*2*cos(2*x)", -1.57, 4.71);
    TF1 *quadrupole3    = new TF1("Quadrupole", "[0] + [1]*2*cos(2*x)", -1.57, 4.71);
    TF1 *pileupLine1    = new TF1("PileupLine1", "[0] + [1]*x", -2.0, -1.0);
    TF1 *pileupLine2    = new TF1("PileupLine2", "[0] + [1]*x", -1.0, 0);
    TF1 *pileupLine3    = new TF1("PileupLine3", "[0] + [1]*x", 0, 1.0);
    TF1 *pileupLine4    = new TF1("PileupLine4", "[0] + [1]*x", 1.0, 2.0);
    
    TH1D *testPileupFunc = new TH1D("test", "test", 100, -2, 2);
    
    //testPileupFunc->SetMaximum(1);
    //testPileupFunc->SetMinimum(-1);
    
    pileupLine1->FixParameter(0, -2);
    pileupLine1->FixParameter(1, -1);
    pileupLine2->FixParameter(0, 1);
    pileupLine2->FixParameter(1, 2);
    pileupLine3->FixParameter(0, 1);
    pileupLine3->FixParameter(1, -2);
    pileupLine4->FixParameter(0, -2);
    pileupLine4->FixParameter(1, 1);
    
    str1 = path + outputFolders[5] + outputFolders[8] + fileType;
    testPileupFunc->Draw();
    pileupLine1->Draw("SAME");
    pileupLine2->Draw("SAME");
    pileupLine3->Draw("SAME");
    pileupLine4->Draw("SAME");
    
    c->SaveAs(str1); 
 
    TH2D* fullySubtractedCorrCustomCentralitySideBandFit[4][3];
    TH2D* fullySubtractedCorrCustomCentralitySideBandRes[4][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandFitPhiProj[4][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandResPhiProj[4][3];
    
    TH2D* fullySubtractedCorrCustomCentralitySideBandPileupCorrected[4][3];
    //TH2D* fullySubtractedCorrCustomCentralitySideBandRes[4][3];
    //TH1D* fullySubtractedCorrCustomCentralitySideBandFitPhiProj[4][3];
    //TH1D* fullySubtractedCorrCustomCentralitySideBandResPhiProj[4][3];
    
    TString fitLabel = "_Fit";
    TString resLabel = "_Residual";
    TString pileupLabel = "PileupCorrected";
    
    ofstream fitParameterFile2D;
    TString fitParameterFileName2D = "Extracted_Parameters_from_fitting_2D.txt";
    TString fitParamterOutputFile2D = path + fitParameterFileName2D;
    fitParameterFile2D.open (fitParamterOutputFile2D);
    
    ofstream fitParameterFile1D;
    TString fitParameterFileName1D = "Extracted_Parameters_from_fitting_1D.txt";
    TString fitParamterOutputFile1D = path + fitParameterFileName1D;
    fitParameterFile1D.open (fitParamterOutputFile1D);
    
    double v2_2D;
    double v2_2D_Error;
    double v2_1D;
    double v2_1D_Error;
    
    int bin1;
    int bin2;
    int bin3;
    int bin4;
    
    double GaussIntegralParameter;
 
	for(int k = 0; k < 4; k++){
        for(int i = 0; i < 3; i++){ 
 
            str1 = fitFunctionLabel + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            myfit[k][i] = new TF2(str1,"[0] + [1]*cos(y) + [2]*2*cos(2*y) + [3]*(exp(-0.5*(y*y)/([4]*[4]))+exp(-0.5*((y-6.28)*(y-6.28))/([4]*[4])))*exp(-0.5*(x*x)/([5]*[5]))", -ETA_RANGE, ETA_RANGE,-1.57,4.71);
            //myfit[k][i] = new TF2(str1,"[0] + [1]*cos(y) + [2]*2*cos(2*y) + (1/(4*[4]))*(TMath::Power(1/TMath::CosH((x-[3])/(2*[4])), 2))*(exp(-0.5*(y*y)/([5]*[5]))+exp(-0.5*((y-6.28)*(y-6.28))/([5]*[5])))", -ETA_RANGE, ETA_RANGE,-1.57,4.71);
            
            //(exp(-0.5*(y*y)/([4]*[4]))+exp(-0.5*((y-6.28)*(y-6.28))/([4]*[4])))*exp(-0.5*(x*x)/([5]*[5]))", -ETA_RANGE, ETA_RANGE,-1.57,4.71);
		}
	}
    
    /*for(int k = 2; k < 3; k++){
        for(int i = 0; i < 3; i++){ 
 
            str1 = fitFunctionLabel + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            myfit[k][i] = new TF2(str1,"[0] + [1]*cos(y) + [2]*2*cos(2*y) + [3]*(exp(-0.5*(y*y)/([4]*[4]))+exp(-0.5*((y-6.28)*(y-6.28))/([4]*[4])))*(1/TMath::Beta(.5,.5*[5]))*TMath::Sqrt([6]/[5])*TMath::Power((1.+[6]*(x-[7])*(x-[7]))/[5], -0.5-0.5*[5])", -1.0, 1.0,-1.57,4.71);
		}
	}
    
    for(int k = 3; k < 4; k++){
        for(int i = 0; i < 3; i++){ 
 
            str1 = fitFunctionLabel + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            myfit[k][i] =  new TF2(str1,"[0] + [1]*cos(y) + [2]*2*cos(2*y) + [3]*(exp(-0.5*(y*y)/([4]*[4]))+exp(-0.5*((y-6.28)*(y-6.28))/([4]*[4])))*exp(-0.5*(x*x)/([5]*[5]))", -ETA_RANGE, ETA_RANGE,-1.57,4.71);
		}
	}*/
	
	myfit[0][0]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[0][0]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[0][0]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[0][0]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[0][0]->SetParLimits(4, 0.1, 2);     //phi width
    myfit[0][0]->SetParLimits(5, 0.1, 2);    //eta width
	
	myfit[0][1]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[0][1]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[0][1]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[0][1]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[0][1]->SetParLimits(4, 0.1, 2);     //phi width
    myfit[0][1]->SetParLimits(5, 0.1, 2);    //eta width
	
	myfit[0][2]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[0][2]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[0][2]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[0][2]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[0][2]->SetParLimits(4, 0.1, 2);     //phi width
    myfit[0][2]->SetParLimits(5, 0.1, 2);    //eta width
	
	myfit[1][0]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[1][0]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[1][0]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[1][0]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[1][0]->SetParLimits(4, 0.1, 2);     //phi width
    myfit[1][0]->SetParLimits(5, 0.1, 2);    //eta width
	
	myfit[1][1]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[1][1]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[1][1]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[1][1]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[1][1]->SetParLimits(4, 0.1, 2);     //phi width
    myfit[1][1]->SetParLimits(5, 0.1, 2);    //eta width
	
	myfit[1][2]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[1][2]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[1][2]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[1][2]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[1][2]->SetParLimits(4, 0.1, 2);     //phi width
    myfit[1][2]->SetParLimits(5, 0.1, 2);    //eta width
    
    myfit[2][0]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[2][0]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[2][0]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[2][0]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[2][0]->SetParLimits(4, 0.1, 2);     //phi width
    myfit[2][0]->SetParLimits(5, 0.1, 2);    //eta width
	
	myfit[2][1]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[2][1]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[2][1]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[2][1]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[2][1]->SetParLimits(4, 0.1, 2);     //phi width
    myfit[2][1]->SetParLimits(5, 0.1, 2);    //eta width
	
	myfit[2][2]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[2][2]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[2][2]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[2][2]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[2][2]->SetParLimits(4, 0.1, 2);     //phi width
    myfit[2][2]->SetParLimits(5, 0.1, 2);    //eta width
    
    myfit[3][0]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[3][0]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[3][0]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[3][0]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[3][0]->SetParLimits(4, 0.1, 2);     //phi width
    myfit[3][0]->SetParLimits(5, 0.1, 2);    //eta width
	
	myfit[3][1]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[3][1]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[3][1]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[3][1]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[3][1]->SetParLimits(4, 0.1, 2);     //phi width
    myfit[3][1]->SetParLimits(5, 0.1, 2);    //eta width
	
	myfit[3][2]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[3][2]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[3][2]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[3][2]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[3][2]->SetParLimits(4, 0.1, 2);     //phi width
    myfit[3][2]->SetParLimits(5, 0.1, 2);    //eta width
 
    for(int k = 0; k < 4; k++){
        for(int i = 0; i < 3; i++){ 
 
            //str1 = fitFunctionLabel + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            //myfit[k][i] = new TF2(str1,"[0] + [1]*cos(y) + [2]*2*cos(2*y) + [3]*(exp(-0.5*(y*y)/([4]*[4]))+exp(-0.5*((y-6.28)*(y-6.28))/([4]*[4])))*exp(-0.5*(x*x)/([5]*[5]))", -ETA_RANGE, ETA_RANGE,-1.57,4.71);
           
            str1 = fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fitLabel;
            fullySubtractedCorrCustomCentralitySideBandFit[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            str1 = fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + resLabel;
            fullySubtractedCorrCustomCentralitySideBandRes[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            
            str1 = phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fitLabel;
            fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            str1 = phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + resLabel;
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
 
            str1 = fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + pileupLabel;
            fullySubtractedCorrCustomCentralitySideBandPileupCorrected[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
 
            fitParameterFile2D << "_________________________________________________________________________________________________________" << endl;
            fitParameterFile1D << "_________________________________________________________________________________________________________" << endl;
 
           
            fullySubtractedCorrCustomCentralitySideBandPileupCorrected[k][i] = (TH2D*)fullySubtractedCorrCustomCentralitySideBand[k][i]->Clone();

              
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);

            myfit1D->SetParLimits(0, -.5, 0.0); //offset A0
            //myfit1D->SetParLimits(1, 0.0, 0.0); //v1     AD
            myfit1D->FixParameter(1, 0.0); //v1     AD
            myfit1D->SetParLimits(2, 0.0, .1); //v2     AQ
            
            
            
            //////////////////////2d fitting///////////////////////////////////////
            
            fullySubtractedCorrCustomCentralitySideBand[k][i]->Fit(myfit[k][i], "BR0E");
            //fullySubtractedCorrCustomCentralitySideBandFit[k][i] = (TH2D*)fullySubtractedCorrCustomCentralitySideBand[k][i]->Clone("hfit");
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->Eval(myfit[k][i]);
            
            //GaussIntegralParameter = calcScaleFactorforFit(myfit[k][i]->GetParameter(5));
            
            //cout << "Pt bin: " << k << "   Cent Bin: " << i << "   ProjectionNumber: " << GaussIntegralParameter << endl;
            
			
			fitParameterFile2D << endl << "pt bin (3 is integrated): " << k << "    Centrality Bin (0 - per, 1 - mid cent, 2 - cent): " << i << endl;
			fitParameterFile2D << endl << "A0" << "\t\t" << "AD" << "\t\t" <<"AQ" << "\t\t" <<"A1" << "\t\t" << "phi width" << "\t\t" <<"eta width"  << endl;
			
			fitParameterFile2D << endl << myfit[k][i]->GetParameter(0) << "\t" << myfit[k][i]->GetParameter(1) << "\t" 
                                       << myfit[k][i]->GetParameter(2) << "\t" << myfit[k][i]->GetParameter(3) << "\t" 
                                       << myfit[k][i]->GetParameter(4) << "\t" << myfit[k][i]->GetParameter(5) << endl; 
            
						fitParameterFile2D << myfit[k][i]->GetParError(0) << "\t" << myfit[k][i]->GetParError(1) << "\t" 
                               << myfit[k][i]->GetParError(2) << "\t" << myfit[k][i]->GetParError(3) << "\t" 
                               << myfit[k][i]->GetParError(4) << "\t" << myfit[k][i]->GetParError(5) << endl;
			
            
            fullySubtractedCorrCustomCentralitySideBandRes[k][i] = (TH2D*)fullySubtractedCorrCustomCentralitySideBand[k][i]->Clone("hres");
            //resHistUS[i]->Add(fitHistUS[i], angCorrUS[i], 1, -1);
            str1 = fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + resLabel;
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->Add(fullySubtractedCorrCustomCentralitySideBandFit[k][i], -1);
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->SetTitle(str1);
            
            
            str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fitLabel + fileType;
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->Draw("SURF1");
            c->SaveAs(str1);
 
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->Write();
            
            str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + resLabel + fileType;
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->Draw("SURF1");
            c->SaveAs(str1); 
            
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->Write();
            
            /////////////////////////////1d fitting///////////////////////////////////////////
    
            
    
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Fit(myfit1D, "qR0E");
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i] = (TH1D*)fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Clone("hfit");
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->Eval(myfit1D);
 
            fitParameterFile1D << endl << "pt bin (3 is integrated): " << k << "    Centrality Bin (0 - per, 1 - mid cent, 2 - cent): " << i << endl;
            fitParameterFile1D << endl << "A0" << "\t\t" << "AD" << "\t\t" <<"AQ" << endl; // "\t\t" << "Jet Amp A1" << "\t\t" << "SigPhi" << endl;
            
            fitParameterFile1D << endl << myfit1D->GetParameter(0) << "\t" << myfit1D->GetParameter(1) << "\t" << myfit1D->GetParameter(2) << endl; //<< "\t" << myfit1D->GetParameter(3) << "\t" << myfit1D->GetParameter(5) << endl; 
            fitParameterFile1D << myfit1D->GetParError(0) << "\t" << myfit1D->GetParError(1) << "\t" << myfit1D->GetParError(2) << endl;//<< "\t" << myfit1D->GetParError(3) << "\t" << myfit1D->GetParError(5) << endl; 
            
            v2_2D = TMath::Sqrt(myfit[k][i]->GetParameter(2));
            v2_2D_Error = v2_2D*(.5)*(myfit[k][i]->GetParError(2)/myfit[k][i]->GetParameter(2));
            
            v2Array[k][i] = v2_2D;
            v2ArrayStatError[k][i] = v2_2D_Error;
            
            
            v2_1D = TMath::Sqrt(myfit1D->GetParameter(2));
            v2_1D_Error = v2_2D*(.5)*(myfit1D->GetParError(2)/myfit1D->GetParameter(2));
            
            fitParameterFile2D << "V2 = " << v2_2D << " +/- " << v2_2D_Error << endl;
            fitParameterFile1D << "V2 = " << v2_1D << " +/- " << v2_1D_Error << endl;
            
            /*fitProjection->SetParameter(0, myfit[k][i]->GetParameter(0));
            fitProjection->SetParameter(1, myfit[k][i]->GetParameter(1));
            fitProjection->SetParameter(2, myfit[k][i]->GetParameter(2));
            fitProjection->SetParameter(3, myfit[k][i]->GetParameter(3));
            fitProjection->SetParameter(5, myfit[k][i]->GetParameter(4));
            fitProjection->SetParameter(4, GaussIntegralParameter);*/
            
            quadrupole1->SetParameter(0, myfit[k][i]->GetParameter(2));
            quadrupole2->SetParameter(0, myfit1D->GetParameter(2));
            //quadrupole3->SetParameter(0, fitProjection->GetParameter(2));
            
            if(k == 1 && i == 1) { quadrupole3->SetParameter(0, myfit1D->GetParameter(0)); quadrupole3->SetParameter(1, myfit1D->GetParameter(2)); }
            
            //fitProjection->SetLineColor(2); //red
            //quadrupole1->SetLineColor(3); // Green, from 2D fit
            //quadrupole2->SetLineColor(4); //blue, from 1D fit
            myfit1D->SetLineColor(2); // light blue
            //quadrupole3->SetLineColor(41);
            
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i] = (TH1D*)fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Clone("hres");
            //resHistUS[i]->Add(fitHistUS[i], angCorrUS[i], 1, -1);
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->Add(myfit1D, -1);
            
            
            str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fitLabel + fileType;
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->Draw();
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Draw("EX0");
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->SetMarkerStyle(20);
            
            //fitProjection->Draw("SAME");
            //quadrupole1->Draw("SAME");
            //quadrupole2->Draw("SAME");
            //quadrupole3->Draw("SAME");
            myfit1D->Draw("SAME");
            c->SaveAs(str1);
 
            str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + resLabel + fileType;
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->Draw();
            c->SaveAs(str1); 
 
 
        }
    }    
 
    fitParameterFile2D.close();
    fitParameterFile1D.close();
 
    cout << endl << "Fit Parameter Files closed" << endl << endl;
    
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------Histogram formatting and output for talks----------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
 
    TCanvas *talkCanvas = new TCanvas("c3", "Histograms", 3900, 1700);
    talkCanvas->Divide(3,2);
    
    int padNumber = 1;
    TString talkOutput = "Corrected_Correlations";
    TString ptRanges[2] = {" 1 < p_{T} < 4 GeV/c ", " 4 < p_{T} < 20 GeV/c "};
    TString centRanges[3] = {" 50-100% ", " 20-50% ", " 0-20% "};
    
    double xValueforPave[3] = {1000, 2100, 3200};
    double yValueforPave[2] = {300, 1150};
    
    TString firstLine;
    TString secondLine;
    
    TString starPrelim = "STAR Preliminary";
    TString auau200Label = "Au+Au @ 200 GeV";
    
    //PaveText *textBox = new TPaveText(
    
    for(int k = 1; k < 3; k++){
            for(int i = 0; i < 3; i++){
 
                
                formatCorrHist(fullySubtractedCorrCustomCentralitySideBand[k][i]);
                
                talkCanvas->cd(padNumber);
                talkCanvas->SetRightMargin(0.1);
                talkCanvas->SetLeftMargin(.3);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->SetStats(false);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->SetTitle();
                //fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetTitleOffset(
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Draw("SURF1");
                
                TPaveText *textBox = new TPaveText(0.7, 0.7, 1.0, 1.0, "NB NDC");
                
                firstLine = ptRanges[k-1];
                secondLine = centRanges[i];
                
                textBox->AddText(firstLine.Data());
                textBox->AddText(secondLine.Data());
                textBox->GetLine(0)->SetTextSize(.055);
                textBox->GetLine(1)->SetTextSize(.055);
                textBox->Draw();
                
                TPaveText *starPrelimTextBox = new TPaveText(0.2, 0.8, .5, .88, "NB NDC");
                starPrelimTextBox->SetFillColor(0);
                starPrelimTextBox->AddText(starPrelim);
                starPrelimTextBox->GetLine(0)->SetTextSize(.055);
                starPrelimTextBox->GetLine(0)->SetTextColor(2);
                starPrelimTextBox->Draw("SAME");
                
                TPaveText *auau200TextBox = new TPaveText(0.2, 0.92, .5, 1.0, "NB NDC");
                auau200TextBox->SetFillColor(0);
                auau200TextBox->AddText(auau200Label);
                auau200TextBox->GetLine(0)->SetTextSize(.055);
                auau200TextBox->Draw("SAME");
 
                padNumber++;
        }    
    }
 
    str1 = path + talkOutput + fileTypeEps;
    talkCanvas->SaveAs(str1);

    //////////////////////////////////////////////////////////////////////////////
	
	TCanvas *talkCanvas2 = new TCanvas("c3", "Histograms", 1200, 1700);
    talkCanvas2->Divide(1,2);
    
    int padNumber = 1;
    TString talkOutput2 = "Invariant_Mass_Plots";
   
                
				firstLine = ptRanges[0];
                secondLine = centRanges[1];
	             
                talkCanvas2->cd(1);
                oneDUSHistos[1][1]->SetStats(false);
				oneDSubtractedInvMassHistosFunction[1][1]->SetStats(false);
				
				oneDUSHistos[1][1]->SetTitle();
				oneDSubtractedInvMassHistosFunction[1][1]->SetTitle();
                oneDSubtractedInvMassHistosFunction[1][1]->GetXaxis()->SetTitleSize(.04);
                oneDSubtractedInvMassHistosFunction[1][1]->GetYaxis()->SetTitleSize(.04);
                oneDSubtractedInvMassHistosFunction[1][1]->GetXaxis()->SetTitleFont(62);
                oneDSubtractedInvMassHistosFunction[1][1]->GetYaxis()->SetTitleFont(62);
                oneDSubtractedInvMassHistosFunction[1][1]->SetLabelSize(.04, "X");
                oneDSubtractedInvMassHistosFunction[1][1]->SetLabelSize(.04, "Y");
                oneDSubtractedInvMassHistosFunction[1][1]->SetLabelFont(62, "X");
                oneDSubtractedInvMassHistosFunction[1][1]->SetLabelFont(62, "Y");
                
                oneDUSHistos[1][1]->GetXaxis()->SetTitleSize(.04);
                oneDUSHistos[1][1]->GetYaxis()->SetTitleSize(.04);
                oneDUSHistos[1][1]->GetXaxis()->SetTitleFont(62);
                oneDUSHistos[1][1]->GetYaxis()->SetTitleFont(62);
                oneDUSHistos[1][1]->SetLabelSize(.04, "X");
                oneDUSHistos[1][1]->SetLabelSize(.04, "Y");
                oneDUSHistos[1][1]->SetLabelFont(62, "X");
                oneDUSHistos[1][1]->SetLabelFont(62, "Y");
                
                oneDUSHistos[1][1]->Draw();
                
                TPaveText *textBox2 = new TPaveText(0.8, 0.8, 1.0, 1.0, "NB NDC");
				textBox2->AddText(firstLine.Data());
                textBox2->AddText(secondLine.Data());
                textBox2->Draw();
                
                talkCanvas2->cd(2);
                
				oneDSubtractedInvMassHistosFunction[1][1]->Draw();
				
                TPaveText *textBox3 = new TPaveText(0.8, 0.8, 1.0, 1.0, "NB NDC");
				textBox3->AddText(firstLine.Data());
                textBox3->AddText(secondLine.Data());
                textBox3->Draw();
 
 
				str1 = path + talkOutput2 + fileType;
				talkCanvas2->SaveAs(str1);

				//////////////////////////////////////////////////////////////////////////////
                
				TString talkOutput3 = "Data_Fit_Residual";
				
				TCanvas *talkCanvas3 = new TCanvas("c3", "Histograms", 3700, 850);
                talkCanvas3->Divide(3,1, .01, .01);

                

                formatCorrHist(fullySubtractedCorrCustomCentralitySideBandFit[1][1]);
				formatCorrHist(fullySubtractedCorrCustomCentralitySideBandRes[1][1]);
                
                talkCanvas3->cd(1);
                fullySubtractedCorrCustomCentralitySideBand[1][1]->SetStats(false);
                fullySubtractedCorrCustomCentralitySideBand[1][1]->SetTitle();
                fullySubtractedCorrCustomCentralitySideBand[1][1]->Draw("SURF1");
                
                TPaveText *textBox4 = new TPaveText(0.7, 0.7, 1.0, 1.0, "NB NDC");
                
                firstLine = ptRanges[0];
                secondLine = centRanges[1];
                
                textBox4->AddText(firstLine.Data());
                textBox4->AddText(secondLine.Data());
                textBox4->GetLine(0)->SetTextSize(.055);
                textBox4->GetLine(1)->SetTextSize(.055);
                textBox4->Draw();
                
                TPaveText *starPrelimTextBox = new TPaveText(0.2, 0.8, .5, .88, "NB NDC");
                starPrelimTextBox->SetFillColor(0);
                starPrelimTextBox->AddText(starPrelim);
                starPrelimTextBox->GetLine(0)->SetTextSize(.055);
                starPrelimTextBox->GetLine(0)->SetTextColor(2);
                starPrelimTextBox->Draw("SAME");
                
                TPaveText *auau200TextBox = new TPaveText(0.2, 0.9, .5, .98, "NB NDC");
                auau200TextBox->SetFillColor(0);
                auau200TextBox->AddText(auau200Label);
                auau200TextBox->GetLine(0)->SetTextSize(.04);
                auau200TextBox->Draw("SAME");

                talkCanvas3->cd(2);
                fullySubtractedCorrCustomCentralitySideBandFit[1][1]->SetStats(false);
                fullySubtractedCorrCustomCentralitySideBandFit[1][1]->SetTitle();
                fullySubtractedCorrCustomCentralitySideBandFit[1][1]->Draw("SURF1");
                
                TPaveText *textBox5 = new TPaveText(0.7, 0.7, 1.0, 1.0, "NB NDC");
                
                firstLine = ptRanges[0];
                secondLine = centRanges[1];
                
                textBox5->AddText(firstLine.Data());
                textBox5->AddText(secondLine.Data());
                textBox5->GetLine(0)->SetTextSize(.055);
                textBox5->GetLine(1)->SetTextSize(.055);
                textBox5->Draw();
                
                TPaveText *starPrelimTextBox = new TPaveText(0.2, 0.8, .5, .88, "NB NDC");
                starPrelimTextBox->SetFillColor(0);
                starPrelimTextBox->AddText(starPrelim);
                starPrelimTextBox->GetLine(0)->SetTextSize(.055);
                starPrelimTextBox->GetLine(0)->SetTextColor(2);
                starPrelimTextBox->Draw("SAME");

                TPaveText *auau200TextBox = new TPaveText(0.2, 0.9, .5, .98, "NB NDC");
                auau200TextBox->SetFillColor(0);
                auau200TextBox->AddText(auau200Label);
                auau200TextBox->GetLine(0)->SetTextSize(.04);
                auau200TextBox->Draw("SAME");
                
                talkCanvas3->cd(3);
                fullySubtractedCorrCustomCentralitySideBandRes[1][1]->SetStats(false);
                fullySubtractedCorrCustomCentralitySideBandRes[1][1]->SetTitle();
				fullySubtractedCorrCustomCentralitySideBandRes[1][1]->SetAxisRange(-.02,.05,"Z");
				//fullySubtractedCorrCustomCentralitySideBandRes[1][1]->SetMinimum(-.02)
                fullySubtractedCorrCustomCentralitySideBandRes[1][1]->Draw("SURF1");
                
                TPaveText *textBox6 = new TPaveText(0.7, 0.7, 1.0, 1.0, "NB NDC");
                
                firstLine = ptRanges[0];
                secondLine = centRanges[1];
                
                textBox6->AddText(firstLine.Data());
                textBox6->AddText(secondLine.Data());
                textBox6->GetLine(0)->SetTextSize(.055);
                textBox6->GetLine(1)->SetTextSize(.055);
                textBox6->Draw();
                
                TPaveText *starPrelimTextBox = new TPaveText(0.2, 0.8, .5, .88, "NB NDC");
                starPrelimTextBox->SetFillColor(0);
                starPrelimTextBox->AddText(starPrelim);
                starPrelimTextBox->GetLine(0)->SetTextSize(.055);
                starPrelimTextBox->GetLine(0)->SetTextColor(2);
                starPrelimTextBox->Draw("SAME");
                
                TPaveText *auau200TextBox = new TPaveText(0.2, 0.9, .5, .98, "NB NDC");
                auau200TextBox->SetFillColor(0);
                auau200TextBox->AddText(auau200Label);
                auau200TextBox->GetLine(0)->SetTextSize(.04);
                auau200TextBox->Draw("SAME");

                str1 = path + talkOutput3 + fileTypeEps;
				talkCanvas3->SaveAs(str1);

                //////////////////////////////////////////////////////--------------------------------------------------------
				
				TString talkOutput4 = "delta_phi_projections_with_fits";
				
				TCanvas *talkCanvas4 = new TCanvas("c3", "Histograms", 1700, 1150);
                
                talkCanvas4->SetRightMargin(.15);
                talkCanvas4->SetLeftMargin(.15);    
                
                TF1 *offset = new TF1("offset" , "[0]", -1.57, 4.71);
                TF1 *offsetPlusDipole = new TF1("offset+dipole", "[0]+[1]*cos(x)", -1.57, 4.71);
                TF1 *offsetPlusDipolePlusQuad = new TF1("offset+dipole+quad", "[0]+[1]*cos(x)+2*[2]*cos(2*x)", -1.57, 4.71);
                TF1 *offsetPlusQuad = new TF1("offsetquad", "[0]+2*[1]*cos(2*x)", -1.57, 4.71);
                TF1 *dipole = new TF1("dipole", "[0]*cos(x)", -1.57, 4.71);
                //TF1 *offsetPlusDipole = new TF1("offset+dipole", "[0]+[1]*cos(x)", -1.57, 4.71);
                
                phiProjSyst = new TH1D("", "", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                
                
                //USE QUADRUPOLE 3
                //TF1 *gaussian = new TF1("
				
				fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->SetStats(false);
				//fullySubtractedCorrCustomCentralitySideBandPhiProjFit[1][1]->SetStats(false);
                fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetXaxis()->SetTitle("#Delta#phi");
				fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetXaxis()->CenterTitle();
				fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetXaxis()->SetTitleOffset(1.3);
				fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetYaxis()->SetTitle("#frac{#Delta#rho}{#rho_{ref}}");
				fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetYaxis()->SetTitleOffset(1.55);
                fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->SetAxisRange(-.019,.02,"Y");
				fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->SetTitle();
                
                //fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetZaxis()->SetTitleOffset(1.5);
                //fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetZaxis()->SetTitleSize(.04);
				
                double phiValue[12];
                double corrValue[12];
                double ePhi[12] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
                double eCorrSyst[12] = {0.000694801, 0.000293818, 0.000769675, 0.000769675, 0.000293818, 0.000694801,
                                        0.000758256, 0.000476921, 0.000569372, 0.000569372, 0.000476921, 0.000758256};
                
                
                
                
                
                for( int i = 1; i < NUM_PHI_BINS+1; i++){
                
                    phiValue[i-1] = fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetBinCenter(i);
                    corrValue[i-1] = fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetBinContent(i);
                }
                
                
                
                TGraphErrors *PhiDataSyst = new TGraphErrors(12, phiValue, corrValue, ePhi, eCorrSyst);
  
                //PhiDataSyst->SetMarkerColor(16);
                //PhiDataSyst->SetMarkerStyle(21);
               // PhiDataSyst->SetMarkerSize(1.1);
                PhiDataSyst->SetLineColor(1);
                PhiDataSyst->SetLineWidth(2);
                
                fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->SetMarkerStyle(20);
                fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->SetMarkerSize(2);
                fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->SetLineWidth(2);
				fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->Draw("EX0");
                
                for(int i=0;i<NUM_PHI_BINS;i++) {
                    double x1 = phiValue[i]-0.08;
                    double x2 = phiValue[i]+0.08;
                    double y1 = corrValue[i]-eCorrSyst[i];
                    double y2 = corrValue[i]+eCorrSyst[i];

                    TLine *la = new TLine(x1, y1, x1, y1+0.0003);
                    la->SetLineColor(1);
                    la->Draw("same");
                    TLine *lb = new TLine(x2, y1, x2, y1+0.0003);
                    lb->SetLineColor(1);
                    lb->Draw("same");
                    TLine *lc = new TLine(x1, y2, x1, y2-0.0003);
                    lc->SetLineColor(1);
                    lc->Draw("same");
                    TLine *ld = new TLine(x2, y2, x2, y2-0.0003);
                    ld->SetLineColor(1);
                    ld->Draw("same");
                    TLine *le = new TLine(x1, y1, x2, y1);
                    le->SetLineWidth(2);
                    le->SetLineColor(1);
                    le->Draw("same");
                    TLine *lf = new TLine(x1, y2, x2, y2);
                    lf->SetLineWidth(2);
                    lf->SetLineColor(1);
                    lf->Draw("same");
                }
                
                
                //PhiDataSyst->Draw("[] SAME");
            
			    fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->Fit(myfit1D, "qR0E");
				fullySubtractedCorrCustomCentralitySideBand[1][1]->Fit(myfit[1][1], "qR0E");
			
			    GaussIntegralParameter = calcScaleFactorforFit(myfit[1][1]->GetParameter(5));
			
			    fitProjection->SetParameter(0, myfit[1][1]->GetParameter(0));
				fitProjection->SetParameter(1, myfit[1][1]->GetParameter(1));
				fitProjection->SetParameter(2, myfit[1][1]->GetParameter(2));
				fitProjection->SetParameter(3, myfit[1][1]->GetParameter(3));
				fitProjection->SetParameter(5, myfit[1][1]->GetParameter(4));
				fitProjection->SetParameter(4, GaussIntegralParameter);
			
                offset->SetParameter(0, myfit[1][1]->GetParameter(0));
                offsetPlusDipole->SetParameter(0, myfit[1][1]->GetParameter(0));
                offsetPlusDipole->SetParameter(1, myfit[1][1]->GetParameter(1));
                offsetPlusDipolePlusQuad->SetParameter(0, myfit[1][1]->GetParameter(0));
                offsetPlusDipolePlusQuad->SetParameter(1, myfit[1][1]->GetParameter(1));
                offsetPlusDipolePlusQuad->SetParameter(2, myfit[1][1]->GetParameter(2));
                offsetPlusQuad->SetParameter(0, myfit[1][1]->GetParameter(0));
                offsetPlusQuad->SetParameter(1, myfit[1][1]->GetParameter(2));
                dipole->SetParameter(0, myfit[1][1]->GetParameter(1));
                //quadrupole3->SetParameter(0, myfit[1][1]->GetParameter(2));
                //quadrupole3->SetParameter(0, myfit[k][i]->GetParameter(2));
            
                //dipole->SetLineColor(4); //blue
                //offsetPlusDipole->SetLineColor(3); //green
                offsetPlusQuad->SetLineColor(4); //blue
                quadrupole3->SetLineColor(3);  //green
                //quadrupole3->SetLineWidth(2);
                fitProjection->SetLineColor(2); //red
                fitProjection->SetLineWidth(3);
            
                //dipole->Draw("SAME");
                //offsetPlusDipole->Draw("SAME");
                offsetPlusQuad->Draw("SAME");
                //quadrupole3->Draw("SAME");
                fitProjection->Draw("SAME"); 

				TPaveText *textBox7 = new TPaveText(0.75, 0.75, .95, .95, "NB NDC");
                
                firstLine = ptRanges[0];
                secondLine = centRanges[1];
                
                textBox7->AddText(firstLine.Data());
                textBox7->AddText(secondLine.Data());
                textBox7->GetLine(0)->SetTextSize(.041);
                textBox7->GetLine(1)->SetTextSize(.041);
                textBox7->Draw();
                
                TPaveText *starPrelimTextBox = new TPaveText(0.28, 0.15, 0.37, 0.22, "NB NDC");
                starPrelimTextBox->SetFillColor(0);
                starPrelimTextBox->AddText(starPrelim);
                starPrelimTextBox->GetLine(0)->SetTextSize(.04);
                starPrelimTextBox->GetLine(0)->SetTextColor(2);
                starPrelimTextBox->Draw("SAME");
                
                TPaveText *auau200TextBox = new TPaveText(0.2, 0.92, .5, 1.0, "NB NDC");
                auau200TextBox->SetFillColor(0);
                auau200TextBox->AddText(auau200Label);
                auau200TextBox->GetLine(0)->SetTextSize(.04);
                auau200TextBox->Draw("SAME");
                
                TLegend *leg = new TLegend(0.2, 0.8, .5, .88);
                leg->SetFillStyle(0);
                leg->SetLineStyle(4000);
                leg->SetLineColor(10);
                //leg->SetLineWidth(1.0);
                leg->SetBorderSize(0.0);
                leg->SetTextSize(0.045);
                leg->AddEntry(offsetPlusQuad, "Offset+Quadrupole", "l");
                //leg->AddEntry(fitProjection, "Full Fit, projected on #Delta#phi", "l");
                leg->Draw("SAME");

                str1 = path + talkOutput4 + fileTypePDF;
				talkCanvas4->SaveAs(str1);
				

    //------------------------------------------------------v2 plots-----------------------------------------------

   /* TCanvas *v2Canvas = new TCanvas("c3", "Histograms", 1100, 850);
    
    TString v2 = "v2";
    
    TH1D *v2Histogram = new TH1D("V_{2} 20-50%","V_{2} 20-50%", 10, 0, 9);

    v2Histogram->SetBinContent(1, v2Array[0][1]);  
    v2Histogram->SetBinContent(3, v2Array[1][1]);
    v2Histogram->SetBinContent(5, v2Array[2][1]);

    v2Histogram->SetBinError(1, v2ArrayStatError[0][1]);  
    v2Histogram->SetBinError(3, v2ArrayStatError[1][1]);
    v2Histogram->SetBinError(5, v2ArrayStatError[2][1]);    
        
    str1 = path + v2 + fileType;     
    v2Histogram->Draw("E");     
    v2Canvas->SaveAs(str1);    */
    
    
   /* TString christinaCheck = "check";
    TString christinaCheck1 = "check1";
    TString christinaCheck2 = "check2";
    
    
    double mixedIntegralSBL;
    double mixedIntegralSBR;
    double mixedIntegralAVG;
    double siblingIntegral;
    
    double mixedIntegralUS;
    
    TCanvas *mixCheck = new TCanvas("c3", "Histograms", 3700, 850);
    mixCheck->Divide(3,1);
    
    TH2D *mixingCheck = new TH2D("mixing check", "mixing check", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D *mixingCheckSBL = new TH2D("SBL", "SBL", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D *mixingCheckSBR = new TH2D("SBR", "SBR", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D *mixingCheckAVG = new TH2D("avg1", "avg1", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D *SBAVG = new TH2D("avg", "avg", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    //         band    pt      vz      multiplcity
    //         left    int     5          10
    
    //siblingIntegral = sibCorrBin[0][3][5][10]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
    mixedIntegralSBL = mixCorrBin[0][3][5][10]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
    mixedIntegralSBR = mixCorrBin[2][3][5][10]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
    mixedIntegralUS  = mixCorrBin[1][3][5][10]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
    
    SBAVG->Add(mixCorrBin[0][3][5][10], mixCorrBin[2][3][5][10], .5, .5);
    
    mixedIntegralAVG = SBAVG->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
                
    //sibCorrBin[0][3][5][10]->Scale(1/siblingIntegral);     
    mixCorrBin[0][3][5][10]->Scale(1/mixedIntegralSBL);
    mixCorrBin[2][3][5][10]->Scale(1/mixedIntegralSBR);
    SBAVG->Scale(1/mixedIntegralAVG);
    mixCorrBin[1][3][5][10]->Scale(1/mixedIntegralUS);
    
   // mixingCheck->Divide(sibCorrBin[0][3][5][10], mixCorrBin[0][3][5][10], 1, 1);
    
   // mixingCheck->Draw("SURF1");
    
    mixingCheckSBL->Divide(mixCorrBin[0][3][5][10], mixCorrBin[1][3][5][10], 1, 1);
    mixingCheckSBR->Divide(mixCorrBin[2][3][5][10], mixCorrBin[1][3][5][10], 1, 1);
    mixingCheckAVG->Divide(SBAVG, mixCorrBin[1][3][5][10], 1, 1);
    
    mixCheck->cd(1);
    mixingCheckSBL->Draw("SURF1");
    mixCheck->cd(2);
    mixingCheckSBR->Draw("SURF1");
    mixCheck->cd(3);
    mixingCheckAVG->Draw("SURF1");
    
    str1 = path + christinaCheck + fileType;
    mixCheck->SaveAs(str1);   */


   /////////////////////////////////////

    TCanvas *otherHistograms = new TCanvas("c3", "Histograms", 1700, 1150);
    
    TString strTemp1 = "US_Correlation";
    TString strTemp2 = "SBAverage_Correlation";
     
    otherHistograms->SetRightMargin(.15);
    otherHistograms->SetLeftMargin(.15);    
   
    delRhoOverRhoRefBinVzInt[1][1][11]->SetTitle();  
    sideBandAverage[1][11]->SetTitle(); 
    
    delRhoOverRhoRefBinVzInt[1][1][11]->SetStats(false);  
    sideBandAverage[1][11]->SetStats(false); 

    formatCorrHist(delRhoOverRhoRefBinVzInt[1][1][11]);  
    formatCorrHist(sideBandAverage[1][11]);  

    delRhoOverRhoRefBinVzInt[1][1][11]->GetZaxis()->SetTitleOffset(1.55);  
    sideBandAverage[1][11]->GetZaxis()->SetTitleOffset(1.55);     
    
    delRhoOverRhoRefBinVzInt[1][1][11]->GetZaxis()->SetTitleSize(.04);  
    sideBandAverage[1][11]->GetZaxis()->SetTitleSize(.04);     

    delRhoOverRhoRefBinVzInt[1][1][11]->Draw("SURF1"); 
    TPaveText *starPrelimTextBox = new TPaveText(0.2, 0.8, .5, .88, "NB NDC");
    starPrelimTextBox->SetFillColor(0);
    starPrelimTextBox->AddText(starPrelim);
    starPrelimTextBox->GetLine(0)->SetTextSize(.055);
    starPrelimTextBox->GetLine(0)->SetTextColor(2);
    starPrelimTextBox->Draw("SAME");
    
    TPaveText *auau200TextBox = new TPaveText(0.6, 0.75, .9, .95, "NB NDC");
    auau200TextBox->SetFillColor(0);
    auau200TextBox->AddText(auau200Label);
    auau200TextBox->GetLine(0)->SetTextSize(.045);
    auau200TextBox->AddText("1 GeV/c < p_{T} < 4 GeV/c");
    auau200TextBox->GetLine(1)->SetTextSize(.045);
    auau200TextBox->AddText("~20-30%");
    auau200TextBox->GetLine(2)->SetTextSize(.045);
    auau200TextBox->Draw("SAME");
    
    str1 = path + strTemp1 + fileTypeEps;
    otherHistograms->SaveAs(str1);
    
    sideBandAverage[1][11]->Draw("SURF1");
    TPaveText *starPrelimTextBox = new TPaveText(0.2, 0.8, .5, .88, "NB NDC");
    starPrelimTextBox->SetFillColor(0);
    starPrelimTextBox->AddText(starPrelim);
    starPrelimTextBox->GetLine(0)->SetTextSize(.055);
    starPrelimTextBox->GetLine(0)->SetTextColor(2);
    starPrelimTextBox->Draw("SAME");
    
    TPaveText *auau200TextBox = new TPaveText(0.6, 0.75, .9, .95, "NB NDC");
    auau200TextBox->SetFillColor(0);
    auau200TextBox->AddText(auau200Label);
    auau200TextBox->GetLine(0)->SetTextSize(.045);
    auau200TextBox->AddText("1 GeV/c < p_{T} < 4 GeV/c");
    auau200TextBox->GetLine(1)->SetTextSize(.045);
    auau200TextBox->AddText("~20-30%");
    auau200TextBox->GetLine(2)->SetTextSize(.045);
    auau200TextBox->Draw("SAME");
    str1 = path + strTemp2 + fileTypeEps;
    otherHistograms->SaveAs(str1);
    
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    
   /* TCanvas *invMassWithFilledBands = new TCanvas("c3", "Histograms", 1700, 1150);
    
    strTemp1 = "InvMassWithSB";
    
    //xAxisUS = oneDUSHistos[k][i]->GetXaxis();
    
    binMassLowUS  = oneDUSHistos[1][1]->GetXaxis()->FindBin(sideBandLow);            
    binMassHighUS = oneDUSHistos[1][1]->GetXaxis()->FindBin(sideBandHigh);
    
    TH1D *fillHisto = new TH1D("", "", 50, 1.6, 2.1);
    TH1D *fillHistoSignal = new TH1D("", "", 50, 1.6, 2.1);
    
    for(int bin = binMassLowUS; bin < binMassHighUS+1; bin++){
    
        fillHisto->SetBinContent(bin,oneDUSHistos[1][1]->GetBinContent(bin));
    }
 
    binMassLowUS  = oneDUSHistos[1][1]->GetXaxis()->FindBin(1.62);            
    binMassHighUS = oneDUSHistos[1][1]->GetXaxis()->FindBin(1.7);
 
    for(int bin = binMassLowUS ; bin < binMassHighUS+1; bin++){
    
        fillHisto->SetBinContent(bin,oneDUSHistos[1][1]->GetBinContent(bin));
    }
    
    binMassLowUS  = oneDUSHistos[1][1]->GetXaxis()->FindBin(1.82);            
    binMassHighUS = oneDUSHistos[1][1]->GetXaxis()->FindBin(1.9);
 
    for(int bin = binMassLowUS ; bin < binMassHighUS+1; bin++){
    
        fillHistoSignal->SetBinContent(bin,oneDUSHistos[1][1]->GetBinContent(bin));
    }
 
    oneDUSHistos[1][1]->Draw(); 
    fillHisto->SetFillColor(3);
    fillHisto->Draw("SAME");
    fillHistoSignal->SetFillColor(2);
    fillHistoSignal->Draw("SAME");
    
    str1 = path + strTemp1 + fileType;
    
    invMassWithFilledBands->SaveAs(str1);*/
    
    
    file->Close();
    output->Close();
   
}

