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

void formatCorrHist(TH2D* hist, TString title);

using namespace std;


int HadronHadronCorrHistMaker(){

    TFile *file = new TFile("C:/Users/ajentsch/Desktop/HadronHadronCorrHistMaker/Hadron_Hadron_Correlations_Root_Files/DA13ADB9E1B3E6AE4FE0E62C8B6F3A2A_All.root");
    
    //file.open("C:/Users/ajentsch/Desktop/3FF0C11E3F4496C22654B1411E336923_All.root",ifstream::in);
    
    TFile *output = new TFile("C:/Users/ajentsch/Desktop/HadronHadronCorrHistMaker/Hadron_Hadron_Correlations_Root_Files/HadronHadronCorrHistMakerOutput.root", "RECREATE");
    
    
    TH2D* sibCorrPlusPlusBin[10][11];    
    TH2D* sibCorrPlusMinusBin[10][11];
    TH2D* sibCorrMinusPlusBin[10][11];
    TH2D* sibCorrMinusMinusBin[10][11];
    
    TH2D* mixCorrPlusPlusBin[10][11];    
    TH2D* mixCorrPlusMinusBin[10][11];
    TH2D* mixCorrMinusPlusBin[10][11];
    TH2D* mixCorrMinusMinusBin[10][11];
    
    TH2D* scaledMixCorrPlusPlusBin[10][11];
    TH2D* scaledMixCorrPlusMinusBin[10][11]; 
    TH2D* scaledMixCorrMinusPlusBin[10][11];
    TH2D* scaledMixCorrMinusMinusBin[10][11];     
    
    TH2D* deltaRhoCorrPlusPlusBin[10][11];
    TH2D* deltaRhoCorrPlusMinusBin[10][11]; 
    TH2D* deltaRhoCorrMinusPlusBin[10][11];
    TH2D* deltaRhoCorrMinusMinusBin[10][11];  
    
    TH2D* deltaRhoOverRhoCorrPlusPlusBin[10][11];
    TH2D* deltaRhoOverRhoCorrPlusMinusBin[10][11]; 
    TH2D* deltaRhoOverRhoCorrMinusPlusBin[10][11];
    TH2D* deltaRhoOverRhoCorrMinusMinusBin[10][11]; 
    
    TH2D* deltaRhoOverRhoCorrPlusPlusBinVzInt[11];
    TH2D* deltaRhoOverRhoCorrPlusMinusBinVzInt[11]; 
    TH2D* deltaRhoOverRhoCorrMinusPlusBinVzInt[11];
    TH2D* deltaRhoOverRhoCorrMinusMinusBinVzInt[11];
    
    TH2D* chargeIndependent[11];
    
    TH1D* chargeIndependentPhiProj[11];
    
    TString path = "C:/Users/ajentsch/Desktop/";
    TString fileType = ".png";
    
    TString binLabelVz[9] = {"0", "1", "2", "3", "4", "5", "6", "7", "8","9"};
    TString binLabelCent[11] = {"0", "1", "2", "3", "4", "5", "6", "7", "8","9", "10"};
    
    TString outputFolders[5] = {"SibCorr/", "MixedCorr/", "ScaledMixedCorr/", "FullCorr/", "SubtractedCorr/"};
    
    TString sibCorrPlusPlusLabel1   = "Sibling_++_correlation_VzBin_";
    TString sibCorrPlusMinusLabel1  = "Sibling_+-_correlation_VzBin_";
    TString sibCorrMinusPlusLabel1  = "Sibling_-+_correlation_VzBin_";
    TString sibCorrMinusMinusLabel1 = "Sibling_--_correlation_VzBin_";
    
    TString mixCorrPlusPlusLabel1   = "Mixed_++_correlation_VzBin_";
    TString mixCorrPlusMinusLabel1  = "Mixed_+-_correlation_VzBin_";
    TString mixCorrMinusPlusLabel1  = "Mixed_-+_correlation_VzBin_";
    TString mixCorrMinusMinusLabel1 = "Mixed_--_correlation_VzBin_";
    
    TString scaledMixCorrPlusPlusLabel1   = "Scaled_Mixed_++_correlation_VzBin_";
    TString scaledMixCorrPlusMinusLabel1  = "Scaled_Mixed_+-_correlation_VzBin_";
    TString scaledMixCorrMinusPlusLabel1  = "Scaled_Mixed_-+_correlation_VzBin_";
    TString scaledMixCorrMinusMinusLabel1 = "Scaled_Mixed_--_correlation_VzBin_";
    
    TString deltaRhoCorrPlusPlusLabel1   = "delta_Rho_++_correlation_VzBin_";
    TString deltaRhoCorrPlusMinusLabel1  = "delta_Rho_+-_correlation_VzBin_";
    TString deltaRhoCorrMinusPlusLabel1  = "delta_Rho_-+_correlation_VzBin_";
    TString deltaRhoCorrMinusMinusLabel1 = "delta_Rho_--_correlation_VzBin_";
    
    TString deltaRhoOverRhoCorrPlusPlusLabel1   = "delta_Rho_Over_Rho_++_correlation_VzBin_";
    TString deltaRhoOverRhoCorrPlusMinusLabel1  = "delta_Rho_Over_Rho_+-_correlation_VzBin_";
    TString deltaRhoOverRhoCorrMinusPlusLabel1  = "delta_Rho_Over_Rho_-+_correlation_VzBin_";
    TString deltaRhoOverRhoCorrMinusMinusLabel1 = "delta_Rho_Over_Rho_--_correlation_VzBin_";
    
    TString deltaRhoOverRhoCorrPlusPlusVzIntLabel1   = "delta_Rho_Over_Rho_++_correlation_Vz_Int";
    TString deltaRhoOverRhoCorrPlusMinusVzIntLabel1  = "delta_Rho_Over_Rho_+-_correlation_Vz_Int";
    TString deltaRhoOverRhoCorrMinusPlusVzIntLabel1  = "delta_Rho_Over_Rho_-+_correlation_Vz_Int";
    TString deltaRhoOverRhoCorrMinusMinusVzIntLabel1 = "delta_Rho_Over_Rho_--_correlation_Vz_Int";
    
    TString chargeIndependentLabel1 = "Charge_Independent_Corr";
    
    TString centBinLabel = "_CentBin_";
    
    TString phiProjLabel = "_PhiProj_";
    
    TString outputFolder = "HadronHadronCorrHistMaker/Output/";
    
    TString str1;
    TString str2;
    TString str3;
    TString str4;
    
    TH2D* hist1;
    TH2D* hist2;
    TH2D* hist3;
    TH2D* hist4;
    
    double numSib = 0;
    double numMix = 0;
    
    double NUM_ETA_BINS = 25;
    double NUM_PHI_BINS = 24;
    
    double weightFactorPlusPlusMixDenominator[11] = {0,0,0,0,0,0,0,0,0,0};
    double weightFactorPlusPlusMixNumerator[10][11];
    double weightFactorPlusMinusMixDenominator[11] = {0,0,0,0,0,0,0,0,0,0};
    double weightFactorPlusMinusMixNumerator[10][11];
    double weightFactorMinusPlusMixDenominator[11] = {0,0,0,0,0,0,0,0,0,0};
    double weightFactorMinusPlusMixNumerator[10][11];
    double weightFactorMinusMinusMixDenominator[11] = {0,0,0,0,0,0,0,0,0,0};
    double weightFactorMinusMinusMixNumerator[10][11];
    
    double scaleFactorPlusPlus = 0;
    double scaleFactorPlusMinus = 0;
    double scaleFactorMinusPlus = 0;
    double scaleFactorMinusMinus = 0;
    
    
    //fullCorrUSBin[i][j]        = new TH2D(str1, str1, 100, -2, 2, 64, -TMath::PiOver2(), 3*TMath::PiOver2());
    //fullCorrLSBin[i][j]        = new TH2D(str2, str2, 25, -2, 2, 25, -TMath::PiOver2(), 3*TMath::PiOver2());
    
    
    for(int i = 0; i < 10; i++){ 
        for(int j = 0; j < 11; j++){
         
            
             
            //sibling histograms
         
            str1 = sibCorrPlusPlusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str2 = sibCorrPlusMinusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str3 = sibCorrMinusPlusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str4 = sibCorrMinusMinusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            
            hist1 = (TH2D*)file->Get(str1);
            hist2 = (TH2D*)file->Get(str2);
            hist3 = (TH2D*)file->Get(str3);
            hist4 = (TH2D*)file->Get(str4);
            
            sibCorrPlusPlusBin[i][j]   = (TH2D*)hist1->Clone();    
            sibCorrPlusMinusBin[i][j]  = (TH2D*)hist2->Clone();
            sibCorrMinusPlusBin[i][j]  = (TH2D*)hist3->Clone();
            sibCorrMinusMinusBin[i][j] = (TH2D*)hist4->Clone();
            
            formatCorrHist(sibCorrPlusPlusBin[i][j], str1);
            formatCorrHist(sibCorrPlusMinusBin[i][j], str2);
            formatCorrHist(sibCorrMinusPlusBin[i][j], str3);
            formatCorrHist(sibCorrMinusMinusBin[i][j], str4);
            
            //mixed histograms
            
            str1 = mixCorrPlusPlusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str2 = mixCorrPlusMinusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str3 = mixCorrMinusPlusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str4 = mixCorrMinusMinusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            
            hist1 = (TH2D*)file->Get(str1);
            hist2 = (TH2D*)file->Get(str2);
            hist3 = (TH2D*)file->Get(str3);
            hist4 = (TH2D*)file->Get(str4);
            
            mixCorrPlusPlusBin[i][j]   = (TH2D*)hist1->Clone();    
            mixCorrPlusMinusBin[i][j]  = (TH2D*)hist2->Clone();
            mixCorrMinusPlusBin[i][j]  = (TH2D*)hist3->Clone();
            mixCorrMinusMinusBin[i][j] = (TH2D*)hist4->Clone();
            
            formatCorrHist(mixCorrPlusPlusBin[i][j], str1);
            formatCorrHist(mixCorrPlusMinusBin[i][j], str2);
            formatCorrHist(mixCorrMinusPlusBin[i][j], str3);
            formatCorrHist(mixCorrMinusMinusBin[i][j], str4);
            
         }
    }    
    
    for(int i = 0; i < 10; i++){ 
        for(int j = 0; j < 11; j++){
         
            //sibling histograms
         
            str1 = scaledMixCorrPlusPlusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str2 = scaledMixCorrPlusMinusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str3 = scaledMixCorrMinusPlusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str4 = scaledMixCorrMinusMinusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            
            scaledMixCorrPlusPlusBin[i][j] = (TH2D*)mixCorrPlusPlusBin[i][j]->Clone();    
            scaledMixCorrPlusMinusBin[i][j]  = (TH2D*)mixCorrPlusMinusBin[i][j]->Clone();
            scaledMixCorrMinusPlusBin[i][j]  = (TH2D*)mixCorrMinusPlusBin[i][j]->Clone();
            scaledMixCorrMinusMinusBin[i][j] = (TH2D*)mixCorrMinusMinusBin[i][j]->Clone();
            
            formatCorrHist(scaledMixCorrPlusPlusBin[i][j], str1);
            formatCorrHist(scaledMixCorrPlusMinusBin[i][j], str2);
            formatCorrHist(scaledMixCorrMinusPlusBin[i][j], str3);
            formatCorrHist(scaledMixCorrMinusMinusBin[i][j], str4);
            
            numSib = sibCorrPlusPlusBin[i][j]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
            numMix = mixCorrPlusPlusBin[i][j]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
            if(numMix == 0) { weightFactorPlusPlusMixNumerator[i][j] = 0;
                              numMix =1;                                }
            else weightFactorPlusPlusMixNumerator[i][j] = numMix;
            scaledMixCorrPlusPlusBin[i][j]->Scale(numSib/numMix);
            
            numSib = sibCorrPlusMinusBin[i][j]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
            numMix = mixCorrPlusMinusBin[i][j]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
            if(numMix == 0) { weightFactorPlusMinusMixNumerator[i][j] = 0;
                              numMix =1;                                 }
            else weightFactorPlusMinusMixNumerator[i][j] = numMix;
            scaledMixCorrPlusMinusBin[i][j]->Scale(numSib/numMix);
            
            numSib = sibCorrMinusPlusBin[i][j]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
            numMix = mixCorrMinusPlusBin[i][j]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
            if(numMix == 0) { weightFactorMinusPlusMixNumerator[i][j] = 0; 
                              numMix =1;                                 }
            else weightFactorMinusPlusMixNumerator[i][j] = numMix;
            scaledMixCorrMinusPlusBin[i][j]->Scale(numSib/numMix);
            
            numSib = sibCorrMinusMinusBin[i][j]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
            numMix = mixCorrMinusMinusBin[i][j]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
            if(numMix == 0) { weightFactorMinusMinusMixNumerator[i][j] = 0; 
                              numMix =1;                                  }
            else weightFactorMinusMinusMixNumerator[i][j] = numMix;
            scaledMixCorrMinusMinusBin[i][j]->Scale(numSib/numMix);
            
            str1 = deltaRhoCorrPlusPlusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str2 = deltaRhoCorrPlusMinusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str3 = deltaRhoCorrMinusPlusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str4 = deltaRhoCorrMinusMinusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            
            deltaRhoCorrPlusPlusBin[i][j]   = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            deltaRhoCorrPlusMinusBin[i][j]  = new TH2D(str2, str2, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2()); 
            deltaRhoCorrMinusPlusBin[i][j]  = new TH2D(str3, str3, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            deltaRhoCorrMinusMinusBin[i][j] = new TH2D(str4, str4, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            
            
            deltaRhoCorrPlusPlusBin[i][j]->Add(sibCorrPlusPlusBin[i][j], scaledMixCorrPlusPlusBin[i][j], 1, -1);
            deltaRhoCorrPlusMinusBin[i][j]->Add(sibCorrPlusMinusBin[i][j], scaledMixCorrPlusMinusBin[i][j], 1, -1); 
            deltaRhoCorrMinusPlusBin[i][j]->Add(sibCorrMinusPlusBin[i][j], scaledMixCorrMinusPlusBin[i][j], 1, -1);
            deltaRhoCorrMinusMinusBin[i][j]->Add(sibCorrMinusMinusBin[i][j], scaledMixCorrMinusMinusBin[i][j], 1, -1);
            
            formatCorrHist(deltaRhoCorrPlusPlusBin[i][j], str1);
            formatCorrHist(deltaRhoCorrPlusMinusBin[i][j], str2);
            formatCorrHist(deltaRhoCorrMinusPlusBin[i][j], str3);
            formatCorrHist(deltaRhoCorrMinusMinusBin[i][j], str4);
            
            
         }
         
    }    

    for(j = 0; j < 11; j++){
        for(i = 0; i < 10; i++){
        
            weightFactorPlusPlusMixDenominator[j] = weightFactorPlusPlusMixNumerator[i][j] + weightFactorPlusPlusMixDenominator[j] ;
            weightFactorPlusMinusMixDenominator[j] = weightFactorPlusMinusMixNumerator[i][j] + weightFactorPlusMinusMixDenominator[j];
            weightFactorMinusPlusMixDenominator[j] = weightFactorMinusPlusMixNumerator[i][j] + weightFactorMinusPlusMixDenominator[j];
            weightFactorMinusMinusMixDenominator[j] = weightFactorMinusMinusMixNumerator[i][j] + weightFactorMinusMinusMixDenominator[j];
        }
    }
    
    
    for(int i = 0; i < 10; i++){ 
        for(int j = 0; j < 11; j++){
    
            str1 = deltaRhoOverRhoCorrPlusPlusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str2 = deltaRhoOverRhoCorrPlusMinusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str3 = deltaRhoOverRhoCorrMinusPlusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            str4 = deltaRhoOverRhoCorrMinusMinusLabel1 + binLabelVz[i] + centBinLabel + binLabelCent[j];
            
            deltaRhoOverRhoCorrPlusPlusBin[i][j]   = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            deltaRhoOverRhoCorrPlusMinusBin[i][j]  = new TH2D(str2, str2, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2()); 
            deltaRhoOverRhoCorrMinusPlusBin[i][j]  = new TH2D(str3, str3, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            deltaRhoOverRhoCorrMinusMinusBin[i][j] = new TH2D(str4, str4, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    
            deltaRhoOverRhoCorrPlusPlusBin[i][j]->Divide(deltaRhoCorrPlusPlusBin[i][j], scaledMixCorrPlusPlusBin[i][j], 1, 1);
            deltaRhoOverRhoCorrPlusMinusBin[i][j]->Divide(deltaRhoCorrPlusMinusBin[i][j], scaledMixCorrPlusMinusBin[i][j], 1, 1);
            deltaRhoOverRhoCorrMinusPlusBin[i][j]->Divide(deltaRhoCorrMinusPlusBin[i][j], scaledMixCorrMinusPlusBin[i][j], 1, 1);
            deltaRhoOverRhoCorrMinusMinusBin[i][j]->Divide(deltaRhoCorrMinusMinusBin[i][j], scaledMixCorrMinusMinusBin[i][j], 1, 1);
        }
    }    
    
    TCanvas *c = new TCanvas("c2", "Histograms", 1100, 850);
    
    for(int j = 0; j < 11; j++){
    
        str1 = deltaRhoOverRhoCorrPlusPlusVzIntLabel1 + centBinLabel + binLabelCent[j];
        str2 = deltaRhoOverRhoCorrPlusMinusVzIntLabel1 + centBinLabel + binLabelCent[j];
        str3 = deltaRhoOverRhoCorrMinusPlusVzIntLabel1 + centBinLabel + binLabelCent[j];
        str4 = deltaRhoOverRhoCorrMinusMinusVzIntLabel1 + centBinLabel + binLabelCent[j];
        
        //cout << "Test 5" << endl;
        
        deltaRhoOverRhoCorrPlusPlusBinVzInt[j]   = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
        deltaRhoOverRhoCorrPlusMinusBinVzInt[j]  = new TH2D(str2, str2, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2()); 
        deltaRhoOverRhoCorrMinusPlusBinVzInt[j]  = new TH2D(str3, str3, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
        deltaRhoOverRhoCorrMinusMinusBinVzInt[j] = new TH2D(str4, str4, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
        
        cout << endl;
        
        for(int i = 0; i < 10; i++){

            //cout << "Test 3" << endl;
            
            scaleFactorPlusPlus = weightFactorPlusPlusMixNumerator[i][j]/weightFactorPlusPlusMixDenominator[j];
            scaleFactorPlusMinus = weightFactorPlusMinusMixNumerator[i][j]/weightFactorPlusMinusMixDenominator[j];
            scaleFactorMinusPlus = weightFactorMinusPlusMixNumerator[i][j]/weightFactorMinusPlusMixDenominator[j];
            scaleFactorMinusMinus = weightFactorMinusMinusMixNumerator[i][j]/weightFactorMinusMinusMixDenominator[j];
            
            
            //cout << scaleFactorPlusPlus << endl;
            
            deltaRhoOverRhoCorrPlusPlusBinVzInt[j]->Add(deltaRhoOverRhoCorrPlusPlusBin[i][j], scaleFactorPlusPlus);
            deltaRhoOverRhoCorrPlusMinusBinVzInt[j]->Add(deltaRhoOverRhoCorrPlusMinusBin[i][j], scaleFactorPlusMinus); 
            deltaRhoOverRhoCorrMinusPlusBinVzInt[j]->Add(deltaRhoOverRhoCorrMinusPlusBin[i][j], scaleFactorMinusPlus);
            deltaRhoOverRhoCorrMinusMinusBinVzInt[j]->Add(deltaRhoOverRhoCorrMinusMinusBin[i][j],scaleFactorMinusMinus);
            
            //deltaRhoOverRhoCorrPlusPlusBinVzInt[j]->Add(deltaRhoOverRhoCorrPlusPlusBin[i][j]);
            //deltaRhoOverRhoCorrPlusMinusBinVzInt[j]->Add(deltaRhoOverRhoCorrPlusMinusBin[i][j]); 
            //deltaRhoOverRhoCorrMinusPlusBinVzInt[j]->Add(deltaRhoOverRhoCorrMinusPlusBin[i][j]);
            //deltaRhoOverRhoCorrMinusMinusBinVzInt[j]->Add(deltaRhoOverRhoCorrMinusMinusBin[i][j]);
            
        }
        
        cout << endl;
        
        str1 = path + outputFolder + deltaRhoOverRhoCorrPlusPlusVzIntLabel1 + centBinLabel + binLabelCent[j] + fileType;
        str2 = path + outputFolder + deltaRhoOverRhoCorrPlusMinusVzIntLabel1 + centBinLabel + binLabelCent[j] + fileType;
        str3 = path + outputFolder + deltaRhoOverRhoCorrMinusPlusVzIntLabel1 + centBinLabel + binLabelCent[j] + fileType;
        str4 = path + outputFolder + deltaRhoOverRhoCorrMinusMinusVzIntLabel1 + centBinLabel + binLabelCent[j] + fileType;
        
        deltaRhoOverRhoCorrPlusPlusBinVzInt[j]->Draw("SURF1");
        //deltaRhoOverRhoCorrPlusPlusBinVzInt[j]->Rebin2D(4,4);
        deltaRhoOverRhoCorrPlusPlusBinVzInt[j]->GetXaxis()->SetRangeUser(-2.0, 2.0); 
        c->SaveAs(str1);
        
        deltaRhoOverRhoCorrPlusMinusBinVzInt[j]->Draw("SURF1");
       // deltaRhoOverRhoCorrPlusMinusBinVzInt[j]->Rebin2D(4,4);
        deltaRhoOverRhoCorrPlusMinusBinVzInt[j]->GetXaxis()->SetRangeUser(-2.0, 2.0); 
        c->SaveAs(str2);
        
        deltaRhoOverRhoCorrMinusPlusBinVzInt[j]->Draw("SURF1");
        //deltaRhoOverRhoCorrMinusPlusBinVzInt[j]->Rebin2D(4,4);
        deltaRhoOverRhoCorrMinusPlusBinVzInt[j]->GetXaxis()->SetRangeUser(-2.0, 2.0); 
        c->SaveAs(str3);
        
        deltaRhoOverRhoCorrMinusMinusBinVzInt[j]->Draw("SURF1");
        //deltaRhoOverRhoCorrMinusMinusBinVzInt[j]->Rebin2D(4,4);
        deltaRhoOverRhoCorrMinusMinusBinVzInt[j]->GetXaxis()->SetRangeUser(-2.0, 2.0); 
        c->SaveAs(str4);
        
    }
        
    for(int i = 0; i < 11; i++){
        
        str1 = chargeIndependentLabel1 + centBinLabel + binLabelCent[i];
        
        chargeIndependent[i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
        
        str1 = chargeIndependentLabel1 + phiProjLabel + centBinLabel + binLabelCent[i];
        
        chargeIndependentPhiProj[i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
        
        str1 = path + outputFolder + chargeIndependentLabel1 + centBinLabel + binLabelCent[i] + fileType;
        
        chargeIndependent[i]->Add(deltaRhoOverRhoCorrPlusPlusBinVzInt[i]);
        chargeIndependent[i]->Add(deltaRhoOverRhoCorrPlusMinusBinVzInt[i]);
        chargeIndependent[i]->Add(deltaRhoOverRhoCorrMinusPlusBinVzInt[i]);
        chargeIndependent[i]->Add(deltaRhoOverRhoCorrMinusMinusBinVzInt[i]);  
        chargeIndependent[i]->Scale(1.0/4.0);
        //chargeIndependent[i]->Rebin2D(4,4);
        chargeIndependent[i]->GetXaxis()->SetRangeUser(-1.8, 1.8);
        chargeIndependent[i]->Draw("SURF1");
        c->SaveAs(str1);        
        
        str1 = path + outputFolder + chargeIndependentLabel1 + phiProjLabel + centBinLabel + binLabelCent[i] + fileType;
        chargeIndependentPhiProj[i] = (TH1D*)chargeIndependent[i]->ProjectionY();
        chargeIndependentPhiProj[i]->Draw();
        c->SaveAs(str1);   
        
    }    
    
    file->Close();
    
    for(int i = 0; i < 10; i++){ 
        for(int j = 0; j < 11; j++){
    
            sibCorrPlusPlusBin[i][j]->Write();    
            sibCorrPlusMinusBin[i][j]->Write();
            sibCorrMinusPlusBin[i][j]->Write();
            sibCorrMinusMinusBin[i][j]->Write();
            mixCorrPlusPlusBin[i][j]->Write();    
            mixCorrPlusMinusBin[i][j]->Write();
            mixCorrMinusPlusBin[i][j]->Write();
            mixCorrMinusMinusBin[i][j]->Write();
            scaledMixCorrPlusPlusBin[i][j]->Write();
            scaledMixCorrPlusMinusBin[i][j]->Write();
            scaledMixCorrMinusPlusBin[i][j]->Write();
            scaledMixCorrMinusMinusBin[i][j]->Write();
            deltaRhoCorrPlusPlusBin[i][j]->Write();
            deltaRhoCorrPlusMinusBin[i][j]->Write(); 
            deltaRhoCorrMinusPlusBin[i][j]->Write();
            deltaRhoCorrMinusMinusBin[i][j]->Write();
            deltaRhoOverRhoCorrPlusPlusBin[i][j]->Write();
            deltaRhoOverRhoCorrPlusMinusBin[i][j]->Write(); 
            deltaRhoOverRhoCorrMinusPlusBin[i][j]->Write();
            deltaRhoOverRhoCorrMinusMinusBin[i][j]->Write();
        }
    }
    
    //cout << "Test 1" << endl;
    
    for(int j = 0; j < 11; j++){    
        
        deltaRhoOverRhoCorrPlusPlusBinVzInt[j]->Write();
        deltaRhoOverRhoCorrPlusMinusBinVzInt[j]->Write(); 
        deltaRhoOverRhoCorrMinusPlusBinVzInt[j]->Write();
        deltaRhoOverRhoCorrMinusMinusBinVzInt[j]->Write();
        
    }    
    
    
    output->Close();
    
}

void formatCorrHist(TH2D* hist, TString title) {

    hist->GetXaxis()->SetTitle("#Delta#eta");
    hist->GetYaxis()->SetTitle("#Delta#phi");
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleOffset(1.3); 
    hist->GetXaxis()->SetTitle("#Delta#eta");
    hist->GetYaxis()->SetTitle("#Delta#phi");
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleOffset(1.3); 
    hist->SetNameTitle(title, title);
    
    return;
    
}

 /*double numSib = angCorr2DAll->Integral(1,25,1,25);
   double numMix = mixedEventAngCorr2DAll->Integral(1,25,1,25);

   angCorr2DAll->Scale(numMix/numSib);

   fullCorr->Add(angCorr2DAll, mixedEventAngCorr2DAll, 1, -1);   
   
   fullCorr->Divide(mixedEventAngCorr2DAll);
   
   
   for(int i = 0; i < 9; i++){
        for( int j = 0; j < 11; j++){
        
            numSib = scaledSibCorrBin[i][j]->Integral(1,25,1,25); 
            numMix = mixCorrBin[i][j]->Integral(1,25,1,25);
           
            scaledSibCorrBin[i][j]->Scale(numMix/numSib);  
             
            fullCorrBin[i][j]->Add(scaledSibCorrBin[i][j], mixCorrBin[i][j], 1, -1);   
   
            fullCorrBin[i][j]->Divide(mixCorrBin[i][j]); 
            
               
            delete eventBuffer[i][j]; 
        
        }
   }    
   
   for(int i = 0; i < 9; i++){
        for(int j = 0; j < 11; j++){
        
            sibCorrBin[i][j]->Write();
            mixCorrBin[i][j]->Write();
            scaledSibCorrBin[i][j]->Write();
            fullCorrBin[i][j]->Write();
            
        }
   }    */
    