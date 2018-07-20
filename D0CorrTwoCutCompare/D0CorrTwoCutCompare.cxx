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

//---------------------------------------------------------------------------------------------
//This code is a utility used to compare two D0-hadron correlation results using two different 
//sets of cuts for systematic comparsion, via the difference of the results and analysis of the
//residual stucture.
//
//Author: Alex Jentsch
//Date of last update: 6-26-2017
//---------------------------------------------------------------------------------------------

using namespace std;


int D0CorrTwoCutCompare(){

    //Input Files to compare
    //TString inputRootFile1Path = "d0HadronCorrMaker_OUTPUT_2_4_GeV_with_eff_corr.root";   ------- GENERAL NOMINAL MAIN DATA GOLD STANDARD RHIC/AGS MEETING PRELIM---------
    TString inputRootFile1Path = "d0HadronCorrMaker_OUTPUT_12_phi_13_eta_new_dStar_range.root";
    TString inputRootFile2Path = "d0HadronCorrMaker_OUTPUT_systematic_test_13.root";
    
    TFile * file1 = new TFile(inputRootFile1Path);
    TFile * file2 = new TFile(inputRootFile2Path);
    
    int NUM_PT_BINS = 5;
    int NUM_CENT_BINS = 3;
    
    //TString corrHistLabel       = "FullSubtractedCorr_SideBand__PtBin_";
    TString corrHistLabel       = "Full_DStar_Correction_FullSubtractedCorr_SideBand__PtBin_";
    TString ptBinLabel[5]       = {"0", "1", "2", "3", "4"};
    TString centBinLabel[3]     = {"Peripheral", "MidCentral", "Central"};
    TString diffLabel           = "Difference_PtBin_";
    TString pngLabel            = ".png";
    
    TString nSigmaLabel = "nSigma_difference_PtBin_";
    
    double phiBinShift = (TMath::Pi()/12.0);
    
    TString str1;
    TString str2;
    
    TH2D* corr1Hist;
    TH2D* corr2Hist;
    TH2D* difference[3][5];
    TH2D* nSigmaPlots[3][5];

    for(int i = 0; i < 3; i++){
        for(int k = 0; k < 5; k++){
        
            difference[i][k] = new TH2D("", "", 13, -2, 2, 12, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            str1 = nSigmaLabel + ptBinLabel[k] + centBinLabel[i];
            //cout << str1 << endl;
            nSigmaPlots[i][k] = new TH2D(str1, str1, 13, -2, 2, 12, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
        }
    }        
   
    //TCanvas * canvas = new TCanvas("","", 800, 800);
    
    TCanvas * finalCorrGridCanvas = new TCanvas("", "", 3300, 2200);
    finalCorrGridCanvas->Divide(3,2);
    
    TString ptBinTitle[5] = {"0 GeV/c < p_{t, D^{0}} < 1 GeV/c", "1 GeV/c < p_{t, D^{0}} < 2 GeV/c", "2 GeV/c < p_{t, D^{0}} < 10 GeV/c", 
                             "3 GeV/c < p_{t, D^{0}} < 5 GeV/c", "5 GeV/c < p_{t, D^{0}} < 10 GeV/c"}; 
    int padCount = 1;
    
    double tmp, nSigmaNum, nSigmaDem, nSigmaFinal, dataNom, dataTest, error;
    
    for(int centBin = 0; centBin < NUM_CENT_BINS; centBin++){
        for(int ptBin = 2; ptBin < 3; ptBin++){
        
    
            finalCorrGridCanvas->cd(padCount);
            
            str1 = corrHistLabel + ptBinLabel[ptBin] + centBinLabel[centBin];
            str2 = diffLabel + ptBinLabel[ptBin] + centBinLabel[centBin];
            
            cout << str1 << endl;
        
            corr1Hist = (TH2D*)file1->Get(str1);
            corr1Hist->SetDirectory(0);
            corr2Hist = (TH2D*)file2->Get(str1);
            corr2Hist->SetDirectory(0);
        
            if(centBin == 0) { difference[centBin][ptBin]->SetTitle(ptBinTitle[ptBin]); }
            else difference[centBin][ptBin]->SetTitle("");
            
            difference[centBin][ptBin]->SetStats(0);
            difference[centBin][ptBin]->GetXaxis()->SetTitle("#Delta#eta");
            difference[centBin][ptBin]->GetXaxis()->CenterTitle();
            difference[centBin][ptBin]->GetXaxis()->SetTitleOffset(1.5);
            difference[centBin][ptBin]->GetXaxis()->SetTitleSize(.065);
            difference[centBin][ptBin]->GetYaxis()->SetTitle("#Delta#phi");
            difference[centBin][ptBin]->GetYaxis()->CenterTitle();
            difference[centBin][ptBin]->GetYaxis()->SetTitleOffset(1.5);
            difference[centBin][ptBin]->GetYaxis()->SetTitleSize(.065);
            //difference->Draw("SURF1");
        
            //corr1Hist
        
            for(int etaBin = 1; etaBin < 14; etaBin++){
                for(int phiBin = 1; phiBin < 13; phiBin++){
                
                   dataNom  = corr1Hist->GetBinContent(etaBin, phiBin);
                   dataTest = corr2Hist->GetBinContent(etaBin, phiBin);
                   error    = corr1Hist->GetBinError(etaBin, phiBin);
                   
                   tmp = dataNom - dataTest;
                   nSigmaNum = TMath::Abs(tmp);
                   
                   nSigmaFinal = (nSigmaNum/error);
                   
                   nSigmaPlots[centBin][ptBin]->SetBinContent(etaBin, phiBin, nSigmaFinal);
                   nSigmaPlots[centBin][ptBin]->SetBinError(etaBin, phiBin, error);
                }
            }
        
            difference[centBin][ptBin]->Add(corr1Hist, corr2Hist, 1, -1);
            //difference->SetNameTitle(str2, str2);
            difference[centBin][ptBin]->GetXaxis()->SetRangeUser(-1.38, 1.38); 
            difference[centBin][ptBin]->SetMaximum(corr1Hist->GetMaximum());
            difference[centBin][ptBin]->Draw("SURF1");
            
            //str2 = str2 + pngLabel;
            //canvas->SaveAs(str2);
            
            finalCorrGridCanvas->cd(padCount+3);
            nSigmaPlots[centBin][ptBin]->SetStats(0);
            nSigmaPlots[centBin][ptBin]->Draw("COLZ");
            
            padCount++;
        
        }
    }
    
    //finalCorrGridCanvas->SaveAs("Diff_between_nominal_systematic_test_13.png");
    finalCorrGridCanvas->SaveAs("Diff_between_nominal_systematic_test_13_WITH_DSTAR_CORRECTION.png");
    
    file1->Close();
    file2->Close();

}