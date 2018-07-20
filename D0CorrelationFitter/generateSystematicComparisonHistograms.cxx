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

int generateSystematicComparisonHistograms(){


    TString path       = "C:/Users/ajentsch/desktop/D0_Analysis_Code_Local_Copies/D0CorrelationFitter/";
    TString file1Label = "fittingOutput_best_nominal.root";
    TString file2Label = "fittingOutput_double_mis_pid.root";
    TString str1;

    str1 = path + file1Label;
    TFile * inputFile1 = new TFile(str1);
    str1 = path + file2Label;
    TFile * inputFile2 = new TFile(str1);


    TString centralityBin[3] = {"Peripheral", "MidCentral", "Central"};
    TString fileType = ".png";
    
    TString offsetLabel = "Offset_A0_";
    TString dipoleLabel = "Dipole_AD_";
    TString quadLabel   = "Quadrupole_AQ_";
    TString jetAmpLabel = "JetAmp_A1_";
    TString sigEtaLabel = "SigmaEta_";
    TString sigPhiLabel = "SigmaPhi_";
    TString comparisonLabel = "comparison_";
	
	
   
	TH1D* fitParameterOffsetPlot1[3];
    TH1D* fitParameterDipolePlot1[3];
    TH1D* fitParameterQuadrupolePlot1[3];
    TH1D* fitParameterJetAmpPlot1[3];
    TH1D* fitParameterSigmaEtaPlot1[3];
    TH1D* fitParameterSigmePhiPlot1[3];
    
    TH1D* fitParameterOffsetPlot2[3];
    TH1D* fitParameterDipolePlot2[3];
    TH1D* fitParameterQuadrupolePlot2[3];
    TH1D* fitParameterJetAmpPlot2[3];
    TH1D* fitParameterSigmaEtaPlot2[3];
    TH1D* fitParameterSigmePhiPlot2[3];
    
    float Lower[6];
    Lower[0] = 0;
    Lower[1] = 1;
    Lower[2] = 2;
    Lower[3] = 3;
    Lower[4] = 5;
    Lower[5] = 10;
    
    for(int i = 0; i < 3; i++){
    
        str1 = offsetLabel + centralityBin[i];
        fitParameterOffsetPlot1[i] = (TH1D*) inputFile1->Get(str1);
        fitParameterOffsetPlot2[i] = (TH1D*) inputFile2->Get(str1);
        fitParameterOffsetPlot1[i]->SetDirectory(0);
        fitParameterOffsetPlot2[i]->SetDirectory(0);
        
        str1 = dipoleLabel + centralityBin[i];
        fitParameterDipolePlot1[i] = (TH1D*) inputFile1->Get(str1);
        fitParameterDipolePlot2[i] = (TH1D*) inputFile2->Get(str1);
        fitParameterDipolePlot1[i]->SetDirectory(0);
        fitParameterDipolePlot2[i]->SetDirectory(0);
        
        str1 = quadLabel + centralityBin[i];
        fitParameterQuadrupolePlot1[i] = (TH1D*) inputFile1->Get(str1);
        fitParameterQuadrupolePlot2[i] = (TH1D*) inputFile2->Get(str1);
        fitParameterQuadrupolePlot1[i]->SetDirectory(0);
        fitParameterQuadrupolePlot2[i]->SetDirectory(0);
        
        str1 = jetAmpLabel + centralityBin[i];
        fitParameterJetAmpPlot1[i] = (TH1D*) inputFile1->Get(str1);
        fitParameterJetAmpPlot2[i] = (TH1D*) inputFile2->Get(str1);
        fitParameterJetAmpPlot1[i]->SetDirectory(0);
        fitParameterJetAmpPlot2[i]->SetDirectory(0);
        
        str1 = sigEtaLabel + centralityBin[i];
        fitParameterSigmaEtaPlot1[i] = (TH1D*) inputFile1->Get(str1);
        fitParameterSigmaEtaPlot2[i] = (TH1D*) inputFile2->Get(str1);
        fitParameterSigmaEtaPlot1[i]->SetDirectory(0);
        fitParameterSigmaEtaPlot2[i]->SetDirectory(0);
        
        str1 = sigPhiLabel + centralityBin[i];
        fitParameterSigmePhiPlot1[i] = (TH1D*) inputFile1->Get(str1);
        fitParameterSigmePhiPlot2[i] = (TH1D*) inputFile2->Get(str1);
        fitParameterSigmePhiPlot1[i]->SetDirectory(0);
        fitParameterSigmePhiPlot2[i]->SetDirectory(0);
        
       
     
    }
    
    TCanvas * canvas = new TCanvas("", "", 800, 800);
    
    for(int i = 0; i < 3; i++){
    
        str1 = path + comparisonLabel + offsetLabel + centralityBin[i] + fileType;
        fitParameterOffsetPlot1[i]->Draw();
        fitParameterOffsetPlot2[i]->SetMarkerColor(3);
        fitParameterOffsetPlot2[i]->SetMarkerStyle(22);
        fitParameterOffsetPlot2[i]->Draw("SAME");
        canvas->SaveAs(str1);
       
        str1 = path + comparisonLabel + dipoleLabel + centralityBin[i] + fileType;
        fitParameterDipolePlot1[i]->Draw();
        fitParameterDipolePlot2[i]->SetMarkerColor(3);
        fitParameterDipolePlot2[i]->SetMarkerStyle(22);
        fitParameterDipolePlot2[i]->Draw("SAME");
        canvas->SaveAs(str1);
        
        str1 = path + comparisonLabel + quadLabel + centralityBin[i] + fileType;
        fitParameterQuadrupolePlot1[i]->Draw();
        fitParameterQuadrupolePlot2[i]->SetMarkerColor(3);
        fitParameterQuadrupolePlot2[i]->SetMarkerStyle(22);
        fitParameterQuadrupolePlot2[i]->Draw("SAME");
        canvas->SaveAs(str1);
        
        str1 = path + comparisonLabel + jetAmpLabel + centralityBin[i] + fileType;
        fitParameterJetAmpPlot1[i]->Draw();
        fitParameterJetAmpPlot2[i]->SetMarkerColor(3);
        fitParameterJetAmpPlot2[i]->SetMarkerStyle(22);
        fitParameterJetAmpPlot2[i]->Draw("SAME");
        canvas->SaveAs(str1);
        
        str1 = path + comparisonLabel + sigEtaLabel + centralityBin[i] + fileType;
        fitParameterSigmaEtaPlot1[i]->Draw();
        fitParameterSigmaEtaPlot2[i]->SetMarkerColor(3);
        fitParameterSigmaEtaPlot2[i]->SetMarkerStyle(22);
        fitParameterSigmaEtaPlot2[i]->Draw("SAME");
        canvas->SaveAs(str1);
        
        str1 = path + comparisonLabel + sigPhiLabel + centralityBin[i] + fileType;
        fitParameterSigmePhiPlot1[i]->Draw();
        fitParameterSigmePhiPlot2[i]->SetMarkerColor(3);
        fitParameterSigmePhiPlot2[i]->SetMarkerStyle(22);
        fitParameterSigmePhiPlot2[i]->Draw("SAME");
        canvas->SaveAs(str1);

    }

    inputFile1->Close();
    inputFile2->Close();
    

}