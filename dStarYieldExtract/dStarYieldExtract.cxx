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



int dStarYieldExtract(){


    TFile * inputFile = new TFile("systematic_test_13.root");
    
    TString str1;
    TString DStarLabel = "KPiPi_minus_KPi_";
    TString siblingLabel      = "Sibling_";
    TString mixLabel          = "Mixed_";
    TString diffLabel         = "Diff_";
    TString sameLabel         = "SAME_";
    TString angularLabel = "_Angular_Distribution";
    TString DStarLabelAng = "DStarHistograms_";
    TString PtBinLabel   = "_PtBin_";
    TString binLabelPt[6]              = {"0", "1", "2", "3", "4", "5"};
    TString binLabelCentralityClass[3] = {"_Peripheral", "_MidCentral", "_Central"};
    double scaleFactors[5][3];
    
    
    int NUM_PHI_BINS       = 12;                //number of delPhi bins to use for correlations
    int NUM_ETA_BINS       = 13;                 //number of delEta bins to use for correlations
    double phiBinShift        = (TMath::Pi()/12.0);     //This number shifts the phi bins to ensure that 0 and pi are at the center of a bin
    

    TH1D* dStarHistogramsSib[5][3];
    TH1D* dStarHistogramsMix[5][3];
    TH1D* dStarHistogramsDiff[5][3];
    TH2D* dStarAngularDistSib[5][3];
    TH2D* dStarAngularDistMix[5][3];
    TH2D* dStarAngularDistDiff[5][3];
    
    for(int k = 2; k < 3; k++){
        for(int j = 0; j < 3; j++){
    
            str1 = siblingLabel + angularLabel + DStarLabelAng + PtBinLabel + binLabelPt[k] + binLabelCentralityClass[j];
            dStarAngularDistSib[k][j] = (TH2D*) inputFile->Get(str1);
            str1 = mixLabel + angularLabel + DStarLabelAng + PtBinLabel + binLabelPt[k] + binLabelCentralityClass[j];
            dStarAngularDistMix[k][j] = (TH2D*) inputFile->Get(str1);
            str1 = diffLabel + angularLabel + DStarLabelAng + PtBinLabel + binLabelPt[k] + binLabelCentralityClass[j];
            dStarAngularDistDiff[k][j] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, (3*TMath::PiOver2())+phiBinShift);
    
        }
    }
    
    for(int k = 2; k < 3; k++){
        for(int j = 0; j < 3; j++){
            
            
            str1 = siblingLabel + DStarLabel + PtBinLabel + binLabelPt[k] + binLabelCentralityClass[j];
            dStarHistogramsSib[k][j] = (TH1D*) inputFile->Get(str1);
            str1 = mixLabel + DStarLabel + PtBinLabel + binLabelPt[k] + binLabelCentralityClass[j];
            dStarHistogramsMix[k][j] = (TH1D*) inputFile->Get(str1);
            
            dStarHistogramsSib[k][j]->Sumw2();
            dStarHistogramsMix[k][j]->Sumw2();
        }
    }
    
   
    
    TCanvas * can = new TCanvas("", "", 1100, 1000);
    
    double NSib;
    double NMix;
    double ratio;
    int bin1, bin2;
    double S, B, SPlusB, SOverSPlusB;
    
    for(int k = 2; k < 3; k++){
        for(int j = 0; j < 3; j++){
            
            cout << "--------------------------------------------------" << endl;
            
            NSib = dStarHistogramsSib[k][j]->Integral(50, 100);
            NMix = dStarHistogramsMix[k][j]->Integral(50, 100);
            ratio = NSib/NMix;
            //ratio = .2;
            //cout << ratio << endl;
            
            scaleFactors[k][j] = ratio;
            
            str1 = diffLabel + DStarLabel + PtBinLabel + binLabelPt[k] + binLabelCentralityClass[j];
            dStarHistogramsMix[k][j]->Scale(ratio);
            dStarHistogramsDiff[k][j] = new TH1D(str1, str1, 400, 0.135, .3);
            dStarHistogramsDiff[k][j]->Add(dStarHistogramsSib[k][j], dStarHistogramsMix[k][j], 1, -1);
            dStarHistogramsDiff[k][j]->GetXaxis()->SetRangeUser(.13, .165);
            dStarHistogramsDiff[k][j]->SetMarkerStyle(20);
            dStarHistogramsDiff[k][j]->Draw("P E");
            str1 = str1 + ".png";
            can->SaveAs(str1);
            
            bin1 = dStarHistogramsDiff[k][j]->FindBin(.143);
            bin2 = dStarHistogramsDiff[k][j]->FindBin(.147);
            
            B = dStarHistogramsMix[k][j]->Integral(bin1, bin2);
            SPlusB = dStarHistogramsSib[k][j]->Integral(bin1, bin2);

            S = SPlusB - B;
            SOverSPlusB = S/SPlusB;
            
            cout << "PtBin: " << k << "   CentBin: " << j << endl;
            cout << "S + B = " << SPlusB << ",  S = " << S << "  , B = " << B << endl; 
            cout << "S/S+B = " << SOverSPlusB << endl;
            cout << "Integral of D* peak region: " <<  dStarHistogramsDiff[k][j]->Integral(bin1, bin2) << "   " << bin1 << " , " << bin2 << endl;
            cout << endl;
            
            str1 = sameLabel + DStarLabel + PtBinLabel + binLabelPt[k] + binLabelCentralityClass[j];
            
            //dStarHistogramsSib[k][j]->GetXaxis()->SetRangeUser(.13, .2);
            //dStarHistogramsMix[k][j]->GetXaxis()->SetRangeUser(.13, .2);
            
            dStarHistogramsSib[k][j]->GetXaxis()->SetRangeUser(.13, .165);
            dStarHistogramsMix[k][j]->GetXaxis()->SetRangeUser(.13, .165);
            
            TLegend * leg = new TLegend( .5, .3, .8, .45);
            leg->AddEntry(dStarHistogramsSib[k][j], "Same-Event", "P");
            leg->AddEntry(dStarHistogramsMix[k][j], "Mixed-Event", "P");
            leg->AddEntry(dStarHistogramsDiff[k][j], "SE-ME", "P");
            gStyle->SetLegendBorderSize(0);
            
            dStarHistogramsSib[k][j]->SetMarkerStyle(20);
            dStarHistogramsSib[k][j]->SetMarkerColor(2);
            dStarHistogramsSib[k][j]->SetStats(0);
            dStarHistogramsSib[k][j]->SetTitle("");
            dStarHistogramsSib[k][j]->GetXaxis()->SetTitle("Inv. Mass Difference (M_{K#pi#pi}-M_{k#pi}) GeV/c^{2}");
            dStarHistogramsSib[k][j]->GetXaxis()->SetTitleSize(.035);
            dStarHistogramsSib[k][j]->GetXaxis()->SetLabelSize(.03);
            dStarHistogramsSib[k][j]->GetYaxis()->SetLabelSize(.03);
            dStarHistogramsSib[k][j]->Draw("P E");
            dStarHistogramsMix[k][j]->SetMarkerStyle(20);
            dStarHistogramsMix[k][j]->Draw("P SAME");
            dStarHistogramsDiff[k][j]->SetMarkerColor(4);
            dStarHistogramsDiff[k][j]->Draw("P SAME");
            dStarHistogramsDiff[k][j]->SetStats(0);
            leg->Draw("SAME");
            //dStarHistogramsDiff[k][j]->GetXaxis()->SetTitle("Inv. Mass Difference (M_{K#pi#pi}-M_{k#pi}) GeV/c");
           
            str1 = str1 + ".png";
            can->SaveAs(str1);
            
            dStarAngularDistDiff[k][j]->Add(dStarAngularDistSib[k][j], dStarAngularDistMix[k][j], 1, -ratio);
            str1 = diffLabel + angularLabel + DStarLabel + PtBinLabel + binLabelPt[k] + binLabelCentralityClass[j];
            dStarAngularDistDiff[k][j]->Draw("SURF1");
            str1 = str1 + ".png";
            can->SaveAs(str1);
            
        }
    }
    
    for(int k = 2; k < 3; k++){
        for(int j = 0; j < 3; j++){
        
        cout << scaleFactors[k][j] << ", ";

        
        
        }
    }
    
    inputFile->Close();

}