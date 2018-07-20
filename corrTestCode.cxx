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


int corrTestCode(){


    TString inputRootFilePath = "C:/Users/ajentsch/desktop/D0_Hadron_Correlation_Root_Files/optimized_cuts.root";
    TString siblingHistLabel       = "Sibling_US_correlation_PtBin_1_VzBin_";
    TString mixedHistLabel         = "Mixed_US_correlation_PtBin_1_VzBin_";
    
    TString vzBinLabel[10]         = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
    TString centBinLabel           = "_CentBin_8";
    
    double numSibVz[10];
    double numSibTotal = 0;
    double weightFactor;
    
    TString deltaRhoLabel          = "deltaRho";
    TString deltaRhoOverRhoLabel   = "deltaRhoOverRho";
    
    double NSibling;
    double NMixed;
    double normFactor;
    
    double phiBinShift = (TMath::Pi()/12.0);
    
    TString str1;
    TString str2;
    
    TH2D* siblingHist;
    TH2D* mixedHist;
    TH2D* deltaRho            = new TH2D(deltaRhoLabel, deltaRhoLabel, 9, -2, 2, 12, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D* deltaRhoOverRhoTemp = new TH2D("","", 9, -2, 2, 12, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D* deltaRhoOverRho     = new TH2D(deltaRhoLabel, deltaRhoLabel, 9, -2, 2, 12, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    TFile * inputRootFile = new TFile(inputRootFilePath);
    
    TCanvas * canvas = new TCanvas("","", 800, 800);
    
    for(int i = 0; i < 10; i++){
    
        str1 = mixedHistLabel + vzBinLabel[i] + centBinLabel;
        
        cout << str1 << endl;
    
        mixedHist = (TH2D*)inputRootFile->Get(str1);
        mixedHist->SetDirectory(0);
        
        NSibling = mixedHist->Integral(1, 9, 1, 12);
        
        numSibVz[i] = NSibling;
        numSibTotal = numSibTotal + NSibling;
    }    
    
    for(int i = 0; i < 10; i++){
    
        str1 = siblingHistLabel + vzBinLabel[i] + centBinLabel;
        str2 = mixedHistLabel + vzBinLabel[i] + centBinLabel;
    
        siblingHist = (TH2D*)inputRootFile->Get(str1);
        siblingHist->SetDirectory(0);
        mixedHist = (TH2D*)inputRootFile->Get(str2);
        mixedHist->SetDirectory(0);

        NSibling = siblingHist->Integral(1,9, 1, 12);
        NMixed   = mixedHist->Integral(1, 9, 1, 12);
        normFactor = NSibling/NMixed;
        
        weightFactor = numSibVz[i]/numSibTotal;
    
        deltaRho->Add(siblingHist, mixedHist, 1, -normFactor);
        deltaRhoOverRhoTemp->Divide(deltaRho, mixedHist, 1, normFactor);
        deltaRhoOverRho->Add(deltaRhoOverRhoTemp, weightFactor);
    
        //siblingHist->Draw("SURF1");
        //canvas->SaveAs("C:/Users/ajentsch/desktop/corrTestCode/sibling.png");
        // mixedHist->Draw("SURF1");
        //canvas->SaveAs("C:/Users/ajentsch/desktop/corrTestCode/mixed.png");
        //deltaRho->Draw("SURF1");
        //canvas->SaveAs("C:/Users/ajentsch/desktop/corrTestCode/deltaRho.png");
        //deltaRhoOverRho->Draw("SURF1");
        //canvas->SaveAs("C:/Users/ajentsch/desktop/corrTestCode/deltaRhoOverRho.png");
    
    }
    
        deltaRhoOverRho->GetXaxis()->SetRangeUser(-1.55, 1.55);        
        deltaRhoOverRho->Draw("SURF1");
        canvas->SaveAs("C:/Users/ajentsch/desktop/corrTestCode/deltaRhoOverRho.png");
    
    inputRootFile->Close();



}