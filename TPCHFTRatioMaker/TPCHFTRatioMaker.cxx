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
        

int TPCHFTRatioMaker(){

    double TITLE_OFFSET=1.3;
	double integralFactor;
	double norm;

	TFile * inputRootFile = new TFile("TPC_HFT_ratio.root");
	
	TH1D * ratioPeripheral = new TH1D("Ratio_Peripheral", "Ratio_Peripheral", 100, 0, 10);
	TH1D * ratioMidCentral = new TH1D("Ratio_MidCentral", "Ratio_MidCentral", 100, 0, 10);
	TH1D * ratioCentral    = new TH1D("Ratio_Central", "Ratio_Central", 100, 0, 10);
    
    TH1D * ptTPCTotal         = new TH1D("Pt_integrated_TPC", "Pt_integrated_TPC", 100, 0, 10);
    TH1D * ptHFTTotal         = new TH1D("Pt_integrated_HFT", "Pt_integrated_HFT", 100, 0, 10);     
    TH1D * ratioTotal      = new TH1D("Ratio_total", "Ratio_total", 100, 0, 10);
	
	ratioPeripheral->GetXaxis()->SetTitle("P_{t}");
	ratioPeripheral->GetYaxis()->SetTitle("HFT Ratio");
	ratioPeripheral->SetTitleOffset(TITLE_OFFSET, "Y");
	ratioMidCentral->GetXaxis()->SetTitle("P_{t}");
	ratioMidCentral->GetYaxis()->SetTitle("HFT Ratio");
	ratioMidCentral->SetTitleOffset(TITLE_OFFSET ,"Y");
	ratioCentral->GetXaxis()->SetTitle("P_{t}");
	ratioCentral->GetYaxis()->SetTitle("HFT Ratio");
	ratioCentral->SetTitleOffset(TITLE_OFFSET, "Y");
	
	ratioPeripheral->SetMarkerStyle(20);
	
	TH1D * TPCPer;
	TH1D * TPCMidCent;
	TH1D * TPCCent;
	TH1D * HFTPer;
	TH1D * HFTMidCent;
	TH1D * HFTCent;
	
	TPCPer = (TH1D*) inputRootFile->Get("TPC_Hadron_Pt_Peripheral");
	TPCPer->SetDirectory(0);
	TPCMidCent = (TH1D*) inputRootFile->Get("TPC_Hadron_Pt_MidCentral");
	TPCMidCent->SetDirectory(0);
	TPCCent = (TH1D*) inputRootFile->Get("TPC_Hadron_Pt_Central");
	TPCCent->SetDirectory(0);
	
    ptTPCTotal->Add(TPCPer, 1);
    ptTPCTotal->Add(TPCMidCent, 1);
    ptTPCTotal->Add(TPCCent, 1);
    
    
	HFTPer = (TH1D*) inputRootFile->Get("HFT_Hadron_Pt_Peripheral");
	HFTPer->SetDirectory(0);
	HFTMidCent = (TH1D*) inputRootFile->Get("HFT_Hadron_Pt_MidCentral");
	HFTMidCent->SetDirectory(0);
	HFTCent = (TH1D*) inputRootFile->Get("HFT_Hadron_Pt_Central");
	HFTCent->SetDirectory(0);
    
    ptHFTTotal->Add(HFTPer, 1);
    ptHFTTotal->Add(HFTMidCent, 1);
    ptHFTTotal->Add(HFTCent, 1);
	
	ratioPeripheral->Divide(HFTPer, TPCPer, 1, 1);
	ratioMidCentral->Divide(HFTMidCent, TPCMidCent, 1, 1);
	ratioCentral->Divide(HFTCent, TPCCent, 1, 1);
    
    ratioPeripheral->GetXaxis()->SetRangeUser(.15,6);
	ratioMidCentral->GetXaxis()->SetRangeUser(.15,6);
	ratioCentral->GetXaxis()->SetRangeUser(.15,6);
    
    ratioTotal->Divide(ptHFTTotal, ptTPCTotal, 1, 1);
    ratioTotal->GetXaxis()->SetRangeUser(.15,6);
	
	TF1 * expoPer = new TF1("per","pol6", .15, 6);
	TF1 * expoMidCent = new TF1("midcent","pol6", .15, 6);
	TF1 * expoCent = new TF1("cent","pol6", .15, 6);
	
	ratioPeripheral->Fit(expoPer);
	ratioMidCentral->Fit(expoMidCent);
	ratioCentral->Fit(expoCent);
    
    //for(int i = 0; i < 3; i++){
    for(int j = 0; j < 7; j++){
        
        cout << "HFTRatioFunc[0]->FixParameter(" << j << ", " << expoPer->GetParameter(j) << ");" << endl;
	}
    for(int j = 0; j < 7; j++){
        
        cout << "HFTRatioFunc[1]->FixParameter(" << j << ", " << expoMidCent->GetParameter(j) << ");" << endl;
	}
    for(int j = 0; j < 7; j++){
        
        cout << "HFTRatioFunc[2]->FixParameter(" << j << ", " << expoCent->GetParameter(j) << ");" << endl;
	}
    
    for(int j = 0; j < 7; j++){
        
        cout << "HFTRatioFunc[0]->SetParError(" << j << ", " << expoPer->GetParError(j) << ");" << endl;
	}
    for(int j = 0; j < 7; j++){
        
        cout << "HFTRatioFunc[1]->SetParError(" << j << ", " << expoMidCent->GetParError(j) << ");" << endl;
	}
    for(int j = 0; j < 7; j++){
        
        cout << "HFTRatioFunc[2]->SetParError(" << j << ", " << expoCent->GetParError(j) << ");" << endl;
	}
    
    
    
	TCanvas * c1 = new TCanvas("c1", "c1", 2400, 1600);
	c1->Divide(3,2);
	
	c1->cd(1);
	ratioPeripheral->Draw();
	c1->cd(2);
	ratioMidCentral->Draw();
	c1->cd(3);
	ratioCentral->Draw();
    c1->cd(5);
	ratioTotal->Draw();
	
	c1->SaveAs("ratios.png");
	
	inputRootFile->Close();




}