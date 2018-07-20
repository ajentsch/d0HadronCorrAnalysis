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

int extractEfficiencyAverages(){


    TFile * inputFile = new TFile("C:/Users/ajentsch/desktop/average_efficiency_products.root");
    ofstream outputfile; 
    outputfile.open("C:/Users/ajentsch/desktop/output_efficiency_averages.txt");
    
    
    
    TString PtBinLabel   = "_PtBin_";
    TString siblingLabel = "Sibling_";
    TString mixLabel     = "Mixed_";
    TString effCorrInfoLabel  = "eff_corr_info_";
    TString binLabelPt[6]              = {"0", "1", "2", "3", "4", "5"};
    TString binLabelCentralityClass[3] = {"_Peripheral", "_MidCentral", "_Central"};
    TString str1;
    TString str2;
    TString str3;
    TString str4;
    
   
    
    
    /*  for(int k = 0; k < 5; k++){
        
        //eff_corr_info_Sibling__Peripheral_PtBin_0
        
        str1 = effCorrInfoLabel + siblingLabel + binLabelCentralityClass[0] + PtBinLabel + binLabelPt[k];
        str2 = effCorrInfoLabel + siblingLabel + binLabelCentralityClass[1] + PtBinLabel + binLabelPt[k];
        str3 = effCorrInfoLabel + siblingLabel + binLabelCentralityClass[2] + PtBinLabel + binLabelPt[k];
   
        efficiencyCorrInfoSib[k][0]  = new TH1D(str1, str1, 3, 0, 3);   //d0-h product, hadron-kaon-pion product, #[kpi]-h pairs
        efficiencyCorrInfoSib[k][1]  = new TH1D(str2, str2, 3, 0, 3);
        efficiencyCorrInfoSib[k][2]  = new TH1D(str3, str3, 3, 0, 3);
   
        str1 = effCorrInfoLabel + mixLabel + binLabelCentralityClass[0] + PtBinLabel + binLabelPt[k];
        str2 = effCorrInfoLabel + mixLabel + binLabelCentralityClass[1] + PtBinLabel + binLabelPt[k];
        str3 = effCorrInfoLabel + mixLabel + binLabelCentralityClass[2] + PtBinLabel + binLabelPt[k];
   
        efficiencyCorrInfoMix[k][0]  = new TH1D(str1, str1, 3, 0, 3);   //d0-h product, hadron-kaon-pion product, #[kpi]-h pairs
        efficiencyCorrInfoMix[k][1]  = new TH1D(str2, str2, 3, 0, 3);
        efficiencyCorrInfoMix[k][2]  = new TH1D(str3, str3, 3, 0, 3);
   
    } */
    
   
    
    TH1D* sibPer[5];
    TH1D* sibMid[5];
    TH1D* sibCent[5];
    TH1D* mixPer[5];
    TH1D* mixMid[5];
    TH1D* mixCent[5];
    
    for(int k = 0; k < 5; k++){
    
        str1 = effCorrInfoLabel + siblingLabel + binLabelCentralityClass[0] + PtBinLabel + binLabelPt[k];
        sibPer[k] = (TH1D*)inputFile->Get(str1);
        str1 = effCorrInfoLabel + siblingLabel + binLabelCentralityClass[1] + PtBinLabel + binLabelPt[k];
        sibMid[k] = (TH1D*)inputFile->Get(str1);
        str1 = effCorrInfoLabel + siblingLabel + binLabelCentralityClass[2] + PtBinLabel + binLabelPt[k];
        sibCent[k] = (TH1D*)inputFile->Get(str1);
        str1 = effCorrInfoLabel + mixLabel + binLabelCentralityClass[0] + PtBinLabel + binLabelPt[k];
        mixPer[k] = (TH1D*)inputFile->Get(str1);
        str1 = effCorrInfoLabel + mixLabel + binLabelCentralityClass[1] + PtBinLabel + binLabelPt[k];
        mixMid[k] = (TH1D*)inputFile->Get(str1);
        str1 = effCorrInfoLabel + mixLabel + binLabelCentralityClass[2] + PtBinLabel + binLabelPt[k];
        mixCent[k] = (TH1D*)inputFile->Get(str1);
    }
    
    double tripleHadronSibPer = 1;
    double tripleHadronSibMid = 1;
    double tripleHadronSibCent = 1;
    double tripleHadronMixPer = 1;
    double tripleHadronMixMid = 1;
    double tripleHadronMixCent = 1;
    
    double doubleHadronSibPer = 1;
    double doubleHadronSibMid = 1;
    double doubleHadronSibCent = 1;
    double doubleHadronMixPer = 1;
    double doubleHadronMixMid = 1;
    double doubleHadronMixCent = 1;
    
    
    
    for(int k = 0; k < 5; k++){
    
        tripleHadronSibPer = sibPer[k]->GetBinContent(1)/sibPer[k]->GetBinContent(3);
        tripleHadronSibMid = sibMid[k]->GetBinContent(1)/sibMid[k]->GetBinContent(3);
        tripleHadronSibCent = sibCent[k]->GetBinContent(1)/sibCent[k]->GetBinContent(3);
        
        outputfile << "siblingTripleEff[" << k << "][0] = " << tripleHadronSibPer << ";" << endl;
        outputfile << "siblingTripleEff[" << k << "][1] = " << tripleHadronSibMid << ";" << endl;
        outputfile << "siblingTripleEff[" << k << "][2] = " << tripleHadronSibCent << ";" << endl;
        
        
    }
    
    outputfile << endl;
    
    for(int k = 0; k < 5; k++){
    
        tripleHadronMixPer = mixPer[k]->GetBinContent(1)/mixPer[k]->GetBinContent(3);
        tripleHadronMixMid = mixMid[k]->GetBinContent(1)/mixMid[k]->GetBinContent(3);
        tripleHadronMixCent = mixCent[k]->GetBinContent(1)/mixCent[k]->GetBinContent(3);
    
        outputfile << "mixedTripleEff[" << k << "][0] = " << tripleHadronMixPer << ";" << endl;
        outputfile << "mixedTripleEff[" << k << "][1] = " << tripleHadronMixMid << ";" << endl;
        outputfile << "mixedTripleEff[" << k << "][2] = " << tripleHadronMixCent << ";" << endl;
    
    }
    
    outputfile << endl;
    
    for(int k = 0; k < 5; k++){
    
        doubleHadronSibPer = sibPer[k]->GetBinContent(2)/sibPer[k]->GetBinContent(3);
        doubleHadronSibMid = sibMid[k]->GetBinContent(2)/sibMid[k]->GetBinContent(3);
        doubleHadronSibCent = sibCent[k]->GetBinContent(2)/sibCent[k]->GetBinContent(3);
   
        outputfile << "siblingDoubleEff[" << k << "][0] = " << doubleHadronSibPer << ";" << endl;
        outputfile << "siblingDoubleEff[" << k << "][1] = " << doubleHadronSibMid << ";" << endl;
        outputfile << "siblingDoubleEff[" << k << "][2] = " << doubleHadronSibCent << ";" << endl;
    }    
    
    for(int k = 0; k < 5; k++){
    
        doubleHadronMixPer = mixPer[k]->GetBinContent(2)/mixPer[k]->GetBinContent(3);
        doubleHadronMixMid = mixMid[k]->GetBinContent(2)/mixMid[k]->GetBinContent(3);
        doubleHadronMixCent = mixCent[k]->GetBinContent(2)/mixCent[k]->GetBinContent(3);
        
        outputfile << "mixedDoubleEff[" << k << "][0] = " << doubleHadronMixPer << ";" << endl;
        outputfile << "mixedDoubleEff[" << k << "][1] = " << doubleHadronMixMid << ";" << endl;
        outputfile << "mixedDoubleEff[" << k << "][2] = " << doubleHadronMixCent << ";" << endl;

    }        
    
      
    cout << "Output file created." << endl;
    
    inputFile->Close();
    outputfile.close();
    
    
    





}