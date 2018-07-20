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


int quickRootFileReadIn(){

    TFile *file1 = new TFile("C:/Users/ajentsch/Desktop/F67341F52F967C4F2E39026DF61ACF60_422.root");
    TFile *file2 = new TFile("C:/Users/ajentsch/Desktop/F67341F52F967C4F2E39026DF61ACF60_422.root");
    
    //file.open("C:/Users/ajentsch/Desktop/3FF0C11E3F4496C22654B1411E336923_All.root",ifstream::in);
    
    TFile *output = new TFile("C:/Users/ajentsch/Desktop/ptDistCompare.root", "RECREATE");
    
    TH1D* hadronPtDistTPCOnly = new TH1D("Inclusive Hadron pt - TPC Only", "Inclusive Hadron pt - TPC Only", 1000, 0, 10);
    TH1D* hadronPtDistHFTOnly = new TH1D("Inclusive Hadron pt - HFT Only", "Inclusive Hadron pt - HFT Only", 1000, 0, 10);
    TH1D* tmp = new TH1D("Inclusive Hadron pt", "Inclusive Hadron pt", 1000, 0, 10);
    TH1D* ratio = new TH1D("Ratio of HFT/TPC", "Ratio of HFT/TPC", 1000, 0, 10);
    TH1D* eventCounter = new TH1D("Event Counter", "Event Counter", 6, 0, 6);
    
    tmp = (TH1D*)file1->Get("Inclusive Hadron pt");
    
    hadronPtDistTPCOnly = (TH1D*)tmp->Clone();
    
    tmp = (TH1D*)file2->Get("Inclusive Hadron pt");
    
    hadronPtDistHFTOnly = (TH1D*)tmp->Clone();
    
    ratio->Divide(hadronPtDistHFTOnly, hadronPtDistTPCOnly, 1, 1);
    
    double HFTLowPt = 0;
    double HFTHighPt = 0;
    double TPCLowPt = 0;
    double TPCHighPt = 0;
    
    HFTLowPt = hadronPtDistHFTOnly->Integral(1,75);
    HFTHighPt = hadronPtDistHFTOnly->Integral(76,1000);
    TPCLowPt = hadronPtDistTPCOnly->Integral(1,75);
    TPCHighPt = hadronPtDistTPCOnly->Integral(76,1000);
    
    cout << endl;
    cout << "HFT pt counts low: " << HFTLowPt << "  " << "High: " << HFTHighPt << endl;
    cout << "TPC pt counts low: " << TPCLowPt << "  " << "High: " << TPCHighPt << endl;
    
    
    file1->Close();
    file2->Close();
    
    hadronPtDistTPCOnly->Write();
    hadronPtDistHFTOnly->Write();
    ratio->Write();
    
    output->Close();
    
    
    TFile *file3 = new TFile("C:/Users/ajentsch/Desktop/F67341F52F967C4F2E39026DF61ACF60_422.root");
    
    double numD0Events =0;
    double numLeftSBEvents = 0;
    double numRightSBEvents = 0;
    
    tmp = (TH1D*)file3->Get("number of events used");
    eventCounter = (TH1D*)tmp->Clone();
    
    
    numD0Events = eventCounter->GetBinContent(1);
    numLeftSBEvents = eventCounter->GetBinContent(5);
    numRightSBEvents = eventCounter->GetBinContent(6);
    
    cout << "Num D0 Events: " << numD0Events << endl; 
    cout << "Num Left SB Events: " << numLeftSBEvents << endl;
    cout << "Num Right Events: " << numRightSBEvents << endl;
    
    file3->Close();
    
    
    
    
}    