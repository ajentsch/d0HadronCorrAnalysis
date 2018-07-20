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


TString convertDoubleToString(double number){
    TString str;
    str.Form("%.5f", number);
    return str;   
}


int histogramExtractAndSave(TFile * inputRootFile, TString outputFolder  = "C:/Users/ajentsch/Desktop/Histograms/"){

    //TString inputFolder   = "C:/Users/ajentsch/Desktop/";
    //TString fileName      = "all";
    //TString rootString    = ".root";
    //TString inputFileFinal = inputFolder + fileName + rootString;
    //TString outputFolder  = "C:/Users/ajentsch/Desktop/Histograms/";
    
    //TFile * inputRootFile = new TFile(inputFileFinal);
    
    double numEvents = 0;
    
    int NUM_PT_BINS = 5;    
        
    TString D0EtaLabel     = "D0_Eta_Dist";
    TString D0PhiLabel     = "D0_Phi_Dist";
    TString D0KaonPtLabel  = "D0_Kaon_Candidate_Pt_Dist";
    TString D0KaonEtaLabel = "D0_Kaon_Eta_Dist";
    TString D0KaonPhiLabel = "D0_Kaon_Phi_Dist";
    TString D0PionPtLabel  = "D0_Candidate_Pion_Pt_Dist";
    TString D0PionEtaLabel = "D0_Pion_Eta_Dist";
    TString D0PionPhiLabel = "D0_Pion_Phi_Dist";
    TString D0DecayLengthLabel = "D0_Decay_Length";
    TString D0DCAToPVLabel = "D0_DCA_to_PV";
    TString D0KaonPVLabel = "D0_Daughter_Kaon_DCA_to_PV";
    TString D0PionPVLabel = "D0_Daughter_Pion_DCA_to_PV";
    TString D0DaughterDCALabel = "D0_Daughter_Pair_DCA";
    
    TString outputLabelDaughterKaons = "Daughter_Kaons_";
    TString outputLabelDaughterPions = "Daughter_Pions_";
    
    
    if(NUM_PT_BINS == 3){ 
    
        TString numbers[3] = {"0_1", "1_2", "2_10"};
        TString D0CutPtBinsLabel[4] = {"0-1", "1-2", "2-10", "Integrated"};
    }
    
    if(NUM_PT_BINS == 5){
    
        TString numbers[5] = {"0_1", "1_2", "2_3", "3_4", "4_10"}; //changed
        TString D0CutPtBinsLabel[6] = {"0-1", "1-2", "2-3", "3-4", "4-10", "Integrated"};
    
    }
    //TString numbers[3] = {"0_1", "1_2", "2_10"};
    //TString numbers[5] = {"0_1", "1_2", "2_3","3_5", "5_10"};
    //TString numbers[4] = {"0_2", "2_3","3_5", "5_10"};
    
    TString PtBinLabel   = "_PtBin_";
    
    //TString D0CutPtBinsLabel[4] = {"0-1", "1-2", "2-10", "Integrated"};
    //TString D0CutPtBinsLabel[6] = {"0-1", "1-2", "2-3", "3-5", "5-10", "Integrated"};
    //TString D0CutPtBinsLabel[5] = {"0-2", "2-3", "3-5", "5-10", "Integrated"};
    
    TString str1;
    TString str2;
    TString fileNumber;
    TString outputPNG     = ".png";
    
    TH1D* D0KaonPtDist[5];
    TH1D* D0KaonEtaDist[5]; 
    TH1D* D0KaonPhiDist[5]; 
    TH1D* D0PionPtDist[5];
    TH1D* D0PionEtaDist[5]; 
    TH1D* D0PionPhiDist[5]; 
    TH1D* D0EtaDist[5];
    TH1D* D0PhiDist[5];
    TH1D* D0DecayLengthDist[5];
    TH1D* D0DCAToPVDist[5];
    TH1D* D0KaonPVDist[5];
    TH1D* D0PionPVDist[5];
    TH1D* D0DaughterDCADist[5];
    
    TH1D* D0EtaDistInt;
    TH1D* D0PhiDistInt;
    TH1D* D0KaonPtDistInt;
    TH1D* D0KaonEtaDistInt;
    TH1D* D0KaonPhiDistInt;
    TH1D* D0PionPtDistInt;
    TH1D* D0PionEtaDistInt;
    TH1D* D0PionPhiDistInt;
    TH1D* D0DecayLengthDistInt;
    TH1D* D0DCAToPVDistInt;
    TH1D* D0KaonPVDistInt;
    TH1D* D0PionPVDistInt;
    TH1D* D0DaughterDCADistInt;
    
    
   

    for(int band = 0; band < NUM_PT_BINS; band++){
    
        str1 = D0EtaLabel + PtBinLabel + D0CutPtBinsLabel[band];
        D0EtaDist[band] = (TH1D*) inputRootFile->Get(str1);
        D0EtaDist[band]->SetDirectory(0);
        str1 = D0PhiLabel + PtBinLabel + D0CutPtBinsLabel[band];
        D0PhiDist[band] = (TH1D*) inputRootFile->Get(str1);
        D0PhiDist[band]->SetDirectory(0);
        str1 = D0KaonPtLabel + PtBinLabel + D0CutPtBinsLabel[band];
        D0KaonPtDist[band]    = (TH1D*) inputRootFile->Get(str1);
        D0KaonPtDist[band]->SetDirectory(0);
        str1 = D0KaonEtaLabel + PtBinLabel + D0CutPtBinsLabel[band];
        D0KaonEtaDist[band]   = (TH1D*) inputRootFile->Get(str1);
        D0KaonEtaDist[band]->SetDirectory(0);
        str1 = D0KaonPhiLabel + PtBinLabel + D0CutPtBinsLabel[band];
        D0KaonPhiDist[band]   = (TH1D*) inputRootFile->Get(str1);
        D0KaonPhiDist[band]->SetDirectory(0);
        str1 = D0PionPtLabel + PtBinLabel + D0CutPtBinsLabel[band];
        D0PionPtDist[band]    = (TH1D*) inputRootFile->Get(str1);
        D0PionPtDist[band]->SetDirectory(0);
        str1 = D0PionEtaLabel + PtBinLabel + D0CutPtBinsLabel[band];
        D0PionEtaDist[band]   = (TH1D*) inputRootFile->Get(str1);
        D0PionEtaDist[band]->SetDirectory(0);
        str1 = D0PionPhiLabel + PtBinLabel + D0CutPtBinsLabel[band];
        D0PionPhiDist[band]   = (TH1D*) inputRootFile->Get(str1);
        D0PionPhiDist[band]->SetDirectory(0);
        
        str1 = D0DecayLengthLabel + PtBinLabel + D0CutPtBinsLabel[band];
        D0DecayLengthDist[band] = (TH1D*) inputRootFile->Get(str1);
        D0DecayLengthDist[band]->SetDirectory(0);
        str1 = D0DCAToPVLabel + PtBinLabel + D0CutPtBinsLabel[band];
        D0DCAToPVDist[band] = (TH1D*) inputRootFile->Get(str1);
        D0DCAToPVDist[band]->SetDirectory(0);
        str1 = D0KaonPVLabel + PtBinLabel + D0CutPtBinsLabel[band];
        D0KaonPVDist[band] = (TH1D*) inputRootFile->Get(str1);
        D0KaonPVDist[band]->SetDirectory(0);
        str1 = D0PionPVLabel + PtBinLabel + D0CutPtBinsLabel[band];
        D0PionPVDist[band] = (TH1D*) inputRootFile->Get(str1);
        D0PionPVDist[band]->SetDirectory(0);
        str1 = D0DaughterDCALabel + PtBinLabel + D0CutPtBinsLabel[band];
        D0DaughterDCADist[band] = (TH1D*) inputRootFile->Get(str1);
        D0DaughterDCADist[band]->SetDirectory(0);
        
    }
    
    cout << "HERE" << endl;
    
    TH1D* kaonEtaDist = (TH1D*) inputRootFile->Get("Kaon_Eta_Distribution"); 
    kaonEtaDist->SetDirectory(0);
    TH1D* pionEtaDist = (TH1D*) inputRootFile->Get("Pion_Eta_Distribution");
    pionEtaDist->SetDirectory(0);
    TH1D* kaonPhiDist = (TH1D*) inputRootFile->Get("Kaon_Phi_Distribution");    
    kaonPhiDist->SetDirectory(0);  
    TH1D* pionPhiDist = (TH1D*) inputRootFile->Get("Pion_Phi_Distribution"); 
    pionPhiDist->SetDirectory(0);  
    TH1D* DCAtoPrimaryVertex = (TH1D*) inputRootFile->Get("track_DCA_to_PV");   
    DCAtoPrimaryVertex->SetDirectory(0);
    TH1D* DCAtoPrimaryVertexCut= (TH1D*) inputRootFile->Get("track_DCA_to_PV_cut_check");
    DCAtoPrimaryVertexCut->SetDirectory(0);
    TH1D* hadronPtDist = (TH1D*) inputRootFile->Get("Inclusive_Hadron_pt");          
    hadronPtDist->SetDirectory(0);
    TH1D* hadronPhiDist = (TH1D*) inputRootFile->Get("Inclusive_Hadron_Phi");         
    hadronPhiDist->SetDirectory(0);
    TH1D* hadronEtaDist  = (TH1D*) inputRootFile->Get("Inclusive_Hadron_Eta");       
    hadronEtaDist->SetDirectory(0);
    TH1D* hadronChi2 = (TH1D*) inputRootFile->Get("Chi2_for_hadron_tracks");            
    hadronChi2->SetDirectory(0);
    TH1D* hadronNHitsFit = (TH1D*) inputRootFile->Get("nHitsFit");       
    hadronNHitsFit->SetDirectory(0);
    TH1D* hadronNHitsFitOverMax = (TH1D*) inputRootFile->Get("nHitsFitOverMax");
    hadronNHitsFitOverMax->SetDirectory(0);
    TH2D* dEdxVsPt = (TH2D*) inputRootFile->Get("dEdx_vs_P");             
    dEdxVsPt->SetDirectory(0); 
    TH2D* invBetaVsPt = (TH2D*) inputRootFile->Get("inv_Beta_Vs_P");         
    invBetaVsPt->SetDirectory(0);
    
    
    TH1D* D0ptDist = (TH1D*) inputRootFile->Get("D0_Candidate_Pt_Dist_PtBin_Integrated");           
    TH1D* D0EtaDistInt = (TH1D*) inputRootFile->Get("D0_Eta_Dist_PtBin_Integrated");         
    TH1D* D0PhiDistInt = (TH1D*) inputRootFile->Get("D0_Phi_Dist_PtBin_Integrated");       
    TH1D* D0KaonPtDistInt = (TH1D*) inputRootFile->Get("D0_Kaon_Candidate_Pt_Dist_PtBin_Integrated");    
    TH1D* D0KaonEtaDistInt = (TH1D*) inputRootFile->Get("D0_Kaon_Eta_Dist_PtBin_Integrated");    
    TH1D* D0KaonPhiDistInt = (TH1D*) inputRootFile->Get("D0_Kaon_Phi_Dist_PtBin_Integrated");    
    TH1D* D0PionPtDistInt = (TH1D*) inputRootFile->Get("D0_Candidate_Pion_Pt_Dist_PtBin_Integrated");   
    TH1D* D0PionEtaDistInt = (TH1D*) inputRootFile->Get("D0_Pion_Eta_Dist_PtBin_Integrated"); 
    TH1D* D0PionPhiDistInt = (TH1D*) inputRootFile->Get("D0_Pion_Phi_Dist_PtBin_Integrated");    
    TH1D* D0DecayLengthDistInt= (TH1D*) inputRootFile->Get("D0_Decay_Length_PtBin_Integrated"); 
    TH1D* D0DCAToPVDistInt = (TH1D*) inputRootFile->Get("D0_DCA_to_PV_PtBin_Integrated");   
    TH1D* D0KaonPVDistInt = (TH1D*) inputRootFile->Get("D0_Daughter_Kaon_DCA_to_PV_PtBin_Integrated");     
    TH1D* D0PionPVDistInt = (TH1D*) inputRootFile->Get("D0_Daughter_Pion_DCA_to_PV_PtBin_Integrated");     
    TH1D* D0DaughterDCADistInt = (TH1D*) inputRootFile->Get("D0_Daughter_Pair_DCA_PtBin_Integrated");
        
    /*TH1D* D0ptDist;
    TH1D* D0EtaDistInt;
    TH1D* D0PhiDistInt;
    TH1D* D0KaonPtDistInt;
    TH1D* D0KaonEtaDistInt;
    TH1D* D0KaonPhiDistInt;
    TH1D* D0PionPtDistInt;
    TH1D* D0PionEtaDistInt;
    TH1D* D0PionPhiDistInt;
    TH1D* D0DecayLengthDistInt;
    TH1D* D0DCAToPVDistInt;
    TH1D* D0KaonPVDistInt;
    TH1D* D0PionPVDistInt;
    TH1D* D0DaughterDCADistInt;*/
        
    TCanvas *canvas = new TCanvas("c2","c2", 2400, 2400);
    
    TString pdfFile = "QA_Histograms.pdf[";
    TString pdfOpen = outputFolder + pdfFile;
    pdfFile = "QA_Histograms.pdf";
    TString pdfPage = outputFolder + pdfFile;
    pdfFile = "QA_Histograms.pdf]";
    TString pdfClose = outputFolder + pdfFile;
    
    
    
    canvas->Print(pdfOpen); //opens pdf file
    
    //TPaveText *cutTextBox = new TPaveText(0.0, 0.0, 1.0, 1.0, "NB NDC");
    
    
    //str1 = "D0 InvMass Signal Band -\t Low: ";
    //str2 = convertDoubleToString(inputRootFile->Get("HistOfCuts"))->GetBinContent(1));
    
    
    canvas->Divide(2,2);
    
        
    canvas->cd(1);
    kaonEtaDist->Draw();
    canvas->cd(2);
    kaonPhiDist->Draw();
    canvas->Print(pdfPage);
        
    canvas->cd(1);
    pionEtaDist->Draw();
    canvas->cd(2);
    pionPhiDist->Draw();
    canvas->Print(pdfPage);
    
    //TCanvas *canvas2 = new TCanvas("c3","c3", 1600, 1600);
    //canvas2->Divide(2,2);
    
    canvas->cd(1);
    hadronPtDist->Draw();
    canvas->cd(2);
    hadronPhiDist->Draw();
    canvas->cd(3);
    hadronEtaDist->Draw();
    canvas->Print(pdfPage);        
    
    canvas->Clear();
    
    canvas->Divide(2,3);
    
    for(int band = 0; band < NUM_PT_BINS; band++){
    
        canvas->cd(1);
        D0EtaDist[band]->Draw();
        canvas->cd(2);
        D0PhiDist[band]->Draw();
        canvas->cd(3);
        D0DecayLengthDist[band]->Draw();
        canvas->cd(4);
        D0DCAToPVDist[band]->Draw();
        canvas->cd(5);
        D0DaughterDCADist[band]->Draw();
        
        canvas->Print(pdfPage);
        
    }

    canvas->Clear();
    
    canvas->Divide(2,2);
    
    for(int band = 0; band < NUM_PT_BINS; band++){
    
        //str1 = outputLabelDaughterKaons + numbers[band];
        //canvas->cd(1);
        //D0KaonPtDist[band]->Draw();
        canvas->cd(1);
        D0KaonPtDist[band]->Draw();
        canvas->cd(2);
        D0KaonEtaDist[band]->Draw();
        canvas->cd(3);
        D0KaonPhiDist[band]->Draw();
        canvas->cd(4);
        D0KaonPVDist[band]->Draw();
        
        canvas->Print(pdfPage);
        
    }    
        
    for(int band = 0; band < NUM_PT_BINS; band++){    
        
        //str1 = outputLabelDaughterPions + numbers[band];
        //canvas->cd(1);
        //D0PionPtDist[band]->Draw();
        canvas->cd(1);
        D0PionPtDist[band]->Draw();
        canvas->cd(2);
        D0PionEtaDist[band]->Draw();
        canvas->cd(3);
        D0PionPhiDist[band]->Draw();
        canvas->cd(4);
        D0PionPVDist[band]->Draw();
        
        canvas->Print(pdfPage);
    }
    
    
    canvas->Clear();
    canvas->Divide(2,3);
    
    canvas->cd(1);
    D0ptDist->Draw();
    canvas->cd(2);
    D0EtaDistInt->Draw();
    canvas->cd(3);
    D0PhiDistInt->Draw();
    canvas->cd(4);
    D0DecayLengthDistInt->Draw();
    canvas->cd(5);
    D0DCAToPVDistInt->Draw();
    canvas->cd(6);
    D0DaughterDCADistInt->Draw();
    
    canvas->Print(pdfPage);
    
    canvas->Clear();
    canvas->Divide(2,2);
    
    canvas->cd(1);
    D0KaonPtDistInt->Draw();
    canvas->cd(2);
    D0KaonEtaDistInt->Draw();
    canvas->cd(3);
    D0KaonPhiDistInt->Draw();
    canvas->cd(4);
    D0KaonPVDistInt->Draw();
    
    canvas->Print(pdfPage);
    
    canvas->cd(1);
    D0PionPtDistInt->Draw();
    canvas->cd(2);
    D0PionEtaDistInt->Draw();
    canvas->cd(3);
    D0PionPhiDistInt->Draw();
    canvas->cd(4);
    D0PionPVDistInt->Draw();
    
    canvas->Print(pdfPage);
    
    canvas->Print(pdfClose); //closes pdf file
    
    //inputRootFile->Close();
    
}