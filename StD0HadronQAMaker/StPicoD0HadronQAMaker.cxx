#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TF1.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"

#include "StPicoCharmContainers/StPicoD0Event.h"
#include "StPicoCharmContainers/StKaonPion.h"



#include "StPicoD0HadronQAMaker.h"
#include "StPicoHFMaker/StHFCuts.h"


/****

Author: Alex Jentsch



****/

ClassImp(StPicoD0HadronQAMaker)

StPicoD0HadronQAMaker::StPicoD0HadronQAMaker(char const * name, char const * inputFilesList, 
                                   char const * outName, StPicoDstMaker* picoDstMaker): 
                            StMaker(name),mPicoDstMaker(picoDstMaker),mPicoD0Event(NULL),
                            mOutFileName(outName), mInputFileList(inputFilesList),
                            mOutputFile(NULL), mChain(NULL), mEventCounter(0), mHFCuts(NULL)
{}

//---------------------------------------Initialization-----------------------------------------
Int_t StPicoD0HadronQAMaker::Init()
{
   mPicoD0Event = new StPicoD0Event();

   mChain = new TChain("T");
   std::ifstream listOfFiles(mInputFileList.Data());
   if (listOfFiles.is_open())
   {
      std::string file;
      while (getline(listOfFiles, file))
      {
         LOG_INFO << "StPicoD0HadronQAMaker - Adding :" << file << endm;
         mChain->Add(file.c_str());
      }
   }
   else
   {
      LOG_ERROR << "StPicoD0HadronQAMaker - Could not open list of files. ABORT!" << endm;
      return kStErr;
   }

   mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
   mChain->SetBranchAddress("dEvent", &mPicoD0Event);

   mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
   mOutputFile->cd();

   if (!mHFCuts)
    mHFCuts = new StHFCuts;   
   mHFCuts->init();
   
   
   //-----------------------FLAGS--------------------------------//
   DEBUG               = false;
   DEBUG_MIX_BUFFER    = false;
   USE_CENT_BINS       = true;
   USE_VZ_BINS         = true;
   USE_PT_BINS         = true;
   USE_FINE_PT_BINS    = false;
   SINGLES_DISTS       = false;
   
   //----------------------------------------------------------//
   
   //---------------------Important constants-----------------//
   BUFFER_SIZE   = 5;
   NUM_PHI_BINS  = 12;
   NUM_ETA_BINS  = 9;
   NUM_VZ_BINS   = 10;
   NUM_CENT_BINS = 16; 
   NUM_PT_BINS   = 3;
  
   
   
   // --------------------Event Mixer Buffer-------------------------------------
   
    if(USE_VZ_BINS) { nVzBins = NUM_VZ_BINS; }  //flag to set binning on Vz
    else nVzBins = 1;
    
    if(USE_CENT_BINS) { nCentBins = NUM_CENT_BINS; }
    else nCentBins = 1;
    
   
   
    //--------------------CUTS------------------------------------------
   
    // D0 Cuts
   
    kaonPtCut           = .15;  //daughter low pt acceptance cut
    pionPtCut           = .15;  //daughter low pt acceptance cut
    
    D0InvMassLow        = 1.82; //signal region US invariant mass
    D0InvMassHigh       = 1.90; //signal region US invariant mass
    USSideBandLeftLow   = 1.62;
    USSideBandLeftHigh  = 1.70;
    USSideBandRightLow  = 2.0;
    USSideBandRightHigh = 2.1;
    
    d0PtLow             = 0.0;    //nominal cuts
    d0PtHigh            = 10.0;
    //d0DecayLengthMin    = .0180;
    d0DecayLengthMax    = 999999.0;
    //daughterDCA         = .0055;
    d0DaughterPionPtMin = .15;
    d0DaughterKaonPtMin = .15;
    //kaonDCA             = .008;
    //pionDCA             = .008;
    //d0DCAtoPV           = .0065;
    
    
    //d0TopologicalCutArray[5][5]; //first index is pt bin, second index is cuts -- decayLength, DCADaughters, d0DCAPV, PiDCAPV, KDCAPV
    
    //pt = 0-1 GeV/c
    d0TopologicalCutArray[0][0] = .0145; //.0180;//minimum decay length 
    d0TopologicalCutArray[0][1] = .0084; //.0055; //maximum DCA daughters 
    d0TopologicalCutArray[0][2] = .0065; //.0061; //maximum DCA between d0 vector and PV 
    d0TopologicalCutArray[0][3] = .008;  //.0110; //minimum DCA between Pi and PV 
    d0TopologicalCutArray[0][4] = .008;  //.0103; //minimum DCA between K and PV
    
    //pt = 1-2 GeV/c
    d0TopologicalCutArray[1][0] = .0181; //minimum decay length
    d0TopologicalCutArray[1][1] = .0066; //maximum DCA daughters
    d0TopologicalCutArray[1][2] = .0049; //maximum DCA between d0 vector and PV 
    d0TopologicalCutArray[1][3] = .0111; //minimum DCA between Pi and PV
    d0TopologicalCutArray[1][4] = .0091; //minimum DCA between K and PV
    
    //pt = 2-3 GeV/c
    d0TopologicalCutArray[2][0] = .0212; //minimum decay length
    d0TopologicalCutArray[2][1] = .0057; //maximum DCA daughters
    d0TopologicalCutArray[2][2] = .0038; //maximum DCA between d0 vector and PV 
    d0TopologicalCutArray[2][3] = .0086; //minimum DCA between Pi and PV
    d0TopologicalCutArray[2][4] = .0095; //minimum DCA between K and PV
    
    //pt = 3-5 GeV/c
    d0TopologicalCutArray[3][0] = .0247; //minimum decay length
    d0TopologicalCutArray[3][1] = .0050; //maximum DCA daughters
    d0TopologicalCutArray[3][2] = .0038; //maximum DCA between d0 vector and PV 
    d0TopologicalCutArray[3][3] = .0081; //minimum DCA between Pi and PV
    d0TopologicalCutArray[3][4] = .0079; //minimum DCA between K and PV
    
    //pt = 5-10 GeV/c
    d0TopologicalCutArray[4][0] = .0259; //minimum decay length
    d0TopologicalCutArray[4][1] = .0060; //maximum DCA daughters
    d0TopologicalCutArray[4][2] = .0040; //maximum DCA between d0 vector and PV 
    d0TopologicalCutArray[4][3] = .0062; //minimum DCA between Pi and PV
    d0TopologicalCutArray[4][4] = .0058; //minimum DCA between K and PV
    
    
    
    //hadron cuts
    
    hadronPtMin         = .15;
    hadronPtMax         = 15.45;
    
    trackChi2max        = 3.0;
    trackDCAtoPvtx      = 3.0;
    
	//Run 14 MinBias triggers
    trigger[0] = 450050;
    trigger[1] = 450060;
    trigger[2] = 450005;
    trigger[3] = 450015;
    trigger[4] = 450025;
    //-----------------------------------------------------------------
    
    eventNumber = 1;
    
    TString VzBinLabel   = "_VzBin_";
    TString CentBinLabel = "_CentBin_";
    TString PtBinLabel   = "_PtBin_";
    
    //Labels for D0-hadron US and LS corr histograms
    
   
    
    //other labels
   
    TString eventCounterLabel = "Event Count Vz ";
    TString etaLabel          = "Inclusive Hadron Eta";
    TString phiLabel          = "Inclusive Hadron Phi";
    TString etaPhiLabel       = "Inclusive 2D Hadron Eta/Phi";
    
    TString phiD0vsPhiHLabel  = "phiD0_vs_phiH";
    TString etaD0vsEtaHLabel  = "etaD0_vs_etaH";
    TString phiD0vsEtaD0Label = "phiD0_vs_etaD0";
    TString phiHvsEtaHLabel   = "phiH_vs_etaH";
    
    
    TString invMassD0PtBin    = "D0_US_invMass_Pt_Bin_";
    TString invMassLSPtBin    = "LS_invMass_Pt_Bin_";
    
    TString fine = "FINE_";
    
    TString str1;
    TString str2;
    TString str3;
    TString str4;
    
    TString binLabelVz[10]             = {"0", "1", "2", "3", "4", "5", "6", "7", "8","9"};
    TString binLabelCent[16]           = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"};
    TString binLabelPt[6]              = {"0", "1", "2", "3", "4", "5"};
    TString binLabelFinePt[11]         = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
    TString binLabelCentralityClass[3] = {"_Peripheral", "_MidCentral", "_Central"};
   
   // --------------------Begin User Variables-----------------------------------
   
   //ptDist          = new TH1D("Pt Distribution", "Pt Distribution", 1000, 0, 10);              
   invMass         = new TH1D("unlikeSign", "unlikeSign", 50, 1.6, 2.1);
   likeSignBG      = new TH1D("LikeSignBG", "LikeSignBG", 50, 1.6, 2.1);
   
   invMassPer      = new TH1D("unlikeSign_Peripheral", "unlikeSign_Peripheral", 50, 1.6, 2.1);
   invMassMidCent  = new TH1D("unlikeSign_MidCentral", "unlikeSign_MidCentral", 50, 1.6, 2.1);
   invMassCent     = new TH1D("unlikeSign_Central", "unlikeSign_Central", 50, 1.6, 2.1);
   
   likeSignBGPer     = new TH1D("LikeSignBG_Peripheral", "LikeSignBG_Peripheral", 50, 1.6, 2.1);
   likeSignBGMidCent = new TH1D("LikeSignBG_MidCentral", "LikeSignBG_MidCentral", 50, 1.6, 2.1);
   likeSignBGCent    = new TH1D("LikeSignBG_Central", "LikeSignBG_Central", 50, 1.6, 2.1);
   
   
   if(USE_PT_BINS){
        for(int k = 0; k < NUM_PT_BINS; k++){
            for(int i = 0; i < 3; i++){
        
                str1 =  invMassD0PtBin + binLabelPt[k] + binLabelCentralityClass[i];       
                D0InvMassPtBin[k][i] = new TH1D(str1, str1, 50, 1.6, 2.1);
                str1 =  invMassLSPtBin + binLabelPt[k] + binLabelCentralityClass[i];       
                LSInvMassPtBin[k][i] = new TH1D(str1, str1, 50, 1.6, 2.1);
            }
            
        }
    }        
   
   //Fine pt binning
   if(USE_FINE_PT_BINS){
        for(int k = 0; k < 11; k++){
            for(int i = 0; i < 3; i++){
        
                str1 =  fine + invMassD0PtBin + binLabelFinePt[k] + binLabelCentralityClass[i];       
                D0InvMassFinePtBin[k][i] = new TH1D(str1, str1, 50, 1.6, 2.1);
                str1 =  fine + invMassLSPtBin + binLabelFinePt[k] + binLabelCentralityClass[i];       
                LSInvMassFinePtBin[k][i] = new TH1D(str1, str1, 50, 1.6, 2.1);
            }
            
        }
   }       
   
  
   //QA Histograms
   eventCounter    = new TH1I("number of events used", "number of events used", 6, 0, 6);
   trackCounter    = new TH1I("number of tracks per event", "number of tracks per event", 1500, 0, 1499);
   usedTracks      = new TH1I("Used tracks", "Used tracks", 1500, 0, 1499);
   kaonEtaDist     = new TH1D("Kaon Eta Distribution", "Kaon Eta Distribution", 250, -1, 1);
   pionEtaDist     = new TH1D("Pion Eta Distribution", "Pion Eta Distribution", 250, -1, 1);
   kaonPhiDist     = new TH1D("Kaon Phi Distribution", "Kaon Phi Distribution", 250, -2*TMath::Pi(), 2*TMath::Pi());
   pionPhiDist     = new TH1D("Pion Phi Distribution", "Pion Phi Distribution", 250, -2*TMath::Pi(), 2*TMath::Pi());
   DCAtoPrimaryVertex    = new TH1D("track DCA to PV", "track DCA to PV", 500, 0.0, 10.0);
   DCAtoPrimaryVertexCut = new TH1D("track DCA to PV cut check", "track DCA to PV cut check", 500, 0.0, 10.0);
   hadronPtDist    = new TH1D("Inclusive Hadron pt", "Inclusive Hadron pt", 250, 0, 10);
   hadronPhiDist   = new TH1D("Inclusive Hadron Phi", "Inclusive Hadron Phi", 250, -TMath::Pi(), TMath::Pi());
   hadronEtaDist   = new TH1D("Inclusvie Hadron Eta", "Inclusive Hadron Eta", 250, -1, 1);
   hadronChi2      = new TH1D("Chi2 for hadron tracks", "Chi2 for hadron tracks", 500, 0, 10);
   pVtxX           = new TH1D("X position of pVtx", "X position of pVtx", 500, -10, 10);
   pVtxY           = new TH1D("Y position of pVtx", "Y position of pVtx", 500, -10, 10);
   pVtxZ           = new TH1D("Z position of pVtx", "Z position of pVtx", 500, -7, 7);
   dEdxVsPt        = new TH2D("dEdx_vs_P", "dEdx_vs_P", 250, 0, 10, 250, 0, 10);
   invBetaVsPt     = new TH2D("#Beta^{-1} Vs. P", "#Beta^{-1} Vs. P", 250, 0, 10, 250, 0, 4);
   vZandCentBinPerEvent = new TH2I("event counts per Vz/Cent bin", "event counts per Vz/Cent bin", 16, 0, 16, 10, 0, 10);
   vZandCentBinPerEvent->GetXaxis()->SetTitle("Centrality Bin");
   vZandCentBinPerEvent->GetYaxis()->SetTitle("Vz Bin");
   //QA for mass-cut D0  
   D0ptDist        = new TH1D("D0 Candidate pt Dist (mass cut)", "D0 Candidate pt Dist (mass cut)", 250, 0, 10);
   D0EtaDist       = new TH1D("D0 Eta Dist", "D0 #eta Dist. (mass cut)", 250, -1, 1);
   D0PhiDist       = new TH1D("D0 Phi Dist", "D0 #phi Dist. (mass cut)", 250, -2*TMath::Pi(), 2*TMath::Pi());
   D0KaonPtDist    = new TH1D("D0 Kaon Candidate pt Dist (mass cut)", "D0  Kaon Candidate pt Dist (mass cut)", 250, 0, 10);
   D0KaonEtaDist   = new TH1D("D0 Kaon Eta Dist", "D0 Kaon #eta Dist. (mass cut)", 250, -1, 1);
   D0KaonPhiDist   = new TH1D("D0 Kaon Phi Dist", "D0 Kaon #phi Dist. (mass cut)", 250, -2*TMath::Pi(), 2*TMath::Pi());
   D0PionPtDist    = new TH1D("D0 Candidate Pion pt Dist (mass cut)", "D0 Pion Candidate pt Dist (mass cut)", 250, 0, 10);
   D0PionEtaDist   = new TH1D("D0 Pion Eta Dist", "D0 Pion #eta Dist. (mass cut)", 250, -1, 1);
   D0PionPhiDist   = new TH1D("D0 Pion Phi Dist", "D0 Pion #phi Dist. (mass cut)", 250, -2*TMath::Pi(), 2*TMath::Pi());
   d0CountPerEvent = new TH1I("number of D0 candidates per event", "number of D0 candidates per event", 50, 0, 50);
   histOfCuts      = new TH1D("HistOfCuts", "HistOfCuts", 55, 1, 55); 

   //single particle distributions for errors
   if(SINGLES_DISTS){
        for(int i = 0; i < nVzBins; i++){
            for(int j = 0; j < nCentBins; j++){
        
                str1 = phiD0vsPhiHLabel + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                str2 = etaD0vsEtaHLabel + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                str3 = phiD0vsEtaD0Label + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                str4 = phiHvsEtaHLabel + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
        
                phiD0vsPhiH[i][j]     = new TH2D(str1, str1, 24, -TMath::Pi(), TMath::Pi(), 24, -TMath::Pi(), TMath::Pi());  //d0 on y axis, h on x axis     
                etaD0vsEtaH[i][j]     = new TH2D(str2, str2, 9, -1, 1, 9, -1, 1);  //d0 on y axis, h on x axis
                phiD0vsEtaD0[i][j]    = new TH2D(str3, str3, 9, -1, 1, 24, -TMath::Pi(), TMath::Pi());
                phiHvsEtaH[i][j]      = new TH2D(str4, str4, 9, -1, 1, 24, -TMath::Pi(), TMath::Pi());
            
            }
        }        
   }
   
   //Histogram formatting

   eventCounter->GetXaxis()->SetBinLabel(1,"D0 Cand. events");
   eventCounter->GetXaxis()->SetBinLabel(2,"minBias");
   eventCounter->GetXaxis()->SetBinLabel(3,"total");
   eventCounter->GetXaxis()->SetBinLabel(4,"events from bad runs");
   eventCounter->GetXaxis()->SetBinLabel(5,"Left SB events");
   eventCounter->GetXaxis()->SetBinLabel(6,"Right SB Events");
   
//----------------------End User Variables------------------------------------
  
   //Fill Cuts histogram here for producing output text file
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
    histOfCuts->SetBinContent(49, 1); //Scale factor counter to produce text file
    
    
   
   return kStOK;
}

//-----------------------------------------------------------------------------

//------------------------------------Destructor------------------------------
StPicoD0HadronQAMaker::~StPicoD0HadronQAMaker()
{
   /*  */
}
//-----------------------------------------------------------------------------

//------------------------------------Finish-----------------------------------
Int_t StPicoD0HadronQAMaker::Finish()
{
   LOG_INFO << " StPicoD0HadronQAMaker - writing data and closing output file " <<endm;
   mOutputFile->cd();
   // save user variables here
   
  
   if(USE_PT_BINS){
        for(int k = 0; k < NUM_PT_BINS; k++){
            for(int i = 0; i < 3; i++){
            
                D0InvMassPtBin[k][i]->Write();
                LSInvMassPtBin[k][i]->Write();
            }
        }
   }        
    
    if(USE_FINE_PT_BINS){
        for(int k = 0; k < 11; k++){
            for(int i = 0; i < 3; i++){
            
                D0InvMassFinePtBin[k][i]->Write();
                LSInvMassFinePtBin[k][i]->Write();
            }
        }
   }
    
   if(SINGLES_DISTS){    
        for(int i = 0; i < nVzBins; i++){
            for(int j = 0; j < nCentBins ; j++){
  
                phiD0vsPhiH[i][j]->Write();  
                etaD0vsEtaH[i][j]->Write();
                phiD0vsEtaD0[i][j]->Write();
                phiHvsEtaH[i][j]->Write();
            }
        }
   }        
   
   invMass->Write();
   likeSignBG->Write();
   invMassPer->Write();
   invMassMidCent->Write();
   invMassCent->Write();
   likeSignBGPer->Write();
   likeSignBGMidCent->Write();
   likeSignBGCent->Write();
   
   D0EtaDist->Write();
   D0PhiDist->Write();
   D0ptDist->Write(); 
   D0KaonPtDist->Write();
   D0KaonEtaDist->Write(); 
   D0KaonPhiDist->Write(); 
   D0PionPtDist->Write();
   D0PionEtaDist->Write(); 
   D0PionPhiDist->Write(); 
   eventCounter->Write();
   kaonEtaDist->Write();
   pionEtaDist->Write();
   kaonPhiDist->Write();
   pionPhiDist->Write();
   hadronPtDist->Write();
   hadronPhiDist->Write();
   hadronEtaDist->Write();
   trackCounter->Write();
   dEdxVsPt->Write();
   invBetaVsPt->Write();
   DCAtoPrimaryVertex->Write();
   DCAtoPrimaryVertexCut->Write();
   d0CountPerEvent->Write();
   vZandCentBinPerEvent->Write();
   usedTracks->Write();
   histOfCuts->Write();
   pVtxX->Write();
   pVtxY->Write();
   pVtxZ->Write();
   hadronChi2->Write();
   
   mOutputFile->Close();
  

   return kStOK;
}
//-----------------------------------------------------------------------------

//----------------------------------Make---------------------------------------
Int_t StPicoD0HadronQAMaker::Make(){ //begin Make member function
   readNextEvent();

   if (!mPicoDstMaker)
   {
      LOG_WARN << " StPicoD0HadronQAMaker - No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   StPicoDst const* picoDst = mPicoDstMaker->picoDst();


   if (!picoDst)
   {
      LOG_WARN << "StPicoD0HadronQAMaker - No PicoDst! Skip! " << endm;
      return kStWarn;
   }

   if(mPicoD0Event->runId() != picoDst->event()->runId() ||
       mPicoD0Event->eventId() != picoDst->event()->eventId())
   {
     LOG_ERROR <<" StPicoD0HadronQAMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
     LOG_ERROR <<" StPicoD0HadronQAMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
     exit(1);
   }

   //-------------------Begin User Analysis---------------------------
 
/****************************************************************************************************************/
/****************************GLOBAL AND LOCAL VARIABLES SET AND INITIALIZED BEGIN*********************************/  
/***************************************************************************************************************/ 
    
    
    //-------------- Various global cuts and pt ranges---------------- 
  
    trackCount = 0;
    centralityBin = 0;  //This number will be between 0 and 15 -- 16 bins total
    
    VzBin = 0;
    centralityClass = -1;
    topologicalCutPtBin = -1;
	
    //--------------Local Variables used in the code-------------------
    
    double pt              = 0;
    double phi             = 0;
    double eta             = 0;
    double dEdx            = 0;
    double beta            = 0;
    int    PIDflag         = 0;
    int    numD0s          = 0;
    bool   minBiasFlag     = false;
    int    realTracks      = 0;
    
    double daughterPionPt = 0;
    double daughterKaonPt = 0;
    
    //----------------------------------------------------------------------------------------------------------------
    
    StPicoTrack* trk;
    //StPicoTrack* trk2;
    StThreeVectorF trackMom;
    StThreeVectorF trackMom2;
    StThreeVectorF daughterKaonMom;
    StThreeVectorF daughterPionMom;
    double bField         = picoDst->event()->bField();
    StThreeVectorF pVtx   = picoDst->event()->primaryVertex();
    StThreeVectorF kaonPionMom;
    
    int d0Counter = 0;
    
    bool storedD0Event    = false;
    bool storedLSBEvent   = false;
    bool storedRSBEvent   = false;
    
    std::vector <StThreeVectorF> mAssociatedHadronList;
    
/****************************************************************************************************************/
/****************************GLOBAL AND LOCAL VARIABLES SET AND INITIALIZED END********************************/  
/***************************************************************************************************************/   
    
 
/****************************************************************************************************************/
/****************************PRELIMINARY EVENT CUTS BEGIN********************************************************/  
/***************************************************************************************************************/    
 
    if(!mHFCuts->isGoodRun(picoDst->event())){
      
        eventCounter->Fill(2);
        eventCounter->Fill(3);
        return kStOk; //makes sure the event comes from good run
    }

    eventCounter->Fill(2);
    
    for(unsigned int i = 0; i < 5; i++){
    
        if(picoDst->event()->isTrigger(trigger[i])){
        
            minBiasFlag = true;
            break;
        }
    }        
    
    if(!minBiasFlag) {return kStOK;}
    
    eventCounter->Fill(1);

    TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();    //The Kaon-Pion list is generated here

    trackCount = picoDst->numberOfTracks();            
    trackCounter->Fill(picoDst->numberOfTracks());
   
/****************************************************************************************************************/
/****************************PRELIMINARY EVENT CUTS END********************************************************/  
/***************************************************************************************************************/ 
  
   
/****************************************************************************************************************/  
/********************QA Loop for storing all of the basic information about the events and tracks BEGIN**********/
/****************************************************************************************************************/   
   
    for(unsigned int i = 0; i < picoDst->numberOfTracks(); ++i){ // Begin loop to fill basic track information and make nega/posi list
                                                                 //gets pt, eta, phi, dEdx and beta from TOF
        trk = picoDst->track(i);
        trackMom = trk->gMom(pVtx, bField);
        
        if(!mHFCuts->isGoodTrack(trk)) { continue; }    //Checks the HFT requirement for a track
        
        trackDCA = ((trk->helix().origin())-pVtx).mag();
        DCAtoPrimaryVertex->Fill(trackDCA);            //QA plot to see full DCA distribution
        
        if(trk->chi2() > trackChi2max) { continue; }  //check chi2 cut
        if(!checkDCAtoPV(trackDCA))  { continue; }   // track quality cut for DCA to PV
        
        pt   = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y()));
        if(pt < hadronPtMin || pt > hadronPtMax){continue;}              //check pt cut
        
        if(trackMom.mag() > .20 && trackMom.mag() < .45 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }  //remove electron contamination
        if(trackMom.mag() > .70 && trackMom.mag() < .80 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }
        
        eta  = trackMom.pseudoRapidity();
		if(eta > 1 || eta < -1) { continue; }
		
        realTracks++;                        //after a track passes all cuts, it is considered a "real track"
        
        dEdx = trk->dEdx();
        
        phi = trackMom.phi();
        
        dEdxVsPt->Fill(trackMom.mag(), dEdx);
        DCAtoPrimaryVertexCut->Fill(trackDCA);
        hadronChi2->Fill(trk->chi2());
        
        if(mHFCuts->hasTofPid(trk)){    
         
            beta = mHFCuts->getTofBeta(trk);                                                
            invBetaVsPt->Fill(trackMom.mag(), (1/beta));        
        }
        
        hadronPtDist->Fill(pt); 
        hadronPhiDist->Fill(phi);        
        hadronEtaDist->Fill(eta);
        
        if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaKaon()) < 2.0)){       //ONLY check nSigma for TPC track -- used to be mHFCuts->isTPCKaon(trk), but this include pt cut
                                                                                  //need to fix this. Need to figure out how to access the nSigma from the cuts
              kaonEtaDist->Fill(eta);
              kaonPhiDist->Fill(phi);
             
        }

        else if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaPion()) < 3.0)){     //ONLY check nSigma for TPC track
               
              pionEtaDist->Fill(eta);
              pionPhiDist->Fill(phi);
              
        }
        
        mAssociatedHadronList.push_back(trackMom); //store the associated tracks to a list for quicker pairing later
        
    } //End loop to fill basic track information and make 
    
    //Fill track distributions and get centrality bins////////////////////
    
        usedTracks->Fill(realTracks);   
    
        centralityBin = getCentralityBin(realTracks);  //get multiplicity bin
        centralityClass = getCentralityClass(realTracks);
        
        if(USE_VZ_BINS){
            Vz = picoDst->event()->primaryVertex().z();
            VzBin = getVzBin(Vz);                  //get Vz bin
        }
   
        else VzBin = 0;
    
        if(centralityBin == -1 || VzBin == -1) { return kStOk; }
        
       
    
        pVtxX->Fill(picoDst->event()->primaryVertex().x());
        pVtxY->Fill(picoDst->event()->primaryVertex().y());
        pVtxZ->Fill(picoDst->event()->primaryVertex().z()); 
   
        vZandCentBinPerEvent->Fill(centralityBin, VzBin); 
    
    ///////////////////////////////////////////////////////////////////
    
    if(DEBUG){ 
        
                cout << endl << endl;
                cout << "*********************EVENT START******************" << endl;
                cout << "We are on event # " << eventNumber << endl;
                cout << "This event has " << realTracks << " tracks." << endl;
                if(USE_VZ_BINS){ cout << "Vz: " << Vz << "   VzBin: " << VzBin << "    Centrality Bin: " << centralityBin << endl; }
                else cout << "Centrality Bin: " << centralityBin << endl;
                cout << endl;
                cout << endl;
                eventNumber++;
    }         
           
        
/****************************************************************************************************************/  
/********************QA Loop for storing all of the basic information about the events and hadron tracks END**********/
/****************************************************************************************************************/         
 

//***************************************************************************************
//BEGIN D0 BLOCK OF CODE HERE -- MAKE SURE THE NON_IDENTIFIED FLAG IS SET TO FALSE!!!!
//***************************************************************************************    

//Still need to add checks for TOF information-----   

 
   
/****************************************************************************************************************/   
/*****************************************BEGIN MAIN D0 LOOP************************************************/
/****************************************************************************************************************/
    if(DEBUG) { cout << "   Begin D0 QA code block" << endl << endl; }
   
    for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx){ // begin main D0 loop
     
        StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
    
        StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
        StPicoTrack const* pion = picoDst->track(kp->pionIdx());
      
        topologicalCutPtBin = getTopologicalCutPtBin(kp->pt());
      
        if(topologicalCutPtBin == -1 ) {continue;}
      
       
        if(!cutCheck(kp, d0PtLow,  d0PtHigh,  d0TopologicalCutArray[topologicalCutPtBin][0],  d0DecayLengthMax,  d0TopologicalCutArray[topologicalCutPtBin][1],
                     d0DaughterPionPtMin, d0DaughterKaonPtMin, d0TopologicalCutArray[topologicalCutPtBin][4], d0TopologicalCutArray[topologicalCutPtBin][3],
                     d0TopologicalCutArray[topologicalCutPtBin][2])) { continue; }    
      
        ptBin     = getPtBin(kp->pt());
        finePtBin = getFinePtBin(kp->pt());
         
        
        if(kaon->charge()*pion->charge() < 0){// begin Unlike-sign conditional 
	      
            invMass->Fill(kp->m()); 
            if (centralityClass == 0){ 
                invMassPer->Fill(kp->m());
                if(USE_FINE_PT_BINS) { D0InvMassFinePtBin[finePtBin][centralityClass]->Fill(kp->m()); }
            }
                                              
            else if(centralityClass == 1){ 
                invMassMidCent->Fill(kp->m());
                if(USE_FINE_PT_BINS) { D0InvMassFinePtBin[finePtBin][centralityClass]->Fill(kp->m()); }
            }
            
            else if(centralityClass == 2){ 
                invMassCent->Fill(kp->m());
                if(USE_FINE_PT_BINS) { D0InvMassFinePtBin[finePtBin][centralityClass]->Fill(kp->m()); }
            }              
            if(ptBin > -1) { 
            
                D0InvMassPtBin[ptBin][centralityClass]->Fill(kp->m()); 
                    
            }           
            
	    }//end Unlike-sign conditional 
      
	    if(kaon->charge()*pion->charge() > 0){//begin Like-sign conditional 
          
            likeSignBG->Fill(kp->m()); 
            if(centralityClass == 0) {
                likeSignBGPer->Fill(kp->m()); 
                if(USE_FINE_PT_BINS) { LSInvMassFinePtBin[finePtBin][centralityClass]->Fill(kp->m()); }
            }  
            else if(centralityClass == 1){ 
                likeSignBGMidCent->Fill(kp->m());
                if(USE_FINE_PT_BINS) { LSInvMassFinePtBin[finePtBin][centralityClass]->Fill(kp->m()); }
            }
            else if(centralityClass == 2){ 
                likeSignBGCent->Fill(kp->m()); 
                if(USE_FINE_PT_BINS) { LSInvMassFinePtBin[finePtBin][centralityClass]->Fill(kp->m()); }  
            }                
            if(ptBin > -1) { 
            
                LSInvMassPtBin[ptBin][centralityClass]->Fill(kp->m());

            }
          
        }//end Like-sign conditional 
 
        if(kaon->charge()*pion->charge() < 0){// begin US conditional
                
            if(kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh){ // begin loop over both unlike-sign AND LIKE SIGN D0 candidates 

                if(DEBUG) { cout << "   D0 Candidate Mass: " << kp->m() << "    PtBin : " << ptBin << endl; }
            
                if(!storedD0Event){
                
                    eventCounter->Fill(0);
                    storedD0Event = true;
                }
   
                d0Counter = d0Counter + 1;
                D0ptDist->Fill(kp->pt()); 
                D0EtaDist->Fill(kp->eta());
                D0PhiDist->Fill(kp->phi());
                
                daughterKaonMom = kaon->gMom(pVtx, bField);
                daughterPionMom = pion->gMom(pVtx, bField);
                
                D0KaonPtDist->Fill(daughterKaonMom.perp());
                D0KaonEtaDist->Fill(daughterKaonMom.pseudoRapidity()); 
                D0KaonPhiDist->Fill(daughterKaonMom.phi()); 
                D0PionPtDist->Fill(daughterPionMom.perp());
                D0PionEtaDist->Fill(daughterPionMom.pseudoRapidity()); 
                D0PionPhiDist->Fill(daughterPionMom.phi()); 
                
                if(SINGLES_DISTS){
                    for(unsigned int i = 0; i < mAssociatedHadronList.size(); i++){ // begin loop for single part distributions on D0/h+/-
           
                        if(i == kp->kaonIdx() || i == kp->pionIdx()) { continue; }   
                    
                        trackMom = mAssociatedHadronList[i];
        
                        phi  = TMath::ATan2(trackMom.y(),trackMom.x());
                        eta = trackMom.pseudoRapidity();
                    
                        phiD0vsPhiH[VzBin][centralityBin]->Fill(phi, kp->phi());
                        etaD0vsEtaH[VzBin][centralityBin]->Fill(eta, kp->eta());     
                        phiD0vsEtaD0[VzBin][centralityBin]->Fill(kp->eta(), kp->phi());    
                    
                    }
                }
            }
                
            else if(kp->m() > USSideBandLeftLow &&  kp->m() < USSideBandLeftHigh){ //Side Band Left
                
                if(!storedLSBEvent){
                
                    eventCounter->Fill(4);
                    storedLSBEvent = true;
                }
    
                if(DEBUG) { cout << "   SideBandLeft: " << kp->m() << "    PtBin : " << ptBin << endl;}    
                
            }
            
            else if(kp->m() > USSideBandRightLow &&  kp->m() < USSideBandRightHigh){ //Side Band Right
                 
                if(!storedRSBEvent){
                
                    eventCounter->Fill(5);
                    storedRSBEvent = true;
                }
                
                if(DEBUG) { cout << "   SideBandRight: " << kp->m() << "    PtBin : " << ptBin << endl;}    
                
            }
        }  //end US conditional  
            
        else if(kaon->charge()*pion->charge() > 0 && kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh) {// like sign conditional 

            if(DEBUG) { cout << "   LS Mass: " << kp->m() << "    PtBin : " << ptBin << endl; }
             
        }          
    }// End main D0 loop

    
    d0CountPerEvent->Fill(d0Counter);
       
    //-------------------End Current Event Analysis--------------------------------
   
   return kStOK;

}//end Make member function

//---------------------User Functions-------------------------------------

bool StPicoD0HadronQAMaker::isGoodPair(StKaonPion const* const kp) const
{
  if(!kp) return false;

  StPicoTrack const* kaon = mPicoDstMaker->picoDst()->track(kp->kaonIdx());
  StPicoTrack const* pion = mPicoDstMaker->picoDst()->track(kp->pionIdx());

  //  To be replaced by mHFCuts->isGoodSecondaryVertexPair(kp))
  bool pairCuts = kp->m() > mHFCuts->cutSecondaryPairMassMin() && 
    kp->m() < mHFCuts->cutSecondaryPairMassMax() &&
    std::cos(kp->pointingAngle()) > mHFCuts->cutSecondaryPairCosThetaMin() &&
    kp->decayLength()  > mHFCuts->cutSecondaryPairDecayLengthMin() && 
    kp->decayLength()  < mHFCuts->cutSecondaryPairDecayLengthMax() &&
    kp->dcaDaughters() < mHFCuts->cutSecondaryPairDcaDaughtersMax();

  return (mHFCuts->isGoodTrack(kaon) && mHFCuts->isGoodTrack(pion) &&
	  mHFCuts->isTPCKaon(kaon) && mHFCuts->isTPCPion(pion) && 
	  pairCuts);

}


bool StPicoD0HadronQAMaker::cutCheck(StKaonPion const* const kp, double ptMin, double ptMax, double decayLengthMin, double decayLengthMax, 
                                                            double dcaDaughters, double kaonPtCut, double pionPtCut, 
                                                            double dcaKaon, double dcaPion, double dcaV0toPV) const
{
  if(!kp) return false;

  StPicoTrack const* kaon = mPicoDstMaker->picoDst()->track(kp->kaonIdx());
  StPicoTrack const* pion = mPicoDstMaker->picoDst()->track(kp->pionIdx());

  
      bool truthCuts = kp->pt() > ptMin && kp->pt() < ptMax &&
                       kp->decayLength() > decayLengthMin && 
                       kp->decayLength() < decayLengthMax &&
                       kp->dcaDaughters() < dcaDaughters  &&
                       kaon->gPt() > kaonPtCut && pion->gPt() > pionPtCut &&
                       kp->kaonDca() > dcaKaon && kp->pionDca() > dcaPion &&
                       kp->perpDcaToVtx() < dcaV0toPV;
  
      return truthCuts;

}


int StPicoD0HadronQAMaker::getCentralityBin(int nTracks){

    /*if(nTracks >= 2   && nTracks < 16)  { return 0;  }
    if(nTracks >= 16  && nTracks < 31)  { return 1;  }
    if(nTracks >= 31  && nTracks < 51)  { return 2;  }
    if(nTracks >= 51  && nTracks < 79)  { return 3;  }
    if(nTracks >= 79  && nTracks < 116) { return 4;  }
    if(nTracks >= 116 && nTracks < 164) { return 5;  }
    if(nTracks >= 164 && nTracks < 224) { return 6;  }
    if(nTracks >= 224 && nTracks < 295) { return 7;  }
    if(nTracks >= 295 && nTracks < 335) { return 8;  }
    if(nTracks >= 335 && nTracks < 376) { return 9;  }
    if(nTracks >= 376 && nTracks < 420) { return 10; }
    if(nTracks >= 420 && nTracks < 470) { return 11; }
    if(nTracks >= 470 && nTracks < 550) { return 12; }
    if(nTracks >= 550 && nTracks < 620) { return 13; }
    if(nTracks >= 620)                  { return 14; }*/
    
    if(nTracks >= 1   && nTracks < 42)  { return 0;  }
    if(nTracks >= 42  && nTracks < 86)  { return 1;  }
    if(nTracks >= 86  && nTracks < 131)  { return 2;  }
    if(nTracks >= 131  && nTracks < 183) { return 3;  }
    if(nTracks >= 183 && nTracks < 235) { return 4;  }
    if(nTracks >= 235 && nTracks < 288) { return 5;  }
    if(nTracks >= 288 && nTracks < 340) { return 6;  }
    if(nTracks >= 340 && nTracks < 392) { return 7;  }
    if(nTracks >= 392 && nTracks < 440) { return 8;  }
    if(nTracks >= 440 && nTracks < 491) { return 9;  }
    if(nTracks >= 491 && nTracks < 542) { return 10;  }
    if(nTracks >= 542 && nTracks < 593) { return 11;  }
    if(nTracks >= 593 && nTracks < 644) { return 12;  }
    if(nTracks >= 644 && nTracks < 695) { return 13;  }
    if(nTracks >= 695 && nTracks < 746) { return 14;  }
    if(nTracks >= 746)                  { return 15; }

    else return -1;

}    

int StPicoD0HadronQAMaker::getVzBin(double Vz){

    if(Vz >= -6.0  && Vz < -4.8)  { return 0;  }   //bin 0: -6 to -4.8
    if(Vz >= -4.8  && Vz < -3.6)  { return 1;  }      //bin 1: -4.8 to -3.6
    if(Vz >= -3.6  && Vz < -2.4)  { return 2;  }          //bin 2: -3.6 to -2.4
    if(Vz >= -2.4  && Vz < -1.2)  { return 3;  }      //bin 3: -2.4 to -1.2
    if(Vz >= -1.2  && Vz < 0)     { return 4;  }      //bin 4: -1.2 to 0
    if(Vz >= 0     && Vz < 1.2)   { return 5;  }
    if(Vz >= 1.2   && Vz < 2.4)   { return 6;  }
    if(Vz >= 2.4   && Vz < 3.6)   { return 7;  }
    if(Vz >= 3.6   && Vz < 4.8)   { return 8;  }
    if(Vz >= 4.8   && Vz < 6.0)   { return 9;  }
    

    else return -1;

}    

int StPicoD0HadronQAMaker::getPtBin(double pt){

    if(pt >=  0.0  && pt <  1.0)    { return 0;  }   
    if(pt >=  1.0  && pt <  4.0)    { return 1;  }
    if(pt >=  4.0  && pt <  10.0)   { return 2;  }        
    //if(pt >=  3.0  && pt <  4.0)   { return 3;  }
    //if(pt >=  4.0  && pt <  5.0)   { return 4;  }
    //if(pt >=  5.0  && pt <  10.0)  { return 5;  }
    
    else return -1;

}    

int StPicoD0HadronQAMaker::getFinePtBin(double pt){

    if(pt >=  0.0  && pt <  1.0)    { return 0;  }   
    if(pt >=  1.0  && pt <  2.0)    { return 1;  }
    if(pt >=  2.0  && pt <  3.0)    { return 2;  }
    if(pt >=  3.0  && pt <  4.0)    { return 3;  }
    if(pt >=  4.0  && pt <  5.0)    { return 4;  }
    if(pt >=  5.0  && pt <  6.0)    { return 5;  }
    if(pt >=  6.0  && pt <  7.0)    { return 6;  }
    if(pt >=  7.0  && pt <  8.0)    { return 7;  }
    if(pt >=  8.0  && pt <  9.0)    { return 8;  }
    if(pt >=  9.0  && pt < 10.0)    { return 9;  }
    if(pt >=  10.0)                 { return 10;  }
    
    
    else return -1;

}    

int StPicoD0HadronQAMaker::getTopologicalCutPtBin(double pt){

    if(pt >=  0.0  && pt <  1.0)    { return 0;  }   
    if(pt >=  1.0  && pt <  2.0)    { return 1;  }
    if(pt >=  2.0  && pt <  3.0)    { return 2;  }
    if(pt >=  3.0  && pt <  5.0)    { return 3;  }
    if(pt >=  5.0  && pt <  10.0)   { return 4;  }
    
    else return -1;

}    

int StPicoD0HadronQAMaker::getCentralityClass(double nTracks){

    if(nTracks >= 1   && nTracks < 131)  { return 0;  }  //peripheral
    if(nTracks >= 131  && nTracks < 440) { return 1;  }  //midcentral
    if(nTracks >= 440)                   { return 2;  }  //central
    
    else return -1;

}    


bool StPicoD0HadronQAMaker::checkDCAtoPV(float trackDCA){

     return (trackDCA <= trackDCAtoPvtx);
     
}



