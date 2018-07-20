#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"

#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "StPicoD0EventMaker/StKaonPion.h"

#include "StMixedEventBuffer/StMixedEventBuffer.h"
#include "StMixedEventBuffer/StMixerEvent.h"
#include "StMixedEventBuffer/StMixerTrack.h"

#include "StPicoD0AnaMaker.h"
#include "StPicoHFMaker/StHFCuts.h"


/****

Author: Alex Jentsch

Current version of AnaMaker (by date): 2/22/2016

Description of current functionality:

1. Can read in both Trees and PicoDst -> Stores all basic QA info for PicoDsts
2. Can make invMass Histos for various pt bins
3. Can produce sibling histograms for various pt bins


Still Needs:

1. Various BG production methods (side band, mixing, etc.)


Update (3/10/2016)

Currently adding event mixing to code. Using the StPicoMixedEventMaker as a baseline, but changing its implementation significantly.
I want to be able to handle all of the mixing in my AnaMaker with all of the mixing functions being used as if from any other normal class.


Update (3/31/2016)

Buffering code built and seems to function. Lots of testing left to do. Need to make sure the stored event information can be retrieved 
and can be binned into histograms.


Update (4/19/2016)

Event mixing is working and the buffering is not causing any problems. May be changing the event mixing algorithm to better the statistics.

Update (6/3/2016)

Code is fully functional for both non-identified and D0-Hadron correlations. Still working on getting all of the cuts right.
Code has been updated for the new code structure with the updated production.

****/

ClassImp(StPicoD0AnaMaker)

StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name, char const * inputFilesList, 
                                   char const * outName, StPicoDstMaker* picoDstMaker): 
                            StMaker(name),mPicoDstMaker(picoDstMaker),mPicoD0Event(NULL),
                            mOutFileName(outName), mInputFileList(inputFilesList),
                            mOutputFile(NULL), mChain(NULL), mEventCounter(0), mHFCuts(NULL)
{}

//---------------------------------------Initialization-----------------------------------------
Int_t StPicoD0AnaMaker::Init()
{
   mPicoD0Event = new StPicoD0Event();

   mChain = new TChain("T");
   std::ifstream listOfFiles(mInputFileList.Data());
   if (listOfFiles.is_open())
   {
      std::string file;
      while (getline(listOfFiles, file))
      {
         LOG_INFO << "StPicoD0AnaMaker - Adding :" << file << endm;
         mChain->Add(file.c_str());
      }
   }
   else
   {
      LOG_ERROR << "StPicoD0AnaMaker - Could not open list of files. ABORT!" << endm;
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
   DEBUG = false;
   DEBUG_MIX_BUFFER = false;
   USE_CENT_BINS = true;
   USE_VZ_BINS = true;
   D0_HADRON_CORR = false;
   NON_IDENTIFIED_CORR = false;
   //----------------------------------------------------------//
   
   //---------------------Important constants-----------------//
   BUFFER_SIZE   = 5;
   NUM_PHI_BINS  = 24;
   NUM_ETA_BINS  = 25;
   NUM_VZ_BINS   = 10;
   NUM_CENT_BINS = 11; 
   
   
   // --------------------Event Mixer Buffer-------------------------------------
   
    if(USE_VZ_BINS) { nVzBins = NUM_VZ_BINS; }  //flag to set binning on Vz
    else nVzBins = 1;
    
    if(USE_CENT_BINS) { nCentBins = NUM_CENT_BINS; }
    else nCentBins = 1;
    
    //********************D0_Hadron section************************//
    
    if(D0_HADRON_CORR){
    
        for(int i = 0; i < nCentBins; i++){
            for( int j = 0; j < nVzBins; j++){
        
                eventBufferUS[j][i] = new StMixedEventBuffer();
                eventBufferUS[j][i]->setBufferSize(BUFFER_SIZE);            //set buffer size here -- the amount of events in each 2d bin
                eventBufferUS[j][i]->setBufferCounter(0);
            }
        }
    }        
   
    if(D0_HADRON_CORR){
   
        for(int i = 0; i < nCentBins; i++){
            for( int j = 0; j < nVzBins; j++){
        
                eventBufferLS[j][i] = new StMixedEventBuffer();
                eventBufferLS[j][i]->setBufferSize(BUFFER_SIZE);            //set buffer size here -- the amount of events in each 2d bin
                eventBufferLS[j][i]->setBufferCounter(0);
            }
        }
    }       
   
    //********************Unidentified Hadron section************************//
   
    if(NON_IDENTIFIED_CORR){
    
        for(int i = 0; i < nCentBins; i++){
            for( int j = 0; j < nVzBins; j++){
        
                eventBuffer[j][i] = new StMixedEventBuffer();
                eventBuffer[j][i]->setBufferSize(BUFFER_SIZE);            //set buffer size here -- the amount of events in each 2d bin
                eventBuffer[j][i]->setBufferCounter(0);
            }
        }
    }  
    
    //--------------------CUTS------------------------------------------
   
    kaonPtCut      = .15;
    pionPtCut      = .15;
    
    D0InvMassLow   = 1.82;
    D0InvMassHigh  = 1.90;
    
    d0PtLow             = .15;
    d0PtHigh            = 20.0;
    d0DecayLengthMin    = .0200;
    d0DecayLengthMax    = 999999.0;
    daughterDCA         = .0055;
    d0DaughterPionPtMin = 1.2;
    d0DaughterKaonPtMin = 1.2;
    kaonDCA             = .008;
    pionDCA             = .008;
    d0DCAtoPV           = .0065;
    
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
    
    
    
    //Labels for D0-hadron US and LS corr histograms
    
    TString sibCorrUSLabel1 = "Sibling_US_correlation";
    TString mixCorrUSLabel1 = "Mixed_US_correlation";
    TString sibCorrLSLabel1 = "Sibling_LS_correlation";
    TString mixCorrLSLabel1 = "Mixed_LS_correlation";
    
    //Labels for unidentified hadron corr histograms
    
    TString sibCorrPlusPlusLabel1   = "Sibling_++_correlation";
    TString sibCorrPlusMinusLabel1  = "Sibling_+-_correlation";
    TString sibCorrMinusPlusLabel1  = "Sibling_-+_correlation";
    TString sibCorrMinusMinusLabel1 = "Sibling_--_correlation";
    
    TString mixCorrPlusPlusLabel1   = "Mixed_++_correlation";
    TString mixCorrPlusMinusLabel1  = "Mixed_+-_correlation";
    TString mixCorrMinusPlusLabel1  = "Mixed_-+_correlation";
    TString mixCorrMinusMinusLabel1 = "Mixed_--_correlation";
    
    //other labels
   
    TString eventCounterLabel = "Event Count Vz ";
    TString etaLabel          = "Inclusive Hadron Eta";
    TString phiLabel          = "Inclusive Hadron Phi";
    
    TString str1;
    TString str2;
    TString str3;
    TString str4;
    
    TString binLabelVz[10] = {"0", "1", "2", "3", "4", "5", "6", "7", "8","9"};
    TString binLabelCent[11] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
   
   // --------------------Begin User Variables-----------------------------------
   
   //ptDist          = new TH1D("Pt Distribution", "Pt Distribution", 1000, 0, 10);              
   invMass         = new TH1D("unlikeSign", "unlikeSign", 50, 1.6, 2.1);
   kaonDist        = new TH1D("Kaon Distribution", "Kaon Distribution", 500, 0 , 500);
   pionDist        = new TH1D("Pion Distribution", "Pion Distribution", 2000, 0 , 2000);
   likeSignBG      = new TH1D("LikeSignBG", "LikeSignBG", 50, 1.6, 2.1);
   invMassMinusBG  = new TH1D("D0MinusLSBG", "D0MinusLSBG", 50, 1.6, 2.1);
   
   
   
   if(D0_HADRON_CORR){//begin D0-Hadron Correlation Histograms
        for(int i = 0; i < nVzBins; i++){ //Initialize all of the histograms for storing the sibling and mixed information
            for(int j = 0; j < nCentBins; j++){
        
                if(USE_VZ_BINS){
                
                    str1 = sibCorrUSLabel1 + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                    str2 = mixCorrUSLabel1 + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                    str3 = sibCorrLSLabel1 + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                    str4 = mixCorrLSLabel1 + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                } 
             
                else { 
                   
                    str1 = sibCorrUSLabel1 + CentBinLabel + binLabelCent[j];
                    str2 = mixCorrUSLabel1 + CentBinLabel + binLabelCent[j];
                    str3 = sibCorrLSLabel1 + CentBinLabel + binLabelCent[j];
                    str4 = mixCorrLSLabel1 + CentBinLabel + binLabelCent[j]; 
                } 
             
                sibCorrBinUS[i][j]       = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                mixCorrBinUS[i][j]       = new TH2D(str2, str2, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                sibCorrBinLS[i][j]       = new TH2D(str3, str3, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                mixCorrBinLS[i][j]       = new TH2D(str4, str4, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            }
        }    
   
        for(int i = 0; i < nVzBins; i++){
    
            if(USE_VZ_BINS) { str1 = eventCounterLabel + binLabelVz[i]; }
        
            else str1 = "event count per centrality bin (for events containing D0s)";
        
            eventCategoryCounter[i] = new TH1D(str1, str1, 11, 0, 11);
            eventCategoryCounter[i]->GetXaxis()->SetTitle("CentralityBin");
        }        
    }//end D0-Hadron Correlation Histograms   
   
    
    
    if(NON_IDENTIFIED_CORR){//begin unidentified Correlation Histograms
        for(int i = 0; i < nVzBins; i++){ //Initialize all of the histograms for storing the sibling and mixed information
            for(int j = 0; j < nCentBins; j++){
        
                str1 = sibCorrPlusPlusLabel1 + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                str2 = sibCorrPlusMinusLabel1 + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                str3 = sibCorrMinusPlusLabel1 + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                str4 = sibCorrMinusMinusLabel1 + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                
                sibCorrPlusPlus[i][j]    = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                sibCorrPlusMinus[i][j]   = new TH2D(str2, str2, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                sibCorrMinusPlus[i][j]   = new TH2D(str3, str3, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                sibCorrMinusMinus[i][j]  = new TH2D(str4, str4, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                
                str1 = mixCorrPlusPlusLabel1 + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                str2 = mixCorrPlusMinusLabel1 + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                str3 = mixCorrMinusPlusLabel1 + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                str4 = mixCorrMinusMinusLabel1 + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                
                mixCorrPlusPlus[i][j]    = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                mixCorrPlusMinus[i][j]   = new TH2D(str2, str2, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                mixCorrMinusPlus[i][j]   = new TH2D(str3, str3, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                mixCorrMinusMinus[i][j]  = new TH2D(str4, str4, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            }
        }    
    }//end unidentified Correlation Histograms
 
 
    if(USE_VZ_BINS){
        for(int i = 0; i < nVzBins; i++){ 
            
            str1 = phiLabel + VzBinLabel + binLabelVz[i];
            str2 = etaLabel + VzBinLabel + binLabelVz[i];
            phiDistVz[i]   = new TH1D(str1, str1, 1000, -2*TMath::Pi(), 2*TMath::Pi());
            etaDistVz[i]   = new TH1D(str2, str2, 1000, -2, 2);
               
        }    
    }
 
   
   //QA Histograms
   eventCounter    = new TH1D("number of events used", "number of events used", 4, 0, 4);
   trackCounter    = new TH1D("number of tracks per event", "number of tracks per event", 2000, 0, 1999);
   kaonPtDist      = new TH1D("Kaon Pt Distribution", "Kaon Pt Distribution", 1000, 0, 5);
   pionPtDist      = new TH1D("Pion Pt Distribution", "Pion Pt Distribution", 1000, 0, 5);
   kaonEtaDist     = new TH1D("Kaon Eta Distribution", "Kaon Eta Distribution", 1000, -1, 1);
   pionEtaDist     = new TH1D("Pion Eta Distribution", "Pion Eta Distribution", 1000, -1, 1);
   kaonPhiDist     = new TH1D("Kaon Phi Distribution", "Kaon Phi Distribution", 1000, -2*TMath::Pi(), 2*TMath::Pi());
   pionPhiDist     = new TH1D("Pion Phi Distribution", "Pion Phi Distribution", 1000, -2*TMath::Pi(), 2*TMath::Pi());
   kaonDCAprimary  = new TH1D("DCA kaons from primary", "DCA kaons from primary", 500, 0.0, 0.5);
   pionDCAprimary  = new TH1D("DCA pions from primary", "DCA pions from primary", 500, 0.0, 0.5);
   
   hadronPtDist    = new TH1D("Inclusive Hadron pt", "Inclusive Hadron pt", 1000, 0, 10);
   hadronPhiDist   = new TH1D("Inclusive Hadron Phi", "Inclusive Hadron Phi", 1000, -2*TMath::Pi(), 2*TMath::Pi());
   hadronEtaDist   = new TH1D("Inclusvie Hadron Eta", "Inclusive Hadron Eta", 1000, -2, 2);
   
   kaonDCAfromD0   = new TH1D("DCA for kaons from D0", "DCA for kaons from D0", 500, 0.0, 0.5);
   pionDCAfromD0   = new TH1D("DCA for pions from D0", "DCA for pions from D0", 500, 0.0, 0.65);
   decayLengthQA   = new TH1D("D0 Candidate Decay Length (no mass cut)", "D0 Candidate Decay Length (no mass cut)", 500, 0.0, 1.5);
   pointingAngleQA = new TH1D("D0 Candidate Pointing Angle(no mass cut)", "D0 Candidate Pointing Angle (no mass cut)", 500, 0.0, 1.7);
   daughterDCAQA   = new TH1D("D0 Daughter DCA", "D0 Daughter DCA (no mass cut)", 500, 0.0, .01);

   dEdxVsPt        = new TH2D("dEdx_vs_Pt", "dEdx_vs_Pt", 500, 0, 10, 500, 0, 15);
   invBetaVsPt     = new TH2D("#Beta^{-1} Vs. Pt", "#Beta^{-1} Vs. Pt", 500, 0, 10, 500, 0, 4);
   
   //QA for mass-cut D0  
   D0ptDist        = new TH1D("D0 Candidate pt Dist (mass cut)", "D0 Candidate pt Dist (mass cut)", 1000, 0, 10);
   D0EtaDist       = new TH1D("D0 Eta Dist", "D0 #eta Dist. (mass cut)", 500, -1, 1);
   D0PhiDist       = new TH1D("D0 Phi Dist", "D0 #phi Dist. (mass cut)", 500, -2*TMath::Pi(), 2*TMath::Pi());
   D0PeakPlusBG    = new TH1D("D0PeakPlusBG", "D0PeakPlusBG", 50, 1.6, 2.1);
   D0LikeSignBG    = new TH1D("LikeSign peak range", "LikeSign peak range", 50, 1.6, 2.1);
   D0PeakMinusBG   = new TH1D("D0Peak", "D0Peak", 50, 1.6, 2.1);
   d0CountPerEvent = new TH1I("number of D0 candidates per event", "number of D0 candidates per event", 50, 0, 50);
   
   //Histogram formatting

   eventCounter->GetXaxis()->SetBinLabel(1,"minBias+D0 Cand.");
   eventCounter->GetXaxis()->SetBinLabel(2,"minBias");
   eventCounter->GetXaxis()->SetBinLabel(3,"total");
   eventCounter->GetXaxis()->SetBinLabel(4,"events from bad runs");
   
//----------------------End User Variables------------------------------------
   
   
   return kStOK;
}

//-----------------------------------------------------------------------------

//------------------------------------Destructor------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker()
{
   /*  */
}
//-----------------------------------------------------------------------------

//------------------------------------Finish-----------------------------------
Int_t StPicoD0AnaMaker::Finish()
{
   LOG_INFO << " StPicoD0AnaMaker - writing data and closing output file " <<endm;
   mOutputFile->cd();
   // save user variables here
   
   //StMixerEvent* event;
   //StMixerTrack mixedTrack;
   //StMixerTrack kaonPionTrack;
   //int nTracks = 0;
    
    if(D0_HADRON_CORR){
        for(int i = 0; i < nVzBins; i++){
            for(int j = 0; j < nCentBins; j++){
        
            sibCorrBinUS[i][j]->Write();
            mixCorrBinUS[i][j]->Write();
            sibCorrBinLS[i][j]->Write();
            mixCorrBinLS[i][j]->Write();
            
            if(DEBUG){ cout << "deleting buffer[" << i << "][" << j << "]" << endl;}
            
            delete eventBufferUS[i][j];
            delete eventBufferLS[i][j];
            
            }
        }
    }        
    
    if(NON_IDENTIFIED_CORR){
        for(int i = 0; i < nVzBins; i++){
            for(int j = 0; j < nCentBins; j++){
        
                sibCorrPlusPlus[i][j]->Write();
                sibCorrPlusMinus[i][j]->Write();
                sibCorrMinusPlus[i][j]->Write();
                sibCorrMinusMinus[i][j]->Write();
                
                mixCorrPlusPlus[i][j]->Write();
                mixCorrPlusMinus[i][j]->Write();
                mixCorrMinusPlus[i][j]->Write();
                mixCorrMinusMinus[i][j]->Write();
            
            if(DEBUG){ cout << "deleting buffer[" << i << "][" << j << "]" << endl;}
            
            delete eventBuffer[i][j];
            
            }
        }
    }        
    
    if(USE_VZ_BINS){
        for(int i = 0; i < nVzBins; i++){
           
                etaDistVz[i]->Write();
                phiDistVz[i]->Write();
            
        }
    }        
  
   //ptDist->Write();
   invMass->Write();
   kaonDist->Write();
   pionDist->Write();
   likeSignBG->Write();
   
   invMassMinusBG->Write();
   D0EtaDist->Write();
   D0PhiDist->Write();
   D0PeakMinusBG->Write();
   eventCounter->Write();
   kaonDCAfromD0->Write();
   pionDCAfromD0->Write();
   decayLengthQA->Write();
   pointingAngleQA->Write();
   daughterDCAQA->Write();
   D0ptDist->Write();   
   kaonPtDist->Write();
   pionPtDist->Write();
   kaonEtaDist->Write();
   pionEtaDist->Write();
   kaonPhiDist->Write();
   pionPhiDist->Write();
   kaonDCAprimary->Write();
   pionDCAprimary->Write();
     
   
   hadronPtDist->Write();
   hadronPhiDist->Write();
   hadronEtaDist->Write();
  
   D0PeakPlusBG->Write();
   D0LikeSignBG->Write();
   trackCounter->Write();
   dEdxVsPt->Write();
   invBetaVsPt->Write();
   
   d0CountPerEvent->Write();
   
   //for(int i = 0; i < nVzBins; i++){
   
        //eventCategoryCounter[i]->Write();
        
   //}     
   

   mOutputFile->Close();
  

   return kStOK;
}
//-----------------------------------------------------------------------------

//----------------------------------Make---------------------------------------
Int_t StPicoD0AnaMaker::Make(){ //begin Make member function
   readNextEvent();

   if (!mPicoDstMaker)
   {
      LOG_WARN << " StPicoD0AnaMaker - No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   StPicoDst const* picoDst = mPicoDstMaker->picoDst();


   if (!picoDst)
   {
      LOG_WARN << "StPicoD0AnaMaker - No PicoDst! Skip! " << endm;
      return kStWarn;
   }

   if(mPicoD0Event->runId() != picoDst->event()->runId() ||
       mPicoD0Event->eventId() != picoDst->event()->eventId())
   {
     LOG_ERROR <<" StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
     LOG_ERROR <<" StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
     exit(1);
   }

   //-------------------Begin User Analysis---------------------------
 
/****************************************************************************************************************/
/****************************GLOBAL AND LOCAL VARIABLES SET AND INITIALIZED BEGIN*********************************/  
/***************************************************************************************************************/ 
    
    
    //-------------- Various global cuts and pt ranges---------------- 
  
    trackCount = 0;
    centralityBin = 0;  //This number will be between 0 and 8 -- 9 bins total
    
    VzBin = 0;
	
    //--------------Local Variables used in the code-------------------
    double delPhi         = 0;
    double delEta         = 0;
    double pt             = 0;
    double phi            = 0;
    double eta            = 0;
    double dEdx           = 0;
    double beta           = 0;
    int    pairCharge     = 0;
    int    PIDflag        = 0;
    bool   minBiasFlag    = false;
    
    StPicoTrack* trk;
    //StPicoTrack* trk2;
    StThreeVectorF trackMom;
    StThreeVectorF trackMom2;
    double bField         = picoDst->event()->bField();
    StThreeVectorF pVtx   = picoDst->event()->primaryVertex();
    StThreeVectorF kaonPionMom;
    
    std::vector <StThreeVectorF> mPosiList;
    std::vector <StThreeVectorF> mNegaList;
    
    
    bool eventStoredInBuffer    = false;
    bool USEventStoredInBuffer  = false;
    bool LSEventStoredInBuffer  = false;
    int d0Counter = 0;
    
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

    if(mPicoD0Event->nKaons() > 0){ kaonDist->Fill(mPicoD0Event->nKaons());}
    if(mPicoD0Event->nPions() > 0){ pionDist->Fill(mPicoD0Event->nPions());}

    TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();    //The Kaon-Pion list is generated here

    trackCount = picoDst->numberOfTracks();            
   
    //cout << "picoDst->numTracks: " << trackCount << endl;
    //cout << "picoDst->numGlobalTracks: " << picoDst->event()->numberOfGlobalTracks() << endl;
   
    trackCounter->Fill(picoDst->numberOfTracks());
   
    centralityBin = getCentralityBin(trackCount);  //get centrality bin -- need to use the Run14 trackCount to get the correct bins
   
    if(USE_VZ_BINS){
        Vz = picoDst->event()->primaryVertex().z();
        VzBin         = getVzBin(Vz);                  //get Vz bin
    }
   
    else VzBin = 0;
   
    if(centralityBin == -1 || VzBin == -1) { return kStOk; }
   
      

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
        if(!mHFCuts->isGoodTrack(trk)) { continue; }
        
        pt   = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y()));
        if(pt<.15){continue;}  
        
        if(trackMom.mag() > .20 && trackMom.mag() < .45 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }
        if(trackMom.mag() > .70 && trackMom.mag() < .80 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }
        eta  = trackMom.pseudoRapidity();
		if(eta > 1 || eta < -1) { continue; }
		
        if(trk->charge() > 0) { mPosiList.push_back(trackMom); }
        else if(trk->charge() < 0) { mNegaList.push_back(trackMom); } 
        else { continue; }
        
        phi  = TMath::ATan2(trackMom.y(),trackMom.x());  
        dEdx = trk->dEdx();
            
        if(mHFCuts->hasTofPid(trk)){    
         
            beta = mHFCuts->getTofBeta(trk);                                                //Basic cut to ensure the tracks are in the TPC acceptance.
            invBetaVsPt->Fill(pt, (1/beta));        
        }
        
        dEdxVsPt->Fill(pt, dEdx);
        
        hadronPtDist->Fill(pt);                                                 //Fill pt dist. for all hadrons
        hadronPhiDist->Fill(phi);                                               //fill hists with phi and eta of hadrons for QA
        hadronEtaDist->Fill(eta);
        
        etaDistVz[VzBin]->Fill(eta);
        phiDistVz[VzBin]->Fill(phi);
        
   
        if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaKaon()) < 2.0)){       //ONLY check nSigma for TPC track -- used to be mHFCuts->isTPCKaon(trk), but this include pt cut
                                                                                  //need to fix this. Need to figure out how to access the nSigma from the cuts
              kaonPtDist->Fill(pt);
              kaonEtaDist->Fill(eta);
              kaonPhiDist->Fill(phi);
              continue;
        }

        if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaPion()) < 3.0)){     //ONLY check nSigma for TPC track
               
              pionPtDist->Fill(pt);
              pionEtaDist->Fill(eta);
              pionPhiDist->Fill(phi);
              continue;
        }
        
    } //End loop to fill basic track information and make nega/posi list   
       
/****************************************************************************************************************/  
/********************QA Loop for storing all of the basic information about the events and tracks END**********/
/****************************************************************************************************************/         
 
if(!NON_IDENTIFIED_CORR && !D0_HADRON_CORR) { return kStOk; }
/****************************************************************************************************************/  
/********************BEGIN NON-IDENTIFIED HADRON CODE BLOCK ---- MAKE SURE THE CORRECT FLAG IS SELECTED**********/
/****************************************************************************************************************/          
       
    //*********************************************************************************************************************** 
    //START SIBLING NONIDENTIFIED PARTICLE CORR HERE*************************************************************************
    //***********************************************************************************************************************
        if(NON_IDENTIFIED_CORR){//begin NON_IDENTIFIED_CORR conditional
        
            if(!eventStoredInBuffer){
                
                eventBuffer[VzBin][centralityBin]->addEvent(picoDst);
                eventStoredInBuffer = true;
            }
            
                
            
            for(unsigned int i = 0; i < picoDst->numberOfTracks(); ++i){//store tracks in buffer
                
                trk = picoDst->track(i);
                trackMom = trk->gMom(pVtx, bField);
                
                if(mHFCuts->isGoodTrack(trk)){
                
                    pt   = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y()));
                    if(pt<.15){continue;}  
        
                    if(trackMom.mag() > .20 && trackMom.mag() < .45 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }
                    if(trackMom.mag() > .70 && trackMom.mag() < .80 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }
                    eta  = trackMom.pseudoRapidity();
					if(eta > 1 || eta < -1) { continue; }                    
					
                    eventBuffer[VzBin][centralityBin]->addTrackToEvent(eventBuffer[VzBin][centralityBin]->getBufferIndex()-1, trackMom, trk->charge(), 0);
                
                }
                
            }
                 
               
            for(unsigned int i = 0; i < mPosiList.size(); i++){ // plus plus calc
               
                
                trackMom = mPosiList[i];
                phi  = trackMom.phi();  
                eta  = trackMom.pseudoRapidity();
                
                for(unsigned int j = 0; j < mPosiList.size(); j++){  //associated loop              
                    if(j == i) { continue; }
                   
                    trackMom2 = mPosiList[j];           
                    delEta = eta - trackMom2.pseudoRapidity();
            
                    delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                    if(delPhi < -TMath::PiOver2()){ delPhi = delPhi + 2*TMath::Pi(); }
                        else if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }
                
                    sibCorrPlusPlus[VzBin][centralityBin]->Fill(delEta, delPhi);
                }    
            }            
            
            for(unsigned int i = 0; i < mPosiList.size(); i++){ // plus minus calc
               
              
                trackMom = mPosiList[i];
                phi  = trackMom.phi();  
                eta  = trackMom.pseudoRapidity();
                
                for(unsigned int j = 0; j < mNegaList.size(); j++){  //associated loop              
                    //if(j == i) { continue; }
                   
                    trackMom2 = mNegaList[j];           
                    delEta = eta - trackMom2.pseudoRapidity();
            
                    delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                    if(delPhi < -TMath::PiOver2()){ delPhi = delPhi + 2*TMath::Pi(); }
                        else if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }
                
                    sibCorrPlusMinus[VzBin][centralityBin]->Fill(delEta, delPhi);
                }    
            }         

             for(unsigned int i = 0; i < mNegaList.size(); i++){ // minus plus calc
               
                
                trackMom =  mNegaList[i];
                phi  = trackMom.phi();  
                eta  = trackMom.pseudoRapidity();
                
                for(unsigned int j = 0; j < mPosiList.size(); j++){  //associated loop              
                    //if(j == i) { continue; }
                   
                    trackMom2 = mPosiList[j];           
                    delEta = eta - trackMom2.pseudoRapidity();
            
                    delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                    if(delPhi < -TMath::PiOver2()){ delPhi = delPhi + 2*TMath::Pi(); }
                        else if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }
                
                    sibCorrMinusPlus[VzBin][centralityBin]->Fill(delEta, delPhi);
                }    
            }                        
            
            for(unsigned int i = 0; i < mNegaList.size(); i++){ // minus minus calc
               
                
                trackMom = mNegaList[i];
                phi  = trackMom.phi();  
                eta  = trackMom.pseudoRapidity();
                
                for(unsigned int j = 0; j < mNegaList.size(); j++){  //associated loop              
                    if(j == i) { continue; }
                    
                    trackMom2 = mNegaList[j];           
                    delEta = eta - trackMom2.pseudoRapidity();
            
                    delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                    if(delPhi < -TMath::PiOver2()){ delPhi = delPhi + 2*TMath::Pi(); }
                        else if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }
                
                    sibCorrMinusMinus[VzBin][centralityBin]->Fill(delEta, delPhi);
                }    
            }                        
        
        //************************************************************************************ 
        //END SIBLING NONIDENTIFIED PARTICLE CORR HERE****************************************
        //************************************************************************************         

          
            
        //************************************************************************************       
        //BEGIN NON-IDENTIFIED EVENT MIXING HERE**********************************************
        //************************************************************************************
        
        if(eventBuffer[VzBin][centralityBin]->getBufferMaxSize() == eventBuffer[VzBin][centralityBin]->getBufferSize()){ // buffer-full mixing conditional for non-identified
    
            if(DEBUG){cout << "Mixing events in Vz/Centrality Bin " << VzBin<< "/" << centralityBin << endl << endl;}
    
            StMixerEvent* event1;
            StMixerEvent* event2;
        
            StMixerTrack mixedTrack1;
            StMixerTrack mixedTrack2;
        
            //int nTracks1 = 0;
            //int nTracks2 = 0;
        
            event1 = eventBuffer[VzBin][centralityBin]->getEvent(0); //The first event in the buffer will be the trigger D0 event
            //nTracks1 = event1->getNoTracks();
        
            for(int i = 1; i < eventBuffer[VzBin][centralityBin]->getBufferSize(); i++){ //begin associated event buffer loop
            
                event2 = eventBuffer[VzBin][centralityBin]->getEvent(i); //This is the event for the associated hadrons
            
                if(DEBUG) { 
        
                    //cout << endl << endl;
                    cout << "Mixing Event 0 with Event " << i << "." << endl;
                    cout << endl;
                }     
            
                for(int i = 0; i < event1->getNoPosiTracks(); i++){ // plus plus calc
               
                    mixedTrack1 = event1->getPosiTrack(i);
                    trackMom = mixedTrack1.gMom();
                    phi  = trackMom.phi();  
                    eta  = trackMom.pseudoRapidity();
                
                    for(int j = 0; j < event2->getNoPosiTracks(); j++){  //associated loop              
                        
                   
                        mixedTrack2 = event2->getPosiTrack(j);
                        trackMom2 = mixedTrack2.gMom();          
                        delEta = eta - trackMom2.pseudoRapidity();
            
                        delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                        if(delPhi < -TMath::PiOver2()){ delPhi = delPhi + 2*TMath::Pi(); }
                            else if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }
                
                        mixCorrPlusPlus[VzBin][centralityBin]->Fill(delEta, delPhi);
                    }    
                }            
            
                for(int i = 0; i < event1->getNoPosiTracks(); i++){ // plus minus calc
               
                    mixedTrack1 = event1->getPosiTrack(i);
                    trackMom = mixedTrack1.gMom();
                    phi  = trackMom.phi();  
                    eta  = trackMom.pseudoRapidity();
                
                    for(int j = 0; j < event2->getNoNegaTracks(); j++){  //associated loop              
                        
                   
                        mixedTrack2 = event2->getNegaTrack(j);
                        trackMom2 = mixedTrack2.gMom();               
                        delEta = eta - trackMom2.pseudoRapidity();
            
                        delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                        if(delPhi < -TMath::PiOver2()){ delPhi = delPhi + 2*TMath::Pi(); }
                            else if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }
                
                        mixCorrPlusMinus[VzBin][centralityBin]->Fill(delEta, delPhi);
                    }    
                }         

                for(int i = 0; i < event1->getNoNegaTracks(); i++){ // minus plus calc
               
                    mixedTrack1 = event1->getNegaTrack(i);
                    trackMom = mixedTrack1.gMom();
                    phi  = trackMom.phi();  
                    eta  = trackMom.pseudoRapidity();
                
                    for(int j = 0; j < event2->getNoPosiTracks(); j++){  //associated loop              
                        
                   
                        mixedTrack2 = event2->getPosiTrack(j);
                        trackMom2 = mixedTrack2.gMom();               
                        delEta = eta - trackMom2.pseudoRapidity();
            
                        delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                        if(delPhi < -TMath::PiOver2()){ delPhi = delPhi + 2*TMath::Pi(); }
                            else if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }
                
                        mixCorrMinusPlus[VzBin][centralityBin]->Fill(delEta, delPhi);
                    }    
                }                        
            
                for(int i = 0; i < event1->getNoNegaTracks(); i++){ // plus plus calc
               
                    mixedTrack1 = event1->getNegaTrack(i);
                    trackMom = mixedTrack1.gMom();
                    phi  = trackMom.phi();  
                    eta  = trackMom.pseudoRapidity();
                
                    for(int j = 0; j < event2->getNoNegaTracks(); j++){  //associated loop              
                        
                   
                        mixedTrack2 = event2->getNegaTrack(j);
                        trackMom2 = mixedTrack2.gMom();              
                        delEta = eta - trackMom2.pseudoRapidity();
            
                        delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                        if(delPhi < -TMath::PiOver2()){ delPhi = delPhi + 2*TMath::Pi(); }
                            else if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }
                
                        mixCorrMinusMinus[VzBin][centralityBin]->Fill(delEta, delPhi);
                    }    
                }                            
            
            }// end associated event buffer loop  
        
            eventBuffer[VzBin][centralityBin]->removeFirstEvent();
      
            if(DEBUG) { 
        
                cout << endl << endl;
                cout << "First Event in non-identified buffer for Vz/centrality bin " << VzBin << "/"<< centralityBin<< " cleared" << endl;
                cout << endl;
           }     
        
        }// end buffer-full mixing conditional for non-identified         
    
    }// end NON_IDENTIFIED_CORR conditional
       
    
     if(NON_IDENTIFIED_CORR){
        mPosiList.clear();
        mNegaList.clear();
    }
    
    if(NON_IDENTIFIED_CORR){

        if(DEBUG) {cout << "event complete" << endl;}
        return kStOk;
    } //Stops the analysis code before doing any D0 stuff if I am doing the non-identified correlations
    
    //mPosiList.clear();
    //mNegaList.clear(); 
    
    //************************************************************************************       
    //END NON-IDENTIFIED EVENT MIXING HERE************************************************
    //************************************************************************************

/****************************************************************************************************************/  
/********************BEGIN NON-IDENTIFIED HADRON CODE BLOCK ---- MAKE SURE THE CORRECT FLAG IS SELECTED**********/
/****************************************************************************************************************/     
    
//***************************************************************************************
//BEGIN D0 BLOCK OF CODE HERE -- MAKE SURE THE NON_IDENTIFIED FLAG IS SET TO FALSE!!!!
//***************************************************************************************    

//Still need to add checks for TOF information-----   

/****************************************************************************************************************/  
/***************************************LOOP TO STORE EVENTS IN MIXER BEGIN***************************************/
/****************************************************************************************************************/   
   
    for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx){//begin loop to check if event can be used (check if a D0 exists)
   
        StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
        StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
        StPicoTrack const* pion = picoDst->track(kp->pionIdx());
                    
                    //  ptmin ptmax   decayLenMin&Max   daughterDCA kaon/pion pt kaon/pion DCA  DCA to PV
        if(!isGoodPair(kp)) continue;             
        if(!cutCheck(kp, d0PtLow,  d0PtHigh,  d0DecayLengthMin,  d0DecayLengthMax,  daughterDCA,
                     d0DaughterPionPtMin, d0DaughterKaonPtMin, kaonDCA,  pionDCA,  d0DCAtoPV)) { continue; }
        
        if(kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh){ 
            
            if(pion->charge()*kaon->charge() < 0 && !USEventStoredInBuffer){ //If an US pair exists and the event was not previously stored, store it
            
                eventBufferUS[VzBin][centralityBin]->addEvent(picoDst);
                USEventStoredInBuffer = true;
            }
            
            else if(pion->charge()*kaon->charge() > 0 && !LSEventStoredInBuffer){ //If a LS pair exists and the event was not previously stored, store it
            
                eventBufferLS[VzBin][centralityBin]->addEvent(picoDst);
                LSEventStoredInBuffer = true;
            }
             
        }
   
    }//end loop to check if event can be used (check if a D0 exists)
   
    if(USEventStoredInBuffer || LSEventStoredInBuffer){ //begin conditional for the event being stored in the mixer
   
        if(DEBUG){ //This conditional would only trip if we made it past the previous conditional -- this means the event was stored
        
            cout << endl << endl;
            cout << "*********************EVENT START******************" << endl;
            cout << "We are on event # " << eventNumber << endl;
            cout << "This event has " << trackCount << " tracks." << endl;
            if(USE_VZ_BINS){ cout << "Vz: " << Vz << "   VzBin: " << VzBin << "    Centrality Bin: " << centralityBin << endl; }
            else cout << "Centrality Bin: " << centralityBin << endl;
            cout << endl;
            cout << "event " << eventNumber << " stored." << endl;
            cout << endl;
            //cout << "D0 candidate stats that tripped the event mixer: " << endl;
            //cout << "Mass: "<< kp->m() << "       US (-1) or LS (1): " << kaon->charge()*pion->charge() << endl;
            //cout << endl;
            eventNumber++;
        }
        
        eventCategoryCounter[VzBin]->Fill(centralityBin);
        eventCounter->Fill(0);
   
   //} //end conditional for the event being stored in the mixer
    
   
   //This block stores the tracks in the buffer and then stores the kaonPion pairs in the buffer in a separate list
    //if(USEventStoredInBuffer || LSEventStoredInBuffer){
    
        if(USEventStoredInBuffer){// begin conditional for US
        
            for(unsigned int i = 0; i < picoDst->numberOfTracks(); ++i){//begin loop to add picoDST tracks to US event
    
                trk = picoDst->track(i);                                                             
                trackMom = trk->gMom(pVtx, bField); 
                //add 
                if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaKaon()) < 2.0)){ PIDflag = 1; }
                else if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaPion()) < 3.0)) { PIDflag = 2; }
                else PIDflag = 0;
            
                eventBufferUS[VzBin][centralityBin]->addTrackToEvent(eventBufferUS[VzBin][centralityBin]->getBufferIndex()-1, trackMom, trk->charge(), PIDflag);   //0 -- any hadron, 1 -- kaon, 2 -- pion
                  
            }//end loop to add picoDST tracks to US event
        }// end conditional for US
        
        if(LSEventStoredInBuffer){//begin conditional for LS
        
            for(unsigned int i = 0; i < picoDst->numberOfTracks(); ++i){//begin loop to add picoDST tracks to LS event
    
                trk = picoDst->track(i);                                                             
                trackMom = trk->gMom(pVtx, bField); 
                //add 
                if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaKaon()) < 2.0)){ PIDflag = 1; }
                else if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaPion()) < 3.0)) { PIDflag = 2; }
                else PIDflag = 0;
            
                eventBufferLS[VzBin][centralityBin]->addTrackToEvent(eventBufferLS[VzBin][centralityBin]->getBufferIndex()-1, trackMom, trk->charge(), PIDflag);   //0 -- any hadron, 1 -- kaon, 2 -- pion
                  
            }//end loop to add picoDST tracks to LS event
        }//end conditional for LS
        
        int i = 0;
        //int buffIdx = 0;
        
        for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx){//begin loop to add kaonPions to event
   
            StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
            StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
            StPicoTrack const* pion = picoDst->track(kp->pionIdx());
   
            if(!isGoodPair(kp)) continue; 
            if(!cutCheck(kp, d0PtLow,  d0PtHigh,  d0DecayLengthMin,  d0DecayLengthMax,  daughterDCA,
                     d0DaughterPionPtMin, d0DaughterKaonPtMin, kaonDCA,  pionDCA,  d0DCAtoPV)) { continue; }
            
            if(kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh){ 
                
                kaonPionMom.set(kp->lorentzVector().px(),kp->lorentzVector().py(),kp->lorentzVector().pz());
                pairCharge = kaon->charge()*pion->charge(); //THIS IS NOT THE ACTUAL CHARGE -- SIMPLY THE US/LS VALUE
                if(pairCharge < 0 ){
                    
                    eventBufferUS[VzBin][centralityBin]->addKaonPionToEvent(eventBufferUS[VzBin][centralityBin]->getBufferIndex()-1, kaonPionMom, kp->m(), kp->kaonIdx(), kp->pionIdx(), pairCharge);
                }
                
                else if(pairCharge > 0){
                
                    eventBufferLS[VzBin][centralityBin]->addKaonPionToEvent(eventBufferLS[VzBin][centralityBin]->getBufferIndex()-1, kaonPionMom, kp->m(), kp->kaonIdx(), kp->pionIdx(), pairCharge);
                
                }
                
                if(DEBUG_MIX_BUFFER){ 
                    
                    cout << "Actual eta from picoDst: " << kp->eta() << endl;
                    
                    cout << "KaonPion charge from event info(LS or US): " << kaon->charge()*pion->charge() << endl;
                    cout << "KaonPion mass: " << kp->m() << endl;
                    cout << "_______________________________________________" << endl;
                    cout << endl;
                    
                    
                }    
                
                i++;               
            }
   
        }//end loop to add kaonPions to event
        
    } // end conditional for the event being stored in the mixer
   
   
   
/****************************************************************************************************************/  
/***************************************LOOP TO STORE EVENTS IN MIXER END***************************************/
/****************************************************************************************************************/   
   
/****************************************************************************************************************/   
/*****************************************BEGIN MAIN SIBLING LOOP************************************************/
/****************************************************************************************************************/
   
    for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx){ // begin main sibling loop
     
        StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
    
        StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
        StPicoTrack const* pion = picoDst->track(kp->pionIdx());
      
        if(!isGoodPair(kp)) continue;   
        if(!cutCheck(kp, d0PtLow,  d0PtHigh,  d0DecayLengthMin,  d0DecayLengthMax,  daughterDCA,
                     d0DaughterPionPtMin, d0DaughterKaonPtMin, kaonDCA,  pionDCA,  d0DCAtoPV)) { continue; }    
      
      
        if(kaon->charge()*pion->charge() < 0){// begin Unlike-sign conditional 
	      
            invMass->Fill(kp->m());     
            if(kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh) { D0PeakPlusBG->Fill(kp->m()); } 
             
	    }//end Unlike-sign conditional 
      
	  
        if(kaon->charge()*pion->charge() > 0){//begin Like-sign conditional 
          
            likeSignBG->Fill(kp->m());   
            if(kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh) { D0LikeSignBG->Fill(kp->m()); }
          
        }//end Like-sign conditional 
 
            
        
     ////////////////Fill QA Histograms from pair trees //////////////////////

        kaonDCAfromD0->Fill(kp->kaonDca());
        pionDCAfromD0->Fill(kp->pionDca());
        decayLengthQA->Fill(kp->decayLength());
        pointingAngleQA->Fill(kp->pointingAngle());
        daughterDCAQA->Fill(kp->dcaDaughters());

     /////////////////////////////////////////////////////////////////////////
     
        
    /****************************************************************************************************************/
    /****************************SIBLING EVENT PAIRS FORMATIION BEGINS HERE*****************************************/  
    /***************************************************************************************************************/     
        
	    if(kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh){ // begin loop over both unlike-sign AND LIKE SIGN D0 candidates 

            if(kaon->charge()*pion->charge() < 0){
                d0Counter = d0Counter + 1;
                D0ptDist->Fill(kp->pt()); 
                D0EtaDist->Fill(kp->eta());
                D0PhiDist->Fill(kp->phi());
            }
           
            for(unsigned int i = 0; i < picoDst->numberOfTracks(); ++i){ // begin picoDST loop for d0-hadron correlations
           
         	    if(i == kp->kaonIdx() || i == kp->pionIdx()) { continue; }    // Need to check this -- should avoid doing correlations with a D0 candidate daughter
             
                    ////////NEED TO ADD SOMETHING TO REJECT ELECTRONS AND MUONS///////////
                trk = picoDst->track(i);                                   //extract track from picoDst and store as StPicoTrack
           	    trackMom = trk->gMom(pVtx, bField);
               
                pt = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y())); //calculate tranverse momentum of track
                    
                    //if(pt < 1.0 ) { continue; }                           //Ensure track is within TPC acceptance and choose a specific sample of hadrons.
                                                                                
                phi = TMath::ATan2(trackMom.y(),trackMom.x());  
                eta = trackMom.pseudoRapidity();
                 
                delPhi = kp->phi()-phi;
                if(delPhi < -TMath::PiOver2()){ delPhi = delPhi + 2*TMath::Pi(); }
             	    else if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }
                
                 delEta =kp->eta()-eta;
                   
                 if(kaon->charge()*pion->charge() < 0){//Unlike Sign histograms
                    
                        sibCorrBinUS[VzBin][centralityBin]->Fill(delEta, delPhi);
                 }     
                 else if(kaon->charge()*pion->charge() > 0){// like sign

                        sibCorrBinLS[VzBin][centralityBin]->Fill(delEta, delPhi);
                 } 
            
            } // end picoDST loop for d0-hadron correlations
        }// End d0 sibling conditional
    }// End main sibling loop

    /****************************************************************************************************************/
    /****************************SIBLING EVENT PAIRS FORMATIION ENDS HERE********************************************/  
    /****************************************************************************************************************/
    
/****************************************************************************************************************/   
/*****************************************END MAIN SIBLING LOOP**************************************************/
/****************************************************************************************************************/    
    
    
//------------------------------------------------------------------------------------------------------------------------//   
    
/***************************************************************************************************************/
/*************************************EVENT MIXING BEINGS HERE**************************************************/
/***************************************************************************************************************/
    
    /************************************************BEGIN UNLIKE SIGN EVENT MIXER BLOCK*****************************************/

    if(eventBufferUS[VzBin][centralityBin]->getBufferMaxSize() == eventBufferUS[VzBin][centralityBin]->getBufferSize()){ //begin buffer-full US mixing conditional
    
        if(DEBUG){cout << "Mixing US events in Vz/Centrality Bin " << VzBin<< "/" << centralityBin << endl << endl;}
    
        StMixerEvent* event1;
        StMixerEvent* event2;
        
        StMixerTrack kaonPionTrack1;
        StMixerTrack mixedTrack2;
        
        //int nTracks1 = 0;
        int nTracks2 = 0;
        
        event1 = eventBufferUS[VzBin][centralityBin]->getEvent(0); //The first event in the buffer will be the trigger D0 event
        
        for(int i = 1; i < eventBufferUS[VzBin][centralityBin]->getBufferSize(); i++){ //begin associated event buffer loop
            
            if(DEBUG) { 
        
                //cout << endl << endl;
                cout << "Mixing US Event 0 with US Event " << i << "." << endl;
                cout << endl;
            
            }     
            
            event2 = eventBufferUS[VzBin][centralityBin]->getEvent(i); //This is the event for the associated hadrons
            nTracks2 = event2->getNoTracks();
            
            for (int idx = 0; idx < event1->getKaonPionListSize(); idx++){//begin loop over event 1 kaon-pion list
                
                kaonPionTrack1 = event1->getKaonPionAt(idx);
                StMixerTrack mixerKaon = event1->getTrack(kaonPionTrack1.kaonIdx());
                StMixerTrack mixerPion = event1->getTrack(kaonPionTrack1.pionIdx());
                
                if(mixerKaon.charge()*mixerPion.charge() > 0){ continue; }
                
                if(kaonPionTrack1.mass() > D0InvMassLow && kaonPionTrack1.mass() < D0InvMassHigh){//begin D0 trigger conditional 
                 
                    for(int j = 0; j < nTracks2; j++){//begin loop over event 2 tracks
                        
                        //if(j == kaonPionTrack1.kaonIdx() || j == kaonPionTrack1.pionIdx()) { continue; }                                
                        
                        mixedTrack2 = event2->getTrack(j);                                                            
           	            trackMom = mixedTrack2.gMom();
               
                        pt = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y())); 
                        phi = TMath::ATan2(trackMom.y(),trackMom.x());  
                        eta = trackMom.pseudoRapidity();
                 
                        delPhi = kaonPionTrack1.gMom().phi()-phi;
                        if(delPhi < -TMath::PiOver2()){ delPhi = delPhi + 2*TMath::Pi(); }
                            else if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }
                
                        delEta =kaonPionTrack1.gMom().pseudoRapidity()-eta;
                
                        mixCorrBinUS[VzBin][centralityBin]->Fill(delEta, delPhi); //delEta & delPhi stored in US mixed correlation histogram
            
                    }// end loop over event 2 tracks
                }// end D0 trigger conditional
            }// end loop over event 1 kaon-pion list
        }// end associated event buffer loop
        
        eventBufferUS[VzBin][centralityBin]->removeFirstEvent();
      
        if(DEBUG) { 
        
            cout << endl << endl;
            cout << "First Event in US buffer for Vz/centrality bin "<< VzBin << "/"<< centralityBin<< " cleared" << endl;
            cout << endl;
       }     
        
    }// end buffer-full US mixing conditional        

    /************************************************END UNLIKE SIGN EVENT MIXER BLOCK*****************************************/    
      
   
   /***********************************************BEGIN LIKE SIGN EVENT MIXER BLOCK*******************************************/  
      
    if(eventBufferLS[VzBin][centralityBin]->getBufferMaxSize() == eventBufferLS[VzBin][centralityBin]->getBufferSize()){ //begin buffer-full US mixing conditional
    
        if(DEBUG){cout << "Mixing LS events in Vz/Centrality Bin " << VzBin<< "/" << centralityBin << endl << endl;}
    
        StMixerEvent* event1;
        StMixerEvent* event2;
        
        StMixerTrack kaonPionTrack1;
        StMixerTrack mixedTrack2;
        
        //int nTracks1 = 0;
        int nTracks2 = 0;
        
        event1 = eventBufferLS[VzBin][centralityBin]->getEvent(0); //The first event in the buffer will be the trigger D0 event
        
        for(int i = 1; i < eventBufferLS[VzBin][centralityBin]->getBufferSize(); i++){ //begin associated event buffer loop
            
            if(DEBUG) { 
        
                //cout << endl << endl;
                cout << "Mixing LS Event 0 with LS Event " << i << "." << endl;
                cout << endl;
            
            }     
            
            event2 = eventBufferLS[VzBin][centralityBin]->getEvent(i); //This is the event for the associated hadrons
            nTracks2 = event2->getNoTracks();
            
            for (int idx = 0; idx < event1->getKaonPionListSize(); idx++){//begin loop over event 1 kaon-pion list
                
                kaonPionTrack1 = event1->getKaonPionAt(idx);
                StMixerTrack mixerKaon = event1->getTrack(kaonPionTrack1.kaonIdx());
                StMixerTrack mixerPion = event1->getTrack(kaonPionTrack1.pionIdx());
                
                if(mixerKaon.charge()*mixerPion.charge() < 0){ continue; }
                
                if(kaonPionTrack1.mass() > D0InvMassLow && kaonPionTrack1.mass() < D0InvMassHigh){//begin D0 trigger conditional 
                 
                    for(int j = 0; j < nTracks2; j++){//begin loop over event 2 tracks
                        
                        //if(j == kaonPionTrack1.kaonIdx() || j == kaonPionTrack1.pionIdx()) { continue; }                                
                        
                        mixedTrack2 = event2->getTrack(j);                                                            
           	            trackMom = mixedTrack2.gMom();
               
                        pt = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y())); 
                        phi = TMath::ATan2(trackMom.y(),trackMom.x());  
                        eta = trackMom.pseudoRapidity();
                 
                        delPhi = kaonPionTrack1.gMom().phi()-phi;
                        if(delPhi < -TMath::PiOver2()){ delPhi = delPhi + 2*TMath::Pi(); }
                            else if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }
                
                        delEta = kaonPionTrack1.gMom().pseudoRapidity()-eta;
                
                        mixCorrBinLS[VzBin][centralityBin]->Fill(delEta, delPhi); //delEta & delPhi stored in US mixed correlation histogram
            
                    }// end loop over event 2 tracks
                }// end D0 trigger conditional
            }// end loop over event 1 kaon-pion list
        }// end associated event buffer loop
        
        eventBufferLS[VzBin][centralityBin]->removeFirstEvent();
      
        if(DEBUG) { 
        
            cout << endl << endl;
            cout << "First Event in LS buffer for Vz/centrality bin "<< VzBin << "/"<< centralityBin<< " cleared" << endl;
            cout << endl;
            cout << "buffer size after removal of first event: " << eventBufferLS[VzBin][centralityBin]->getBufferSize() << endl;
       }     
        
    }// end buffer-full LS mixing conditional    
    
    /***********************************************END LIKE SIGN EVENT MIXER BLOCK*******************************************/
    
    
/***************************************************************************************************************/
/*************************************EVENT MIXING ENDS HERE****************************************************/
/***************************************************************************************************************/   
       
       
       
       d0CountPerEvent->Fill(d0Counter);
       invMassMinusBG->Add(invMass, likeSignBG, 1, -1);
       D0PeakMinusBG->Add(D0PeakPlusBG, D0LikeSignBG, 1, -1);
    
   //if(DEBUG) { cout << endl << "***************EVENT END****************************" << endl; }
    
   //-------------------End Current Event Analysis--------------------------------
   
   return kStOK;

}//end Make member function







//---------------------User Functions-------------------------------------

bool StPicoD0AnaMaker::isGoodPair(StKaonPion const* const kp) const
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


bool StPicoD0AnaMaker::cutCheck(StKaonPion const* const kp, double ptMin, double ptMax, double decayLengthMin, double decayLengthMax, 
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


int StPicoD0AnaMaker::getCentralityBin(int nTracks){

    if(nTracks >= 2   && nTracks < 14)  { return 0;  }
    if(nTracks >= 14  && nTracks < 32)  { return 1;  }
    if(nTracks >= 34  && nTracks < 67)  { return 2;  }
    if(nTracks >= 67  && nTracks < 115) { return 3;  }
    if(nTracks >= 115 && nTracks < 183) { return 4;  }
    if(nTracks >= 183 && nTracks < 275) { return 5;  }
    if(nTracks >= 275 && nTracks < 392) { return 6;  }
    if(nTracks >= 392 && nTracks < 537) { return 7;  }
    if(nTracks >= 537 && nTracks < 720) { return 8;  }
    if(nTracks >= 720 && nTracks < 829) { return 9;  }
    if(nTracks >= 829)                  { return 10; }

    else return -1;

}    

int StPicoD0AnaMaker::getVzBin(double Vz){

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






