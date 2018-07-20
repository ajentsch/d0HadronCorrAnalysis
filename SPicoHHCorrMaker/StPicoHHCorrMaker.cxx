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

#include "StPicoCharmContainers/StPicoD0Event.h"
#include "StPicoCharmContainers/StKaonPion.h"

#include "StMixedEventBuffer/StMixedEventBuffer.h"
#include "StMixedEventBuffer/StMixerEvent.h"
#include "StMixedEventBuffer/StMixerTrack.h"

#include "StPicoHHCorrMaker.h"
#include "StPicoHFMaker/StHFCuts.h"

/***********************************

    This version of the code is used to just calculate h-h correlations, with charge separation. (update: 10/6/2016)

***********************************/


ClassImp(StPicoHHCorrMaker)

StPicoHHCorrMaker::StPicoHHCorrMaker(char const * name, char const * inputFilesList, 
                                   char const * outName, StPicoDstMaker* picoDstMaker): 
                            StMaker(name),mPicoDstMaker(picoDstMaker),mPicoD0Event(NULL),
                            mOutFileName(outName), mInputFileList(inputFilesList),
                            mOutputFile(NULL), mChain(NULL), mEventCounter(0), mHFCuts(NULL)
{}

//---------------------------------------Initialization-----------------------------------------
Int_t StPicoHHCorrMaker::Init()
{
   mPicoD0Event = new StPicoD0Event();

   mChain = new TChain("T");
   std::ifstream listOfFiles(mInputFileList.Data());
   if (listOfFiles.is_open())
   {
      std::string file;
      while (getline(listOfFiles, file))
      {
         LOG_INFO << "StPicoHHCorrMaker - Adding :" << file << endm;
         mChain->Add(file.c_str());
      }
   }
   else
   {
      LOG_ERROR << "StPicoHHCorrMaker - Could not open list of files. ABORT!" << endm;
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
   NON_IDENTIFIED_CORR = true;
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
   
    kaonPtCut         = .15;    // in GeV
    pionPtCut         = .15;
    
    hadronPtCutLow    = 0.15;
    hadronPtCutHigh   = 15.45;
    
    trackDCAmax       = 3.0;    //in centimeters
    trackChi2max      = 3.0;
   
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
   
   for(int i = 0; i < nVzBins; i++){
    
            if(USE_VZ_BINS) { str1 = eventCounterLabel + binLabelVz[i]; }
        
            else str1 = "event count per centrality bin (for events containing D0s)";
        
            eventCategoryCounter[i] = new TH1D(str1, str1, 11, 0, 11);
            eventCategoryCounter[i]->GetXaxis()->SetTitle("CentralityBin");
        }        
      
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
   eventCounter        = new TH1D("number of events used", "number of events used", 4, 0, 4);
   trackCounter        = new TH1D("number of tracks per event", "number of tracks per event", 2000, 0, 1999);
   kaonPtDist          = new TH1D("Kaon Pt Distribution", "Kaon Pt Distribution", 1000, 0, 5);
   pionPtDist          = new TH1D("Pion Pt Distribution", "Pion Pt Distribution", 1000, 0, 5);
   kaonEtaDist         = new TH1D("Kaon Eta Distribution", "Kaon Eta Distribution", 1000, -1, 1);
   pionEtaDist         = new TH1D("Pion Eta Distribution", "Pion Eta Distribution", 1000, -1, 1);
   kaonPhiDist         = new TH1D("Kaon Phi Distribution", "Kaon Phi Distribution", 1000, -2*TMath::Pi(), 2*TMath::Pi());
   pionPhiDist         = new TH1D("Pion Phi Distribution", "Pion Phi Distribution", 1000, -2*TMath::Pi(), 2*TMath::Pi());
   DCAtoPrimaryVertex  = new TH1D("track DCA to PV", "track DCA to PV", 500, 0.0, 10.0);
   DCAtoPrimaryVertexCut = new TH1D("track DCA to PV cut check", "track DCA to PV cut check", 500, 0.0, 10.0);
   //pionDCAprimary  = new TH1D("DCA pions from primary", "DCA pions from primary", 500, 0.0, 0.5);
   
   hadronPtDist    = new TH1D("Inclusive Hadron pt", "Inclusive Hadron pt", 1000, 0, 10);
   hadronPhiDist   = new TH1D("Inclusive Hadron Phi", "Inclusive Hadron Phi", 1000, -2*TMath::Pi(), 2*TMath::Pi());
   hadronEtaDist   = new TH1D("Inclusvie Hadron Eta", "Inclusive Hadron Eta", 1000, -2, 2);
   hadronChi2      = new TH1D("Chi2 for hadron tracks", "Chi2 for hadron tracks", 500, 0, 10);
   
   pVtxX           = new TH1D("X position of pVtx", "X position of pVtx", 500, -10, 10);
   pVtxY           = new TH1D("Y position of pVtx", "Y position of pVtx", 500, -10, 10);
   pVtxZ           = new TH1D("Z position of pVtx", "Z position of pVtx", 500, -7, 7);
   
   
   dEdxVsPt        = new TH2D("dEdx_vs_P", "dEdx_vs_P", 500, 0, 10, 500, 0, 15);
   invBetaVsPt     = new TH2D("#Beta^{-1} Vs. P", "#Beta^{-1} Vs. P", 500, 0, 10, 500, 0, 4);
   
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
StPicoHHCorrMaker::~StPicoHHCorrMaker()
{
   /*  */
}
//-----------------------------------------------------------------------------

//------------------------------------Finish-----------------------------------
Int_t StPicoHHCorrMaker::Finish()
{
   LOG_INFO << " StPicoHHCorrMaker - writing data and closing output file " <<endm;
   mOutputFile->cd();
   // save user variables here
   
   
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
   eventCounter->Write();
   kaonPtDist->Write();
   pionPtDist->Write();
   kaonEtaDist->Write();
   pionEtaDist->Write();
   kaonPhiDist->Write();
   pionPhiDist->Write();
   DCAtoPrimaryVertex->Write();
   DCAtoPrimaryVertexCut->Write();
   //pionDCAprimary->Write();
   hadronPtDist->Write();
   hadronPhiDist->Write();
   hadronEtaDist->Write();
   trackCounter->Write();
   dEdxVsPt->Write();
   invBetaVsPt->Write();
   pVtxX->Write();
   pVtxY->Write();
   pVtxZ->Write();
   hadronChi2->Write();
   
   mOutputFile->Close();
  

   return kStOK;
}
//-----------------------------------------------------------------------------

//----------------------------------Make---------------------------------------
Int_t StPicoHHCorrMaker::Make(){ //begin Make member function
   readNextEvent();

   if (!mPicoDstMaker)
   {
      LOG_WARN << " StPicoHHCorrMaker - No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   StPicoDst const* picoDst = mPicoDstMaker->picoDst();


   if (!picoDst)
   {
      LOG_WARN << "StPicoHHCorrMaker - No PicoDst! Skip! " << endm;
      return kStWarn;
   }

   if(mPicoD0Event->runId() != picoDst->event()->runId() ||
       mPicoD0Event->eventId() != picoDst->event()->eventId())
   {
     LOG_ERROR <<" StPicoHHCorrMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
     LOG_ERROR <<" StPicoHHCorrMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
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
    double delPhiCp       = 0;
    double delEta         = 0;
    double pt             = 0;
    double phi            = 0;
    double eta            = 0;
    double dEdx           = 0;
    double beta           = 0;
    float  trackDCA       = 0;
    //int    pairCharge     = 0;
    //int    PIDflag        = 0;
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

    
    if(USE_VZ_BINS){
        Vz = picoDst->event()->primaryVertex().z();
        VzBin         = getVzBin(Vz);                  //get Vz bin
    }
   
    else VzBin = 0;
    
    if(VzBin == -1) { return kStOk; }    
      
    pVtxX->Fill(picoDst->event()->primaryVertex().x());
    pVtxY->Fill(picoDst->event()->primaryVertex().y());
    pVtxZ->Fill(picoDst->event()->primaryVertex().z()); 
      
/****************************************************************************************************************/
/****************************PRELIMINARY EVENT CUTS END********************************************************/  
/***************************************************************************************************************/ 
  
  
/****************************************************************************************************************/  
/********************QA Loop for storing all of the basic information about the events and tracks BEGIN**********/
/****************************************************************************************************************/   
   
    for(unsigned int i = 0; i < picoDst->numberOfTracks(); ++i){ // Begin loop to fill basic track information and make nega/posi list
                                                                 //gets pt, eta, phi, dEdx and beta from TOF
        trk = picoDst->track(i);
        
        trackDCA = ((trk->helix().origin())-pVtx).mag();
        
        DCAtoPrimaryVertex->Fill(trackDCA);
        
        if(trk->chi2() > trackChi2max) { continue; }
        if(!checkDCAtoPV(trackDCA))  { continue; }   // track quality cut for DCA to PV
        //if(!(trk->pMom().mag() > 0)) { continue; }   //only use primary tracks in the calculations
        
        trackMom = trk->gMom(pVtx, bField);
        
        if(!mHFCuts->isGoodTrack(trk)) { continue; }
        
        
        pt   = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y()));
        if(pt<hadronPtCutLow || pt>hadronPtCutHigh){continue;}  
        
        if(trackMom.mag() > .20 && trackMom.mag() < .45 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }
        if(trackMom.mag() > .70 && trackMom.mag() < .80 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }
        eta  = trackMom.pseudoRapidity();
		if(eta > 1 || eta < -1) { continue; }
        
        DCAtoPrimaryVertexCut->Fill(trackDCA);
		
        if(trk->charge() > 0) { mPosiList.push_back(trackMom); }
        else if(trk->charge() < 0) { mNegaList.push_back(trackMom); } 
        else { continue; }
        
        trackCount++;
        
        phi  = TMath::ATan2(trackMom.y(),trackMom.x());  
        dEdx = trk->dEdx();
        hadronChi2->Fill(trk->chi2());
            
        if(mHFCuts->hasTofPid(trk)){    
         
            beta = mHFCuts->getTofBeta(trk);                                                //Basic cut to ensure the tracks are in the TPC acceptance.
            invBetaVsPt->Fill(trackMom.mag(), (1/beta));        
        }
        
        dEdxVsPt->Fill(trackMom.mag(), dEdx);
        
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
       
       
    centralityBin = getCentralityBin(trackCount);  //get centrality bin -- need to use the Run14 trackCount to get the correct bins
   
   
    if(centralityBin == -1) { return kStOk; }    
    
    
    trackCounter->Fill(trackCount);    
       
       
/****************************************************************************************************************/  
/********************QA Loop for storing all of the basic information about the events and tracks END**********/
/****************************************************************************************************************/         
 
if(!NON_IDENTIFIED_CORR) { return kStOk; }
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
                trackDCA = ((trk->helix().origin())-pVtx).mag();
                
                if(mHFCuts->isGoodTrack(trk)){
                
                    //if(!(trk->pMom().mag() > 0)) { continue; }
                 
        
                    if(trk->chi2() > trackChi2max) { continue; }
                    if(!checkDCAtoPV(trackDCA))    { continue; }  
                    pt   = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y()));
                    if(pt<hadronPtCutLow || pt>hadronPtCutHigh){continue;}  
        
                    if(trackMom.mag() > .20 && trackMom.mag() < .45 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }
                    if(trackMom.mag() > .70 && trackMom.mag() < .80 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }
                    eta  = trackMom.pseudoRapidity();
					if(eta > 1 || eta < -1) { continue; }                    
					if(trk->charge() == 0)  { continue; }
       
                    
                    eventBuffer[VzBin][centralityBin]->addTrackToEvent(eventBuffer[VzBin][centralityBin]->getBufferIndex()-1, trackMom, trk->charge(), 0);
                
                }
                
            }
                 
               
            for(unsigned int i = 0; i < mPosiList.size(); i++){ // plus plus calc
               
                
                trackMom = mPosiList[i];
                phi  = trackMom.phi();  
                eta  = trackMom.pseudoRapidity();
                
                for(unsigned int j = i+1; j < mPosiList.size(); j++){  //associated loop              
                    if(j == i) { continue; }
                   
                    trackMom2 = mPosiList[j];           
                    delEta = TMath::Abs(eta - trackMom2.pseudoRapidity());
            
                    delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                    
                    if(delPhi < -TMath::Pi()) { delPhi = delPhi + 2*TMath::Pi(); }     //shift [-2Pi, 2Pi] -> [-Pi,Pi]
                    else if(delPhi > TMath::Pi()) { delPhi = delPhi - 2*TMath::Pi(); }
                    
                    delPhi = TMath::Abs(delPhi);//Gets absolute value of angle 
                    if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }//....and shifts it
                    
                    delPhiCp = -delPhi;//Gets negative copy of absolute value of angle.....
                    if(delPhiCp < -TMath::PiOver2()){ delPhiCp = delPhiCp + 2*TMath::Pi(); }//....and shifts it
                
                    sibCorrPlusPlus[VzBin][centralityBin]->Fill(delEta, delPhi);
                    sibCorrPlusPlus[VzBin][centralityBin]->Fill(-delEta, delPhiCp);
                    sibCorrPlusPlus[VzBin][centralityBin]->Fill(-delEta, delPhi);
                    sibCorrPlusPlus[VzBin][centralityBin]->Fill(delEta, delPhiCp);
                }    
            }            
            
            for(unsigned int i = 0; i < mPosiList.size(); i++){ // plus minus calc
               
              
                trackMom = mPosiList[i];
                phi  = trackMom.phi();  
                eta  = trackMom.pseudoRapidity();
                
                for(unsigned int j = 0; j < mNegaList.size(); j++){  //associated loop              
                    //if(j == i) { continue; }
                   
                    trackMom2 = mNegaList[j];           
                    delEta = TMath::Abs(eta - trackMom2.pseudoRapidity());
            
                    delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                    
                    if(delPhi < -TMath::Pi()) { delPhi = delPhi + 2*TMath::Pi(); }     //shift [-2Pi, 2Pi] -> [-Pi,Pi]
                    else if(delPhi > TMath::Pi()) { delPhi = delPhi - 2*TMath::Pi(); }
                    
                    delPhi = TMath::Abs(delPhi);//Gets absolute value of angle 
                    if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }//....and shifts it
                    
                    delPhiCp = -delPhi;//Gets negative copy of absolute value of angle.....
                    if(delPhiCp < -TMath::PiOver2()){ delPhiCp = delPhiCp + 2*TMath::Pi(); }//....and shifts it
                
                    sibCorrPlusMinus[VzBin][centralityBin]->Fill(delEta, delPhi);
                    sibCorrPlusMinus[VzBin][centralityBin]->Fill(-delEta, delPhiCp);
                    sibCorrPlusMinus[VzBin][centralityBin]->Fill(-delEta, delPhi);
                    sibCorrPlusMinus[VzBin][centralityBin]->Fill(delEta, delPhiCp);
                }    
            }         

             for(unsigned int i = 0; i < mNegaList.size(); i++){ // minus plus calc
               
                
                trackMom =  mNegaList[i];
                phi  = trackMom.phi();  
                eta  = trackMom.pseudoRapidity();
                
                for(unsigned int j = 0; j < mPosiList.size(); j++){  //associated loop              
                    //if(j == i) { continue; }
                   
                    trackMom2 = mPosiList[j];           
                    delEta = TMath::Abs(eta - trackMom2.pseudoRapidity());
            
                    delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                    
                    if(delPhi < -TMath::Pi()) { delPhi = delPhi + 2*TMath::Pi(); }     //shift [-2Pi, 2Pi] -> [-Pi,Pi]
                    else if(delPhi > TMath::Pi()) { delPhi = delPhi - 2*TMath::Pi(); }
                    
                    delPhi = TMath::Abs(delPhi);//Gets absolute value of angle 
                    if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }//....and shifts it
                    
                    delPhiCp = -delPhi;//Gets negative copy of absolute value of angle.....
                    if(delPhiCp < -TMath::PiOver2()){ delPhiCp = delPhiCp + 2*TMath::Pi(); }//....and shifts it
                
                    sibCorrMinusPlus[VzBin][centralityBin]->Fill(delEta, delPhi);
                    sibCorrMinusPlus[VzBin][centralityBin]->Fill(-delEta, delPhiCp);
                    sibCorrMinusPlus[VzBin][centralityBin]->Fill(-delEta, delPhi);
                    sibCorrMinusPlus[VzBin][centralityBin]->Fill(delEta, delPhiCp);
                }    
            }                        
            
            for(unsigned int i = 0; i < mNegaList.size(); i++){ // minus minus calc
               
                
                trackMom = mNegaList[i];
                phi  = trackMom.phi();  
                eta  = trackMom.pseudoRapidity();
                
                for(unsigned int j = i+1; j < mNegaList.size(); j++){  //associated loop              
                    if(j == i) { continue; }
                    
                    trackMom2 = mNegaList[j];           
                    delEta = TMath::Abs(eta - trackMom2.pseudoRapidity());
            
                    delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                    
                    if(delPhi < -TMath::Pi()) { delPhi = delPhi + 2*TMath::Pi(); }     //shift [-2Pi, 2Pi] -> [-Pi,Pi]
                    else if(delPhi > TMath::Pi()) { delPhi = delPhi - 2*TMath::Pi(); }
                    
                    delPhi = TMath::Abs(delPhi);//Gets absolute value of angle 
                    if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }//....and shifts it
                    
                    delPhiCp = -delPhi;//Gets negative copy of absolute value of angle.....
                    if(delPhiCp < -TMath::PiOver2()){ delPhiCp = delPhiCp + 2*TMath::Pi(); }//....and shifts it
                
                    sibCorrMinusMinus[VzBin][centralityBin]->Fill(delEta, delPhi);
                    sibCorrMinusMinus[VzBin][centralityBin]->Fill(-delEta, delPhiCp);
                    sibCorrMinusMinus[VzBin][centralityBin]->Fill(-delEta, delPhi);
                    sibCorrMinusMinus[VzBin][centralityBin]->Fill(delEta, delPhiCp);
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
                        delEta = TMath::Abs(eta - trackMom2.pseudoRapidity());
            
                        delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                    
                        if(delPhi < -TMath::Pi()) { delPhi = delPhi + 2*TMath::Pi(); }     //shift [-2Pi, 2Pi] -> [-Pi,Pi]
                        else if(delPhi > TMath::Pi()) { delPhi = delPhi - 2*TMath::Pi(); }
                    
                        delPhi = TMath::Abs(delPhi);//Gets absolute value of angle 
                        if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }//....and shifts it
                    
                        delPhiCp = -delPhi;//Gets negative copy of absolute value of angle.....
                        if(delPhiCp < -TMath::PiOver2()){ delPhiCp = delPhiCp + 2*TMath::Pi(); }//....and shifts it
                
                        mixCorrPlusPlus[VzBin][centralityBin]->Fill(delEta, delPhi);
                        mixCorrPlusPlus[VzBin][centralityBin]->Fill(-delEta, delPhiCp);
                        mixCorrPlusPlus[VzBin][centralityBin]->Fill(-delEta, delPhi);
                        mixCorrPlusPlus[VzBin][centralityBin]->Fill(delEta, delPhiCp);
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
                        delEta = TMath::Abs(eta - trackMom2.pseudoRapidity());
            
                        delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                    
                        if(delPhi < -TMath::Pi()) { delPhi = delPhi + 2*TMath::Pi(); }     //shift [-2Pi, 2Pi] -> [-Pi,Pi]
                        else if(delPhi > TMath::Pi()) { delPhi = delPhi - 2*TMath::Pi(); }
                    
                        delPhi = TMath::Abs(delPhi);//Gets absolute value of angle 
                        if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }//....and shifts it
                    
                        delPhiCp = -delPhi;//Gets negative copy of absolute value of angle.....
                        if(delPhiCp < -TMath::PiOver2()){ delPhiCp = delPhiCp + 2*TMath::Pi(); }//....and shifts it
                
                        mixCorrPlusMinus[VzBin][centralityBin]->Fill(delEta, delPhi);
                        mixCorrPlusMinus[VzBin][centralityBin]->Fill(-delEta, delPhiCp);
                        mixCorrPlusMinus[VzBin][centralityBin]->Fill(-delEta, delPhi);
                        mixCorrPlusMinus[VzBin][centralityBin]->Fill(delEta, delPhiCp);
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
                        delEta = TMath::Abs(eta - trackMom2.pseudoRapidity());
            
                        delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                    
                        if(delPhi < -TMath::Pi()) { delPhi = delPhi + 2*TMath::Pi(); }     //shift [-2Pi, 2Pi] -> [-Pi,Pi]
                        else if(delPhi > TMath::Pi()) { delPhi = delPhi - 2*TMath::Pi(); }
                    
                        delPhi = TMath::Abs(delPhi);//Gets absolute value of angle 
                        if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }//....and shifts it
                    
                        delPhiCp = -delPhi;//Gets negative copy of absolute value of angle.....
                        if(delPhiCp < -TMath::PiOver2()){ delPhiCp = delPhiCp + 2*TMath::Pi(); }//....and shifts it
                
                        mixCorrMinusPlus[VzBin][centralityBin]->Fill(delEta, delPhi);
                        mixCorrMinusPlus[VzBin][centralityBin]->Fill(-delEta, delPhiCp);
                        mixCorrMinusPlus[VzBin][centralityBin]->Fill(-delEta, delPhi);
                        mixCorrMinusPlus[VzBin][centralityBin]->Fill(delEta, delPhiCp);
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
                        delEta = TMath::Abs(eta - trackMom2.pseudoRapidity());
            
                        delPhi = phi - TMath::ATan2(trackMom2.y(), trackMom2.x());
                    
                        if(delPhi < -TMath::Pi()) { delPhi = delPhi + 2*TMath::Pi(); }     //shift [-2Pi, 2Pi] -> [-Pi,Pi]
                        else if(delPhi > TMath::Pi()) { delPhi = delPhi - 2*TMath::Pi(); }
                    
                        delPhi = TMath::Abs(delPhi);//Gets absolute value of angle 
                        if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }//....and shifts it
                    
                        delPhiCp = -delPhi;//Gets negative copy of absolute value of angle.....
                        if(delPhiCp < -TMath::PiOver2()){ delPhiCp = delPhiCp + 2*TMath::Pi(); }//....and shifts it
                
                        mixCorrMinusMinus[VzBin][centralityBin]->Fill(delEta, delPhi);
                        mixCorrMinusMinus[VzBin][centralityBin]->Fill(-delEta, delPhiCp);
                        mixCorrMinusMinus[VzBin][centralityBin]->Fill(-delEta, delPhi);
                        mixCorrMinusMinus[VzBin][centralityBin]->Fill(delEta, delPhiCp);
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


    
   //-------------------End Current Event Analysis--------------------------------
   
   return kStOK;

}//end Make member function







//---------------------User Functions-------------------------------------




int StPicoHHCorrMaker::getCentralityBin(int nTracks){

   /* if(nTracks >= 2   && nTracks < 14)  { return 0;  }
    if(nTracks >= 14  && nTracks < 32)  { return 1;  }
    if(nTracks >= 34  && nTracks < 67)  { return 2;  }
    if(nTracks >= 67  && nTracks < 115) { return 3;  }
    if(nTracks >= 115 && nTracks < 183) { return 4;  }
    if(nTracks >= 183 && nTracks < 275) { return 5;  }
    if(nTracks >= 275 && nTracks < 392) { return 6;  }
    if(nTracks >= 392 && nTracks < 537) { return 7;  }
    if(nTracks >= 537 && nTracks < 720) { return 8;  }
    if(nTracks >= 720 && nTracks < 829) { return 9;  }
    if(nTracks >= 829)                  { return 10; }*/
    
    if(nTracks >= 2   && nTracks < 16)  { return 0;  }
    if(nTracks >= 16  && nTracks < 31)  { return 1;  }
    if(nTracks >= 31  && nTracks < 65)  { return 2;  }
    if(nTracks >= 65  && nTracks < 116) { return 3;  }
    if(nTracks >= 116 && nTracks < 187) { return 4;  }
    if(nTracks >= 187 && nTracks < 237) { return 5;  }
    if(nTracks >= 237 && nTracks < 300) { return 6;  }
    if(nTracks >= 300 && nTracks < 370) { return 7;  }
    if(nTracks >= 370 && nTracks < 440) { return 8;  }
    if(nTracks >= 440 && nTracks < 520) { return 9;  }
    if(nTracks >= 520)                  { return 10; }
    

    else return -1;

}    

int StPicoHHCorrMaker::getVzBin(double Vz){

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

bool StPicoHHCorrMaker::checkDCAtoPV(float trackDCA){

     return (trackDCA <= trackDCAmax);
     
}




