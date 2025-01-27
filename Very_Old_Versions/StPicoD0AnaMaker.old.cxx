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

****/

ClassImp(StPicoD0AnaMaker)

StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name,char const * inputFilesList, 
    char const * outName,StPicoDstMaker* picoDstMaker): 
  StMaker(name),mPicoDstMaker(picoDstMaker),mPicoD0Event(NULL), mOutFileName(outName), mInputFileList(inputFilesList),
  mOutputFile(NULL), mChain(NULL), mEventCounter(0), mHFCuts(NULL)
{}

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

   // --------------------Begin User Variables-----------------------------------
   
   ptDist          = new TH1D("Pt Distribution", "Pt Distribution", 1000, 0, 10);              
   invMass         = new TH1D("unlikeSign", "unlikeSign", 50, 1.6, 2.1);
   angCorrPhi      = new TH1D("#Delta Phi", "#Delta#phi", 500, -TMath::PiOver2(), 3*TMath::PiOver2());
   angCorrEta      = new TH1D("#Delta#Eta", "#Delta#eta", 500, -2 , 2);  
   kaonDist        = new TH1D("Kaon Distribution", "Kaon Distribution", 500, 0 , 500);
   pionDist        = new TH1D("Pion Distribution", "Pion Distribution", 2000, 0 , 2000);
   likeSignBG      = new TH1D("Like Sign BG", "Like Sign BG", 50, 1.6, 2.1);
   invMassMinusBG  = new TH1D("D0 minus LS BG", "D0 minus LS BG", 50, 1.6, 2.1);
   angCorr2DAll    = new TH2D("2D Ang Corr", "2D Ang Corr", 100, -2, 2, 100, -TMath::PiOver2(), 3*TMath::PiOver2());
   
   invMassBin1      = new TH1D("unlikeSign_bin1", "unlikeSign_bin1", 50, 1.6, 2.1);
   invMassBin2      = new TH1D("unlikeSign_bin2", "unlikeSign_bin2", 50, 1.6, 2.1);
   invMassBin3      = new TH1D("unlikeSign_bin3", "unlikeSign_bin3", 50, 1.6, 2.1);
   invMassBin4      = new TH1D("unlikeSign_bin4", "unlikeSign_bin4", 50, 1.6, 2.1);
   invMassBin5      = new TH1D("unlikeSign_bin5", "unlikeSign_bin5", 50, 1.6, 2.1);

   likeSignBin1     = new TH1D("likeSign_bin1", "likeSign_bin1", 50, 1.6, 2.1);
   likeSignBin2     = new TH1D("likeSign_bin2", "likeSign_bin2", 50, 1.6, 2.1);
   likeSignBin3     = new TH1D("likeSign_bin3", "likeSign_bin3", 50, 1.6, 2.1);
   likeSignBin4     = new TH1D("likeSign_bin4", "likeSign_bin4", 50, 1.6, 2.1);
   likeSignBin5     = new TH1D("likeSign_bin5", "likeSign_bin5", 50, 1.6, 2.1);

   USminusLSBin1    = new TH1D("US-LS_bin1", "US-LS_bin1", 50, 1.6, 2.1);
   USminusLSBin2    = new TH1D("US-LS_bin2", "US-LS_bin2", 50, 1.6, 2.1);
   USminusLSBin3    = new TH1D("US-LS_bin3", "US-LS_bin3", 50, 1.6, 2.1);
   USminusLSBin4    = new TH1D("US-LS_bin4", "US-LS_bin4", 50, 1.6, 2.1);
   USminusLSBin5    = new TH1D("US-LS_bin5", "US-LS_bin5", 50, 1.6, 2.1);


   angCorrPhiBin1   = new TH1D("#Delta Phi bin 1", "#Delta#phi bin 1", 100, -TMath::PiOver2(), 3*TMath::PiOver2());
   angCorrEtaBin1   = new TH1D("#Delta#Eta bin 1", "#Delta#eta bin 1", 100, -2 , 2);
   angCorr2DBin1    = new TH2D("2D Ang Corr bin 1", "2D Ang Corr bin 1", 100, -2, 2, 100, -TMath::PiOver2(), 3*TMath::PiOver2());
   angCorrPhiBin2   = new TH1D("#Delta Phi bin 2", "#Delta#phi bin 2", 100, -TMath::PiOver2(), 3*TMath::PiOver2());
   angCorrEtaBin2   = new TH1D("#Delta#Eta bin 2", "#Delta#eta bin 2", 100, -2 , 2);
   angCorr2DBin2    = new TH2D("2D Ang Corr bin 2", "2D Ang Corr bin 2", 100, -2, 2, 100, -TMath::PiOver2(), 3*TMath::PiOver2());
   angCorrPhiBin3   = new TH1D("#Delta Phi bin 3", "#Delta#phi bin 3", 100, -TMath::PiOver2(), 3*TMath::PiOver2());
   angCorrEtaBin3   = new TH1D("#Delta#Eta bin 3", "#Delta#eta bin 3", 100, -2 , 2);
   angCorr2DBin3    = new TH2D("2D Ang Corr bin 3", "2D Ang Corr bin 3", 100, -2, 2, 100, -TMath::PiOver2(), 3*TMath::PiOver2());
   angCorrPhiBin4   = new TH1D("#Delta Phi bin 4", "#Delta#phi bin 4", 100, -TMath::PiOver2(), 3*TMath::PiOver2());
   angCorrEtaBin4   = new TH1D("#Delta#Eta bin 4", "#Delta#eta bin 4", 100, -2 , 2);
   angCorr2DBin4    = new TH2D("2D Ang Corr bin 4", "2D Ang Corr bin 4", 100, -2, 2, 100, -TMath::PiOver2(), 3*TMath::PiOver2());
   angCorrPhiBin5   = new TH1D("#Delta Phi bin 5", "#Delta#phi bin 5", 100, -TMath::PiOver2(), 3*TMath::PiOver2());
   angCorrEtaBin5   = new TH1D("#Delta#Eta bin 5", "#Delta#eta bin 5", 100, -2 , 2);
   angCorr2DBin5    = new TH2D("2D Ang Corr bin 5", "2D Ang Corr bin 5", 100, -2, 2, 100, -TMath::PiOver2(), 3*TMath::PiOver2());
   
   //QA Histograms
   
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
   hadronEtaDist   = new TH1D("Inclusvie Hadron Eta", "Inclusive Hadron Eta", 1000, -1, 1);

   kaonDCAfromD0   = new TH1D("DCA for kaons from D0", "DCA for kaons from D0", 500, 0.0, 0.5);
   pionDCAfromD0   = new TH1D("DCA for pions from D0", "DCA for pions from D0", 500, 0.0, 0.65);
   decayLengthQA   = new TH1D("D0 Candidate Decay Length (no mass cut)", "D0 Candidate Decay Length (no mass cut)", 500, 0.0, 1.5);
   pointingAngleQA = new TH1D("D0 Candidate Pointing Angle(no mass cut)", "D0 Candidate Pointing Angle (no mass cut)", 500, 0.0, 1.7);
   daughterDCAQA   = new TH1D("D0 Daughter DCA", "D0 Daughter DCA (no mass cut)", 500, 0.0, .01);

   //QA for mass-cut D0  

   D0ptDist        = new TH1D("D0 Candidate pt Dist (mass cut)", "D0 Candidate pt Dist (mass cut)", 1000, 0, 10);
   D0EtaDist       = new TH1D("D0 Eta Dist", "D0 #eta Dist. (mass cut)", 500, -1, 1);
   D0PhiDist       = new TH1D("D0 Phi Dist", "D0 #phi Dist. (mass cut)", 500, -2*TMath::Pi(), 2*TMath::Pi());
   D0PeakPlusBG    = new TH1D("D0 Peak + BG", "D0 Peak + BG", 50, 1.6, 2.1);
   D0LikeSignBG    = new TH1D("LikeSign peak range", "LikeSign peak range", 50, 1.6, 2.1);
   D0PeakMinusBG   = new TH1D("D0 Peak", "D0 Peak", 50, 1.6, 2.1);

   eventCounter    = new TH1D("number of events used", "number of events used", 3, 0, 3);
   trackCounter    = new TH1D("number of tracks per event", "number of tracks per event", 2000, 0, 1999);
   
   //Histogram formatting

   eventCounter->GetXaxis()->SetBinLabel(1,"minBias");
   eventCounter->GetXaxis()->SetBinLabel(2,"total");
   eventCounter->GetXaxis()->SetBinLabel(3,"events from bad runs");
   
   //----------------------End User Variables------------------------------------
   
   return kStOK;
}
//-----------------------------------------------------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker()
{
   /*  */
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Finish()
{
   LOG_INFO << " StPicoD0AnaMaker - writing data and closing output file " <<endm;
   mOutputFile->cd();
   // save user variables here
   
   ptDist->Write();
   invMass->Write();
   kaonDist->Write();
   pionDist->Write();
   likeSignBG->Write();
   angCorrEta->Write();
   angCorrPhi->Write();
   angCorr2DAll->Write();
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
   invMassBin1->Write();
   invMassBin2->Write();
   invMassBin3->Write();
   invMassBin4->Write();
   invMassBin5->Write();   
   angCorrPhiBin1->Write();
   angCorrEtaBin1->Write();
   angCorr2DBin1->Write();
   angCorrPhiBin2->Write();
   angCorrEtaBin2->Write();
   angCorr2DBin2->Write();
   angCorrPhiBin3->Write();
   angCorrEtaBin3->Write();
   angCorr2DBin3->Write();
   angCorrPhiBin4->Write();
   angCorrEtaBin4->Write();
   angCorr2DBin4->Write();
   angCorrPhiBin5->Write();
   angCorrEtaBin5->Write();
   angCorr2DBin5->Write();
   hadronPtDist->Write();
   hadronPhiDist->Write();
   hadronEtaDist->Write(); 
   likeSignBin1->Write();
   likeSignBin2->Write();
   likeSignBin3->Write();
   likeSignBin4->Write();
   likeSignBin5->Write();
   USminusLSBin1->Write();
   USminusLSBin2->Write();
   USminusLSBin3->Write();
   USminusLSBin4->Write();
   USminusLSBin5->Write();
   D0PeakPlusBG->Write();
   D0LikeSignBG->Write();
   trackCounter->Write();

   mOutputFile->Close();
  

   return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Make()
{
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
    
    //-------------- Various cuts and pt ranges----------------
    double pt1            = 0.0;
    double pt2            = 1.0;
    double pt3            = 2.0;
    double pt4            = 3.0;
    double pt5            = 5.0;
    double pt6            = 10.0;
    double decayLength1   = .0145;
    double decayLength2   = .0181;
    double decayLength3   = .0212;
    double decayLength4   = .0247;
    double decayLength5   = .0259;
    double dcaDaughters1  = .0084;
    double dcaDaughters2  = .0066;
    double dcaDaughters3  = .0057;
    double dcaDaughters4  = .0050;
    double dcaDaughters5  = .0060;
    double dcaKaonPV1     = .0103;
    double dcaKaonPV2     = .0091;
    double dcaKaonPV3     = .0095;
    double dcaKaonPV4     = .0079;
    double dcaKaonPV5     = .0058;
    double dcaPionPV1     = .0110;
    double dcaPionPV2     = .0111;
    double dcaPionPV3     = .0086;
    double dcaPionPV4     = .0081;
    double dcaPionPV5     = .0062;
    double dcaV0toPV1     = .0061;
    double dcaV0toPV2     = .0049;
    double dcaV0toPV3     = .0038;
    double dcaV0toPV4     = .0038;
    double dcaV0toPV5     = .0040;

    double kaonPtCut      = .15;
    double pionPtCut      = .15;
    
    double D0InvMassLow   = 1.82;
    double D0InvMassHigh  = 1.90;
	//--------------------------------------------
    
    //Variables used in the code
    double delPhi         = 0;
    double delEta         = 0;
    double pt             = 0;
    double phi            = 0;
    double eta            = 0;
    StPicoTrack* trk;
    /////////////////////////////
    
 // if(!mHFCuts->isGoodEvent(const_cast<const StPicoDst*>(picoDst), NULL)) return kStOk; //makes sure the event comes from good run
  if(!mHFCuts->isGoodRun(picoDst->event())){
      
      eventCounter->Fill(1);
      eventCounter->Fill(2);
      return kStOk; //makes sure the event comes from good run
  }

  eventCounter->Fill(1);
  if(!picoDst->event()->isMinBias()) {return kStOK;}          //min bias event flag
  //if(!picoDst->event()->isCentral()) {return kStOK;}          //central events flag ---- What are the criteria for central events
  eventCounter->Fill(0);

  if(mPicoD0Event->nKaons() > 0){ kaonDist->Fill(mPicoD0Event->nKaons());}
  if(mPicoD0Event->nPions() > 0){ pionDist->Fill(mPicoD0Event->nPions());}

  //Fill kaon and pion pt dist
  

   TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();

   //------------------------BEGIN LOOP THROUGH THE CURRENT EVENT'S TRACKS ------------------------------------------------------//////
   
   /********************ALL GENERAL EVENT STATS SHOULD BE COLLECTED IN THIS BLOCK*************************************************************/
   
   //Still need to add checks for TOF information-----
   
   trackCounter->Fill(picoDst->numberOfTracks());
   
   for(unsigned int i = 0; i < picoDst->numberOfTracks(); ++i){

          trk = picoDst->track(i);
          pt = TMath::Sqrt((trk->pMom().x()*trk->pMom().x())+(trk->pMom().y()*trk->pMom().y()));
          phi = trk->pMom().phi();  
          eta = trk->pMom().pseudoRapidity();
          
          if(pt<.15){continue;}                                                   //Basic cut to ensure the tracks are in the TPC acceptance.
          

          hadronPtDist->Fill(pt);                                                 //Fill pt dist. for all hadrons
          hadronPhiDist->Fill(phi);                                               //fill hists with phi and eta of hadrons for QA
          hadronEtaDist->Fill(trk->pMom().pseudoRapidity());
   


          if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaKaon()) < 2.0)){       //ONLY check nSigma for TPC track -- used to be mHFCuts->isTPCKaon(trk), but this include pt cut
                                                                                  //need to fix this. Need to figure out how to access the nSigma from the cuts
                kaonPtDist->Fill(pt);
                kaonEtaDist->Fill(eta);
                kaonPhiDist->Fill(phi);
               // kaonDCAprimary->Fill(      
                continue;
          }

          if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaPion()) < 3.0)){     //ONLY check nSigma for TPC track
                
                pionPtDist->Fill(pt);
                pionEtaDist->Fill(eta);
                pionPhiDist->Fill(phi);
               // pionDCAprimary->Fill(
                continue;
          }
       }
   /******************************************************************************************************************/
   
   /****************************BEGIN BLOCK USING TREE INFORMATION****************************************************/
   for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx){
     
	 // this is an example of how to get the kaonPion pairs and their corresponding tracks
      StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
    
      StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
      StPicoTrack const* pion = picoDst->track(kp->pionIdx());
      
      if(!isGoodPair(kp)) continue;                                 //THIS STOPS THE ANALYSIS. MAKE SURE ANY GENERAL EVENT STATS ARE COLLECTED BEFORE THIS!!!!
      
      ////////////////Fill QA Histograms from pair trees //////////////////////

     
     kaonDCAfromD0->Fill(kp->kaonDca());
     pionDCAfromD0->Fill(kp->pionDca());
     decayLengthQA->Fill(kp->decayLength());
     pointingAngleQA->Fill(kp->pointingAngle());
     daughterDCAQA->Fill(kp->dcaDaughters());

     /////////////////////////////////////////////////////////////////////////
     
      //cutCheck(StKaonPion const* const kp, double ptMin, double ptMax, double decayLengthMin, double decayLengthMax, double dcaDaughters, 
      //         double kaonPtCut, double pionPtCut, double dcaKaontoPV, double dcaPiontoPV, double dcaV0toPV)

      //ptDist->Fill(kp->pt());
      if(kaon->charge()*pion->charge() < 0){
	      if     (cutCheck(kp, pt1, pt2, decayLength1, 999999.0, dcaDaughters1, kaonPtCut, pionPtCut, dcaKaonPV1, dcaPionPV1, dcaV0toPV1)) { invMassBin1->Fill(kp->m()); } 
	      else if(cutCheck(kp, pt2, pt3, decayLength2, 999999.0, dcaDaughters2, kaonPtCut, pionPtCut, dcaKaonPV2, dcaPionPV2, dcaV0toPV2)) { invMassBin2->Fill(kp->m()); }
              else if(cutCheck(kp, pt3, pt4, decayLength3, 999999.0, dcaDaughters3, kaonPtCut, pionPtCut, dcaKaonPV3, dcaPionPV3, dcaV0toPV3)) { invMassBin3->Fill(kp->m()); }     //Fill invmass based on pt of D0
	      else if(cutCheck(kp, pt4, pt5, decayLength4, 999999.0, dcaDaughters4, kaonPtCut, pionPtCut, dcaKaonPV4, dcaPionPV4, dcaV0toPV4)) { invMassBin4->Fill(kp->m()); }
              else if(cutCheck(kp, pt5, pt6, decayLength5, 999999.0, dcaDaughters5, kaonPtCut, pionPtCut, dcaKaonPV5, dcaPionPV5, dcaV0toPV5)) { invMassBin5->Fill(kp->m()); }
	  
	  
	      if(kaon->gPt() > 1.2 && pion->gPt() > 1.2 && kp->decayLength() > .0200 && 
                 kp->dcaDaughters() < .0055 && kp->kaonDca() > .008 && 
                 kp->pionDca() > .008 && kp->perpDcaToVtx() < .0065)
                { 
                   invMass->Fill(kp->m());     
                   if(kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh) { D0PeakPlusBG->Fill(kp->m()); } 
                }     //Pt integrated invmass
	  }
      
	  
      if(kaon->charge()*pion->charge() > 0){
          if     (cutCheck(kp, pt1, pt2, decayLength1, 999999.0, dcaDaughters1, kaonPtCut, pionPtCut, dcaKaonPV1, dcaPionPV1, dcaV0toPV1)) { likeSignBin1->Fill(kp->m()); }
          else if(cutCheck(kp, pt2, pt3, decayLength2, 999999.0, dcaDaughters2, kaonPtCut, pionPtCut, dcaKaonPV2, dcaPionPV2, dcaV0toPV2)) { likeSignBin2->Fill(kp->m()); }
          else if(cutCheck(kp, pt3, pt4, decayLength3, 999999.0, dcaDaughters3, kaonPtCut, pionPtCut, dcaKaonPV3, dcaPionPV3, dcaV0toPV3)) { likeSignBin3->Fill(kp->m()); }     //Fill likesign based on pt of D0
          else if(cutCheck(kp, pt4, pt5, decayLength4, 999999.0, dcaDaughters4, kaonPtCut, pionPtCut, dcaKaonPV4, dcaPionPV4, dcaV0toPV4)) { likeSignBin4->Fill(kp->m()); }
          else if(cutCheck(kp, pt5, pt6, decayLength5, 999999.0, dcaDaughters5, kaonPtCut, pionPtCut, dcaKaonPV5, dcaPionPV5, dcaV0toPV5)) { likeSignBin5->Fill(kp->m()); }

          if(kaon->gPt() > 1.2 && pion->gPt() > 1.2 && kp->decayLength() > .0200 && 
             kp->dcaDaughters() < .0055 &&  kp->kaonDca() > .008 && 
             kp->pionDca() > .008 && kp->perpDcaToVtx() < .0065)
                { 
                   likeSignBG->Fill(kp->m());   
                   if(kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh) { D0LikeSignBG->Fill(kp->m()); }
                }

          }
 
    
       
        //This block is for analyzing anything involving just the D0 peak.***************

	if(kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh && kaon->charge()*pion->charge() < 0){

           //D0Count->Fill(0);
           D0ptDist->Fill(kp->pt()); 
           D0EtaDist->Fill(kp->eta());
           D0PhiDist->Fill(kp->phi());

     	  for(unsigned int i = 0; i < picoDst->numberOfTracks(); ++i){
           
         	 if(i == kp->kaonIdx() || i == kp->pionIdx()) { continue; }                                // Need to check this -- should avoid doing correlations with a D0 candidate daughter
                    
                    ////////NEED TO ADD SOMETHING TO REJECT ELECTRONS AND MUONS///////////
                
          	        trk = picoDst->track(i);                                                             //extract track from picoDst and store as StPicoTrack
           	   
                    //Calculate kinematic variables
               
                    pt = TMath::Sqrt((trk->pMom().x()*trk->pMom().x())+(trk->pMom().y()*trk->pMom().y())); //calculate tranverse momentum of track
                    if(pt < 1.0 ) { continue; }                                                 //Ensure track is within TPC acceptance and choose a specific sample of hadrons.
                                                                                
                    phi = TMath::ATan2(trk->pMom().y(),trk->pMom().x());  
                
                   // hadronPhiDist->Fill(phi);                                                            //fill hists with phi and eta of hadrons for QA
                   // hadronEtaDist->Fill(trk->pMom().pseudoRapidity());
                
                
                    //if(pt < 1.0) { continue; }                                                          //Use a pt cut on associated hadrons to help reduce the background in the delEta/delPhi distributions
                    
                    delPhi = kp->phi()-phi;
                    if(delPhi < -TMath::PiOver2()){ delPhi = delPhi + 2*TMath::Pi(); }
             	       else if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }
                
                   // angCorrPhi->Fill(delPhi);

                    eta = trk->pMom().pseudoRapidity();
                    delEta =kp->eta()-eta;
                   
                 if(kp->decayLength() > .0200 && kp->dcaDaughters() < .0055 && kp->kaonDca() > .008 && kp->pionDca() > .008 && 
                    kp->kaonDca() > .008 && kp->pionDca() > .008 && kp->perpDcaToVtx() < .0065){ //D0 track cuts for pt-inclusive sibling dist.
                    angCorrPhi->Fill(delPhi);
                    angCorrEta->Fill(delEta);
                    angCorr2DAll->Fill(delEta, delPhi);
                 }
                
                 //pt bins for D0 pt -- NEED TO WRITE A FUNCTION TO DO THESE CUT CHECKS
                if(cutCheck(kp, pt1, pt2, decayLength1, 999999, dcaDaughters1, 0.0, 0.0, dcaKaonPV1, dcaPionPV1, dcaV0toPV1)){
                    angCorrPhiBin1->Fill(delPhi);
                    angCorrEtaBin1->Fill(delEta);
                    angCorr2DBin1->Fill(delEta, delPhi); 
                }
                    else if(cutCheck(kp, pt2, pt3, decayLength2, 999999, dcaDaughters2, 0.0, 0.0, dcaKaonPV2, dcaPionPV2, dcaV0toPV2)){
                        angCorrPhiBin2->Fill(delPhi);
                        angCorrEtaBin2->Fill(delEta);
                        angCorr2DBin2->Fill(delEta, delPhi);
                    }
                    else if(cutCheck(kp, pt3, pt4, decayLength3, 999999, dcaDaughters3, 0.0, 0.0, dcaKaonPV3, dcaPionPV3, dcaV0toPV3)){ 
                        angCorrPhiBin3->Fill(delPhi);
                        angCorrEtaBin3->Fill(delEta);
                        angCorr2DBin3->Fill(delEta, delPhi);
                    }
                    else if(cutCheck(kp, pt4, pt5, decayLength4, 999999, dcaDaughters4, 0.0, 0.0, dcaKaonPV4, dcaPionPV4, dcaV0toPV4)){
                        angCorrPhiBin4->Fill(delPhi);
                        angCorrEtaBin4->Fill(delEta);
                        angCorr2DBin4->Fill(delEta, delPhi);
                    }
                    else if(cutCheck(kp, pt5, pt6, decayLength5, 999999, dcaDaughters5, 0.0, 0.0, dcaKaonPV5, dcaPionPV5, dcaV0toPV5)){ 
                        angCorrPhiBin5->Fill(delPhi);
                        angCorrEtaBin5->Fill(delEta);
                        angCorr2DBin5->Fill(delEta, delPhi);
                    }  
         } 
    }

       invMassMinusBG->Add(invMass, likeSignBG, 1, -1);
       USminusLSBin1->Add(invMassBin1, likeSignBin1, 1, -1);
       USminusLSBin2->Add(invMassBin2, likeSignBin2, 1, -1);
       USminusLSBin3->Add(invMassBin3, likeSignBin3, 1, -1);
       USminusLSBin4->Add(invMassBin4, likeSignBin4, 1, -1);
       USminusLSBin5->Add(invMassBin5, likeSignBin5, 1, -1);
       D0PeakMinusBG->Add(D0PeakPlusBG, D0LikeSignBG, 1, -1);
    }


    //START BUILDING SIDE BAND HERE.
    
    
    
   //-------------------End User Analysis--------------------------------
   
   return kStOK;
}

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
  
                                                
 
