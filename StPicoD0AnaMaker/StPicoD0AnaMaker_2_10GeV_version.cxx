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

#include "StMixedEventBuffer/StMixedEventBuffer.h"
#include "StMixedEventBuffer/StMixerEvent.h"
#include "StMixedEventBuffer/StMixerTrack.h"

#include "StPicoD0AnaMaker.h"
#include "StPicoHFMaker/StHFCuts.h"


/***********************************************************************

Author: Alex Jentsch

Most recent update: 8/15/2017

************************************************************************/

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
    DEBUG               = false;         //sets the debug cout statements -- default is FALSE
    DEBUG_MIX_BUFFER    = false;         //sets the debug cout statements for event mixing -- default is FALSE
    USE_CENT_BINS       = true;          //sets the multiplicity binning for the correlations and mixing -- default is TRUE
    USE_VZ_BINS         = true;          //sets the vz-binning for the correlations and mixing -- default is TRUE   
    USE_PT_BINS         = true;          //sets the Pt-binning for the correlations and mixing -- default is TRUE   
    USE_FINE_PT_BINS    = false;         //sets the 10-pt-bin output for calculating efficiency (want narrow pt-bins) -- default is FALSE
    SINGLES_DISTS       = true;         //sets the output for "singles" distributions, which are the 2d-histograms on the eta and phi of the D0 and hadron, rather than delEta/delPhi -- default is FALSE
    D0_HADRON_CORR      = false;          //sets the output for the correlations (turn off for testing other things) -- default is TRUE
    EVENT_MIXING        = false;          //toggle events mixing on and off -- default is TRUE
   
    USE_TWO_TRACK_INEFF = true;          //sets the cut for two-track inefficiency filtration (merging and crossing) -- no default yet
    USE_TOF             = false;         //sets the cut for STRICT use of TOF (requires it for D0 reconstruction)    -- no default yet
    USE_TOF_HYBRID      = false;         //sets the cut for hybrid TOF (use it if you got it)                       -- no default yet
    USE_DOUBLE_MIS_PID_PROTECTION = false; //sets the cut to remove doubley-counted pairs when double mis PID happens with TPC -- no default yet
   
    USE_PAIR_WISE_PT_CUT = false;          //obsolete cut -- will remove later.
    USE_PAIR_WEIGHT      = true;
    //----------------------------------------------------------//
   
    //---------------------Important constants-----------------//
   
    BUFFER_SIZE        = 5;                 //how many events to use in event mixer (per vz/multiplicity bin)      
    NUM_PHI_BINS       = 12;                //number of delPhi bins to use for correlations
    NUM_ETA_BINS       = 9;                 //number of delEta bins to use for correlations
    NUM_VZ_BINS        = 10;                //number of vz bins to use for proper events mixing
    NUM_CENT_BINS      = 16;                //number of multiplicity bins -- I named these "cent" bins early on in the code writing and just never changed it
    NUM_PT_BINS        = 3;                 //number of pt bins used for analysis -- narrow bins down low (more stats), wider bins up high.
    NUM_D0_CUT_PT_BINS = 4;                 //same as above, just includes an "integrated" bin --  this was a screw-up that should have been addressed.
    phiBinShift        = (TMath::Pi()/12.0);     //This number shifts the phi bins to ensure that 0 and pi are at the center of a bin
    phiBinShiftPrimary = (TMath::Pi()/24.0);
    
    
    //------------------efficiency correction constants---DO NOT CHANGE THESE!-----------
    
    siblingTripleEff[0][0] = 0.393047;
    siblingTripleEff[0][1] = 0.394547;
    siblingTripleEff[0][2] = 0.396262;
    siblingTripleEff[1][0] = 0.420809;
    siblingTripleEff[1][1] = 0.422871;
    siblingTripleEff[1][2] = 0.424891;
    siblingTripleEff[2][0] = 0.4424;
    siblingTripleEff[2][1] = 0.439583;
    siblingTripleEff[2][2] = 0.439322;
    siblingTripleEff[3][0] = 0.496483;
    siblingTripleEff[3][1] = 0.498035;
    siblingTripleEff[3][2] = 0.499246;
    siblingTripleEff[4][0] = 0.588417;
    siblingTripleEff[4][1] = 0.605623;
    siblingTripleEff[4][2] = 0.609697;

    mixedTripleEff[0][0] = 0.392859;
    mixedTripleEff[0][1] = 0.394468;
    mixedTripleEff[0][2] = 0.396218;
    mixedTripleEff[1][0] = 0.420582;
    mixedTripleEff[1][1] = 0.422762;
    mixedTripleEff[1][2] = 0.424831;
    mixedTripleEff[2][0] = 0.442251;
    mixedTripleEff[2][1] = 0.439538;
    mixedTripleEff[2][2] = 0.439267;
    mixedTripleEff[3][0] = 0.496691;
    mixedTripleEff[3][1] = 0.497913;
    mixedTripleEff[3][2] = 0.499215;
    mixedTripleEff[4][0] = 0.587526;
    mixedTripleEff[4][1] = 0.605144;
    mixedTripleEff[4][2] = 0.609426;

    siblingDoubleEff[0][0] = 0.00166726;
    siblingDoubleEff[0][1] = 0.00141288;
    siblingDoubleEff[0][2] = 0.00115337;
    siblingDoubleEff[1][0] = 0.00223267;
    siblingDoubleEff[1][1] = 0.00110712;
    siblingDoubleEff[1][2] = 0.00185601;
    siblingDoubleEff[2][0] = 0.0038524;
    siblingDoubleEff[2][1] = 0.00403631;
    siblingDoubleEff[2][2] = 0.00361487;
    siblingDoubleEff[3][0] = 0.00744415;
    siblingDoubleEff[3][1] = 0.0110326;
    siblingDoubleEff[3][2] = 0.00703188;
    siblingDoubleEff[4][0] = 0.0161578;
    siblingDoubleEff[4][1] = 0.0263765;
    siblingDoubleEff[4][2] = 0.0156148;

    mixedDoubleEff[0][0] = 0.00166659;
    mixedDoubleEff[0][1] = 0.00141294;
    mixedDoubleEff[0][2] = 0.00115324;
    mixedDoubleEff[1][0] = 0.00223077;
    mixedDoubleEff[1][1] = 0.00110693;
    mixedDoubleEff[1][2] = 0.00185579;
    mixedDoubleEff[2][0] = 0.0038524;
    mixedDoubleEff[2][1] = 0.00403651;
    mixedDoubleEff[2][2] = 0.00361474;
    mixedDoubleEff[3][0] = 0.00745632;
    mixedDoubleEff[3][1] = 0.0110326;
    mixedDoubleEff[3][2] = 0.00703235;
    mixedDoubleEff[4][0] = 0.0161606;
    mixedDoubleEff[4][1] = 0.026375;
    mixedDoubleEff[4][2] = 0.0156062;
    
    BOverBPlusS[0][0] = .883161;
    BOverBPlusS[0][1] = .961933;
    BOverBPlusS[0][2] = .983925;
    BOverBPlusS[1][0] = .523717;
    BOverBPlusS[1][1] = .771727;
    BOverBPlusS[1][2] = .877849;
    BOverBPlusS[2][0] = .31328;
    BOverBPlusS[2][1] = .433515;
    BOverBPlusS[2][2] = .598426;
    BOverBPlusS[3][0] = .172516;
    BOverBPlusS[3][1] = .246846;
    BOverBPlusS[3][2] = .330899;
    BOverBPlusS[4][0] = .243718;
    BOverBPlusS[4][1] = .191699;
    BOverBPlusS[4][2] = .230888;

    SOverBPlusS[0][0] = .116839;
    SOverBPlusS[0][1] = .038067;
    SOverBPlusS[0][2] = .016075;
    SOverBPlusS[1][0] = .476283;
    SOverBPlusS[1][1] = .228273;
    SOverBPlusS[1][2] = .122151;
    SOverBPlusS[2][0] = .686719;
    SOverBPlusS[2][1] = .566486;
    SOverBPlusS[2][2] = .401574;
    SOverBPlusS[3][0] = .827486;
    SOverBPlusS[3][1] = .753154;
    SOverBPlusS[3][2] = .669102;
    SOverBPlusS[4][0] = .756281;
    SOverBPlusS[4][1] = .808303;
    SOverBPlusS[4][2] = .769112;
    
    
   
    //----------------------------------------------------------------------------
    //--------------------Event Mixer Buffer-------------------------------------
    //----------------------------------------------------------------------------
    
    if(USE_VZ_BINS) { nVzBins = NUM_VZ_BINS; }  //sets the number of vz bins if flag is TRUE
    else nVzBins = 1;
    
    if(USE_CENT_BINS) { nCentBins = NUM_CENT_BINS; } //sets the number of multiplicity bins if flag is TRUE
    else nCentBins = 1;
    
    //if(D0_HADRON_CORR){
    
        for(int i = 0; i < nCentBins; i++){                             //this builds the buffer for the associated hadron events
            for( int j = 0; j < nVzBins; j++){
        
                eventBufferPicoEvent[j][i] = new StMixedEventBuffer();
                eventBufferPicoEvent[j][i]->setBufferSize(BUFFER_SIZE);    //set buffer size here -- the amount of events in each 2d bin
                eventBufferPicoEvent[j][i]->setBufferCounter(0);
            }
        }
    //}        
   
    //if(D0_HADRON_CORR){
   
        for(int i = 0; i < nCentBins; i++){                             //this builds the buffer to store the D0-candidate events (only the D0 candidate itself is stored).
            for( int j = 0; j < nVzBins; j++){
        
                eventBufferD0Candidate[j][i] = new StMixedEventBuffer();
                eventBufferD0Candidate[j][i]->setBufferSize(BUFFER_SIZE);    //set buffer size here -- the amount of events in each 2d bin
                eventBufferD0Candidate[j][i]->setBufferCounter(0);
            }
        }
    //}       
   
   
    //------------------------------------------------------------------   
    //--------------------CUTS------------------------------------------
    //------------------------------------------------------------------
    
    // D0 Cuts
   
    kaonPtCut           = .15;  //daughterKaon low pt acceptance cut
    pionPtCut           = .15;  //daughterPion low pt acceptance cut
    
    D0InvMassLow        = 1.82; //signal region US invariant mass Low cut
    D0InvMassHigh       = 1.90; //signal region US invariant mass High cut
    
    USSideBandLeftLow   = 1.70;  //left sideband low cut
    USSideBandLeftHigh  = 1.80;  //left sideband high cut
    USSideBandRightLow  = 1.92;   //right sideband low cut
    USSideBandRightHigh = 2.1;   //right sideband high cut
    
    d0PtLow             = 0.0;     //reconstructed D0 candidate low pt cut
    d0PtHigh            = 10.0;    //reconstructed D0 candidate high pt cut
    d0DecayLengthMax    = 999999.0;
    
    d0DaughterPionPtMin = .15;
    d0DaughterKaonPtMin = .15;
    
    //Old Average Cuts -- good
    //daughterDCA         = .0055;
    //d0DecayLengthMin    = .0180;
    //kaonDCA             = .008;
    //pionDCA             = .008;
    //d0DCAtoPV           = .0065;
    
    //d0TopologicalCutArray[5][5]; //first index is pt bin, second index is cuts -- decayLength, DCADaughters, d0DCAPV, PiDCAPV, KDCAPV
 
    //CUT SET 2 -- OPTIMIZED HF GROUP V2 CUTS	
	
    /* //pt = 0-1 GeV/c
    d0TopologicalCutArray[0][0] = .0145; //.0180;//minimum decay length 
    d0TopologicalCutArray[0][1] = .0084; //maximum DCA daughters 
    d0TopologicalCutArray[0][2] = .0061; //maximum DCA between d0 vector and PV 
    d0TopologicalCutArray[0][3] = .0110; //minimum DCA between Pi and PV 
    d0TopologicalCutArray[0][4] = .0103; //minimum DCA between K and PV
    
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
    d0TopologicalCutArray[4][4] = .0058; //minimum DCA between K and PV  */
    
    
    //-------------------------------------------------------------------------------------------
    //CUT SET 3 -- AVERAGE WIDE HIGH PT-BIN CUTS	
	
    //pt = 0-1 GeV/c
    d0TopologicalCutArray[0][0] = .0145; //.0180;//minimum decay length 
    d0TopologicalCutArray[0][1] = .0084; //maximum DCA daughters 
    d0TopologicalCutArray[0][2] = .0061; //maximum DCA between d0 vector and PV 
    d0TopologicalCutArray[0][3] = .0110; //minimum DCA between Pi and PV 
    d0TopologicalCutArray[0][4] = .0103; //minimum DCA between K and PV
    
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
    
    /* //pt = 3-5 GeV/c
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
    d0TopologicalCutArray[4][4] = .0058; //minimum DCA between K and PV  */
    
	
    //Associated Hadron Cuts
    
    hadronPtMin         = .15;
    hadronPtMax         = 15.45;
    
    trackChi2max        = 3.0;
    trackDCAtoPvtx      = 3.0;
    nHitsFitMin         = 20;
    nHitsFitMinOverMax  = .52;
    
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
    
    TString SibCorrLabels[4] = {"Sibling_SideBandLeft_correlation", "Sibling_US_correlation", "Sibling_SideBandRight_correlation", "Sibling_LS_correlation"};
    TString MixCorrLabels[4] = {"Mixed_SideBandLeft_correlation", "Mixed_US_correlation", "Mixed_SideBandRight_correlation", "Mixed_LS_correlation"};
    
    //other labels
   
    TString siblingLabel      = "Sibling_";
    TString mixLabel          = "Mixed_";
    
    TString effCorrInfoLabel  = "eff_corr_info_";
   
    TString eventCounterLabel = "Event Count Vz ";
    TString etaLabel          = "Inclusive Hadron Eta";
    TString phiLabel          = "Inclusive Hadron Phi";
    TString etaPhiLabel       = "Inclusive 2D Hadron Eta/Phi";
    
    TString phiD0vsPhiHLabel  = "phiD0_vs_phiH";
    TString etaD0vsEtaHLabel  = "etaD0_vs_etaH";
    TString phiD0vsEtaD0Label = "phiD0_vs_etaD0";
    TString phiHvsEtaHLabel   = "phiH_vs_etaH";
    
    TString bandLabels[4] = {"_SBL_", "_US_Signal_Region_", "_SBR_", "_LS_"};
    
	TString D0ptLabel[3]          = {"SBL_Candidate_Pt_Dist", "D0_Candidate_Pt_Dist", "SBR_Candidate_Pt_Dist"};
    TString D0EtaLabel[3]         = {"SBL_Eta_Dist", "D0_Eta_Dist", "SBR_Eta_Dist"};
    TString D0PhiLabel[3]         = {"SBL_Phi_Dist", "D0_Phi_Dist", "SBR_Phi_Dist"};
    TString D0KaonPtLabel[3]      = {"SBL_Kaon_Candidate_Pt_Dist", "D0_Kaon_Candidate_Pt_Dist", "SBR_Kaon_Candidate_Pt_Dist"};
    TString D0KaonEtaLabel[3]     = {"SBL_Kaon_Eta_Dist", "D0_Kaon_Eta_Dist", "SBR_Kaon_Eta_Dist"};
    TString D0KaonPhiLabel[3]     = {"SBL_Kaon_Phi_Dist", "D0_Kaon_Phi_Dist", "SBR_Kaon_Phi_Dist"};
    TString D0PionPtLabel[3]      = {"SBL_Candidate_Pion_Pt_Dist", "D0_Candidate_Pion_Pt_Dist", "SBR_Candidate_Pion_Pt_Dist"};
    TString D0PionEtaLabel[3]     = {"SBL_Pion_Eta_Dist", "D0_Pion_Eta_Dist", "SBR_Pion_Eta_Dist"};
    TString D0PionPhiLabel[3]     = {"SBL_Pion_Phi_Dist", "D0_Pion_Phi_Dist", "SBR_Pion_Phi_Dist"};
    TString D0DecayLengthLabel[3] = {"SBL_Decay_Length", "D0_Decay_Length", "SBR_Decay_Length"};
    TString D0DCAToPVLabel[3]     = {"SBL_DCA_to_PV", "D0_DCA_to_PV", "SBR_DCA_to_PV"};
    TString D0KaonPVLabel[3]      = {"SBL_Daughter_Kaon_DCA_to_PV", "D0_Daughter_Kaon_DCA_to_PV", "SBR_Daughter_Kaon_DCA_to_PV"};
    TString D0PionPVLabel[3]      = {"SBL_Daughter_Pion_DCA_to_PV", "D0_Daughter_Pion_DCA_to_PV", "SBR_Daughter_Pion_DCA_to_PV"};
    TString D0DaughterDCALabel[3] = {"SBL_Daughter_Pair_DCA", "D0_Daughter_Pair_DCA", "SBR_Daughter_Pair_DCA"};
    
    TString D0PtKaonVsPtPionLabel[3] = {"SBL_Daughter_Kaon_Pt_vs_Pion_Pt","D0_Daughter_Kaon_Pt_vs_Pion_Pt","SBR_Daughter_Kaon_Pt_vs_Pion_Pt"} ;
    TString D0PKaonVsPPionLabel[3] = {"SBL_Daughter_Kaon_P_vs_Pion_P", "D0_Daughter_Kaon_P_vs_Pion_P", "SBR_Daughter_Kaon_P_vs_Pion_P"};
    TString D0RawEtaVsRawPhiLabel[3] = {"SBL_Raw_delEta_Vs_delPhi", "D0_Raw_delEta_Vs_delPhi", "SBR_Raw_delEta_Vs_delPhi"};
    
    TString D0RawEtaVsRawPhiLeftPtBlobLabel[3] = {"SBL_Raw_delEta_Vs_delPhi_LeftPtBlob", "D0_Raw_delEta_Vs_delPhi_LeftPtBlob", "SBR_Raw_delEta_Vs_delPhi_LeftPtBlob"};
    TString D0RawEtaVsRawPhiRightPtBlobLabel[3] = {"SBL_Raw_delEta_Vs_delPhi_RightPtBlob", "D0_Raw_delEta_Vs_delPhi_RightPtBlob", "SBR_Raw_delEta_Vs_delPhi_RightPtBlob"};
    
    // TString D0CutPtBinsLabel[5] = {"0-1", "1-2", "2-3", "3-5", "5-10"};
    
    //TString D0CutPtBinsLabel[6] = {"0-1", "1-2", "2-3", "3-5", "5-10", "Integrated"};   //NOMINAL VALUES -- DON'T DELETE
    TString D0CutPtBinsLabel[4] = {"0-1", "1-2", "2-10", "Integrated"};   //NOMINAL VALUES -- DON'T DELETE
    
    //TString D0CutPtBinsLabel[5] = {"0-2", "2-3", "3-5", "5-10", "Integrated"};
    
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
   
   
   if(D0_HADRON_CORR){//begin D0-Hadron Correlation Histograms
        for(int band = 0; band < 4; band++){
            for(int i = 0; i < nVzBins; i++){ //Initialize all of the histograms for storing the sibling and mixed information
                for(int j = 0; j < nCentBins; j++){
        
                    str1 = SibCorrLabels[band] + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                    str2 = MixCorrLabels[band] + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                
                    sibCorrBin[band][i][j]             = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, (3*TMath::PiOver2())+phiBinShift);
                    mixCorrBin[band][i][j]             = new TH2D(str2, str2, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, (3*TMath::PiOver2())+phiBinShift);
                   
                    sibCorrBin[band][i][j]->Sumw2();
                    mixCorrBin[band][i][j]->Sumw2();
                }
            }
        } 
    }        
    
    if(D0_HADRON_CORR && USE_PT_BINS){//begin D0-Hadron Correlation Histograms -- USE PT BINS HERE
        for(int band = 0; band < 4; band++){
            for(int i = 0; i < nVzBins; i++){ 
                for(int j = 0; j < nCentBins; j++){
                    for(int k = 0; k < NUM_PT_BINS; k++){
                
                
                        str1 = SibCorrLabels[band] + PtBinLabel + binLabelPt[k] + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                        str2 = MixCorrLabels[band] + PtBinLabel + binLabelPt[k] + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                    
                        sibCorrBinPt[band][k][i][j] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, (3*TMath::PiOver2())+phiBinShift);
                        mixCorrBinPt[band][k][i][j] = new TH2D(str2, str2, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, (3*TMath::PiOver2())+phiBinShift);
                        
                        sibCorrBinPt[band][k][i][j]->Sumw2();
                        mixCorrBinPt[band][k][i][j]->Sumw2();
                    }
                }            
            }
        }    
    }//end D0-Hadron Correlation Histograms -- USE PT BINS HERE    
   
    for(int band = 0; band < NUM_D0_CUT_PT_BINS; band++){
		for(int i = 0; i < 3; i++){
        
		    str1 = D0ptLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
			D0ptDist[i][band]       = new TH1D(str1, str1, 250, 0, 10);
			str1 = D0EtaLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
			D0EtaDist[i][band]       = new TH1D(str1, str1, 250, -1, 1);
			str1 = D0PhiLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
			D0PhiDist[i][band]       = new TH1D(str1, str1, 250, -2*TMath::Pi(), 2*TMath::Pi());
			str1 = D0KaonPtLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
			D0KaonPtDist[i][band]    = new TH1D(str1, str1, 250, 0, 10);
			str1 = D0KaonEtaLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
			D0KaonEtaDist[i][band]   = new TH1D(str1, str1, 250, -1, 1);
			str1 = D0KaonPhiLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
			D0KaonPhiDist[i][band]   = new TH1D(str1, str1, 250, -2*TMath::Pi(), 2*TMath::Pi());
			str1 = D0PionPtLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
			D0PionPtDist[i][band]    = new TH1D(str1, str1, 250, 0, 10);
			str1 = D0PionEtaLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
			D0PionEtaDist[i][band]   = new TH1D(str1, str1, 250, -1, 1);
			str1 = D0PionPhiLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
			D0PionPhiDist[i][band]   = new TH1D(str1, str1, 250, -2*TMath::Pi(), 2*TMath::Pi());
			
			str1 = D0DecayLengthLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
			D0DecayLengthDist[i][band]       = new TH1D(str1, str1, 250, 0, .4);
			str1 = D0DCAToPVLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
			D0DCAToPVDist[i][band]       = new TH1D(str1, str1, 500, 0, .02);
			str1 = D0KaonPVLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
			D0KaonPVDist[i][band]    = new TH1D(str1, str1, 500, 0, .1);
			str1 = D0PionPVLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
			D0PionPVDist[i][band]   = new TH1D(str1, str1, 500, 0, .1);
			str1 = D0DaughterDCALabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
			D0DaughterDCADist[i][band]   = new TH1D(str1, str1, 500, 0, .02);
			
			str1 = D0PtKaonVsPtPionLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
            D0PtKaonVsPtPion[i][band] = new TH2D(str1, str1, 250, 0, 6, 250, 0, 6);
            str1 = D0PKaonVsPPionLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
            D0PKaonVsPPion[i][band] = new TH2D(str1, str1, 250, 0, 6, 250, 0, 6);
            str1 = D0RawEtaVsRawPhiLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
            D0RawEtaVsRawPhi[i][band] = new TH2D(str1, str1, 250, 0, 2, 250, 0, TMath::Pi());
            
            str1 = D0RawEtaVsRawPhiLeftPtBlobLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
            D0RawEtaVsRawPhiLeftPtBlob[i][band] = new TH2D(str1, str1, 250, 0, 2, 250, 0, TMath::Pi());
            str1 = D0RawEtaVsRawPhiRightPtBlobLabel[i] + PtBinLabel + D0CutPtBinsLabel[band];
            D0RawEtaVsRawPhiRightPtBlob[i][band] = new TH2D(str1, str1, 250, 0, 2, 250, 0, TMath::Pi());
        }
   }
   
    
    //Event QA Histograms
   
    eventCounter         = new TH1I("number_of_events_used", "number of events used", 6, 0, 6);
    trackCounter         = new TH1I("number_of_tracks_per_event", "number of tracks per event", 1500, 0, 1499);
    usedTracks           = new TH1I("Used_tracks", "Used tracks", 1500, 0, 1499);
    pVtxX                = new TH1D("X_position_of_pVtx", "X position of pVtx", 500, -10, 10);
    pVtxY                = new TH1D("Y_position_of_pVtx", "Y position of pVtx", 500, -10, 10);
    pVtxZ                = new TH1D("Z_position_of_pVtx", "Z position of pVtx", 500, -7, 7);
    
    TString vzCentBinLabel = "event_counts_per_Vz_Cent_bin";
    
    for(int band = 0; band < 4; band++){
        for(int k = 0; k < 5; k++){
            
            str1 = vzCentBinLabel + bandLabels[band] + PtBinLabel + binLabelPt[k];
            
            vZandCentBinPerEvent[band][k] = new TH2D(str1, str1, 16, 0, 16, 10, 0, 10);
            vZandCentBinPerEvent[band][k]->GetXaxis()->SetTitle("Centrality Bin");
            vZandCentBinPerEvent[band][k]->GetYaxis()->SetTitle("Vz Bin");
        }
    }
   
    for(int k = 0; k < 5; k++){
        
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
   
    }
   
   
   //Track QA Histograms
   kaonEtaDist           = new TH1D("Kaon_Eta_Distribution", "Kaon Eta Distribution", 250, -1, 1);
   pionEtaDist           = new TH1D("Pion_Eta_Distribution", "Pion Eta Distribution", 250, -1, 1);
   kaonPhiDist           = new TH1D("Kaon_Phi_Distribution", "Kaon Phi Distribution", 250, -2*TMath::Pi(), 2*TMath::Pi());
   pionPhiDist           = new TH1D("Pion_Phi_Distribution", "Pion Phi Distribution", 250, -2*TMath::Pi(), 2*TMath::Pi());
   DCAtoPrimaryVertex    = new TH1D("track_DCA_to_PV", "track DCA to PV", 500, 0.0, 10.0);
   DCAtoPrimaryVertexCut = new TH1D("track_DCA_to_PV_cut_check", "track DCA to PV cut check", 500, 0.0, 10.0);
   hadronPtDist          = new TH1D("Inclusive_Hadron_pt", "Inclusive Hadron pt", 250, 0, 10);
   hadronPhiDist         = new TH1D("Inclusive_Hadron_Phi", "Inclusive Hadron Phi", 250, -TMath::Pi(), TMath::Pi());
   hadronEtaDist         = new TH1D("Inclusive_Hadron_Eta", "Inclusive Hadron Eta", 250, -1, 1);
   hadronChi2            = new TH1D("Chi2_for_hadron_tracks", "Chi2 for hadron tracks", 500, 0, 10);
   hadronNHitsFit        = new TH1D("nHitsFit", "nHitsFit", 60, 0, 60);
   hadronNHitsFitOverMax = new TH1D("nHitsFitOverMax", "nHitsFitOverMax", 100, 0, 1);
   dEdxVsPt        = new TH2D("dEdx_vs_P", "dEdx_vs_P", 400, 0, 4, 400, 1.5, 4);
   invBetaVsPt     = new TH2D("inv_Beta_Vs_P", "#Beta^{-1} Vs. P", 400, 0, 4, 400, 0, 4);
   
   TPCOnlyAssociatedHadronPtDist[0] = new TH1D("TPC_Hadron_Pt_Peripheral", "TPC_Hadron_Pt_Peripheral", 100, 0, 10);
   TPCOnlyAssociatedHadronPtDist[1] = new TH1D("TPC_Hadron_Pt_MidCentral", "TPC_Hadron_Pt_MidCentral", 100, 0, 10); 
   TPCOnlyAssociatedHadronPtDist[2] = new TH1D("TPC_Hadron_Pt_Central", "TPC_Hadron_Pt_Central", 100, 0, 10);      
   HFTAssociatedHadronPtDist[0]     = new TH1D("HFT_Hadron_Pt_Peripheral", "HFT_Hadron_Pt_Peripheral", 100, 0, 10);     //for hft ratio
   HFTAssociatedHadronPtDist[1]     = new TH1D("HFT_Hadron_Pt_MidCentral", "HFT_Hadron_Pt_MidCentral", 100, 0, 10);
   HFTAssociatedHadronPtDist[2]     = new TH1D("HFT_Hadron_Pt_Central", "HFT_Hadron_Pt_Central", 100, 0, 10);
  
   d0CountPerEvent = new TH1I("number_of_D0_candidates_per_event", "number of D0 candidates per event", 50, 0, 50);
   histOfCuts      = new TH1D("HistOfCuts", "HistOfCuts", 55, 1, 55); 

   
   //----Efficiency Mapping Functions--------------------
   //-----------------------THIS IS THE VERY PRELIMINARY FUNCTION USED TO WEIGHT THE PAIRS BASED ON HADRON EFFICIENCY
   effWeightPions[0]  = new TF1("effWeightPionsPer", "[0]*TMath::Exp(-([1]/x)**[2])", 0, 15.45);
   effWeightPions[1]  = new TF1("effWeightPionsMid", "[0]*TMath::Exp(-([1]/x)**[2])", 0, 15.45);
   effWeightPions[2]  = new TF1("effWeightPionsCent", "[0]*TMath::Exp(-([1]/x)**[2])", 0, 15.45);
  
   effWeightPions[0]->SetParameter(0, .829);
   effWeightPions[0]->SetParameter(1, .12);
   effWeightPions[0]->SetParameter(2, 4.146);
  
   effWeightPions[1]->SetParameter(0, .809);
   effWeightPions[1]->SetParameter(1, .109);
   effWeightPions[1]->SetParameter(2, 3.224);
   
   effWeightPions[2]->SetParameter(0, .77);
   effWeightPions[2]->SetParameter(1, .0835);
   effWeightPions[2]->SetParameter(2, 1.8415);
   
   effWeightKaons[0]  = new TF1("effWeightKaonsPer", "[0]*TMath::Exp(-([1]/x)**[2]) + [3]*x", 0, 15.45);
   effWeightKaons[1]  = new TF1("effWeightKaonsPer", "[0]*TMath::Exp(-([1]/x)**[2]) + [3]*x", 0, 15.45);
   effWeightKaons[2]  = new TF1("effWeightKaonsPer", "[0]*TMath::Exp(-([1]/x)**[2]) + [3]*x", 0, 15.45);
   
   effWeightKaons[0]->SetParameter(0, .5965);
   effWeightKaons[0]->SetParameter(1, .236);
   effWeightKaons[0]->SetParameter(2, 2.7295);
   effWeightKaons[0]->SetParameter(3, .092);
  
   effWeightKaons[1]->SetParameter(0, .503);
   effWeightKaons[1]->SetParameter(1, .231);
   effWeightKaons[1]->SetParameter(2, 3.968);
   effWeightKaons[1]->SetParameter(3, .152);
   
   effWeightKaons[2]->SetParameter(0, .472);
   effWeightKaons[2]->SetParameter(1, .229);
   effWeightKaons[2]->SetParameter(2, 3.7085);
   effWeightKaons[2]->SetParameter(3, .144);
   
   /*effWeightD0     = new TF1("effWeightPD0", "TMath::Exp([0]*([1]-[2]*TMath::Exp(-x/[3])))", 0, 20); 
   effWeightD0->SetParameter(0, 2.30259);
   effWeightD0->SetParameter(1, -1.306);
   effWeightD0->SetParameter(2, 1.9760001);
   effWeightD0->SetParameter(3, 3.3699999);*/ //old efficiency correction as of 6/9/2017, fails at lowest pt
   
   effWeightD0[0] = new TF1("effWeightPD0Per","pol4", .5, 10); //Using this to correct efficiency at low pt 6/9/2017
   effWeightD0[1] = new TF1("effWeightPD0Mid","pol4", .5, 10); //Using this to correct efficiency at low pt 6/9/2017
   effWeightD0[2] = new TF1("effWeightPD0Cent","pol4", .5, 10); //Using this to correct efficiency at low pt 6/9/2017
   
   effWeightD0[0]->FixParameter(0, 0.002258);
   effWeightD0[0]->FixParameter(1, -0.000672);
   effWeightD0[0]->FixParameter(2, 0.000689);
   effWeightD0[0]->FixParameter(3, 0.00003547);
   effWeightD0[0]->FixParameter(4, -0.00000714933);
    
   effWeightD0[1]->FixParameter(0, 0.004414);
   effWeightD0[1]->FixParameter(1, -0.007136);
   effWeightD0[1]->FixParameter(2, 0.003958);
   effWeightD0[1]->FixParameter(3, -0.000377);
   effWeightD0[1]->FixParameter(4, 0.00000954214);
    
   effWeightD0[2]->FixParameter(0, 0.001541);
   effWeightD0[2]->FixParameter(1, -0.000731);
   effWeightD0[2]->FixParameter(2, 0.000957);
   effWeightD0[2]->FixParameter(3, -0.00003644);
   effWeightD0[2]->FixParameter(4, -0.0000022057);   
   
   HFTRatioFunc[0] = new TF1("HFTRatioPer","pol6", 0.0, 15.45); //Using this to correct efficiency at low pt 6/9/2017
   HFTRatioFunc[1] = new TF1("HFTRatioMid","pol6", 0.0, 15.45); //Using this to correct efficiency at low pt 6/9/2017
   HFTRatioFunc[2] = new TF1("HFTRatioCent","pol6", 0.0, 15.45); //Using this to correct efficiency at low pt 6/9/2017
   
   HFTRatioFunc[0]->FixParameter(0, 0.642748);
   HFTRatioFunc[0]->FixParameter(1, -0.294797);
   HFTRatioFunc[0]->FixParameter(2, 0.217882);
   HFTRatioFunc[0]->FixParameter(3, -0.0792603);
   HFTRatioFunc[0]->FixParameter(4, 0.0152961);
   HFTRatioFunc[0]->FixParameter(5, -0.00151916);
   HFTRatioFunc[0]->FixParameter(6, .0000608388);
   
   HFTRatioFunc[1]->FixParameter(0, 0.806221);
   HFTRatioFunc[1]->FixParameter(1, -0.717618);
   HFTRatioFunc[1]->FixParameter(2, 0.615353);
   HFTRatioFunc[1]->FixParameter(3, -0.266492);
   HFTRatioFunc[1]->FixParameter(4, 0.0617126);
   HFTRatioFunc[1]->FixParameter(5, -0.00730829);
   HFTRatioFunc[1]->FixParameter(6, 0.000346782);
   
   HFTRatioFunc[2]->FixParameter(0, 0.945512);
   HFTRatioFunc[2]->FixParameter(1, -0.933711);
   HFTRatioFunc[2]->FixParameter(2, 0.736739);
   HFTRatioFunc[2]->FixParameter(3, -0.29807);
   HFTRatioFunc[2]->FixParameter(4, 0.0649339);
   HFTRatioFunc[2]->FixParameter(5, -0.00729887);
   HFTRatioFunc[2]->FixParameter(6, 0.00033193);
  
   
   
   //levy = new TF1("Levy", "((x*x)/(TMath::Sqrt((x*x)+(1.86*1.86))))*([0]/((1+((TMath::Sqrt((x*x)+(1.86*1.86))-1.86)/([1]*[2])))**[1]))", 0, 10);
   
   
   //single particle distributions for errors
    if(SINGLES_DISTS){
        for(int i = 0; i < nVzBins; i++){
            for(int j = 0; j < nCentBins; j++){
                for(int k = 0; k < NUM_PT_BINS; k++){
                    for(int band = 0; band < 4; band++){
        
                        str1 = siblingLabel + bandLabels[band] +  phiD0vsPhiHLabel  + PtBinLabel + binLabelPt[k] + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                        str2 = siblingLabel + bandLabels[band] +  etaD0vsEtaHLabel  + PtBinLabel + binLabelPt[k] + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                        str3 = siblingLabel + bandLabels[band] +  phiD0vsEtaD0Label + PtBinLabel + binLabelPt[k] + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                        str4 = siblingLabel + bandLabels[band] +  phiHvsEtaHLabel   + PtBinLabel + binLabelPt[k] + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                
                        phiD0vsEtaD0[band][k][i][j]    = new TH2D(str3, str3, 9, -1, 1, 24, -TMath::Pi()+phiBinShiftPrimary, TMath::Pi()+phiBinShiftPrimary);
                        phiHvsEtaH[band][k][i][j]      = new TH2D(str4, str4, 9, -1, 1, 24, -TMath::Pi()+phiBinShiftPrimary, TMath::Pi()+phiBinShiftPrimary);
                        
                        
                    }           
                }
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
    
    /* histOfCuts->SetBinContent(27, d0TopologicalCutArray[3][0]); //d0DecayLengthMin
    histOfCuts->SetBinContent(28, d0TopologicalCutArray[3][1]); //daughterDCA
    histOfCuts->SetBinContent(29, d0TopologicalCutArray[3][2]); //d0DCAtoPV
    histOfCuts->SetBinContent(30, d0TopologicalCutArray[3][3]); //pionDCA
    histOfCuts->SetBinContent(31, d0TopologicalCutArray[3][4]); //kaonDCA
    
    histOfCuts->SetBinContent(32, d0TopologicalCutArray[4][0]); //d0DecayLengthMin
    histOfCuts->SetBinContent(33, d0TopologicalCutArray[4][1]); //daughterDCA
    histOfCuts->SetBinContent(34, d0TopologicalCutArray[4][2]); //d0DCAtoPV
    histOfCuts->SetBinContent(35, d0TopologicalCutArray[4][3]); //pionDCA
    histOfCuts->SetBinContent(36, d0TopologicalCutArray[4][4]); //kaonDCA */
    
    
    histOfCuts->SetBinContent(37, hadronPtMin); //hadronPtMin
    histOfCuts->SetBinContent(38, hadronPtMax); //hadronPtMax
    histOfCuts->SetBinContent(39, BUFFER_SIZE); //BUFFER_SIZE
    histOfCuts->SetBinContent(40, NUM_PHI_BINS); //NUM_PHI_BINS
    histOfCuts->SetBinContent(41, NUM_ETA_BINS); //NUM_ETA_BINS
    histOfCuts->SetBinContent(42, 1); // Require HFT
    histOfCuts->SetBinContent(43, trackChi2max); // chi2 cut
    histOfCuts->SetBinContent(44, trackDCAtoPvtx); // DCA cut
    histOfCuts->SetBinContent(49, 1); //Scale factor counter to produce text file
    histOfCuts->SetBinContent(50, int(USE_TOF)); //USE TOF STRICTLY
    histOfCuts->SetBinContent(51, int(USE_TOF_HYBRID)); //USE TOF HYBRID
    histOfCuts->SetBinContent(52, int(USE_DOUBLE_MIS_PID_PROTECTION)); //USE DOUBLE MIS-PID PROTECTION
    
   
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
   
    if(D0_HADRON_CORR){
        for(int band = 0; band < 4; band++){
            for(int i = 0; i < nVzBins; i++){
                for(int j = 0; j < nCentBins; j++){
        
                    sibCorrBin[band][i][j]->Write();
                    mixCorrBin[band][i][j]->Write();
           
                    for(int k = 0; k < NUM_PT_BINS; k++){
            
                        sibCorrBinPt[band][k][i][j]->Write();
                        mixCorrBinPt[band][k][i][j]->Write();
                    }
                }    
            }
        }
    }        
    
    for(int i = 0; i < nVzBins; i++){
        for(int j = 0; j < nCentBins; j++){
            cout << "deleting buffer[" << i << "][" << j << "]" << endl;
            delete eventBufferPicoEvent[i][j];
            delete eventBufferD0Candidate[i][j];
        }
    }

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
                for(int k = 0; k < NUM_PT_BINS; k++){
                    for(int band = 0; band < 4; band++){ 
                
                        phiD0vsEtaD0[band][k][i][j]->Write();
                        phiHvsEtaH[band][k][i][j]->Write();
                        
                    }                    
                }
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
   
   //D0ptDist->Write();
   
   for(int band = 0; band < NUM_D0_CUT_PT_BINS; band++){
		for(int i = 0; i < 3; i++){
			
			D0ptDist[i][band]->Write();
			D0EtaDist[i][band]->Write();
			D0PhiDist[i][band]->Write();
			D0KaonPtDist[i][band]->Write();
			D0KaonEtaDist[i][band]->Write(); 
			D0KaonPhiDist[i][band]->Write(); 
			D0PionPtDist[i][band]->Write();
			D0PionEtaDist[i][band]->Write(); 
			D0PionPhiDist[i][band]->Write(); 
			D0DecayLengthDist[i][band]->Write();
			D0DCAToPVDist[i][band]->Write();
			D0KaonPVDist[i][band]->Write();
			D0PionPVDist[i][band]->Write();
			D0DaughterDCADist[i][band]->Write();
            D0PtKaonVsPtPion[i][band]->Write();
            D0PKaonVsPPion[i][band]->Write();
            D0RawEtaVsRawPhi[i][band]->Write();
            D0RawEtaVsRawPhiLeftPtBlob[i][band]->Write();
            D0RawEtaVsRawPhiRightPtBlob[i][band]->Write();
		}
	}
       
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
    for(int band = 0; band < 4; band++){
        for(int k = 0; k < 5; k++){
        
            vZandCentBinPerEvent[band][k]->Write();
        }
    }
    usedTracks->Write();
    histOfCuts->Write();
    pVtxX->Write();
    pVtxY->Write();
    pVtxZ->Write();
    hadronChi2->Write();
    hadronNHitsFit->Write();
    hadronNHitsFitOverMax->Write();
    
    for(int k = 0; k < 5; k++){
        for(int i = 0; i < 3; i++){
        
            efficiencyCorrInfoSib[k][i]->Write();
            efficiencyCorrInfoMix[k][i]->Write();
        }
    }
    
    for(int i = 0; i < 3; i++){
    
        TPCOnlyAssociatedHadronPtDist[i]->Write();
        HFTAssociatedHadronPtDist[i]->Write();
    }

   
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
    bandBin = -1; 
    centralityClass = -1;
    topologicalCutPtBin = -1;
	ptBin = -1;
	
    //--------------Local Variables used in the code-------------------
    double delPhi          = 0;
    double delPhiCp        = 0;
    double delEta          = 0;
    double pt              = 0;
    double phi             = 0;
    double kpPhi           = 0;
    double eta             = 0;
    double dEdx            = 0;
    double beta            = 0;
    double nHitsFitRatio   = 0;
    int    PIDflag         = 0;
    int    numD0s          = 0;
    bool   minBiasFlag     = false;
    int    realTracks      = 0;
    int    realTracksFinal = 0;
    
    //Efficiency correction stuff -- this all needs to be updated
    
    double hadronEffWeight = 1;     //this is the efficiency weight for the hadrons, based on the exponential function
    double pionEffWeight = 1;
    double kaonEffWeight = 1;
    
    double tripleProduct;
    double doubleProduct;
    
    double backgroundEffWeight = 1;
    //double D0HadronEffWeight   = 1;
    ////////////////////////////////////////////////////////////////////
    
    double delEtaDaughters = 0;
    double delPhiDaughters = 0;
	
    double D0Eff = 1;
    double HFTRatioKaon = 1; 
    double HFTRatioPion = 1; 
    double HFTRatioHadron = 1;
    //double oneOverD0Eff = 1;
   
    double daughterPionPt = 0;
    double daughterKaonPt = 0;
	int    kPiCharge      = 0;   //US is negative and LS is positive
    double finalPairWeight = 1;
    
    //----------------------------------------------------------------------------------------------------------------
    
    bool storedD0Event    = false;
    bool associatedHadronIsDaughter = false;
   
    StPicoTrack* trk;
    StThreeVectorF trackMom;
    StThreeVectorF trackMom2;
    StThreeVectorF daughterKaonMom;
    StThreeVectorF daughterPionMom;
    double bField         = picoDst->event()->bField();
    StThreeVectorF pVtx   = picoDst->event()->primaryVertex();
    StThreeVectorF kaonPionMom;
    
    std::vector <StThreeVectorF> mAssociatedHadronList;
    std::vector <int> mAssociatedHadronListIndex;
    std::vector <int> mAssociatedHadronListCharge;
    
    bool storeEventToMix        = true;
    bool storeKaonPionEvent     = false;
    int d0Counter = 0;
    
/****************************************************************************************************************/
/****************************GLOBAL AND LOCAL VARIABLES SET AND INITIALIZED END********************************/  
/***************************************************************************************************************/   
    
 
/****************************************************************************************************************/
/****************************PRELIMINARY EVENT CUTS BEGIN********************************************************/  
/***************************************************************************************************************/    
    //if(!mPicoD0Event->isGoodEvent()) { return kStOk; }    //apply basic event cuts here
 
    if(USE_VZ_BINS){
            Vz = picoDst->event()->primaryVertex().z();
            VzBin = getVzBin(Vz);                  //get Vz bin
        }
   
        else VzBin = 0;
    
    if(VzBin == -1) { return kStOk; }
 
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
/****************************PRELIMINARY EVENT CUTS END*********************************************************/  
/***************************************************************************************************************/ 

if(DEBUG){ 
        
        cout << endl << endl;
        cout << "******************************************************************************" << endl;
        cout << "*******************************EVENT START************************************" << endl;
        cout << "We are on event # " << eventNumber << endl << endl;
        eventNumber++;
}        


/****************************************************************************************/
/***************************Initial check for D0 candidate begin*************************/  
/****************************************************************************************/  
  
    for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx){//begin loop to check if event can be used (check if a D0 exists)
   
        StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
        StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
        StPicoTrack const* pion = picoDst->track(kp->pionIdx());
                   
        daughterKaonMom = kaon->gMom(pVtx, bField);
        daughterPionMom = pion->gMom(pVtx, bField);           
        kPiCharge = kaon->charge()*pion->charge();           
       
        //bandBin = getUSInvMassBandBin(kp->m(), kPiCharge); //0 - SBR, 1 - D0 cand. , 2 - SBL, 3 LS   
        bandBin = getUSInvMassBandBin(kp->m(), kPiCharge, USSideBandLeftLow, USSideBandLeftHigh, D0InvMassLow, D0InvMassHigh, USSideBandRightLow,USSideBandRightHigh);		
        topologicalCutPtBin = getTopologicalCutPtBin(kp->pt());
      
	    if(DEBUG) { cout << "bandBin: " << bandBin << endl; }
	  
        if(topologicalCutPtBin == -1 ) {continue;}
        
        if(USE_DOUBLE_MIS_PID_PROTECTION && (((fabs(kaon->nSigmaKaon()) <= 2.0) && (fabs(kaon->nSigmaPion()) <= 2.0)) &&
           ((fabs(pion->nSigmaKaon()) <= 2.0) && (fabs(pion->nSigmaPion()) <= 2.0)))){ continue; }   //TPC nSigma double misPID protection  
           
		if((kaon->nHitsFit() < nHitsFitMin || !kaon->isHFTTrack()) || (pion->nHitsFit() < nHitsFitMin || !pion->isHFTTrack()) ) { continue; }
        if((fabs(kaon->nSigmaKaon()) > 2.0) || (fabs(pion->nSigmaPion()) > 2.0)){ continue; }   //good TPC track   
         
        if(USE_TOF && !mHFCuts->isTOFKaon(kaon, mHFCuts->getTofBeta(kaon))) {continue;}
        if(USE_TOF && !mHFCuts->isTOFPion(pion, mHFCuts->getTofBeta(pion))) {continue;}
        
        if(USE_TOF_HYBRID && (mHFCuts->getTofBeta(kaon) > 0 && !mHFCuts->isTOFKaon(kaon, mHFCuts->getTofBeta(kaon)))) {continue;}
        if(USE_TOF_HYBRID && (mHFCuts->getTofBeta(pion) > 0 && !mHFCuts->isTOFPion(pion, mHFCuts->getTofBeta(pion)))) {continue;}
        
       
		
	    //  ptmin ptmax   decayLenMin&Max   daughterDCA kaon/pion pt kaon/pion DCA  DCA to PV
        if((bandBin == 1 || bandBin == 3) &&                              //Check for candidate in peak range (LS & US)
            cutCheck(kp, d0PtLow,  d0PtHigh,  d0TopologicalCutArray[topologicalCutPtBin][0],  d0DecayLengthMax,  d0TopologicalCutArray[topologicalCutPtBin][1],
                    d0DaughterPionPtMin, d0DaughterKaonPtMin, d0TopologicalCutArray[topologicalCutPtBin][4], d0TopologicalCutArray[topologicalCutPtBin][3],
                    d0TopologicalCutArray[topologicalCutPtBin][2])) {

                    storeEventToMix = false;
                    storeKaonPionEvent = true;                     
                    if(DEBUG) { cout << "D0 Candidate information" << endl << endl; cout << "Mass: "<< kp->m() << "       US (-1) or LS (1): " << kPiCharge << endl << endl;     
                                cout << "Kaon Index: " << kp->kaonIdx() << "    Pion Index: " << kp->pionIdx() << endl << endl;}
                                
                    mAssociatedHadronListIndex.push_back(kp->kaonIdx()); //stores index of candidate daughter so it is not used as an associated track
                    mAssociatedHadronListIndex.push_back(kp->pionIdx()); //stores index of candidate daughter so it is not used as an associated track
                     
                    if(kPiCharge < 0 ) {numD0s++;}
        }     

       else if(bandBin == 0 && kPiCharge < 0 &&             //Check for candidate in left side band (US)
            cutCheck(kp, d0PtLow,  d0PtHigh,  d0TopologicalCutArray[topologicalCutPtBin][0],  d0DecayLengthMax,  d0TopologicalCutArray[topologicalCutPtBin][1],
                     d0DaughterPionPtMin, d0DaughterKaonPtMin, d0TopologicalCutArray[topologicalCutPtBin][4], d0TopologicalCutArray[topologicalCutPtBin][3],
                     d0TopologicalCutArray[topologicalCutPtBin][2])) {

                     if(DEBUG){ cout << "Left side band Candidate information" << endl << endl; cout << "Mass: "<< kp->m() << endl << endl;  
                                cout << "Kaon Index: " << kp->kaonIdx() << "    Pion Index: " << kp->pionIdx() << endl << endl;}
                     storeEventToMix = false;
                     storeKaonPionEvent = true;  

                     mAssociatedHadronListIndex.push_back(kp->kaonIdx()); //stores index of candidate daughter so it is not used as an associated track
                     mAssociatedHadronListIndex.push_back(kp->pionIdx()); //stores index of candidate daughter so it is not used as an associated track                     
            }               
            
      else if(bandBin == 2 && kPiCharge < 0 &&          //Check for candidate in right side band (US)
            cutCheck(kp, d0PtLow,  d0PtHigh,  d0TopologicalCutArray[topologicalCutPtBin][0],  d0DecayLengthMax,  d0TopologicalCutArray[topologicalCutPtBin][1],
                     d0DaughterPionPtMin, d0DaughterKaonPtMin, d0TopologicalCutArray[topologicalCutPtBin][4], d0TopologicalCutArray[topologicalCutPtBin][3],
                     d0TopologicalCutArray[topologicalCutPtBin][2])) {

                    if(DEBUG){ cout << "Right side band Candidate information" << endl << endl; cout << "Mass: "<< kp->m() << endl << endl; 
                               cout << "Kaon Index: " << kp->kaonIdx() << "    Pion Index: " << kp->pionIdx() << endl << endl;} 
            
                    storeEventToMix = false;
                    storeKaonPionEvent = true; 

                    mAssociatedHadronListIndex.push_back(kp->kaonIdx()); //stores index of candidate daughter so it is not used as an associated track
                    mAssociatedHadronListIndex.push_back(kp->pionIdx()); //stores index of candidate daughter so it is not used as an associated track     

                                         
            }                     
        
    }//end loop to check if event can be used (check if a D0 exists)        
    
    if(DEBUG) {cout << "Number of D0 candidates in this event: " << numD0s << endl << endl; }  
    if(DEBUG) {cout << "Number of daughters from kPi pairs (signal or SB): " << mAssociatedHadronListIndex.size() << endl << endl; } 
  
/*****************************************************************************************************************/
/*****************************************Initial check for D0 candidate end**************************************/
/*****************************************************************************************************************/  
  
  
   
/****************************************************************************************************************/  
/********************QA Loop for storing all of the basic information about the events and tracks BEGIN**********/
/****************************************************************************************************************/   
   
    for(unsigned int i = 0; i < picoDst->numberOfTracks(); ++i){ // Begin associated hadron list loop
                         
        associatedHadronIsDaughter = false;
                         
        trk = picoDst->track(i);
        trackMom = trk->gMom(pVtx, bField);
        nHitsFitRatio = (float)trk->nHitsFit()/((float)trk->nHitsMax());
        trackDCA = ((trk->helix().origin())-pVtx).mag();
        dEdx = trk->dEdx();
        phi = trackMom.phi();    
        eta  = trackMom.pseudoRapidity();    
        pt   = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y()));    
            
        if(trk->chi2() > trackChi2max || !checkDCAtoPV(trackDCA) || (trk->nHitsFit() < nHitsFitMin) ) { continue; } //check chi2 cut -- track quality cut for DCA to PV --num fit points check
        if(nHitsFitRatio < nHitsFitMinOverMax)                                                                          {continue;} //num fit points over maximum number check
        if(pt < hadronPtMin || pt > hadronPtMax)                                                                        {continue;}  //check pt cut
        if(trackMom.mag() > .20 && trackMom.mag() < .45 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) {continue;}  //remove electron contamination
        if(trackMom.mag() > .70 && trackMom.mag() < .80 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) {continue;}  //"                            "
        if(eta > 1.0 || eta < -1.0) { continue; }
        //if(picoDst->event()->primaryVertex().z() > 4.0 || picoDst->event()->primaryVertex().z() < -4.0) { continue; } //only for Xin test.
	
        realTracks++;                        //after a track passes all cuts, it is considered a "real track"

        //if(pt < hadronPtMin || pt > 1.0)                                                                        {continue;}  //check pt cut ADDED IN
        if(!mHFCuts->isGoodTrack(trk)) { continue; }    //Checks the HFT requirement for a track
        
        realTracksFinal++;                        //this is for comparing the # of associated hadrons with and without HFT requirement

        DCAtoPrimaryVertex->Fill(trackDCA);  
        dEdxVsPt->Fill(trackMom.mag(), dEdx);
        DCAtoPrimaryVertexCut->Fill(trackDCA);
        hadronChi2->Fill(trk->chi2());
        hadronNHitsFit->Fill(trk->nHitsFit());
        hadronNHitsFitOverMax->Fill(nHitsFitRatio);
        hadronPtDist->Fill(pt); 
        hadronPhiDist->Fill(phi);        
        hadronEtaDist->Fill(eta);
        
        if(mHFCuts->hasTofPid(trk)){    
         
            beta = mHFCuts->getTofBeta(trk);                                                
            invBetaVsPt->Fill(trackMom.mag(), (1/beta));        
        }
        
        if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaKaon()) < 2.0)){       //ONLY check nSigma for TPC track -- used to be mHFCuts->isTPCKaon(trk), but this include pt cut
                                                                                  //need to fix this. Need to figure out how to access the nSigma from the cuts
              kaonEtaDist->Fill(eta);
              kaonPhiDist->Fill(phi);
        }

        else if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaPion()) < 2.0)){     //ONLY check nSigma for TPC track
               
              pionEtaDist->Fill(eta);
              pionPhiDist->Fill(phi);
        }
        
        for(int index = 0; index < mAssociatedHadronListIndex.size(); index++){ //this loop sets the associatedList up to NOT store a kPi candidate daughter
        
            if(mAssociatedHadronListIndex[index] == i) {associatedHadronIsDaughter = true;}
            
        }    
        
        if(associatedHadronIsDaughter) {continue;}
        
        mAssociatedHadronList.push_back(trackMom);
        mAssociatedHadronListCharge.push_back(trk->charge());
    
    } //End associated hadron list loop  
    
    
        
    
    //Fill track distributions and get centrality bins////////////////////
    
        usedTracks->Fill(realTracksFinal);   
    
        centralityBin = getCentralityBin(realTracks);  //get multiplicity bin
        centralityClass = getCentralityClass(realTracks);
        
        if(centralityBin == -1 || VzBin == -1) { return kStOk; }
       
        pVtxX->Fill(picoDst->event()->primaryVertex().x());
        pVtxY->Fill(picoDst->event()->primaryVertex().y());
        pVtxZ->Fill(picoDst->event()->primaryVertex().z()); 
   
        //vZandCentBinPerEvent->Fill(centralityBin, VzBin); //moved to only fill for D0-events
        
        
    for(unsigned int i = 0; i < picoDst->numberOfTracks(); ++i){ // Begin loop to get pt distributions for TPC and HFT tracks to calculate ratio function
                         
        trk = picoDst->track(i);
        trackMom = trk->gMom(pVtx, bField);
        nHitsFitRatio = (float)trk->nHitsFit()/((float)trk->nHitsMax());
        trackDCA = ((trk->helix().origin())-pVtx).mag();
        dEdx = trk->dEdx();
        phi = trackMom.phi();    
        eta  = trackMom.pseudoRapidity();    
        pt   = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y()));    
            
        if(trk->chi2() > trackChi2max || !checkDCAtoPV(trackDCA) || (trk->nHitsFit() < nHitsFitMin) ) { continue; } //check chi2 cut -- track quality cut for DCA to PV --num fit points check
        if(nHitsFitRatio < nHitsFitMinOverMax)                                                                          {continue;} //num fit points over maximum number check
        if(pt < hadronPtMin || pt > hadronPtMax)                                                                        {continue;}  //check pt cut
        if(trackMom.mag() > .20 && trackMom.mag() < .45 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) {continue;}  //remove electron contamination
        if(trackMom.mag() > .70 && trackMom.mag() < .80 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) {continue;}  //"                            "
        if(eta > 1 || eta < -1) { continue; }
	
        TPCOnlyAssociatedHadronPtDist[centralityClass]->Fill(pt);
        
        if(!mHFCuts->isGoodTrack(trk)) { continue; }    //Checks the HFT requirement for a track
        
        HFTAssociatedHadronPtDist[centralityClass]->Fill(pt);     //for hft ratio

    } //End loop to get pt distributions for TPC and HFT tracks to calculate ratio function
    
    
    ///////////////////////////////////////////////////////////////////
    
    if(DEBUG){ 
        
        cout << endl << endl;
        //cout << "******************************************************************************" << endl;
        //cout << "*******************************EVENT START************************************" << endl;
        //cout << "We are on event # " << eventNumber << endl;
        cout << "This event has " << realTracks << " \"good\" quality tracks." << endl;
        cout << endl << "Number of associated tracks: " << mAssociatedHadronList.size() << endl; 
        cout << "Vz: " << Vz << "   VzBin: " << VzBin << "    Centrality Bin: " << centralityBin << endl; 
        cout << endl << endl;
        //eventNumber++;
    }         
           
        
/****************************************************************************************************************/  
/********************QA Loop for storing all of the basic information about the events and tracks END**********/
/****************************************************************************************************************/         
 

//***************************************************************************************
//BEGIN D0 BLOCK OF CODE HERE -- MAKE SURE THE NON_IDENTIFIED FLAG IS SET TO FALSE!!!!
//***************************************************************************************    

//Still need to add checks for TOF information-----   

/****************************************************************************************************************/  
/***************************************LOOP TO STORE EVENTS IN MIXER BEGIN***************************************/
/****************************************************************************************************************/   
   
   
    
    if(EVENT_MIXING){  // begin event mixing on/off switch
        
        if(eventBufferPicoEvent[VzBin][centralityBin]->getBufferSize() == BUFFER_SIZE) { storeEventToMix = false; }
        
        if(storeEventToMix){ //Begin block to store non-d0 containing events 
    
            eventBufferPicoEvent[VzBin][centralityBin]->addEvent(picoDst);
               
            for(unsigned int i = 0; i < picoDst->numberOfTracks(); ++i){//begin loop to add picoDST tracks to buffer
                         
                        trk = picoDst->track(i);                                                             
                        trackMom = trk->gMom(pVtx, bField); 
                        if(!mHFCuts->isGoodTrack(trk)) { continue; }
                        
                        trackDCA = ((trk->helix().origin())-pVtx).mag();
        
                        if(trk->chi2() > trackChi2max) { continue; }  //check chi2 cut
                        if(!checkDCAtoPV(trackDCA))    { continue; }   // track quality cut for DCA to PV
                        if(trk->nHitsFit() < nHitsFitMin) { continue; }      //num fit points check
                        nHitsFitRatio = (float)trk->nHitsFit()/((float)trk->nHitsMax());
        
                        if(nHitsFitRatio < nHitsFitMinOverMax) { continue; } //num fit points over maximum number check
                        
                        pt   = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y()));
                        if(pt < hadronPtMin || pt > hadronPtMax){continue;}   
                        
                        if(trackMom.mag() > .20 && trackMom.mag() < .45 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }
                        if(trackMom.mag() > .70 && trackMom.mag() < .80 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }
                        
                        eta  = trackMom.pseudoRapidity();
                        if(eta > 1 || eta < -1) { continue; }
                        
                        if((fabs(trk->nSigmaKaon()) < 2.0)){ PIDflag = 1; }
                        else if((fabs(trk->nSigmaPion()) < 3.0)) { PIDflag = 2; }
                        else PIDflag = 0;
            
                        eventBufferPicoEvent[VzBin][centralityBin]->addTrackToEvent(eventBufferPicoEvent[VzBin][centralityBin]->getBufferIndex()-1, trackMom, trk->charge(), PIDflag);   //0 -- any hadron, 1 -- kaon, 2 -- pion
                  
                    }//end loop to add picoDST tracks to buffer
                    
                    if(DEBUG){
                         
                        int bufIndex = 0;
                        bufIndex = eventBufferPicoEvent[VzBin][centralityBin]->getBufferIndex();                        
                        cout << "Event stored in buffer." << "Num of tracks = " << eventBufferPicoEvent[VzBin][centralityBin]->getEvent(bufIndex-1)->getNoTracks() << endl;
                   }
                   
        } // end conditional for the event being stored in the mixer
        
        if(storeKaonPionEvent){//begin flag to store D0  
        
            
        
            eventBufferD0Candidate[VzBin][centralityBin]->addEvent(picoDst);
            
            for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx){//begin loop to store D0s to buffer
   
                StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
                StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
                StPicoTrack const* pion = picoDst->track(kp->pionIdx());
				
				daughterKaonMom = kaon->gMom(pVtx, bField);
                daughterPionMom = pion->gMom(pVtx, bField);
                kPiCharge = kaon->charge()*pion->charge();           
       
                bandBin = getUSInvMassBandBin(kp->m(), kPiCharge, USSideBandLeftLow, USSideBandLeftHigh, D0InvMassLow, D0InvMassHigh, USSideBandRightLow,USSideBandRightHigh); //0 - SBR, 1 - D0 cand. , 2 - SBL, 3 LS      
				topologicalCutPtBin = getTopologicalCutPtBin(kp->pt());
                
                if(topologicalCutPtBin == -1 ) {continue;}
                
                if(USE_DOUBLE_MIS_PID_PROTECTION && (((fabs(kaon->nSigmaKaon()) <= 2.0) && (fabs(kaon->nSigmaPion()) <= 2.0)) &&
                    ((fabs(pion->nSigmaKaon()) <= 2.0) && (fabs(pion->nSigmaPion()) <= 2.0)))){ continue; }   //TPC nSigma double misPID protection  
                
                if((kaon->nHitsFit() < nHitsFitMin || !kaon->isHFTTrack()) || (pion->nHitsFit() < nHitsFitMin || !pion->isHFTTrack()) ) { continue; } 
				if((fabs(kaon->nSigmaKaon()) > 2.0) || (fabs(pion->nSigmaPion()) > 2.0)){ continue; }   //good TPC track    

                if(USE_TOF && !mHFCuts->isTOFKaon(kaon, mHFCuts->getTofBeta(kaon))) {continue;}
                if(USE_TOF && !mHFCuts->isTOFPion(pion, mHFCuts->getTofBeta(pion))) {continue;}
                
                if(USE_TOF_HYBRID && (mHFCuts->getTofBeta(kaon) > 0 && !mHFCuts->isTOFKaon(kaon, mHFCuts->getTofBeta(kaon)))) {continue;}
                if(USE_TOF_HYBRID && (mHFCuts->getTofBeta(pion) > 0 && !mHFCuts->isTOFPion(pion, mHFCuts->getTofBeta(pion)))) {continue;}
                
				//if(USE_PAIR_WISE_PT_CUT && !pairWisePtCutCheck(topologicalCutPtBin, daughterKaonMom.perp(), daughterPionMom.perp(), 1)) { continue; }
				
				
                if(kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh  &&                                 //begin conditional to store D0 candidate to buffer
                    cutCheck(kp, d0PtLow,  d0PtHigh,  d0TopologicalCutArray[topologicalCutPtBin][0],  d0DecayLengthMax,  d0TopologicalCutArray[topologicalCutPtBin][1],
                     d0DaughterPionPtMin, d0DaughterKaonPtMin, d0TopologicalCutArray[topologicalCutPtBin][4], d0TopologicalCutArray[topologicalCutPtBin][3],
                     d0TopologicalCutArray[topologicalCutPtBin][2])) {          
                    
                    kaonPionMom.set(kp->lorentzVector().px(),kp->lorentzVector().py(),kp->lorentzVector().pz());
                   
                    eventBufferD0Candidate[VzBin][centralityBin]->addKaonPionToEvent(eventBufferD0Candidate[VzBin][centralityBin]->getBufferIndex()-1, kaonPionMom, kp->m(), 
                                                                                     kp->kaonIdx(), kp->pionIdx(), kaon->charge()*pion->charge(), kaon->charge(), pion->charge(),
                                                                                     kaon->gMom(pVtx, bField), pion->gMom(pVtx, bField));
                    
               }//end conditional to store D0 candidate to buffer
               
               if(kp->m() > USSideBandLeftLow && kp->m() < USSideBandLeftHigh && kaon->charge()*pion->charge() < 0 &&  //begin conditional to store US left sideband
                    cutCheck(kp, d0PtLow,  d0PtHigh,  d0TopologicalCutArray[topologicalCutPtBin][0],  d0DecayLengthMax,  d0TopologicalCutArray[topologicalCutPtBin][1],
                     d0DaughterPionPtMin, d0DaughterKaonPtMin, d0TopologicalCutArray[topologicalCutPtBin][4], d0TopologicalCutArray[topologicalCutPtBin][3],
                     d0TopologicalCutArray[topologicalCutPtBin][2])) {          
                    
                    kaonPionMom.set(kp->lorentzVector().px(),kp->lorentzVector().py(),kp->lorentzVector().pz());
                    
                    eventBufferD0Candidate[VzBin][centralityBin]->addKaonPionToEvent(eventBufferD0Candidate[VzBin][centralityBin]->getBufferIndex()-1, kaonPionMom, kp->m(), 
                                                                                     kp->kaonIdx(), kp->pionIdx(), kaon->charge()*pion->charge(), kaon->charge(), pion->charge(),
                                                                                     kaon->gMom(pVtx, bField), pion->gMom(pVtx, bField));
                    
               }//end conditional to store US left sideband
               
               if(kp->m() > USSideBandRightLow && kp->m() < USSideBandRightHigh && kaon->charge()*pion->charge() < 0 &&  //begin conditional to store US right sideband
                    cutCheck(kp, d0PtLow,  d0PtHigh,  d0TopologicalCutArray[topologicalCutPtBin][0],  d0DecayLengthMax,  d0TopologicalCutArray[topologicalCutPtBin][1],
                     d0DaughterPionPtMin, d0DaughterKaonPtMin, d0TopologicalCutArray[topologicalCutPtBin][4], d0TopologicalCutArray[topologicalCutPtBin][3],
                     d0TopologicalCutArray[topologicalCutPtBin][2])) {          
                    
                    kaonPionMom.set(kp->lorentzVector().px(),kp->lorentzVector().py(),kp->lorentzVector().pz());
                    
                    eventBufferD0Candidate[VzBin][centralityBin]->addKaonPionToEvent(eventBufferD0Candidate[VzBin][centralityBin]->getBufferIndex()-1, kaonPionMom, kp->m(), 
                                                                                     kp->kaonIdx(), kp->pionIdx(), kaon->charge()*pion->charge(), kaon->charge(), pion->charge(),
                                                                                     kaon->gMom(pVtx, bField), pion->gMom(pVtx, bField));
                    
               }//end conditional to store US right sideband
           }//end loop to store D0s to buffer
        }//end flag to store D0   
     }//end event mixing on/off switch    
/****************************************************************************************************************/  
/***************************************LOOP TO STORE EVENTS IN MIXER END***************************************/
/****************************************************************************************************************/   
   
/****************************************************************************************************************/   
/*****************************************BEGIN MAIN SIBLING LOOP************************************************/
/****************************************************************************************************************/
    if(DEBUG) { cout << "   Begin sibling pair code block" << endl << endl; }
   
    for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx){ // begin main sibling loop
     
        StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
    
        StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
        StPicoTrack const* pion = picoDst->track(kp->pionIdx());
		
	    daughterKaonMom = kaon->gMom(pVtx, bField);
        daughterPionMom = pion->gMom(pVtx, bField);
        kPiCharge = kaon->charge()*pion->charge();           
       
        bandBin = getUSInvMassBandBin(kp->m(), kPiCharge, USSideBandLeftLow, USSideBandLeftHigh, D0InvMassLow, D0InvMassHigh, USSideBandRightLow,USSideBandRightHigh); //0 - SBR, 1 - D0 cand. , 2 - SBL, 3 LS            
        topologicalCutPtBin = getTopologicalCutPtBin(kp->pt());
	    ptBin     = getPtBin(kp->pt());
        finePtBin = getFinePtBin(kp->pt());
      
        if(topologicalCutPtBin == -1 ) {continue;}
      
        if(USE_DOUBLE_MIS_PID_PROTECTION && (((fabs(kaon->nSigmaKaon()) <= 2.0) && (fabs(kaon->nSigmaPion()) <= 2.0)) &&
           ((fabs(pion->nSigmaKaon()) <= 2.0) && (fabs(pion->nSigmaPion()) <= 2.0)))){ continue; }   //TPC nSigma double misPID protection  
      
        if((kaon->nHitsFit() < nHitsFitMin || !kaon->isHFTTrack()) || (pion->nHitsFit() < nHitsFitMin || !pion->isHFTTrack()) ) { continue; } 
        if((fabs(kaon->nSigmaKaon()) > 2.0) || (fabs(pion->nSigmaPion()) > 2.0)){ continue; }   //good TPC track   
      
        if(USE_TOF && !mHFCuts->isTOFKaon(kaon, mHFCuts->getTofBeta(kaon))) {continue;}
        if(USE_TOF && !mHFCuts->isTOFPion(pion, mHFCuts->getTofBeta(pion))) {continue;}
        
        if(USE_TOF_HYBRID && (mHFCuts->getTofBeta(kaon) > 0 && !mHFCuts->isTOFKaon(kaon, mHFCuts->getTofBeta(kaon)))) {continue;}
        if(USE_TOF_HYBRID && (mHFCuts->getTofBeta(pion) > 0 && !mHFCuts->isTOFPion(pion, mHFCuts->getTofBeta(pion)))) {continue;}
       
	    //if(USE_PAIR_WISE_PT_CUT && !pairWisePtCutCheck(topologicalCutPtBin, daughterKaonMom.perp(), daughterPionMom.perp(), 1)) { continue; }
	   
	   
        if(!cutCheck(kp, d0PtLow,  d0PtHigh,  d0TopologicalCutArray[topologicalCutPtBin][0],  d0DecayLengthMax,  d0TopologicalCutArray[topologicalCutPtBin][1],
                     d0DaughterPionPtMin, d0DaughterKaonPtMin, d0TopologicalCutArray[topologicalCutPtBin][4], d0TopologicalCutArray[topologicalCutPtBin][3],
                     d0TopologicalCutArray[topologicalCutPtBin][2])) { continue; }    
      
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
        
       
     
    /****************************************************************************************************************/
    /****************************SIBLING EVENT PAIRS FORMATIION BEGINS HERE*****************************************/  
    /***************************************************************************************************************/     
        
        
        
        if(kPiCharge < 0){// begin US conditional
                
				//bandBin = getUSInvMassBandBin(kp->m(), kPiCharge, USSideBandLeftLow, USSideBandLeftHigh, D0InvMassLow, D0InvMassHigh, USSideBandRightLow,USSideBandRightHigh); //0 - SBR, 1 - D0 cand. , 2 - SBL, 3 LS
				
                if(DEBUG && bandBin == 1) { cout << "   D0 Candidate Mass: " << kp->m() << "    PtBin : " << ptBin << endl; }
            
                if(!storedD0Event && bandBin == 1){
                
                    eventCounter->Fill(0);
                    storedD0Event = true;
                }
   
                if(bandBin==1) {d0Counter = d0Counter + 1;}
				D0ptDist[bandBin][NUM_PT_BINS]->Fill(kp->pt());
                D0ptDist[bandBin][ptBin]->Fill(kp->pt()); 
                D0EtaDist[bandBin][ptBin]->Fill(kp->eta());
                D0PhiDist[bandBin][ptBin]->Fill(kp->phi());
                D0KaonPtDist[bandBin][ptBin]->Fill(daughterKaonMom.perp());
                D0KaonEtaDist[bandBin][ptBin]->Fill(daughterKaonMom.pseudoRapidity()); 
                D0KaonPhiDist[bandBin][ptBin]->Fill(daughterKaonMom.phi()); 
                D0PionPtDist[bandBin][ptBin]->Fill(daughterPionMom.perp());
                D0PionEtaDist[bandBin][ptBin]->Fill(daughterPionMom.pseudoRapidity()); 
                D0PionPhiDist[bandBin][ptBin]->Fill(daughterPionMom.phi()); 
                D0DecayLengthDist[bandBin][ptBin]->Fill(kp->decayLength()); 
                D0DCAToPVDist[bandBin][ptBin]->Fill(kp->perpDcaToVtx()); 
                D0KaonPVDist[bandBin][ptBin]->Fill(kp->kaonDca());
                D0KaonPVDist[bandBin][NUM_PT_BINS]->Fill(kp->kaonDca());                 
                D0PionPVDist[bandBin][ptBin]->Fill(kp->pionDca());
                D0PionPVDist[bandBin][NUM_PT_BINS]->Fill(kp->pionDca());                 
                D0DaughterDCADist[bandBin][ptBin]->Fill(kp->dcaDaughters()); 
                D0EtaDist[bandBin][NUM_PT_BINS]->Fill(kp->eta());
                D0PhiDist[bandBin][NUM_PT_BINS]->Fill(kp->phi());
                D0KaonPtDist[bandBin][NUM_PT_BINS]->Fill(daughterKaonMom.perp());
                D0KaonEtaDist[bandBin][NUM_PT_BINS]->Fill(daughterKaonMom.pseudoRapidity());
                D0KaonPhiDist[bandBin][NUM_PT_BINS]->Fill(daughterKaonMom.phi());
                D0PionPtDist[bandBin][NUM_PT_BINS]->Fill(daughterPionMom.perp());
                D0PionEtaDist[bandBin][NUM_PT_BINS]->Fill(daughterPionMom.pseudoRapidity());
                D0PionPhiDist[bandBin][NUM_PT_BINS]->Fill(daughterPionMom.phi());
                D0DecayLengthDist[bandBin][NUM_PT_BINS]->Fill(kp->decayLength());
                D0DCAToPVDist[bandBin][NUM_PT_BINS]->Fill(kp->perpDcaToVtx());
                D0DaughterDCADist[bandBin][NUM_PT_BINS]->Fill(kp->dcaDaughters());
                D0PtKaonVsPtPion[bandBin][ptBin]->Fill(daughterPionMom.perp(), daughterKaonMom.perp());
                D0PKaonVsPPion[bandBin][ptBin]->Fill(daughterPionMom.mag(), daughterKaonMom.mag());
                
                delEtaDaughters = TMath::Abs(daughterKaonMom.pseudoRapidity()-daughterPionMom.pseudoRapidity());
                delPhiDaughters = daughterKaonMom.phi()-daughterPionMom.phi();
               
               // cout << "HERE AGAIN AGAIN AGAIN!" << endl;
 
                if(delPhiDaughters < -TMath::Pi()) { delPhiDaughters = delPhiDaughters + 2*TMath::Pi(); }     //shift [-2Pi, 2Pi] -> [-Pi,Pi]
                else if(delPhiDaughters > TMath::Pi()) { delPhiDaughters = delPhiDaughters - 2*TMath::Pi(); }
                
                delPhiDaughters = TMath::Abs(delPhiDaughters);
               
                if(bandBin != 3) { 
						
						D0RawEtaVsRawPhi[bandBin][ptBin]->Fill(delEtaDaughters, delPhiDaughters); 
						D0RawEtaVsRawPhi[bandBin][NUM_PT_BINS]->Fill(delEtaDaughters, delPhiDaughters);
						
                        if((daughterKaonMom.perp() > 1.2 && daughterKaonMom.perp() < 2.0) && (daughterPionMom.perp() > .35 && daughterPionMom.perp() < 0.6)){
                        
                            D0RawEtaVsRawPhiLeftPtBlob[bandBin][ptBin]->Fill(delEtaDaughters, delPhiDaughters);
                        }
                        
                        if((daughterKaonMom.perp() > .5 && daughterKaonMom.perp() < 1.1) && (daughterPionMom.perp() > 1.0 && daughterPionMom.perp() < 1.3)){
                        
                            D0RawEtaVsRawPhiRightPtBlob[bandBin][ptBin]->Fill(delEtaDaughters, delPhiDaughters);
                        }
					}
               
            }
            
            //bandBin = getUSInvMassBandBin(kp->m(), kPiCharge, USSideBandLeftLow, USSideBandLeftHigh, D0InvMassLow, D0InvMassHigh, USSideBandRightLow,USSideBandRightHigh);	
            
            if(bandBin > -1){// begin sibling pair formation with whatever "band" we are in, if we have a candidate in the band -- NEED TO MAKE A TRACK LIST TO SPEED THIS UP!!!
                
                realTracks = 0;
                
                vZandCentBinPerEvent[bandBin][ptBin]->Fill(centralityBin, VzBin);
                
                if(DEBUG) { cout << endl << "Sibling pair formation with bandBin: " << bandBin << endl; }
                
                kpPhi = kp->phi(); 
                
                if(SINGLES_DISTS){
                        //if(phi < (-TMath::Pi()+phiBinShift)) { phi = phi + 2*TMath::Pi(); }
                        if(kpPhi < (-TMath::Pi()+phiBinShift)) { kpPhi = kpPhi + 2*TMath::Pi(); }
                        //phiD0vsPhiH[bandBin][ptBin][VzBin][centralityBin]->Fill(phi, kpPhi);
                        //etaD0vsEtaH[bandBin][ptBin][VzBin][centralityBin]->Fill(eta, kp->eta());     
                        phiD0vsEtaD0[bandBin][ptBin][VzBin][centralityBin]->Fill(kp->eta(), kpPhi);    
                        //phiHvsEtaH[bandBin][ptBin][VzBin][centralityBin]->Fill(eta,phi);
                    
                    }
                
                for(unsigned int i = 0; i < mAssociatedHadronList.size(); i++){ // begin picoDST loop for d0-hadron correlations
           
                    trackMom = mAssociatedHadronList[i];
                  
                    daughterKaonMom = kaon->gMom(pVtx, bField);
                    daughterPionMom = pion->gMom(pVtx, bField);
                    
                    pt = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y()));
                    
                    realTracks++;
                    
                    phi  = TMath::ATan2(trackMom.y(),trackMom.x());
                    eta = trackMom.pseudoRapidity();
                    
					
				    if(SINGLES_DISTS){
                        
                        if(phi < (-TMath::Pi()+phiBinShift)) { phi = phi + 2*TMath::Pi(); }
                        phiHvsEtaH[bandBin][ptBin][VzBin][centralityBin]->Fill(eta,phi);
                    
                    }
                    
                    delPhi = kp->phi()-phi;
                    delEta = TMath::Abs(kp->eta()-eta);
                    
                    if(delEta < .1 && USE_TWO_TRACK_INEFF && !(passTwoTrackCuts(daughterKaonMom, kaon->charge(), trackMom, mAssociatedHadronListCharge[i], bField) || 
                                                               passTwoTrackCuts(daughterKaonMom, pion->charge(), trackMom, mAssociatedHadronListCharge[i], bField))  ){ continue; }
                    
                    if(delPhi < -TMath::Pi()) { delPhi = delPhi + 2*TMath::Pi(); }     //shift [-2Pi, 2Pi] -> [-Pi,Pi]
                    else if(delPhi > TMath::Pi()) { delPhi = delPhi - 2*TMath::Pi(); }
                    
                    delPhi = TMath::Abs(delPhi);//Gets absolute value of angle 
                    
					if(delPhi >= (3*TMath::PiOver2()) + phiBinShift){ delPhi = delPhi - 2*TMath::Pi(); }//....and shifts it
                    //if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }//....and shifts it
                    
                    delPhiCp = -delPhi;//Gets negative copy of absolute value of angle.....
                    if(delPhiCp < -TMath::PiOver2() + phiBinShift){ delPhiCp = delPhiCp + 2*TMath::Pi(); }//....and shifts it
                    //if(delPhiCp < -TMath::PiOver2()){ delPhiCp = delPhiCp + 2*TMath::Pi(); }//....and shifts it
                
                    //efficiencyCorrInfoSib[ptBin][centralityClass]->Fill(2.5); // increments the pair count
                
                    if(D0_HADRON_CORR){
                    
                        hadronEffWeight = effWeightPions[centralityClass]->Eval(pt, 0, 0, 0);                      //E hadron
                        pionEffWeight   = effWeightPions[centralityClass]->Eval(daughterPionMom.perp(), 0, 0, 0);  //E Pion
                        kaonEffWeight   = effWeightKaons[centralityClass]->Eval(daughterKaonMom.perp(), 0, 0, 0);  //E Kaon
                        D0Eff = effWeightD0[centralityClass]->Eval(kp->pt(), 0, 0, 0);            //E K-Pi pair
                        
                        HFTRatioKaon   = HFTRatioFunc[centralityClass]->Eval(daughterKaonMom.perp(), 0, 0, 0); 
                        HFTRatioPion   = HFTRatioFunc[centralityClass]->Eval(daughterPionMom.perp(), 0, 0, 0); 
                        HFTRatioHadron = HFTRatioFunc[centralityClass]->Eval(pt, 0, 0, 0);
                        
                        
                        tripleProduct = HFTRatioHadron*hadronEffWeight*HFTRatioPion*pionEffWeight*HFTRatioKaon*kaonEffWeight;
                        doubleProduct = D0Eff*HFTRatioHadron*hadronEffWeight;
                        
                        if(tripleProduct <= 0 || doubleProduct <= 0) { continue; }
                        
                        //efficiencyCorrInfoSib[ptBin][centralityClass]->Fill(1.5, D0HadronEffWeight);
                        
                        //oneOverD0Eff = 1/D0Eff;
                        //oneOverEffWeight = 1/hadronEffWeight;
                        
                        if(bandBin == 1) { 
                        
                            finalPairWeight = calculatePairWeightSignalRegion(BOverBPlusS[ptBin][centralityClass], SOverBPlusS[ptBin][centralityClass], 
                                                                              tripleProduct, doubleProduct, siblingTripleEff[ptBin][centralityClass], 
                                                                              siblingDoubleEff[ptBin][centralityClass]);
                        }                                                      
                        
                        else finalPairWeight = siblingTripleEff[ptBin][centralityClass]/tripleProduct;
                    
                        if(!USE_PAIR_WEIGHT) { finalPairWeight = 1; }
                    
                        sibCorrBin[bandBin][VzBin][centralityBin]->Fill(delEta, delPhi, finalPairWeight);
                        sibCorrBin[bandBin][VzBin][centralityBin]->Fill(-delEta, delPhiCp, finalPairWeight);
                        sibCorrBin[bandBin][VzBin][centralityBin]->Fill(-delEta, delPhi, finalPairWeight);
                        sibCorrBin[bandBin][VzBin][centralityBin]->Fill(delEta, delPhiCp, finalPairWeight);
                            if(ptBin > -1) { 
                    
                                sibCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(delEta, delPhi, finalPairWeight);
                                sibCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(-delEta, delPhiCp, finalPairWeight); 
                                sibCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(-delEta, delPhi, finalPairWeight);
                                sibCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(delEta, delPhiCp, finalPairWeight); 
                            }
                    }
                
                }// end picoDST loop for d0-hadron correlations
                
                if(DEBUG) { cout << endl << "   Number of pairs with this candidate: " << realTracks << endl; }
        
                
           }// end sibling pair formation with whatever "band" we are in, if we have a candidate in the band
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
    
    if(EVENT_MIXING){// begin switch to turn on event mixing
        if((eventBufferPicoEvent[VzBin][centralityBin]->getBufferMaxSize() == eventBufferPicoEvent[VzBin][centralityBin]->getBufferSize()) //begin buffer-full mixing conditional
            && eventBufferD0Candidate[VzBin][centralityBin]->getBufferSize() >= 1){ 
    
            if(DEBUG_MIX_BUFFER){cout << "Mixing events in Vz/Centrality Bin " << VzBin<< "/" << centralityBin << endl << endl;}
            
            StMixerEvent* eventWithD0;
            StMixerEvent* hadronEvent;
            StMixerTrack kaonPionTrack1;
            StMixerTrack hadronTrack2;
            int nD0Tracks = 0;
            int nHadronTracks = 0;
            // D0Mom = 0;
            double D0Pt = 0;
            double hadronBField;
        
            if(DEBUG_MIX_BUFFER) { cout << "Code gets here" << endl; }
        
            for(int k = 0; k <  eventBufferD0Candidate[VzBin][centralityBin]->getBufferSize(); k++){ // begin d0 event buffer loop
            
                eventWithD0 = eventBufferD0Candidate[VzBin][centralityBin]->getEvent(k); //This event in the buffer will be the trigger D0 event
                nD0Tracks = eventWithD0->getKaonPionListSize();
        
                for(int i = 0; i < eventBufferPicoEvent[VzBin][centralityBin]->getBufferSize(); i++){ //begin associated event buffer loop
            
                    hadronEvent = eventBufferPicoEvent[VzBin][centralityBin]->getEvent(i); 
                    hadronBField = hadronEvent->field();
                    nHadronTracks = hadronEvent->getNoTracks();
                    if(DEBUG_MIX_BUFFER) { cout << "Number of hadron track: " << nHadronTracks << endl; }   
                    if(DEBUG_MIX_BUFFER) { cout << "buffer position for hadron event: " << i << endl; }
            
                    for (int idx = 0; idx < nD0Tracks; idx++){//begin loop over event D0 kaon-pion list
                
                        kaonPionTrack1 = eventWithD0->getKaonPionAt(idx);
                        if(DEBUG_MIX_BUFFER) { cout << "index position for D0 track: " << idx << endl; }
                        kaonPionMom = kaonPionTrack1.gMom();
                        D0Pt = TMath::Sqrt((kaonPionMom.x()*kaonPionMom.x())+(kaonPionMom.y()*kaonPionMom.y()));
                    
                        daughterPionPt = TMath::Sqrt((kaonPionTrack1.pionGMom().x()*kaonPionTrack1.pionGMom().x())+(kaonPionTrack1.pionGMom().y()*kaonPionTrack1.pionGMom().y()));
                        daughterKaonPt = TMath::Sqrt((kaonPionTrack1.kaonGMom().x()*kaonPionTrack1.kaonGMom().x())+(kaonPionTrack1.kaonGMom().y()*kaonPionTrack1.kaonGMom().y()));
                    
                        ptBin = getPtBin(D0Pt);
                        
                       
                    
                        if(kaonPionTrack1.charge() < 0 && kaonPionTrack1.mass() > D0InvMassLow && kaonPionTrack1.mass() < D0InvMassHigh)                 { bandBin = 1; }    //US center band
                        else if(kaonPionTrack1.charge() < 0 && kaonPionTrack1.mass() > USSideBandLeftLow &&  kaonPionTrack1.mass() < USSideBandLeftHigh) { bandBin = 0; }    //US left band
                        else if(kaonPionTrack1.charge() < 0 && kaonPionTrack1.mass() > USSideBandRightLow && kaonPionTrack1.mass() < USSideBandRightHigh){ bandBin = 2; }    //US right band
                        else if(kaonPionTrack1.charge() > 0 && kaonPionTrack1.mass() > D0InvMassLow && kaonPionTrack1.mass() < D0InvMassHigh)            { bandBin = 3; }    //LS conditional
                    
                            if(DEBUG_MIX_BUFFER) { cout << "Mixing with mother in bin: " << bandBin << endl; }
                    
                            for(int j = 0; j < nHadronTracks; j++){//begin loop over hadron event tracks
                        
                                hadronTrack2 = hadronEvent->getTrack(j);  
                                             
                                trackMom = hadronTrack2.gMom();
                                
                                
                                pt = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y())); 
                                phi = TMath::ATan2(trackMom.y(),trackMom.x());  
                                eta = trackMom.pseudoRapidity();
                 
                                /* if(SINGLES_DISTS){
                                   
                                    phiD0vsPhiHMixed[bandBin][ptBin][VzBin][centralityBin]->Fill(phi, kaonPionTrack1.gMom().phi());
                                    etaD0vsEtaHMixed[bandBin][ptBin][VzBin][centralityBin]->Fill(eta, kaonPionTrack1.gMom().pseudoRapidity());     
                                    phiD0vsEtaD0Mixed[bandBin][ptBin][VzBin][centralityBin]->Fill(kaonPionTrack1.gMom().pseudoRapidity(), kaonPionTrack1.gMom().phi());    
                                    phiHvsEtaHMixed[bandBin][ptBin][VzBin][centralityBin]->Fill(eta,phi);
                                    
                                } */
                 
                                delPhi = kaonPionTrack1.gMom().phi()-phi;
                                if(delPhi <= -TMath::Pi()) { delPhi = delPhi + 2*TMath::Pi(); }     //shift [-2Pi, 2Pi] -> [-Pi,Pi]
                                else if(delPhi >= TMath::Pi()) { delPhi = delPhi - 2*TMath::Pi(); }
                 
                                delPhi = TMath::Abs(delPhi);//Gets absolute value of angle 
                                if(delPhi >= (3*TMath::PiOver2()) + phiBinShift){ delPhi = delPhi - 2*TMath::Pi(); }//....and shifts it
                                //if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }//....and shifts it
                    
                                delPhiCp = -delPhi;//Gets negative copy of absolute value of angle.....
                                if(delPhiCp < -TMath::PiOver2() + phiBinShift){ delPhiCp = delPhiCp + 2*TMath::Pi(); }//....and shifts it
                                //if(delPhiCp < -TMath::PiOver2()){ delPhiCp = delPhiCp + 2*TMath::Pi(); }//....and shifts it
                
                
                                delEta = TMath::Abs(kaonPionTrack1.gMom().pseudoRapidity()-eta);
                                
                                if(delEta < .1 && USE_TWO_TRACK_INEFF && !(passTwoTrackCuts(daughterKaonMom, kaonPionTrack1.kaonCharge(), trackMom, hadronTrack2.charge(), hadronBField) || 
                                                            passTwoTrackCuts(daughterKaonMom, kaonPionTrack1.pionCharge(), trackMom, hadronTrack2.charge(), hadronBField)) ){ continue; }
               
                                hadronEffWeight = effWeightPions[centralityClass]->Eval(pt, 0, 0, 0);                      //E hadron
                                pionEffWeight   = effWeightPions[centralityClass]->Eval(daughterPionMom.perp(), 0, 0, 0);  //E Pion
                                kaonEffWeight   = effWeightKaons[centralityClass]->Eval(daughterKaonMom.perp(), 0, 0, 0);  //E Kaon
                                D0Eff = effWeightD0[centralityClass]->Eval(D0Pt, 0, 0, 0);            //E K-Pi pair
                                
                                HFTRatioKaon   = HFTRatioFunc[centralityClass]->Eval(daughterKaonMom.perp(), 0, 0, 0); 
                                HFTRatioPion   = HFTRatioFunc[centralityClass]->Eval(daughterPionMom.perp(), 0, 0, 0); 
                                HFTRatioHadron = HFTRatioFunc[centralityClass]->Eval(pt, 0, 0, 0);
                                
                                
                                tripleProduct = HFTRatioHadron*hadronEffWeight*HFTRatioPion*pionEffWeight*HFTRatioKaon*kaonEffWeight;
                                doubleProduct = D0Eff*HFTRatioHadron*hadronEffWeight;
                                
                                if(tripleProduct <= 0 || doubleProduct <= 0) { continue; }
                                
                                if(bandBin == 1) { 
                                
                                    finalPairWeight = calculatePairWeightSignalRegion(BOverBPlusS[ptBin][centralityClass], SOverBPlusS[ptBin][centralityClass], 
                                                                                      tripleProduct, doubleProduct, mixedTripleEff[ptBin][centralityClass], 
                                                                                      mixedDoubleEff[ptBin][centralityClass]);
                                }                                                      
                                
                                else finalPairWeight = siblingTripleEff[ptBin][centralityClass]/tripleProduct;
                                
                                if(!USE_PAIR_WEIGHT) { finalPairWeight = 1; }
                               
                                if(D0_HADRON_CORR){
                                    mixCorrBin[bandBin][VzBin][centralityBin]->Fill(delEta, delPhi, finalPairWeight); //delEta & delPhi stored in US mixed correlation histogram
                                    mixCorrBin[bandBin][VzBin][centralityBin]->Fill(-delEta, delPhiCp, finalPairWeight); //delEta & delPhi stored in US mixed correlation histogram
                                    mixCorrBin[bandBin][VzBin][centralityBin]->Fill(-delEta, delPhi, finalPairWeight); //delEta & delPhi stored in US mixed correlation histogram
                                    mixCorrBin[bandBin][VzBin][centralityBin]->Fill(delEta, delPhiCp, finalPairWeight); //delEta & delPhi stored in US mixed correlation histogram
                                    if(ptBin > -1) { 
                                    
                                        mixCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(delEta, delPhi, finalPairWeight);
                                        mixCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(-delEta, delPhiCp, finalPairWeight);
                                        mixCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(-delEta, delPhi, finalPairWeight);
                                        mixCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(delEta, delPhiCp, finalPairWeight);
                                    }
                                }
                            }// end loop over hadron event tracks
                        }// end loop over event D0 kaon-pion list
                    }// end associated event buffer loop
                }// end d0 event buffer loop
            
            eventBufferPicoEvent[VzBin][centralityBin]->removeFirstEvent(); //MAY NEED TO FIX THIS ALGORITHM TO AVOID A BIAS
            
            for(int numEvents = 0; numEvents < eventBufferD0Candidate[VzBin][centralityBin]->getBufferSize(); numEvents++){ //This will remove all of the events in the D0 buffer, since they are now used.
            
                eventBufferD0Candidate[VzBin][centralityBin]->removeFirstEvent(); 
            
            }
      
            if(DEBUG) { 
        
                cout << endl << endl;
                cout << "First Event in US buffer for Vz/centrality bin "<< VzBin << "/"<< centralityBin<< " cleared" << endl;
                cout << endl;
            }     
        
        }// end buffer-full US mixing conditional        
    }// end switch to turn on event mixer  
    
/***************************************************************************************************************/
/*************************************EVENT MIXING ENDS HERE****************************************************/
/***************************************************************************************************************/   
       
       
       
       d0CountPerEvent->Fill(d0Counter);
       //invMassMinusBG->Add(invMass, likeSignBG, 1, -1);
       //D0PeakMinusBG->Add(D0PeakPlusBG, D0LikeSignBG, 1, -1);
    
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


bool StPicoD0AnaMaker::pairWisePtCutCheck(int ptBin, double kpt, double pipt, int band){

    if(ptBin == 3 || ptBin == 4) { return true; }

	                           //0-1,  1-2,  2-3, 3-5, 5-10 GeV/c pt-Bins
			 
	double ptKaonCutsLow[5]  = {.56,  .48,  .75,  0.0, 0.0};
	double ptPionCutsLow[5]  = {1.08, .92,  1.47, 0.0, 0.0};
	double ptKaonCutsHigh[5] = {.67,  1.14, 1.35, 0.0, 0.0};
    double ptPionCutsHigh[5] = {1.24, 1.44, 2.13, 0.0, 0.0};
	
	bool passCuts = (kpt < ptKaonCutsLow[ptBin]  || kpt > ptKaonCutsHigh[ptBin]) &&
	                (pipt < ptPionCutsLow[ptBin] || pipt > ptPionCutsHigh[ptBin]);
					
	return passCuts;
	
}
	
	
int StPicoD0AnaMaker::getCentralityBin(int nTracks){

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

int StPicoD0AnaMaker::getPtBin(double pt){

    if(pt >=  0.0  && pt <  1.0)    { return 0;  }   
    if(pt >=  1.0  && pt <  2.0)    { return 1;  }
    if(pt >=  2.0  && pt <  10.0)   { return 2;  }        
    //if(pt >=  3.0  && pt <  5.0)   { return 3;  }
    //if(pt >=  5.0  && pt <  10.0)   { return 4;  }
    
    
    else return -1;

}    

int StPicoD0AnaMaker::getFinePtBin(double pt){

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

int StPicoD0AnaMaker::getTopologicalCutPtBin(double pt){

    if(pt >=  0.0  && pt <  1.0)    { return 0;  }   
    if(pt >=  1.0  && pt <  2.0)    { return 1;  }
    if(pt >=  2.0  && pt <  10.0)    { return 2;  }
    //if(pt >=  3.0  && pt <  5.0)    { return 3;  }
    //if(pt >=  5.0  && pt <  10.0)   { return 4;  }
    
    else return -1;

}    

int StPicoD0AnaMaker::getCentralityClass(double nTracks){

    if(nTracks >= 1   && nTracks < 131)  { return 0;  }  //peripheral
    if(nTracks >= 131  && nTracks < 440) { return 1;  }  //midcentral
    if(nTracks >= 440)                   { return 2;  }  //central
    
    else return -1;

}    


int StPicoD0AnaMaker::getUSInvMassBandBin(double mass, int charge, double USSideBandLeftLow, double USSideBandLeftHigh, double D0InvMassLow, double D0InvMassHigh,
																	double USSideBandRightLow, double USSideBandRightHigh) const
{

	if(charge > 0 && (mass > D0InvMassLow && mass < D0InvMassHigh)) { return 3; } //LS
    if(charge < 0 && (mass > D0InvMassLow && mass < D0InvMassHigh)) { return 1; }   //US  
    if(charge < 0 && (mass > USSideBandLeftLow && mass < USSideBandLeftHigh)) { return 0; }
    if(charge < 0 && (mass > USSideBandRightLow && mass < USSideBandRightHigh)) { return 2; }
	
	else return -1;
	
}

bool StPicoD0AnaMaker::checkDCAtoPV(float trackDCA){

     return (trackDCA <= trackDCAtoPvtx);
     
}

double StPicoD0AnaMaker::calculatePairWeightSignalRegion(double BOverBPlusS, double SOverBPlusS, double tripleProduct, double doubleProduct, double tripleProductAvg, double doubleProductAvg){

    return ((BOverBPlusS*(tripleProductAvg/tripleProduct)) + (SOverBPlusS*(doubleProductAvg/doubleProduct)));
}

bool StPicoD0AnaMaker::passTwoTrackCuts(StThreeVectorF t1, int t1Charge, StThreeVectorF t2, int t2Charge, double B){


    double delPhi = t1.phi()-t2.phi();
    double pt1 = t1.perp();
    double pt2 = t2.perp();

    if(B > 0){ //full-field
    
        if((t1Charge < 0 && t2Charge > 0) && ( (delPhi > 0) || (delPhi < -TMath::PiOver2())) ) { return true; }
        if((t1Charge > 0 && t2Charge < 0) && ( (delPhi < 0) || (delPhi > TMath::PiOver2())) )  { return true; }
        if((t1Charge > 0 && t2Charge > 0) && ( ((delPhi < 0) && (pt1 < pt2)) || ((delPhi > 0) && (pt1 > pt2)) || (delPhi < -TMath::PiOver4()) || (delPhi > TMath::PiOver4()) ) ) {return true; }
        if((t1Charge < 0 && t2Charge < 0) && ( ((delPhi > 0) && (pt1 < pt2)) || ((delPhi < 0) && (pt1 > pt2)) || (delPhi < -TMath::PiOver4()) || (delPhi > TMath::PiOver4()) ) ) {return true; }
    }
    
    if(B < 0){ //full-field reverse
    
        if((t1Charge < 0 && t2Charge > 0) && ( (delPhi < 0) || (delPhi > -TMath::PiOver2())) ) { return true; }
        if((t1Charge > 0 && t2Charge < 0) && ( (delPhi > 0) || (delPhi < TMath::PiOver2())) )  { return true; }
        if((t1Charge > 0 && t2Charge > 0) && ( ((delPhi > 0) && (pt1 > pt2)) || ((delPhi < 0) && (pt1 < pt2)) || (delPhi > -TMath::PiOver4()) || (delPhi < TMath::PiOver4()) ) ) {return true; }
        if((t1Charge < 0 && t2Charge < 0) && ( ((delPhi < 0) && (pt1 > pt2)) || ((delPhi > 0) && (pt1 < pt2)) || (delPhi > -TMath::PiOver4()) || (delPhi < TMath::PiOver4()) ) ) {return true; }
    }
    
    return false;
    
}
