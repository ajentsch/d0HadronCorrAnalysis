#ifndef StPicoD0AnaMaker_h
#define StPicoD0AnaMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoD0Event
 *  simultaneously and do analysis. 
 *
 *  Please write your analysis in the ::Make() function.
 *
 *  Authors:  **Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  **Code Maintainer
 *
 * **************************************************
 */

#include "TChain.h"
#include "StMaker.h"
#include "TH2.h"
#include "StMixedEventBuffer/StMixedEventBuffer.h"

class TString;
class TFile;
class TNtuple;
class StPicoD0Event;
class StKaonPion;
class StPicoDstMaker;
class StHFCuts;
class StMixedEventBuffer;


class StPicoD0AnaMaker : public StMaker
{
  public:
    StPicoD0AnaMaker(char const * name, char const * inputFilesList, 
                     char const * outName, StPicoDstMaker* picoDstMaker);
    virtual ~StPicoD0AnaMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    int getEntries() const;

    void setHFCuts(StHFCuts* cuts);   

      

  private:
    StPicoD0AnaMaker() {}
    void readNextEvent();

    bool isGoodPair(StKaonPion const*) const;
    bool cutCheck(StKaonPion const*, double, double, double, double, double, double, double, double, double, double) const;
    int  getCentralityBin(int nTracks);
    int  getVzBin(double Vz);
    int  getPtBin(double pt);
    int  getTopologicalCutPtBin(double pt);
    int  getFinePtBin(double pt);
    int  getCentralityClass(double nTracks);
    bool checkDCAtoPV(float trackDCA);
    double calculatePairWeightSignalRegion(double, double, double, double, double, double);
	bool pairWisePtCutCheck(int, double, double , int);
    bool passTwoTrackCuts(StThreeVectorF, int, StThreeVectorF, int, double);
	int getUSInvMassBandBin(double, int, double, double, double, double, double, double) const;
    double calcInvMass(double, double, double, double, double, double, double, double);
    double calcInvMass3Particle(double, double, double, double, double, double, double, double, double, double, double, double);

    int NUM_PHI_BINS;
    int NUM_ETA_BINS;
    int NUM_VZ_BINS;
    int NUM_CENT_BINS;
    int NUM_PT_BINS;
    int NUM_D0_CUT_PT_BINS;
    double phiBinShift;
    double phiBinShiftPrimary;
    
    StPicoDstMaker* mPicoDstMaker;
    StPicoD0Event* mPicoD0Event;
    
    StMixedEventBuffer* eventBufferPicoEvent[10][16];      //To use for the actual PicoEvent
    StMixedEventBuffer* eventBufferD0Candidate[10][16];      //To use for D0 candidate
    //StKaonPion*         d0CandidateBuffer[10][11];  //stores the candidate D0s use for eventual mixing
    
    TString mOutFileName;
    TString mInputFileList;
    TFile* mOutputFile;
    TChain* mChain;
    int mEventCounter;

    StHFCuts* mHFCuts;

    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate
    
    /*double ptRange[6] = {0.0, 1.0, 2.0, 3.0, 5.0, 10.0};
    double decayLengthCuts[5] = {.0145, .0181, .0212, .0247, .0259};
    double daughterDCACuts[5] = {.0084, .0066, .0057, .0050, .0060};
    double dcaKaonPV[5] = {.0103, .0091, .0095, .0079, .0058};
    double dcaPionPV[5] = {.0110, .0111, .0086, .0081, .0062};
    double dcaV0toPV[5] = {.0061, .0049, .0038, .0038, .0040};*/
    
    int eventNumber;
    int trackCount;
    double Vz;
    int centralityBin;
    int VzBin;
    int ptBin;
    int bandBin;
    int finePtBin;
    int centralityClass;
    int topologicalCutPtBin;
    
    int nVzBins;
    int nCentBins;
    int BUFFER_SIZE;
    
             //  ptmin ptmax   decayLenMin&Max   daughterDCA kaon/pion pt kaon/pion DCA  DCA to PV
        //if(!isGoodPair(kp)) continue;             
        //if(!cutCheck(kp, 0.15,  20.0,  .0200,  999999.0,  .0055,  1.2,  1.2,  .008,  .008,  .0065)) { continue; }
    
    double trackDCA;
    double kaonPtCut;
    double pionPtCut;
    double hadronPtMin;
    double hadronPtMax;
    float trackChi2max;
    float trackDCAtoPvtx;
    int   nHitsFitMin;
    float nHitsFitMinOverMax;
    double d0PtLow;
    double d0PtHigh;
    //double d0DecayLengthMin;
    double d0DecayLengthMax;
    //double daughterDCA;
    double d0DaughterPionPtMin;
    double d0DaughterKaonPtMin;
    //double kaonDCA;
    //double pionDCA;
    //double d0DCAtoPV;
    double D0InvMassLow;
    double D0InvMassHigh;
    double USSideBandLeftLow;
    double USSideBandLeftHigh;
    double USSideBandRightLow;
    double USSideBandRightHigh;
    
    double d0TopologicalCutArray[5][5]; //first index is pt bin, second index is cuts -- decayLength, DCADaughters, d0DCAPV, PiDCAPV, KDCAPV
    
    //double d0DecayLengthMin[5];
    //double daughterDCA[5];
    //double kaonDCA[5];
    //double pionDCA[5];
    //double d0DCAtoPV[5];
    
    bool DEBUG;                //important flags for debugging and for switching binning on and off
    bool DEBUG_MIX_BUFFER;
    bool USE_CENT_BINS;
    bool USE_VZ_BINS;
    bool SINGLES_DISTS;
    bool D0_HADRON_CORR;
    bool EVENT_MIXING;
    bool USE_PT_BINS;
    bool USE_FINE_PT_BINS;
	bool USE_TOF;
    bool USE_SYMMETERIZATION;
	bool USE_PAIR_WISE_PT_CUT;
    bool USE_TOF_HYBRID;
    bool USE_DOUBLE_MIS_PID_PROTECTION;
    bool SOFT_MIS_PID_PROTECTION;
    bool USE_TWO_TRACK_INEFF;
    bool USE_PAIR_WEIGHT;
    bool TURN_OFF_LOW_PT;
    
    double siblingTripleEff[5][3];
    double mixedTripleEff[5][3];
    double siblingDoubleEff[5][3];
    double mixedDoubleEff[5][3];
    double BOverBPlusS[5][3];
    double SOverBPlusS[5][3];
                                           // vpdmb-5-p-nobsmd-hlt // vpdmb-5-p-nobsmd-hlt // vpdmb-5-p-nobsmd // vpdmb-5-p-nobsmd // vpdmb-5-p-nobsmd 
                                           

    unsigned int trigger[5]; 

      
    
    //TH2D* sibCorrBin[5][10][16];
    //TH2D* mixCorrBin[5][10][16];
    TH2D* sibCorrBinPt[5][6][10][16];
    TH2D* mixCorrBinPt[5][6][10][16];
    TH1D* etaDistVz[10];
    TH1D* phiDistVz[10];
    TH2D* etaPhiDistVz[10];
    TH1D* D0InvMassPtBin[6][3];
    TH1D* LSInvMassPtBin[6][3];
    TH1D* D0InvMassFinePtBin[11][3];
    TH1D* LSInvMassFinePtBin[11][3];
    TH1D* ptDist;
    TH1D* invMass;
    TH1D* invMassPer;     
    TH1D* invMassMidCent;
    TH1D* invMassCent;
    TH1D* likeSignBGPer;
    TH1D* likeSignBGMidCent;
    TH1D* likeSignBGCent;
    TH1D* likeSignBG;
	TH1I* eventCounter;
    TH1D* D0ptDist[3][6];
	TH1D* D0EtaDist[3][6];
    TH1D* D0PhiDist[3][6];
    TH1D* D0KaonPtDist[3][6];
    TH1D* D0KaonEtaDist[3][6]; 
    TH1D* D0KaonPhiDist[3][6]; 
    TH1D* D0PionPtDist[3][6]; 
    TH1D* D0PionEtaDist[3][6]; 
    TH1D* D0PionPhiDist[3][6];  
    TH1D* D0DecayLengthDist[3][6]; 
    TH1D* D0DCAToPVDist[3][6]; 
    TH1D* D0KaonPVDist[3][6]; 
    TH1D* D0PionPVDist[3][6]; 
    TH1D* D0DaughterDCADist[3][6]; 
    //TH2D* D0PtKaonVsPtPion[3][6];
    //TH2D* D0PKaonVsPPion[3][6];
    //TH2D* D0RawEtaVsRawPhi[3][6];
    //TH2D* D0RawEtaVsRawPhiLeftPtBlob[3][6];
    //TH2D* D0RawEtaVsRawPhiRightPtBlob[3][6]; 
    //TH2D* SBLPtKaonVsPtPion[6];
    //TH2D* SBLPKaonVsPPion[6];
    //TH2D* SBLRawEtaVsRawPhi[6];
    //TH2D* SBRPtKaonVsPtPion[6];
    //TH2D* SBRPKaonVsPPion[6];
    //TH2D* SBRRawEtaVsRawPhi[6];
   
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
   
    TH1D* KPiptDist;
    TH1D* KPiEtaDistInt;
    TH1D* KPiPhiDistInt;
    TH1D* KPiKaonPtDistInt;
    TH1D* KPiKaonEtaDistInt;
    TH1D* KPiKaonPhiDistInt;
    TH1D* KPiPionPtDistInt;
    TH1D* KPiPionEtaDistInt;
    TH1D* KPiPionPhiDistInt;
    TH1D* KPiDecayLengthDistInt;
    TH1D* KPiDCAToPVDistInt;
    TH1D* KPiKaonPVDistInt;
    TH1D* KPiPionPVDistInt;
    TH1D* KPiDaughterDCADistInt;
   
    TH1D* kaonPtDist;
    TH1D* pionPtDist;
    TH1D* kaonEtaDist;
    TH1D* pionEtaDist;
    TH1D* kaonPhiDist;
    TH1D* pionPhiDist;
    TH1I* trackCounter;    
    TH1D* hadronPtDist;
    TH1D* hadronPhiDist;
    TH1D* hadronEtaDist;
    TH2D* dEdxVsPt;
    TH2D* invBetaVsPt;
    TH1I* usedTracks;
    TH1I* d0CountPerEvent;
    TH2D* vZandCentBinPerEvent[4][5];
    TH1D* histOfCuts;
    TH1D* hadronChi2;
    TH1D* hadronNHitsFit;
    TH1D* hadronNHitsFitOverMax;
    TH1D* pVtxX;
    TH1D* pVtxY;
    TH1D* pVtxZ;
    TH1D* DCAtoPrimaryVertex;
    TH1D* DCAtoPrimaryVertexCut;
    TH2D* phiD0vsPhiH[4][5][10][16];     
    TH2D* etaD0vsEtaH[4][5][10][16];
    TH2D* phiD0vsEtaD0[4][5][10][16];
    TH2D* phiHvsEtaH[4][5][10][16];
    TH2D* phiD0vsPhiHMixed[4][5][10][16];     
    TH2D* etaD0vsEtaHMixed[4][5][10][16];
    TH2D* phiD0vsEtaD0Mixed[4][5][10][16];
    TH2D* phiHvsEtaHMixed[4][5][10][16];
    
    TH1D* efficiencyCorrInfoSib[5][3];
    TH1D* efficiencyCorrInfoMix[5][3];
    
    TH1D* TPCOnlyAssociatedHadronPtDist[3]; //for hft ratio
    TH1D* HFTAssociatedHadronPtDist[3];     //for hft ratio
    
    TH1D* dStarHistogramsSib[5][3];    //these are for figuring out D* influence on the correlations
    TH1D* dStarHistogramsMix[5][3];    //these are for figuring out D* influence on the correlations
	
	TH1D* dStarMinusD0HistogramsSib[5][3];    //these are for figuring out D* influence on the correlations
    TH1D* dStarMinusD0HistogramsMix[5][3];    //these are for figuring out D* influence on the correlations
    
    TH2D* dStarAngularDistSib[5][3];
    TH2D* dStarAngularDistMix[5][3];
    
    
    //mathematical functions for efficiency
    
    TF1 *effWeightPions[3];
    TF1 *effWeightKaons[3];
    TF1 *effWeightD0[3];
    TF1 *HFTRatioFunc[3];
    //TF1 *levy;

    ClassDef(StPicoD0AnaMaker, 0)
};

inline int StPicoD0AnaMaker::getEntries() const 
{
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoD0AnaMaker::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}

inline void StPicoD0AnaMaker::setHFCuts(StHFCuts* cuts)   
{ 
  mHFCuts = cuts; 
}

#endif
