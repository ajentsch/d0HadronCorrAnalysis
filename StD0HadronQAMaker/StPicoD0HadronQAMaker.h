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



class StPicoD0HadronQAMaker : public StMaker
{
  public:
    StPicoD0HadronQAMaker(char const * name, char const * inputFilesList, 
                     char const * outName, StPicoDstMaker* picoDstMaker);
    virtual ~StPicoD0HadronQAMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    int getEntries() const;

    void setHFCuts(StHFCuts* cuts);   

      

  private:
    StPicoD0HadronQAMaker() {}
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

    int NUM_PHI_BINS;
    int NUM_ETA_BINS;
    int NUM_VZ_BINS;
    int NUM_CENT_BINS;
    int NUM_PT_BINS;
    
    
    StPicoDstMaker* mPicoDstMaker;
    StPicoD0Event* mPicoD0Event;
    
   
    
    TString mOutFileName;
    TString mInputFileList;
    TFile* mOutputFile;
    TChain* mChain;
    int mEventCounter;

    StHFCuts* mHFCuts;

    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate
   
    
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
    bool USE_PT_BINS;
    bool USE_FINE_PT_BINS;
                                           // vpdmb-5-p-nobsmd-hlt // vpdmb-5-p-nobsmd-hlt // vpdmb-5-p-nobsmd // vpdmb-5-p-nobsmd // vpdmb-5-p-nobsmd 
    unsigned int trigger[5]; 

      
    
    TH2D* sibCorrBin[4][10][16];
    TH2D* mixCorrBin[4][10][16];
    TH2D* sibCorrBinPt[4][6][10][16];
    TH2D* mixCorrBinPt[4][6][10][16];
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
    TH1D* D0EtaDist;
    TH1D* D0PhiDist;
    TH1I* eventCounter;
    TH1D* D0ptDist;
    TH1D* D0KaonPtDist;
    TH1D* D0KaonEtaDist; 
    TH1D* D0KaonPhiDist; 
    TH1D* D0PionPtDist;
    TH1D* D0PionEtaDist; 
    TH1D* D0PionPhiDist; 
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
    TH2I* vZandCentBinPerEvent;
    TH1D* histOfCuts;
    TH1D* hadronChi2;
    TH1D* pVtxX;
    TH1D* pVtxY;
    TH1D* pVtxZ;
    TH1D* DCAtoPrimaryVertex;
    TH1D* DCAtoPrimaryVertexCut;
    TH2D* phiD0vsPhiH[10][16];     
    TH2D* etaD0vsEtaH[10][16];
    TH2D* phiD0vsEtaD0[10][16];
    TH2D* phiHvsEtaH[10][16];
    
    

    ClassDef(StPicoD0HadronQAMaker, 0)
};

inline int StPicoD0HadronQAMaker::getEntries() const 
{
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoD0HadronQAMaker::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}

inline void StPicoD0HadronQAMaker::setHFCuts(StHFCuts* cuts)   
{ 
  mHFCuts = cuts; 
}

#endif
