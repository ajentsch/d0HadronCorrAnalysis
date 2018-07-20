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

class TString;
class TFile;
class TNtuple;
class StPicoD0Event;
class StKaonPion;
class StPicoDstMaker;
class StHFCuts;


class StPicoD0AnaMaker : public StMaker
{
  public:
    StPicoD0AnaMaker(char const * name, char const * inputFilesList, 
        char const * outName,StPicoDstMaker* picoDstMaker);
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
    
    
    TH1D* ptDist;
    TH1D* invMass;
    TH1D* angCorrPhi;
    TH1D* angCorrEta;
    TH1D* kaonDist;
    TH1D* pionDist;
    TH1D* likeSignBG;
    TH1D* invMassMinusBG;
    TH1D* D0EtaDist;
    TH1D* D0PhiDist;
   
    TH1D* eventCounter;
    TH2D* angCorr2DAll;
    TH1D* kaonDCAfromD0;
    TH1D* pionDCAfromD0;
    TH1D* decayLengthQA;
    TH1D* pointingAngleQA;
    TH1D* daughterDCAQA;
    TH1D* D0ptDist;
    TH1D* kaonPtDist;
    TH1D* pionPtDist;
    TH1D* kaonEtaDist;
    TH1D* pionEtaDist;
    TH1D* kaonPhiDist;
    TH1D* pionPhiDist;
    TH1D* kaonDCAprimary;
    TH1D* pionDCAprimary;
    TH1D* invMassBin1;	
    TH1D* invMassBin2;	
    TH1D* invMassBin3;	
    TH1D* invMassBin4;	
    TH1D* invMassBin5;
    TH1D* angCorrPhiBin1;
    TH1D* angCorrEtaBin1;
    TH2D* angCorr2DBin1;
    TH1D* angCorrPhiBin2;
    TH1D* angCorrEtaBin2;
    TH2D* angCorr2DBin2;
    TH1D* angCorrPhiBin3;
    TH1D* angCorrEtaBin3;
    TH2D* angCorr2DBin3;
    TH1D* angCorrPhiBin4;
    TH1D* angCorrEtaBin4;
    TH2D* angCorr2DBin4;
    TH1D* angCorrPhiBin5;
    TH1D* angCorrEtaBin5;
    TH2D* angCorr2DBin5;

    TH1D* likeSignBin1;	
    TH1D* likeSignBin2;	
    TH1D* likeSignBin3;	
    TH1D* likeSignBin4;	
    TH1D* likeSignBin5;	 
    TH1D* USminusLSBin1;
    TH1D* USminusLSBin2;
    TH1D* USminusLSBin3;
    TH1D* USminusLSBin4;
    TH1D* USminusLSBin5;
    
    TH1D* D0PeakPlusBG;
    TH1D* D0LikeSignBG;
    TH1D* D0PeakMinusBG;
    TH1D* trackCounter;    

    TH1D* hadronPtDist;
    TH1D* hadronPhiDist;
    TH1D* hadronEtaDist;

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
