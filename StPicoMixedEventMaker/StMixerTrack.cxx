#include <limits>

#include "StMixerTrack.h"
#include "StPicoDstMaker/StPicoTrack.h"

StMixerTrack::StMixerTrack() : mMom(StThreeVectorF()), mPIDFlag(-1), mCharge(0), mKaonCharge(0), mPionCharge(0), mMass(-1), mKaonIdx(-1), mPionIdx(-1), mKaonMom(StThreeVectorF()), mPionMom(StThreeVectorF())
{
}

StMixerTrack::StMixerTrack(StMixerTrack const * trk) : mMom(trk->mMom), mPIDFlag(trk->mPIDFlag), mCharge(trk->mCharge), mKaonCharge(trk->mKaonCharge), mPionCharge(trk->mPionCharge),
                                                       mMass(trk->mMass), mKaonIdx(trk->mKaonIdx), mPionIdx(trk->mPionIdx)
{
} 

StMixerTrack::StMixerTrack(StThreeVectorF trackMomentum, int charge, int PID){

    mMom = trackMomentum;
    mCharge = charge;
    mPIDFlag = PID;

}    

StMixerTrack::StMixerTrack(StThreeVectorF trackMomentum, double mass, int kaonIdx, int pionIdx, int charge, int kaonCharge, int pionCharge, StThreeVectorF kaonMom, StThreeVectorF pionMom){

    mMom = trackMomentum;
    mMass = mass;
    mKaonIdx = kaonIdx;
    mPionIdx = pionIdx;
    mCharge = charge;  //charge = -1 for D0 and +1 for a LS KPi pair
    mKaonCharge = kaonCharge;
    mPionCharge = pionCharge;
    mKaonMom = kaonMom;
    mPionMom = pionMom;
}    