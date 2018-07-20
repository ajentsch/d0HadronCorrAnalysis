#include <limits>

#include "StMixerEvent.h"
//#include "TClonesArray.h"

StMixerEvent::StMixerEvent() : mVtx(StThreeVectorF()) 
{

}     

StMixerEvent::StMixerEvent(StMixerEvent* event) : mVtx(event->mVtx), mBField(event->mBField),                       
                                                  mTracks(event->mTracks), mKaonPionTracks(event->mKaonPionTracks){ }                                        //Need to make a fucntions


StMixerEvent::StMixerEvent(StThreeVectorF vtx, float b) : mVtx(StThreeVectorF()), mBField(std::numeric_limits<float>::quiet_NaN()){  
                                                                                                                                           //constructor to save event with B Field
    mVtx = vtx;                                                                                                                       //and the primary vertex -- need to start using
    mBField = b;                                                                                                                       //kfVertex
    //mKaonPionList.push_back(kaonPion);
}

void StMixerEvent::addTrack(StMixerTrack t) 
{
    mTracks.push_back(t);
    if(t.charge() > 0){ mPosiTracks.push_back(t); } //h-h
    else if(t.charge() < 0) { mNegaTracks.push_back(t); } //h-h
}

void StMixerEvent::addKaonPion(StMixerTrack kaonPion){

    mKaonPionTracks.push_back(kaonPion);
    
}


/*StMixerEvent::StMixerEvent(StMixerEvent *t) : mVtx(t->mVtx), mBField(t->mBField),
					      mTracks(t->mTracks),
					      mEventKaons(t->mEventKaons), mEventPions(t->mEventPions)
{
}*/

