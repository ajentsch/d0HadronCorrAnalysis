#include <limits>

#include "TTree.h"
#include "TH2F.h"

#include "StMixedEventBuffer.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoDst.h"

#include "StPicoCharmContainers/StPicoD0Event.h"
#include "StPicoCharmContainers/StKaonPion.h"
#include "StMixerEvent.h"

ClassImp(StMixedEventBuffer)

StMixedEventBuffer::StMixedEventBuffer()
{
}


StMixedEventBuffer::~StMixedEventBuffer(){
    
    /*for ( int i = 0; i < mEventsList.size(); i++)
    {
        delete mEventsList[i];
    }

    cout << "buffer destructor called" << endl;*/
    
}

//void addEvent(StPicoDst const* picoDst, TClonesArray const * KaonPionList);
void StMixedEventBuffer::addEvent(StPicoDst const* const picoDst){        

    StThreeVectorF pVertex = picoDst->event()->primaryVertex();
    StMixerEvent* event = new StMixerEvent(pVertex, picoDst->event()->bField());
    
    
    if(mBufferCounter < mBufferSize){
        mEventsList.push_back(event);
        mBufferCounter++;
    } 
     
}


void StMixedEventBuffer::addTrackToEvent(int index, StThreeVectorF trackMom, int charge, int PID){

    StMixerEvent* event = mEventsList.at(index);
    StMixerTrack trk = StMixerTrack(trackMom, charge, PID);
    
    event->addTrack(trk);
    
}    
    
    
void StMixedEventBuffer::addKaonPionToEvent(int index, StThreeVectorF kaonPionMom, double mass, int kaonIdx, int pionIdx, int charge, int kaonCharge, int pionCharge, StThreeVectorF kaonMom, StThreeVectorF pionMom){

    StMixerEvent* event = mEventsList.at(index);
    StMixerTrack trk = StMixerTrack(kaonPionMom, mass, kaonIdx, pionIdx, charge, kaonCharge, pionCharge, kaonMom, pionMom);
    
    event->addKaonPion(trk);
    
}        
    
    
void StMixedEventBuffer::clearBuffer() {

    for (unsigned int i = 0; i < mEventsList.size(); i++)
    {
        delete mEventsList[i];
    }
    mEventsList.erase(mEventsList.begin(), mEventsList.end());   
    //cout << "buffered events deleted" << endl;
    
}
    
void StMixedEventBuffer::removeFirstEvent() {

    
    delete mEventsList[0];
    
    mEventsList.erase(mEventsList.begin());   
    
    mBufferCounter--;
}    
    
    