#ifndef StMixedEventBuffer_h
#define StMixedEventBuffer_h

//This is the class for defining buffers
//used for storing events to do event
//mixing


//Code by Alex Jentsch

#include <vector>
#include "StarClassLibrary/StThreeVectorF.hh"


class TTree;
class TH2F;
class StPicoEvent;
class StPicoTrack;
class StPicoDst;
class StPicoD0Event;
class StKaonPion;
class StMixerEvent;
class StMixerTrack;


class StMixedEventBuffer{
public:
    StMixedEventBuffer();
    virtual ~StMixedEventBuffer();
    //void addEvent(StPicoDst const* picoDst);            // This function should have the event information passed to it
    void addEvent(StPicoDst const* picoDst);
    
    void setBufferSize(int buffer);
    void setBufferCounter(int buffer);
    int  getBufferSize();
    int  getBufferIndex();
    int getBufferMaxSize();
    void addTrackToEvent(int index, StThreeVectorF trackMom, int charge, int PID);
    void addKaonPionToEvent(int, StThreeVectorF kaonPionMom, double mass, int kaonIdx, int pionIdx, int charge, int kaonCharge, int pionCharge, StThreeVectorF kaonMom, StThreeVectorF pionMom);
    StMixerEvent* getEvent(int index);
    void clearBuffer();
    void removeFirstEvent();
    
private:                 
    
    std::vector <StMixerEvent*> mEventsList;
    unsigned short int mBufferSize;    // number of events to store in each buffer category
    unsigned short int mBufferCounter;
   
ClassDef(StMixedEventBuffer, 0)   

};

inline void StMixedEventBuffer::setBufferSize(int buffer){ mBufferSize = buffer;}          //Sets the size of the buffer
inline void StMixedEventBuffer::setBufferCounter(int buffer){ mBufferCounter = buffer;}    //incrementing counter for every event added --- set to zero after a mix is complete
inline int  StMixedEventBuffer::getBufferIndex() { return mBufferCounter; }
inline int  StMixedEventBuffer::getBufferSize()  { return mEventsList.size();    }
inline int  StMixedEventBuffer::getBufferMaxSize() { return mBufferSize; } 
inline StMixerEvent* StMixedEventBuffer::getEvent(int index){ return mEventsList.at(index); }
//inline void StMixedEventBuffer::clearBuffer() { mEventsList.erase(mEventsList.begin(), mEventsList.end()); }

#endif
