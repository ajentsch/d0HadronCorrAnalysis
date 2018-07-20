#ifndef StMixerEvent_hh
#define StMixerEvent_hh

#include <vector>
#include "StarClassLibrary/StThreeVectorF.hh"

#include "StMixerTrack.h"

#include "StPicoCharmContainers/StKaonPion.h"
/* **************************************************
 *
 * Event class used for mixed event buffer.
 * 
 * 
 *
 * ****************************************************/

class StMixerTrack;


class StMixerEvent{
 public:
  StMixerEvent();
  StMixerEvent(StMixerEvent* event);  //copy constructor
  StMixerEvent(StThreeVectorF, float);
  ~StMixerEvent(){;};  //destructor
  
  void addTrack(StMixerTrack);
  void addKaonPion(StMixerTrack);
  void setPos( float const, float const, float const);
  void setField( float const );
  int getNoTracks();
  int getNoPosiTracks(); //h-h
  int getNoNegaTracks(); //h-h
  StMixerTrack getTrack(int index);
  StMixerTrack getKaonPionAt(int index);
  StMixerTrack getPosiTrack(int index);  //h-h
  StMixerTrack getNegaTrack(int index);  //h-h
  int getKaonPionListSize();
  
  StThreeVectorF const & vertex() const;
  double const field() const;
 
 private:
  StThreeVectorF mVtx;
  float mBField;
  std::vector <StMixerTrack> mTracks;
  std::vector <StMixerTrack> mPosiTracks; //h-h
  std::vector <StMixerTrack> mNegaTracks; //h-h
  std::vector <StMixerTrack> mKaonPionTracks;
 
};

inline void StMixerEvent::setPos( float const vx, float const vy, float const vz){
  mVtx = StThreeVectorF(vx, vy, vz);
}
inline void StMixerEvent::setField( float const field ){ mBField = field; }
inline int StMixerEvent::getNoTracks(){ return mTracks.size();}
inline int StMixerEvent::getNoPosiTracks(){ return mPosiTracks.size();} //h-h
inline int StMixerEvent::getNoNegaTracks(){ return mNegaTracks.size();} //h-h
inline StMixerTrack StMixerEvent::getTrack(int index){ return mTracks.at(index); }
inline StThreeVectorF const & StMixerEvent::vertex() const { return mVtx; }
inline double const StMixerEvent::field() const {return mBField; }
inline StMixerTrack StMixerEvent::getKaonPionAt(int index) { return mKaonPionTracks.at(index); }
inline StMixerTrack StMixerEvent::getPosiTrack(int index) { return mPosiTracks.at(index); } //h-h
inline StMixerTrack StMixerEvent::getNegaTrack(int index) { return mNegaTracks.at(index); } //h-h
inline int StMixerEvent::getKaonPionListSize() { return mKaonPionTracks.size(); }
#endif
