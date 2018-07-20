#ifndef StMixerTrack_hh
#define StMixerTrack_hh


/* **************************************************
 *
 * Track class used for mixed event buffer
 *
 *
 * ************************************************/

#include "StarClassLibrary/StThreeVectorF.hh"

class StPicoTrack;

class StMixerTrack{
 public:
  StMixerTrack();
  StMixerTrack(StMixerTrack const * trk);
  ~StMixerTrack(){/*cout << "track destructor called" << endl*/;};
  StMixerTrack(StThreeVectorF trackMomentum, int charge, int PID);
  StMixerTrack(StThreeVectorF trackMomentum, double mass, int kaonIdx, int pionIdx, int charge, int kaonCharge, int pionCharge, StThreeVectorF kaonMom, StThreeVectorF pionMom);
  
  StThreeVectorF const gMom() const;
  int particleType();
  int charge();
  double mass();         
  int kaonIdx();      
  int pionIdx();   
  StThreeVectorF const kaonGMom() const; 
  StThreeVectorF const pionGMom() const;  
  int kaonCharge();
  int pionCharge();
  
  
 private:
  
  StThreeVectorF mMom;
  int mPIDFlag;
  int mCharge;
  double mMass;
  int mKaonIdx;
  int mPionIdx;
  StThreeVectorF mKaonMom;
  StThreeVectorF mPionMom;
  int mKaonCharge;
  int mPionCharge;
  
};

inline StThreeVectorF const StMixerTrack::gMom() const { return mMom; }
inline int StMixerTrack::particleType() { return mPIDFlag; }
inline int StMixerTrack::charge()       { return mCharge;  }
inline double StMixerTrack::mass()         { return mMass;    }
inline int StMixerTrack::kaonIdx()      { return mKaonIdx;    }
inline int StMixerTrack::pionIdx()      { return mPionIdx;    }
inline int StMixerTrack::kaonCharge()      { return mKaonCharge;    }
inline int StMixerTrack::pionCharge()      { return mPionCharge;    }
inline StThreeVectorF const StMixerTrack::kaonGMom() const { return mKaonMom; }
inline StThreeVectorF const StMixerTrack::pionGMom() const { return mPionMom; }

#endif
