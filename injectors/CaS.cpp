#include "NernstChannel.h"



// Non-inactivating calcium channel (mPower = 3, hPower = 0)
// Base formula from Turrigiano et al., 1995 (Ca1 in that paper)
// Cultured spiny lobster STG cells, recorded at room temperature (22-25 C)
// Adjustable with vShift, vScale, tauScale
DECLARE_SOLUTE(Ca, 2);
class CaS : public NernstChannel<3,0,Ca>
{
  public:
    // Define the parameters of the ion channel
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const list<Trace> & traces,
                                  NeuronModel & model);
    // use initialize to effect fudge-factor shifts and scales
    virtual void initialize(void);
  
  protected:
    // model activation dynamics
    //   (assign mInf, mTau)
    virtual void dynamics(const FDouble & v, const double & t,
                          tuple<const FDouble&> & cIntTuple,
                          tuple<const double&> & cExtTuple) const;

    // model parameters
    double mOffset;        // (mV)
    double mScale;         // (mV)
    double mTauA;          // (ms)
    double mTauB;          // (ms)
    double mTauOffset;     // (mV)
    double mTauVScale;     // (mV)
  
    // fudge-factors
    double vShift;         // (mV)
    double mVShift;        // (mV)  
    double vScale;         // (pure)
    double mVScale;        // (pure)
    double tauScale;       // (pure)
    double mTauScale;      // (pure)
    
    // fudged model parameters
    double mOffsetF;       // (mV)
    double mScaleF;        // (mv)
    double mTauAF;         // (ms)
    double mTauBF;         // (ms)
    double mTauOffsetF;    // (mV)
    double mTauVScaleF;    // (mV)
};



void
CaS::defineParameters(ParameterDescriptionList & parametersList,
                      const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 23.5;                    // (C)  (Actually 22-25 C)
  DEFINE_PARAMETER( gBarQ10,    2.0,  UNITLESS ); // wild guess
  DEFINE_PARAMETER( mQ10,       3.0,  UNITLESS ); // wild guess
  
  // parameters of the model
  DEFINE_PARAMETER( gBar,       0.0,  "uS/nF" ); // (uS/nF)

  DEFINE_PARAMETER( mOffset,    21.6, "mV" );  // (mV)
  DEFINE_PARAMETER( mScale,     8.5,  "mV" );  // (mV)
  DEFINE_PARAMETER( mTauA,      16.0, "ms" );  // (ms)
  DEFINE_PARAMETER( mTauB,      13.1, "ms" );  // (ms)
  DEFINE_PARAMETER( mTauOffset, 25.1, "mV" );  // (mV)
  DEFINE_PARAMETER( mTauVScale, 26.4, "mV" );  // (mV)

  // fudge-factors
  DEFINE_PARAMETER( vShift,       0.0, "mV" );
  DEFINE_PARAMETER( mVShift,      0.0, "mV" );
  DEFINE_PARAMETER( vScale,       1.0,  UNITLESS );
  DEFINE_PARAMETER( mVScale,      1.0,  UNITLESS );
  DEFINE_PARAMETER( tauScale,     1.0,  UNITLESS );
  DEFINE_PARAMETER( mTauScale,    1.0,  UNITLESS );
}



void
CaS::initialize(void)
{
  // initialize the NernstChannel base object
  NernstChannel<3,0,Ca>::initialize();

  // use initialize to effect fudge-factor shifts and scales
  
  // effect shifts and scales on fudged model parameters
  mOffsetF = mOffset + vShift + mVShift;
  mTauOffsetF = mTauOffset + vShift + mVShift;
  
  mScaleF = mScale * vScale * mVScale;
  mTauVScaleF = mTauVScale * vScale * mVScale;
  
  mTauAF = mTauA * tauScale * mTauScale;
  mTauBF = mTauB * tauScale * mTauScale;  
}



void
CaS::dynamics(const FDouble & v, const double & t,
              tuple<const FDouble&> & cIntTuple,
              tuple<const double&> & cExtTuple) const
{
  (void)t;
  (void) cIntTuple;
  (void) cExtTuple;
  
  mInf = logistic((v + mOffsetF) / mScaleF);
  mTau = mTauAF - scaleLogistic(mTauBF, (v + mTauOffsetF) / mTauVScaleF);
}


// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, CaS, injectorFactory);
