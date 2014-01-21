#include "OhmicChannel.h"


// Delayed rectifier potassium current (mPower = 4, hPower = 0)
// Base formula from Turrigiano et al., 1995
// Cultured spiny lobster STG cells, recorded at room temperature (22-25 C)
// Adjustable with vShift, vScale, tauScale
class Kd : public OhmicChannel<4,0>
{
  public:
    // Define the parameters of the ion channel
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const list<Trace> & traces,
                                  NeuronModel & model);
    // use initialize to effect fudge-factor shifts and scales
    virtual void initialize(void);

  protected:
    // model activation and inactivation dynamics
    //   (assign mInf, mTau, hInf, and hTau)
    virtual void dynamics(const FDouble & v, const double & t) const;
    
    // model parameters
    double mOffset;        // (mV)
    double mScale;         // (mV)
    double mTauA;          // (ms)
    double mTauB;          // (ms)
    double mTauOffset;     // (mV)
    double mTauVScale;     // (mV)
    
    double eK;             // (mV)

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
Kd::defineParameters(ParameterDescriptionList & parametersList,
                     const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 23.5;                   // (C)  (Actually 22-25 C)
  DEFINE_PARAMETER( gBarQ10,    2.0,  UNITLESS ); // wild guess
  DEFINE_PARAMETER( mQ10,       3.0,  UNITLESS ); // wild guess
  
  // parameters of the model
  DEFINE_PARAMETER( gBar,       0.0,  "uS/nF" );
  DEFINE_PARAMETER( eK,        -80.0, "mV" );

  DEFINE_PARAMETER( mOffset,    12.3, "mV" );
  DEFINE_PARAMETER( mScale,     11.8, "mV" );
  DEFINE_PARAMETER( mTauA,      7.2,  "ms" );
  DEFINE_PARAMETER( mTauB,      6.4,  "ms" );
  DEFINE_PARAMETER( mTauOffset, 28.3, "mV" );
  DEFINE_PARAMETER( mTauVScale, 19.2, "mV" );

  // fudge-factors
  DEFINE_PARAMETER( vShift,     0.0,  "mV" );
  DEFINE_PARAMETER( mVShift,    0.0,  "mV" );
  DEFINE_PARAMETER( vScale,     1.0,  UNITLESS );
  DEFINE_PARAMETER( mVScale,    1.0,  UNITLESS );
  DEFINE_PARAMETER( tauScale,   1.0,  UNITLESS );
  DEFINE_PARAMETER( mTauScale,  1.0,  UNITLESS );
}



void
Kd::initialize(void)
{
  // initialize the OhmicChannel base object
  OhmicChannel<4,0>::initialize();
  
  // use initialize to effect fudge-factor shifts and scales
  
  // set reversal potential
  E = eK;
  
  // effect shifts and scales on fudged model parameters
  mOffsetF = mOffset + vShift + mVShift;
  mTauOffsetF = mTauOffset + vShift + mVShift;
  
  mScaleF = mScale * vScale * mVScale;
  mTauVScaleF = mTauVScale * vScale * mVScale;
  
  mTauAF = mTauA * tauScale * mTauScale;
  mTauBF = mTauB * tauScale * mTauScale;
}



void
Kd::dynamics(const FDouble & v, const double & t) const
{
  (void)t;
  
  mInf = logistic((v + mOffsetF) / mScaleF);
  mTau = mTauAF - scaleLogistic(mTauBF, (v + mTauOffsetF) / mTauVScaleF);
}


// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, Kd, injectorFactory);
