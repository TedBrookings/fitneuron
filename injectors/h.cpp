#include "OhmicChannel.h"



// Hyperpolarization-activated Sodium and Potasium channel
// (mPower = 1, hPower = 0)
// Base formula from Turrigiano et al., 1995
// Cultured spiny lobster STG cells, recorded at room temperature (22-25 C)
// Adjustable with vShift, vScale, tauScale
class h : public OhmicChannel<1,0>
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
h::defineParameters(ParameterDescriptionList & parametersList,
                    const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 23.5;                     // (C)  (Actually 22-25 C)
  DEFINE_PARAMETER( gBarQ10,    3.13,   UNITLESS ); // +- 0.46, Tang et al 2010
  DEFINE_PARAMETER( mQ10,       3.36,   UNITLESS ); // +- 0.18, Tang et al 2010
  
  // parameters of the model
  DEFINE_PARAMETER( gBar,       0.0,    "uS/nF" );
  DEFINE_PARAMETER( E,         -32.3,   "mV" ); // +- 0.6

  DEFINE_PARAMETER( mOffset,    78.3,   "mV" ); // +- 2.0
  DEFINE_PARAMETER( mScale,     8.7,    "mV" ); // +- 0.6
  DEFINE_PARAMETER( mTauA,      272.0,  "ms" );
  DEFINE_PARAMETER( mTauB,     -1499.0, "ms" );
  DEFINE_PARAMETER( mTauOffset, 42.2,   "mV" );
  DEFINE_PARAMETER( mTauVScale, 8.73,   "mV" );

  // fudge-factors
  DEFINE_PARAMETER( vShift,       0.0,  "mV" );
  DEFINE_PARAMETER( mVShift,      0.0,  "mV" );
  DEFINE_PARAMETER( vScale,       1.0,  UNITLESS );
  DEFINE_PARAMETER( mVScale,      1.0,  UNITLESS );
  DEFINE_PARAMETER( tauScale,     1.0,  UNITLESS );
  DEFINE_PARAMETER( mTauScale,    1.0,  UNITLESS );
}



void
h::initialize(void)
{
  // initialize the OhmicChannel base object
  OhmicChannel<1,0>::initialize();

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
h::dynamics(const FDouble & v, const double & t) const
{
  (void)t;
  
  mInf = logistic(-(v + mOffsetF) / mScaleF);
  mTau = mTauAF - scaleLogistic(mTauBF, (v + mTauOffsetF) / mTauVScaleF);
}


// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, h, injectorFactory);
