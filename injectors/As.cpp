#include "OhmicChannel.h"


// Slow component of transient potassium current (mPower = 3, hPower = 1)
// Base formula from Turrigiano et al., 1995
// Cultured spiny lobster STG cells, recorded at room temperature (22-25 C)
// Adjustable with vShift, vScale, tauScale
class As : public OhmicChannel<3,1>
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
  
    double hOffset;        // (mV)
    double hScale;         // (mV)
    double hTauA;          // (ms)
    double hTauB;          // (ms)
    double hTauOffset;     // (mV)
    double hTauVScale;     // (mV)
    
    double eK;             // (mV)

    // fudge-factors
    double vShift;         // (mV)
    double mVShift;        // (mV)  
    double hVShift;        // (mV)  
    double vScale;         // (pure)
    double mVScale;        // (pure)
    double hVScale;        // (pure)
    double tauScale;       // (pure)
    double mTauScale;      // (pure)
    double hTauScale;      // (pure)
    
    // fudged model parameters
    double mOffsetF;       // (mV)
    double mScaleF;        // (mv)
    double mTauAF;         // (ms)
    double mTauBF;         // (ms)
    double mTauOffsetF;    // (mV)
    double mTauVScaleF;    // (mV)
  
    double hOffsetF;       // (mV)
    double hScaleF;        // (mV)
    double hTauAF;         // (ms)
    double hTauBF;         // (ms)
    double hTauOffsetF;    // (mV)
    double hTauVScaleF;    // (mV)
};



void
As::defineParameters(ParameterDescriptionList & parametersList,
                     const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 23.5;                     // (C)  (Actually 22-25 C)
  // A channel Q10s from Tang et al 2010
  DEFINE_PARAMETER( gBarQ10,    1.90,   UNITLESS ); // +-0.44
  DEFINE_PARAMETER( mQ10,       3.00,   UNITLESS ); // +-0.16
  DEFINE_PARAMETER( hQ10,       3.78,   UNITLESS ); // +-0.18
  
  // parameters of the model
  DEFINE_PARAMETER( gBar,       0.0,    "uS/nF" );
  DEFINE_PARAMETER( eK,        -80.0,   "mV" );

  DEFINE_PARAMETER( mOffset,    24.3,   "mV" );
  DEFINE_PARAMETER( mScale,     9.4,    "mV" );
  DEFINE_PARAMETER( mTauA,      13.3,   "ms" );
  DEFINE_PARAMETER( mTauB,      9.0,    "ms" );
  DEFINE_PARAMETER( mTauOffset, 50.3,   "mV" );
  DEFINE_PARAMETER( mTauVScale, 11.8,   "mV" );

  DEFINE_PARAMETER( hOffset,    61.3,   "mV" );
  DEFINE_PARAMETER( hScale,     6.6,    "mV" );
  DEFINE_PARAMETER( hTauA,      9821.0, "ms" );
  DEFINE_PARAMETER( hTauB,      9269.0, "ms" );
  DEFINE_PARAMETER( hTauOffset, 69.9,   "mV" );
  DEFINE_PARAMETER( hTauVScale, 4.6,    "mV" );

  // fudge-factors
  DEFINE_PARAMETER( vShift,     0.0,    "mV" );
  DEFINE_PARAMETER( mVShift,    0.0,    "mV" );
  DEFINE_PARAMETER( hVShift,    0.0,    "mV" );
  DEFINE_PARAMETER( vScale,     1.0,     UNITLESS );
  DEFINE_PARAMETER( mVScale,    1.0,     UNITLESS );
  DEFINE_PARAMETER( hVScale,    1.0,     UNITLESS );
  DEFINE_PARAMETER( tauScale,   1.0,     UNITLESS );
  DEFINE_PARAMETER( mTauScale,  1.0,     UNITLESS );
  DEFINE_PARAMETER( hTauScale,  1.0,     UNITLESS );
}



void
As::initialize(void)
{
  // initialize the OhmicChannel base object
  OhmicChannel<3,1>::initialize();
  
  // use initialize to effect fudge-factor shifts and scales

  // set reversal potential
  E = eK;
  
  // effect shifts and scales on fudged model parameters
  mOffsetF = mOffset + vShift + mVShift;
  mTauOffsetF = mTauOffset + vShift + mVShift;
  hOffsetF = hOffset + vShift + hVShift;
  hTauOffsetF = hTauOffset + vShift + hVShift;
  
  mScaleF = mScale * vScale * mVScale;
  mTauVScaleF = mTauVScale * vScale * mVScale;
  hScaleF = hScale * vScale * hVScale;
  hTauVScaleF = hTauVScale * vScale * hVScale;
  
  mTauAF = mTauA * tauScale * mTauScale;
  mTauBF = mTauB * tauScale * mTauScale;
  hTauAF = hTauA * tauScale * hTauScale;
  hTauBF = hTauB * tauScale * hTauScale;
}



void
As::dynamics(const FDouble & v, const double & t) const
{
  (void)t;
  
  mInf = logistic((v + mOffsetF) / mScaleF);
  mTau = mTauAF - scaleLogistic(mTauBF, (v + mTauOffsetF) / mTauVScaleF);

  hInf = logistic(-(v + hOffsetF) / hScaleF);
  hTau = hTauAF - scaleLogistic(hTauBF, (v + hTauOffsetF) / hTauVScaleF);
}


// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, As, injectorFactory);
