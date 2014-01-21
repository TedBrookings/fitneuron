#include "HodgkinHuxleyChannel.h"



// Connor-Stevens model of delayed rectifier potassium current
// Adjustable with vShift, vScale, tauScale
class KdCS : public HodgkinHuxleyChannel<4,0>
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
    double mAlphaCoef;     // (kHz)
    double mAlphaOffset;   // (mV)
    double mAlphaScale;    // (mV)
    double mBetaCoef;      // (kHz)
    double mBetaOffset;    // (mV)
    double mBetaScale;     // (mV)
    
    double eK;             // (mV)
    
    // fudge-factors
    double vShift;         // (mV)
    double mVShift;        // (mV)  
    double vScale;         // (pure)
    double mVScale;        // (pure)
    double tauScale;       // (pure)
    double mTauScale;      // (pure)
    
    // fudged model parameters
    double mAlphaCoefF;    // (kHz)
    double mAlphaOffsetF;  // (mV)
    double mAlphaScaleF;   // (mV)
    double mBetaCoefF;     // (kHz)
    double mBetaOffsetF;   // (mV)
    double mBetaScaleF;    // (mV)
};



void
KdCS::defineParameters(ParameterDescriptionList & parametersList,
                       const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 6.3;                         // (C)
  DEFINE_PARAMETER( gBarQ10,      2.0,    UNITLESS ); // wild guess
  DEFINE_PARAMETER( mQ10,         3.0,    UNITLESS ); // wild guess
  
  // fitable parameters (from Connor-Stevens, but E changed to match STG)
  DEFINE_PARAMETER( gBar,         20.0,   "uS/nF" );
  DEFINE_PARAMETER( eK,          -80.0,   "mV" );
  
  DEFINE_PARAMETER( mAlphaCoef,   0.05,   "kHz" );
  DEFINE_PARAMETER( mAlphaOffset, 45.7,   "mV" );
  DEFINE_PARAMETER( mAlphaScale,  10.0,   "mV" );
  DEFINE_PARAMETER( mBetaCoef,    0.0625, "kHz" );
  DEFINE_PARAMETER( mBetaOffset,  55.7,   "mV" );
  DEFINE_PARAMETER( mBetaScale,   80.0,   "mV" );

  // fudge-factors
  DEFINE_PARAMETER( vShift,     0.0,      "mV" );
  DEFINE_PARAMETER( mVShift,    0.0,      "mV" );
  DEFINE_PARAMETER( vScale,     1.0,      UNITLESS );
  DEFINE_PARAMETER( mVScale,    1.0,      UNITLESS );
  DEFINE_PARAMETER( tauScale,   1.0,      UNITLESS );
  DEFINE_PARAMETER( mTauScale,  1.0,      UNITLESS );
}



void
KdCS::initialize(void)
{
  // initialize the OhmicChannel base object
  HodgkinHuxleyChannel<4,0>::initialize();

  // use initialize to effect fudge-factor shifts and scales
  
  // set reversal potential
  E = eK;
  
  // effect shifts and scales on fudged model parameters
  mAlphaOffsetF = mAlphaOffset + vShift + mVShift;
  mBetaOffsetF = mBetaOffset + vShift + mVShift;
  
  mAlphaScaleF = mAlphaScale * vScale * mVScale;
  mBetaScaleF = mBetaScale * vScale * mVScale;
  
  mAlphaCoefF = mAlphaCoef / (tauScale * mTauScale);
  mBetaCoefF = mBetaCoef / (tauScale * mTauScale);
}



void
KdCS::dynamics(const FDouble & v, const double & t) const
{
  (void)t;
  
  mAlpha = mAlphaCoefF * linoid((v + mAlphaOffsetF) / mAlphaScaleF);
  mBeta = mBetaCoefF * exp(-(v + mBetaOffsetF) / mBetaScaleF);
}


// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, KdCS, injectorFactory);
