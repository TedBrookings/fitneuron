#include "HodgkinHuxleyChannel.h"


// Connor-Stevens model of Fast Sodium current (mPower = 3, hPower = 1)
// Adjustable with vShift, vScale, tauScale
class NaCS : public HodgkinHuxleyChannel<3,1>
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
      
    double hAlphaCoef;     // (kHz)
    double hAlphaOffset;   // (mV)
    double hAlphaScale;    // (mV)
    double hBetaCoef;      // (kHz)
    double hBetaOffset;    // (mV)
    double hBetaScale;     // (mV)

    double eNa;            // (mV)

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
    double mAlphaCoefF;    // (kHz)
    double mAlphaOffsetF;  // (mV)
    double mAlphaScaleF;   // (mV)
    double mBetaCoefF;     // (kHz)
    double mBetaOffsetF;   // (mV)
    double mBetaScaleF;    // (mV)
      
    double hAlphaCoefF;    // (kHz)
    double hAlphaOffsetF;  // (mV)
    double hAlphaScaleF;   // (mV)
    double hBetaCoefF;     // (kHz)
    double hBetaOffsetF;   // (mV)
    double hBetaScaleF;    // (mV)
};



void
NaCS::defineParameters(ParameterDescriptionList & parametersList,
                       const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 6.3;                        // (C)
  DEFINE_PARAMETER( gBarQ10,      2.0,   UNITLESS ); // wild guess
  DEFINE_PARAMETER( mQ10,         3.0,   UNITLESS ); // wild guess
  DEFINE_PARAMETER( hQ10,         3.0,   UNITLESS ); // wild guess
  
  // fitable parameters
  DEFINE_PARAMETER( gBar,         120.0, "uS/nF" );
  DEFINE_PARAMETER( eNa,          55.0,  "mV" );
  
  DEFINE_PARAMETER( mAlphaCoef,   1.0,   "kHz" );
  DEFINE_PARAMETER( mAlphaOffset, 29.7,  "mV" );
  DEFINE_PARAMETER( mAlphaScale,  10.0,  "mV" );
  DEFINE_PARAMETER( mBetaCoef,    4.0,   "kHz" );
  DEFINE_PARAMETER( mBetaOffset,  54.7,  "mV" );
  DEFINE_PARAMETER( mBetaScale,   18.0,  "mV" );

  DEFINE_PARAMETER( hAlphaCoef,   0.07,  "kHz" );
  DEFINE_PARAMETER( hAlphaOffset, 48.0,  "mV" );
  DEFINE_PARAMETER( hAlphaScale,  20.0,  "mV" );
  DEFINE_PARAMETER( hBetaCoef,    1.0,   "kHz" );
  DEFINE_PARAMETER( hBetaOffset,  18.0,  "mV" );
  DEFINE_PARAMETER( hBetaScale,   10.0,  "mV" );

  // fudge-factors
  DEFINE_PARAMETER( vShift,      0.0,    "mV" );
  DEFINE_PARAMETER( mVShift,     0.0,    "mV" );
  DEFINE_PARAMETER( hVShift,     0.0,    "mV" );
  DEFINE_PARAMETER( vScale,      1.0,    UNITLESS );
  DEFINE_PARAMETER( mVScale,     1.0,    UNITLESS );
  DEFINE_PARAMETER( hVScale,     1.0,    UNITLESS );
  DEFINE_PARAMETER( tauScale,    1.0,    UNITLESS );
  DEFINE_PARAMETER( mTauScale,   1.0,    UNITLESS );
  DEFINE_PARAMETER( hTauScale,   1.0,    UNITLESS );
}



void
NaCS::initialize(void)
{
  // initialize the OhmicChannel base object
  HodgkinHuxleyChannel<3,1>::initialize();
  
  // use initialize to effect fudge-factor shifts and scales

  // set reversal potential
  E = eNa;

  // effect shifts and scales on fudged model parameters
  mAlphaOffsetF = mAlphaOffset + vShift + mVShift;
  mBetaOffsetF = mBetaOffset + vShift + mVShift;
  hAlphaOffsetF = hAlphaOffset + vShift + hVShift;
  hBetaOffsetF = hBetaOffset + vShift + hVShift;
  
  mAlphaScaleF = mAlphaScale * vScale * mVScale;
  mBetaScaleF = mBetaScale * vScale * mVScale;
  hAlphaScaleF = hAlphaScale * vScale * hVScale;
  hBetaScaleF = hBetaScale * vScale * hVScale;
  
  mAlphaCoefF = mAlphaCoef / (tauScale * mTauScale);
  mBetaCoefF = mBetaCoef / (tauScale * mTauScale);
  hAlphaCoefF = hAlphaCoef / (tauScale * hTauScale);
  hBetaCoefF = hBetaCoef / (tauScale * hTauScale);
}



void
NaCS::dynamics(const FDouble & v, const double & t) const
{
  (void)t;
  
  mAlpha = mAlphaCoefF * linoid((v + mAlphaOffsetF) / mAlphaScaleF);
  mBeta = mBetaCoefF * exp(-(v + mBetaOffsetF) / mBetaScaleF);
  
  hAlpha = hAlphaCoefF * exp(-(v + hAlphaOffsetF) / hAlphaScaleF);
  hBeta = scaleLogistic(hBetaCoefF, (v + hBetaOffsetF) / hBetaScaleF);
}


// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, NaCS, injectorFactory);
