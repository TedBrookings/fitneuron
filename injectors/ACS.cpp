#include "OhmicChannel.h"


// Connor-Stevens model of Fast Sodium current (mPower = 3, hPower = 1)
// Adjustable with vShift, vScale, tauScale
class ACS : public OhmicChannel<3,1>
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
    double mCoef;          // (pure)
    double mOffset;        // (mV)
    double mScale;         // (mv)
    double mOffset2;       // (mV)
    double mScale2;        // (mV)
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
    
    double eK;             // (mV);

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
    double mOffset2F;      // (mV)
    double mScale2F;       // (mV)
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
ACS::defineParameters(ParameterDescriptionList & parametersList,
                      const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 6.3;                      // (C)
  // A channel Q10s from Tang et al 2010
  DEFINE_PARAMETER( gBarQ10,    1.90, UNITLESS );  // +-0.44
  DEFINE_PARAMETER( mQ10,       3.00, UNITLESS );  // +-0.16
  DEFINE_PARAMETER( hQ10,       3.78, UNITLESS );  // +-0.18
  
  // fitable parameters (from Connor-Stevens, but E changed to match STG)
  DEFINE_PARAMETER( gBar,       47.7,  "uS/nF" );
  DEFINE_PARAMETER( eK,        -80.0,  "mV" );
  
  DEFINE_PARAMETER( mCoef,      0.0761, UNITLESS );
  DEFINE_PARAMETER( mOffset,    94.22, "mV" );
  DEFINE_PARAMETER( mScale,     31.84, "mV" );
  DEFINE_PARAMETER( mOffset2,   1.17,  "mV" );
  DEFINE_PARAMETER( mScale2,    28.93, "mV" );
  DEFINE_PARAMETER( mTauA,      1.380, "ms" );
  DEFINE_PARAMETER( mTauB,      4.400, "ms" );
  DEFINE_PARAMETER( mTauOffset, 55.96, "mV" );
  DEFINE_PARAMETER( mTauVScale, 20.12, "mV" );

  DEFINE_PARAMETER( hOffset,    53.3,  "mV" );
  DEFINE_PARAMETER( hScale,     14.54, "mV" );
  DEFINE_PARAMETER( hTauA,      4.712, "ms" );
  DEFINE_PARAMETER( hTauB,      10.18, "ms" );
  DEFINE_PARAMETER( hTauOffset, 50.0,  "mV" );
  DEFINE_PARAMETER( hTauVScale, 16.027,"mV" );

  // fudge-factors
  DEFINE_PARAMETER( vShift,     0.0,   "mV" );
  DEFINE_PARAMETER( mVShift,    0.0,   "mV" );
  DEFINE_PARAMETER( hVShift,    0.0,   "mV" );
  DEFINE_PARAMETER( vScale,     1.0,    UNITLESS );
  DEFINE_PARAMETER( mVScale,    1.0,    UNITLESS );
  DEFINE_PARAMETER( hVScale,    1.0,    UNITLESS );
  DEFINE_PARAMETER( tauScale,   1.0,    UNITLESS );
  DEFINE_PARAMETER( mTauScale,  1.0,    UNITLESS );
  DEFINE_PARAMETER( hTauScale,  1.0,    UNITLESS );
}



void
ACS::initialize(void)
{
  // initialize the OhmicChannel base object
  OhmicChannel<3,1>::initialize();

  // use initialize to effect fudge-factor shifts and scales

  // set reversal potential
  E = eK;  
  
  // effect shifts and scales on fudged model parameters
  mOffsetF = mOffset + vShift + mVShift;
  mOffset2F = mOffset2 + vShift + mVShift;
  mTauOffsetF = mTauOffset + vShift + mVShift;
  hOffsetF = hOffset + vShift + hVShift;
  hTauOffsetF = hTauOffset + vShift + hVShift;
  
  mScaleF = mScale * vScale;
  mScale2F = mScale2 * vScale;
  mTauVScaleF = mTauVScale * vScale;
  hScaleF = hScale * vScale;
  hTauVScaleF = hTauVScale * vScale;
  
  mTauAF = mTauA * tauScale * mTauScale;
  mTauBF = mTauB * tauScale * mTauScale;
  hTauAF = hTauA * tauScale * hTauScale;
  hTauBF = hTauB * tauScale * hTauScale;
}



void
ACS::dynamics(const FDouble & v, const double & t) const
{
  (void)t;
  
  mInf = pow(scaleLogistic(mCoef, -(v + mOffset2F)/ mScale2F) *
             exp((v + mOffsetF) / mScaleF), 1.0/3.0);
  if(!(mInf < 1.0))
    mInf = 1.0;
  mTau = mTauAF + scaleLogistic(mTauBF, -(v + mTauOffsetF) / mTauVScaleF);

  hInf = power<4>( logistic(-(v + hOffsetF) / hScaleF) );
  hTau = hTauAF + scaleLogistic(hTauBF, -(v + hTauOffsetF) / hTauVScaleF);
}



// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, ACS, injectorFactory);
