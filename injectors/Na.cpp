#include "OhmicChannel.h"


// Fast Sodium current (mPower = 3, hPower = 1)
// Base formula from Turrigiano et al., 1995
// Cultured spiny lobster STG cells, recorded at room temperature (22-25 C)
// Adjustable with vShift, vScale, tauScale
class Na : public OhmicChannel<3,1>
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

  private:
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
    double hTauOffset1;    // (mV)
    double hTauVScale1;    // (mV)
    double hTauOffset2;    // (mV)
    double hTauVScale2;    // (mV)
    
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
    double hTauOffset1F;   // (mV)
    double hTauVScale1F;   // (mV)
    double hTauOffset2F;   // (mV)
    double hTauVScale2F;   // (mV)
};



void
Na::defineParameters(ParameterDescriptionList & parametersList,
                     const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 23.5;                     // (C)  (Actually 22-25 C)
  DEFINE_PARAMETER( gBarQ10,     2.0,   UNITLESS ); // wild guess
  DEFINE_PARAMETER( mQ10,        3.0,   UNITLESS ); // wild guess
  DEFINE_PARAMETER( hQ10,        3.0,   UNITLESS ); // wild guess
  
  // parameters of the model
  DEFINE_PARAMETER( gBar,        0.0,   "uS/nF" );
  DEFINE_PARAMETER( eNa,         55.0,  "mV" );

  DEFINE_PARAMETER( mOffset,     25.5,  "mV" );
  DEFINE_PARAMETER( mScale,      5.29,  "mV" );
  DEFINE_PARAMETER( mTauA,       1.32,  "ms" );
  DEFINE_PARAMETER( mTauB,       1.26,  "ms" );
  DEFINE_PARAMETER( mTauOffset,  120.0, "mV" );
  DEFINE_PARAMETER( mTauVScale,  25.0,  "mV" );

  DEFINE_PARAMETER( hOffset,     48.9,  "mV" );
  DEFINE_PARAMETER( hScale,      5.18,  "mV" );
  DEFINE_PARAMETER( hTauA,       1.005, "ms" ); // =0.67 * 1.5
  DEFINE_PARAMETER( hTauB,       0.67,  "ms" );
  DEFINE_PARAMETER( hTauOffset1, 62.9,  "mV" );
  DEFINE_PARAMETER( hTauVScale1, 10.0,  "mV" );
  DEFINE_PARAMETER( hTauOffset2, 34.9,  "mV" );
  DEFINE_PARAMETER( hTauVScale2, 3.6,   "mV" );

  // fudge-factors
  DEFINE_PARAMETER( vShift,      0.0,   "mV" );
  DEFINE_PARAMETER( mVShift,     0.0,   "mV" );
  DEFINE_PARAMETER( hVShift,     0.0,   "mV" );
  DEFINE_PARAMETER( vScale,      1.0,   UNITLESS );
  DEFINE_PARAMETER( mVScale,     1.0,   UNITLESS );
  DEFINE_PARAMETER( hVScale,     1.0,   UNITLESS );
  DEFINE_PARAMETER( tauScale,    1.0,   UNITLESS );
  DEFINE_PARAMETER( mTauScale,   1.0,   UNITLESS );
  DEFINE_PARAMETER( hTauScale,   1.0,   UNITLESS );
}



void
Na::initialize(void)
{
  // initialize the OhmicChannel base object
  OhmicChannel<3,1>::initialize();
  
  // use initialize to effect fudge-factor shifts and scales
  
  // set reversal potential
  E = eNa;
  
  // effect shifts and scales on fudged model parameters
  mOffsetF = mOffset + vShift + mVShift;
  mTauOffsetF = mTauOffset + vShift + mVShift;
  hOffsetF = hOffset + vShift + hVShift;
  hTauOffset1F = hTauOffset1 + vShift + hVShift;
  hTauOffset2F = hTauOffset2 + vShift + hVShift;
  
  mScaleF = mScale * vScale * mVScale;
  mTauVScaleF = mTauVScale * vScale * mVScale;
  hScaleF = hScale * vScale * hVScale;
  hTauVScale1F = hTauVScale1 * vScale * hVScale;
  hTauVScale2F = hTauVScale2 * vScale * hVScale;
  
  mTauAF = mTauA * tauScale * mTauScale;
  mTauBF = mTauB * tauScale * mTauScale;
  hTauAF = hTauA * tauScale * hTauScale;
  hTauBF = hTauB * tauScale * hTauScale;
}



void
Na::dynamics(const FDouble & v, const double & t) const
{
  (void)t;
  
  mInf = logistic((v + mOffsetF) / mScaleF);
  mTau = mTauAF - scaleLogistic(mTauBF, (v + mTauOffsetF) / mTauVScaleF);
  
  hInf = logistic(-(v + hOffsetF) / hScaleF);
  hTau = scaleLogistic(
    hTauAF - scaleLogistic(hTauBF, -(v + hTauOffset2F) / hTauVScale2F),
    (v + hTauOffset1F) / hTauVScale1F
  );
}


// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, Na, injectorFactory);
