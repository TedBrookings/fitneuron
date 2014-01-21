#include "OhmicChannel.h"


// Persistant Sodium current
// Base formula from Turrigiano et al., 1995
// Cultured spiny lobster STG cells, recorded at room temperature (22-25 C)
// Adjustable with vShift, vScale, tauScale
class Nap : public OhmicChannel<3,1>
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
    double hTauOffset;     // (mV)
    double hTauVScale;     // (mV)

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
    double hTauOffsetF;    // (mV)
    double hTauVScaleF;    // (mV)
};



void
Nap::defineParameters(ParameterDescriptionList & parametersList,
                      const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 23.5;                    // (C)  (Actually 22-25 C)
  DEFINE_PARAMETER( gBarQ10,    2.0,   UNITLESS ); // wild guess
  DEFINE_PARAMETER( mQ10,       3.0,   UNITLESS ); // wild guess
  DEFINE_PARAMETER( hQ10,       3.0,   UNITLESS ); // wild guess
  
  // parameters of the model
  DEFINE_PARAMETER( gBar,       0.0,   "uS/nF" );
  DEFINE_PARAMETER( eNa,        55.0,  "mV" );

  DEFINE_PARAMETER( mOffset,    26.8,  "mV" );
  DEFINE_PARAMETER( mScale,     8.2,   "mV" );
  DEFINE_PARAMETER( mTauA,      19.8,  "ms" );
  DEFINE_PARAMETER( mTauB,      10.7,  "ms" );
  DEFINE_PARAMETER( mTauOffset, 26.5,  "mV" );
  DEFINE_PARAMETER( mTauVScale, 8.6,   "mV" );

  DEFINE_PARAMETER( hOffset,    48.5,  "mV" );
  DEFINE_PARAMETER( hScale,     4.8,   "mV" );
  DEFINE_PARAMETER( hTauA,      666.0, "ms" );
  DEFINE_PARAMETER( hTauB,      379.0, "ms" );
  DEFINE_PARAMETER( hTauOffset, 33.6,  "mV" );
  DEFINE_PARAMETER( hTauVScale, 11.7,  "mV" );

  // fudge-factors
  DEFINE_PARAMETER( vShift,     0.0,   "mV" );
  DEFINE_PARAMETER( mVShift,    0.0,   "mV" );
  DEFINE_PARAMETER( hVShift,    0.0,   "mV" );
  DEFINE_PARAMETER( vScale,     1.0,   UNITLESS );
  DEFINE_PARAMETER( mVScale,    1.0,   UNITLESS );
  DEFINE_PARAMETER( hVScale,    1.0,   UNITLESS );
  DEFINE_PARAMETER( tauScale,   1.0,   UNITLESS );
  DEFINE_PARAMETER( mTauScale,  1.0,   UNITLESS );
  DEFINE_PARAMETER( hTauScale,  1.0,   UNITLESS );
}



void
Nap::initialize(void)
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
Nap::dynamics(const FDouble & v, const double & t) const
{
  (void)t;
  
  mInf = logistic((v + mOffsetF) / mScaleF);
  mTau = mTauAF - scaleLogistic(mTauBF, (v + mTauOffsetF) / mTauVScaleF);
  
  hInf = logistic(-(v + hOffsetF) / hScaleF);
  hTau = hTauAF - scaleLogistic(hTauBF, (v + hTauOffsetF) / hTauVScaleF);
}


// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, Nap, injectorFactory);
