#include "NernstChannel.h"



// Inactivating calcium channel (mPower = 3, hPower = 1)
// Base formula from Turrigiano et al., 1995 (Ca1 in that paper)
// Cultured spiny lobster STG cells, recorded at room temperature (22-25 C)
// Adjustable with vShift, vScale, tauScale
DECLARE_SOLUTE(Ca, 2);
class CaT : public NernstChannel<3,1,Ca>
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
    //   (assign mInf, mTau, hInf, hTau)
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
  
    double hOffset;        // (mV)
    double hScale;         // (mV)
    double hTauA;          // (ms)
    double hTauB;          // (ms)
    double hTauOffset;     // (mV)
    double hTauVScale;     // (mV)

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
CaT::defineParameters(ParameterDescriptionList & parametersList,
                      const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 23.5;                    // (C)  (Actually 22-25 C)
  DEFINE_PARAMETER( gBarQ10,    2.0,   UNITLESS ); // wild guess
  DEFINE_PARAMETER( mQ10,       3.0,   UNITLESS ); // wild guess
  DEFINE_PARAMETER( hQ10,       3.0,   UNITLESS ); // wild guess
  
  // parameters of the model
  DEFINE_PARAMETER( gBar,       0.0,   "uS/nF" ); // (uS/nF)

  DEFINE_PARAMETER( mOffset,    27.1,  "mV" );  // (mV)
  DEFINE_PARAMETER( mScale,     7.18,  "mV" );  // (mV)
  DEFINE_PARAMETER( mTauA,      21.7,  "ms" );  // (ms)
  DEFINE_PARAMETER( mTauB,      21.3,  "ms" );  // (ms)
  DEFINE_PARAMETER( mTauOffset, 68.1,  "mV" );  // (mV)
  DEFINE_PARAMETER( mTauVScale, 20.5,  "mV" );  // (mV)

  DEFINE_PARAMETER( hOffset,    30.1,  "mV" );  // (mV)
  DEFINE_PARAMETER( hScale,     5.5,   "mV" );  // (mV)
  DEFINE_PARAMETER( hTauA,      105.0, "ms" );  // (ms)
  DEFINE_PARAMETER( hTauB,      89.8,  "ms" );  // (ms)
  DEFINE_PARAMETER( hTauOffset, 55.0,  "mV" );  // (mV)
  DEFINE_PARAMETER( hTauVScale, 16.9,  "mV" );  // (mV)

  // fudge-factors
  DEFINE_PARAMETER( vShift,       0.0, "mV" );
  DEFINE_PARAMETER( mVShift,      0.0, "mV" );
  DEFINE_PARAMETER( hVShift,      0.0, "mV" );
  DEFINE_PARAMETER( vScale,       1.0,  UNITLESS );
  DEFINE_PARAMETER( mVScale,      1.0,  UNITLESS );
  DEFINE_PARAMETER( hVScale,      1.0,  UNITLESS );
  DEFINE_PARAMETER( tauScale,     1.0,  UNITLESS );
  DEFINE_PARAMETER( mTauScale,    1.0,  UNITLESS );
  DEFINE_PARAMETER( hTauScale,    1.0,  UNITLESS );
}



void
CaT::initialize(void)
{
  // initialize the NernstChannel base object
  NernstChannel<3,1,Ca>::initialize();

  // use initialize to effect fudge-factor shifts and scales
  
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
CaT::dynamics(const FDouble & v, const double & t,
              tuple<const FDouble&> & cIntTuple,
              tuple<const double&> & cExtTuple) const
{
  (void)t;
  (void)cIntTuple;
  (void)cExtTuple;
  
  mInf = logistic((v + mOffsetF) / mScaleF);
  mTau = mTauAF - scaleLogistic(mTauBF, (v + mTauOffsetF) / mTauVScaleF);

  hInf = logistic(-(v + hOffsetF) / hScaleF);
  hTau = hTauAF - scaleLogistic(hTauBF, (v + hTauOffsetF) / hTauVScaleF);
}



// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, CaT, injectorFactory);
