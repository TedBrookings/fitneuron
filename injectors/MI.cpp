#include "OhmicChannel.h"


// Modulator-induced Inward current, AKA Proctolin-activated current
// From 1992 J. Neurophys paper, Golowasch, Buchholtz, Epstein, Marder
// Charge carrier is probably Na and Cl
// Adjustable with vShift, vScale, tauScale
class MI : public OhmicChannel<1,0>
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
    double vHalfActive;    // (mV)
    double mScale;         // (mV)
    double mTauA;          // (ms)
  
    // fudge-factors
    double vShift;         // (mV)
    double mVShift;        // (mV)  
    double vScale;         // (pure)
    double mVScale;        // (pure)
    double tauScale;       // (pure)
    double mTauScale;      // (pure)
    
    // fudged model parameters
    double vHalfActiveF;   // (mV)
    double mScaleF;        // (mv)
    double mTauAF;         // (ms)
};



void
MI::defineParameters(ParameterDescriptionList & parametersList,
                     const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 10.0;                      // (C)
  DEFINE_PARAMETER( gBarQ10,      2.0,  UNITLESS ); // wild guess
  DEFINE_PARAMETER( mQ10,         3.0,  UNITLESS ); // wild guess
  
  // parameters of the model
  DEFINE_PARAMETER( gBar,         0.0,  "uS/nF" );
  // E was set to -10 by Adam Taylor
  DEFINE_PARAMETER( E,            0.0,  "mV"); // probably in [-20 , 20]
  // rough guestimate for mOffset:
  //  mOffset = mScale * log((E-vPeak)/mScale - 1) - vPeak
  //    and vPeak in [-15, -5]
  // mOffset was set to 45.0 by Adam Taylor
  DEFINE_PARAMETER( vHalfActive, -30.0, "mV" ); // probably in [-50 , -10]
  DEFINE_PARAMETER( mScale,       5.0,  "mV" );
  DEFINE_PARAMETER( mTauA,        6.0,  "ms" );

  // fudge-factors
  DEFINE_PARAMETER( vShift,       0.0,  "mV" );
  DEFINE_PARAMETER( mVShift,      0.0,  "mV" );
  DEFINE_PARAMETER( vScale,       1.0,  UNITLESS );
  DEFINE_PARAMETER( mVScale,      1.0,  UNITLESS );
  DEFINE_PARAMETER( tauScale,     1.0,  UNITLESS );
  DEFINE_PARAMETER( mTauScale,    1.0,  UNITLESS );
}



void
MI::initialize(void)
{
  // initialize the OhmicChannel base object
  OhmicChannel<1,0>::initialize();

  // use initialize to effect fudge-factor shifts and scales
  
  // effect shifts and scales on fudged model parameters
  vHalfActiveF = vHalfActive - vShift - mVShift;

  mScaleF = mScale * vScale * mVScale;
  mTauAF = mTauA * tauScale * mTauScale;
}



void
MI::dynamics(const FDouble & v, const double & t) const
{
  (void)t;
  
  mInf = logistic((v - vHalfActiveF) / mScaleF);
  mTau = mTauAF;
}


// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, MI, injectorFactory);
