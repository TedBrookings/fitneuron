#include "OhmicChannel.h"


// Calcium-activated Potasium current (mPower = 4, hPower = 0)
// Base formula from Turrigiano et al., 1995
// Cultured spiny lobster STG cells, recorded at room temperature (22-25 C)
// Adjustable with vShift, vScale, tauScale
DECLARE_SOLUTE(Ca, 2);
class KCa : public OhmicChannel<4,0,Ca>
{
  public:
    // Define the parameters of the ion channel
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const list<Trace> & traces,
                                  NeuronModel & model);
    // use initialize to effect fudge-factor shifts and scales
    virtual void initialize(void);
  
  protected:
    // model activation dynamics
    //   (assign mInf, mTau)
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
    
    double eK;             // (mV)

    Solute* ca;            // Solute object for calcium concentrations/current
    double cHalfActive;    // (mM)  

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
KCa::defineParameters(ParameterDescriptionList & parametersList,
                      const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 23.5;                      // (C)  (Actually 22-25 C)
  DEFINE_PARAMETER( gBarQ10,     2.0,    UNITLESS ); // wild guess
  DEFINE_PARAMETER( mQ10,        3.0,    UNITLESS ); // wild guess
  
  // parameters of the model
  DEFINE_PARAMETER( gBar,        0.0,    "uS/nF" );
  DEFINE_PARAMETER( eK,         -80.0,   "mV" );

  DEFINE_PARAMETER( mOffset,     28.3,   "mV" );
  DEFINE_PARAMETER( mScale,      12.6,   "mV" );
  DEFINE_PARAMETER( mTauA,       90.3,   "ms" );
  DEFINE_PARAMETER( mTauB,       75.1,   "ms" );
  DEFINE_PARAMETER( mTauOffset,  46.0,   "mV" );
  DEFINE_PARAMETER( mTauVScale,  22.7,   "mV" );

  DEFINE_PARAMETER( cHalfActive, 3.0e-3, "mm" );

  // fudge-factors
  DEFINE_PARAMETER( vShift,       0.0,   "mV" );
  DEFINE_PARAMETER( mVShift,      0.0,   "mV" );
  DEFINE_PARAMETER( vScale,       1.0,   UNITLESS );
  DEFINE_PARAMETER( mVScale,      1.0,   UNITLESS );
  DEFINE_PARAMETER( tauScale,     1.0,   UNITLESS );
  DEFINE_PARAMETER( mTauScale,    1.0,   UNITLESS );
}



void
KCa::initialize(void)
{
  // initialize the OhmicChannel base object
  OhmicChannel<4,0,Ca>::initialize();

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
KCa::dynamics(const FDouble & v, const double & t,
              tuple<const FDouble&> & cIntTuple,
              tuple<const double&> & cExtTuple) const
{
  (void)t;
  (void)cExtTuple;
  const FDouble & caInt = get<0>(cIntTuple);
  
  mInf = scaleLogistic(caInt / (cHalfActive + caInt),
                       (v + mOffsetF) / mScaleF);
  mTau = mTauAF - scaleLogistic(mTauBF, (v + mTauOffsetF) / mTauVScaleF);
}


// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, KCa, injectorFactory);
