#include "OhmicChannel.h"



/*
Calcium-Activated Nonspecific current (mPower = 1, hPower = 0)
Model based on
Zhang et al., "Calcium-Dependent Plateau Potentials in a Crab Stomatogastric
Ganglion Motor Neuron. II. Calcium-Activated Slow Inward Current",
J. Neurophys 74(5), 1995

Notes
Activated by caffeine ---caffeine is known to evoke Ca++ release from
intracellular stores in many cell types
*/

DECLARE_SOLUTE(Ca, 2);
class CAN : public OhmicChannel<1,0,Ca>
{
  public:
    // Define the parameters of the ion channel
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const list<Trace> & traces,
                                  NeuronModel & model);
  
  protected:
    // model activation dynamics
    //   (assign mInf, mTau)
    virtual void dynamics(const FDouble & v, const double & t,
                          tuple<const FDouble&> & cIntTuple,
                          tuple<const double&> & cExtTuple) const;

    // model parameters
    double tauUnbind;      // (ms)
    double tauBind;        // (ms) 
    double cBind;          // (mM)
};



void
CAN::defineParameters(ParameterDescriptionList & parametersList,
                      const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 21.0;                    // (C)  (Actually 20-22 C)
  DEFINE_PARAMETER( gBarQ10,   2.0,    UNITLESS ); // wild guess
  DEFINE_PARAMETER( mQ10,      3.0,    UNITLESS ); // wild guess
  
  // parameters of the model
  DEFINE_PARAMETER( gBar,      0.0,    "uS/nF" );
  DEFINE_PARAMETER( E,        -27.0,   "mV" ); // +- 3.5

  DEFINE_PARAMETER( tauUnbind, 627.0,  "ms" ); // +- 34.0
  DEFINE_PARAMETER( tauBind,   300.0,  "ms" ); // try 100 - 500 ms
  DEFINE_PARAMETER( cBind,     3.0e-3, "mM" ); // wild guess
}



void
CAN::dynamics(const FDouble & v, const double & t,
              tuple<const FDouble&> & cIntTuple,
              tuple<const double&> & cExtTuple) const
{
  (void) v; // mark v as unused
  (void) t;
  (void) cExtTuple; // mark cExtTuple as unused
  const FDouble & caInt = get<0>(cIntTuple);
  const auto bindRate = caInt / (tauBind * (caInt + cBind));
  
  mTau = tauUnbind / (tauUnbind * bindRate + 1.0);
  mInf = bindRate * mTau;
}



// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, CAN, injectorFactory);
