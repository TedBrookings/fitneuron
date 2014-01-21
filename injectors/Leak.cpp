#include "OhmicChannel.h"


class Leak : public OhmicChannel<0,0>
{
  public:
    // Define the parameters of the ion channel
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const list<Trace> & traces,
                                  NeuronModel & model);
};



void
Leak::defineParameters(ParameterDescriptionList & parametersList,
                       const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 11.0; // (C)
  DEFINE_PARAMETER( gBarQ10, 2.0,  UNITLESS );  // wild guess

  // model parameters, from Connor-Stevens model
  DEFINE_PARAMETER( gBar,    0.3,  "uS/nF" );
  DEFINE_PARAMETER( E,      -17.0, "mV" );
}


// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, Leak, injectorFactory);
