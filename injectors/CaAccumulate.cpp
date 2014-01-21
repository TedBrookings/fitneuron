#include "Accumulator.h"



DECLARE_SOLUTE(Ca, 2);
class CaAccumulate : public Accumulator<Ca>
{
  public:
    // Define the parameters of the accumulator
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const list<Trace> & traces,
                                  NeuronModel & model);
};



void
CaAccumulate::defineParameters(ParameterDescriptionList & parametersList,
                               const list<Trace> & traces, NeuronModel & model)
{
  celsiusBase = 18.0;                      // (C)
  
  DEFINE_PARAMETER( cEquilibrium, 0.5e-3, "mM");
  DEFINE_PARAMETER( tau,          300.0,  "ms");  // PD: 300.0, AB: 303.0
  // domainDepth based on desire to match coefficients from Soto-Trevino et al
  //   2005:
  //  PD iCoef = 0.515e-3 mM/nA, C = 9 nF
  //  AB iCoef = 0.418e-3 mM/nA, C = 12 nF
  DEFINE_PARAMETER( domainDepth,  2.520,  "um");  // PD: 2.520, AB: 4.140
  DEFINE_PARAMETER( q10,          2.0,    UNITLESS);    // (wild guess)
}



// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, CaAccumulate, injectorFactory);
