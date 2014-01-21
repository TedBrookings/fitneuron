#include "../include/CurrentInjector.h"


using namespace std;


class CompartmentConnection : public TwoPointInjector<>
{
  public:
   // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const list<Trace> & traces,
                                  NeuronModel & model) {
      (void) parametersList;
      (void) traces;
      (void) model;
    }
    // Define the state values of the Injector
    virtual void defineStateValues(NeuronModel & model) { (void) model; };
    // Initialize the CurrentInjector after parameters are set
    virtual void initialize(void);
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const vector<FDouble> & states,
                                    const double & time,
                                    vector<FDouble> & dStates) const {
      (void) states;
      (void) time;
      (void) dStates;
    };
    // Compute injected current
    virtual void computeCurrent(const vector<FDouble> & states);
  private:
    double g;
    
    double capacitancePre;
    double gPerC_Pre;
    size_t vIndexPre;

    double capacitancePost;
    double gPerC_Post;
    size_t vIndexPost;
};


void
CompartmentConnection::initialize(void)
{
  // Initialize the CurrentInjector after parameters are set
  
  g = 2.0 / (preCompPtr->axialResistance + postCompPtr->axialResistance);

  capacitancePre = preCompPtr->capacitance;
  gPerC_Pre = g / capacitancePre;
  vIndexPre = preCompPtr->voltageInd;
  
  capacitancePost = postCompPtr->capacitance;
  gPerC_Post = g / capacitancePost;
  vIndexPost = postCompPtr->voltageInd;
}


void
CompartmentConnection::computeCurrent(const vector<FDouble> & states)
{
  // Compute injected current
  const FDouble & vPre = states[vIndexPre];
  const FDouble & vPost = states[vIndexPost];
  
  const auto deltaV = vPre - vPost;
  this->iPre = (-gPerC_Pre) * deltaV;
  this->iPost = gPerC_Post * deltaV;
}


// register channel so it can be loaded dynamically by injectorFactory
REGISTER_CLASS(CurrentInjector, CompartmentConnection, injectorFactory);
