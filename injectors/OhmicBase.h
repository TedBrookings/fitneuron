#ifndef OhmicBase_h
#define OhmicBase_h


#include "../include/CurrentInjector.h"


using namespace std;


// OhmicBase forms the base class for OhmicChannel and HodgkinHuxleyChannel
// It implements much of the logic of a channel with activation and/or
// inactivation gates (m and h) and dependence upon a variable number of
// SoluteTypes
// Derived classes must implement defineParameters() and computeDerivatives()

typedef size_t Power_t;

template <Power_t mPower, Power_t hPower, typename... SoluteTypes>
class OhmicBase
  :  public OnePointInjector<SoluteTypes...>    
{
  public:
    // Return a FDouble* that points to the current injected into a given
    // Solute (return NULL if none)
    virtual FDouble const* getSoluteCurrentPtr(const Solute & solute) const {
      (void)solute;
      return NULL;
    }
    // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const std::list<Trace> & traces,
                                  NeuronModel & model) = 0;
    // Define the state values of the Injector
    virtual void defineStateValues(NeuronModel & model);
    // Initialize the CurrentInjector after parameters are set
    virtual void initialize(void);
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const vector<FDouble> & states,
                                    const double & time,
                                    vector<FDouble> & dStates) const = 0;
    // Compute injected current
    virtual void computeCurrent(const vector<FDouble> & states);
    
  protected:
    size_t mIndex;
    double mQ10;
    double mTempFact;

    size_t hIndex;
    double hQ10;
    double hTempFact;
    
    double gBar;  // uS/nF
    double gBarQ10;
    double E;     // mV
    
    double celsiusBase; // C (base temperature of model)
  
  private:
    double gBar_t; // uS/nF, instantaneous gBar
};


template <Power_t mPower, typename... SoluteTypes>
class OhmicBase<mPower,0,SoluteTypes...>
  : public OnePointInjector<SoluteTypes...>
{
  public:
    // Return a FDouble* that points to the current injected into a given
    // Solute (return NULL if none)
    virtual FDouble const* getSoluteCurrentPtr(const Solute & solute) const {
      (void)solute;
      return NULL;
    }
    // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const std::list<Trace> & traces,
                                  NeuronModel & model) = 0;
    // Define the state values of the Injector
    virtual void defineStateValues(NeuronModel & model);
    // Initialize the CurrentInjector after parameters are set
    virtual void initialize(void);
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const vector<FDouble> & states,
                                    const double & time,
                                    vector<FDouble> & dStates) const = 0;
    // Compute injected current
    virtual void computeCurrent(const vector<FDouble> & states);
    
  protected:
    size_t mIndex;
    double mQ10;
    double mTempFact;

    double gBar;  // uS/nF
    double gBarQ10;
    double E;     // mV
    
    double celsiusBase; // C (base temperature of model)
  
  private:
    double gBar_t; // uS/nF, instantaneous gBar
};


template <Power_t hPower, typename... SoluteTypes>
class OhmicBase<0,hPower,SoluteTypes...>
    : public OnePointInjector<SoluteTypes...>
{
  public:
    // Return a FDouble* that points to the current injected into a given
    // Solute (return NULL if none)
    virtual FDouble const* getSoluteCurrentPtr(const Solute & solute) const {
      (void)solute;
      return NULL;
    }
    // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const std::list<Trace> & traces,
                                  NeuronModel & model) = 0;
    // Define the state values of the Injector
    virtual void defineStateValues(NeuronModel & model);
    // Initialize the CurrentInjector after parameters are set
    virtual void initialize(void);
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const vector<FDouble> & states,
                                    const double & time,
                                    vector<FDouble> & dStates) const = 0;
    // Compute injected current
    virtual void computeCurrent(const vector<FDouble> & states);
    
  protected:
    size_t hIndex;
    double hQ10;
    double hTempFact;

    double gBar;  // uS/nF
    double gBarQ10;
    double E;     // mV
    
    double celsiusBase; // C (base temperature of model)
  
  private:
    double gBar_t; // uS/nF, instantaneous gBar
};

template <typename... SoluteTypes>
class OhmicBase<0, 0, SoluteTypes...>
    : public OnePointInjector<SoluteTypes...>
{
  public:
    // Return a FDouble* that points to the current injected into a given
    // Solute (return NULL if none)
    virtual FDouble const* getSoluteCurrentPtr(const Solute & solute) const {
      (void)solute;
      return NULL;
    }
    // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const std::list<Trace> & traces,
                                  NeuronModel & model) = 0;
    // Define the state values of the Injector
    virtual void defineStateValues(NeuronModel & model);
    // Initialize the CurrentInjector after parameters are set
    virtual void initialize(void);
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const vector<FDouble> & states,
                                    const double & time,
                                    vector<FDouble> & dStates) const = 0;
    // Compute injected current
    virtual void computeCurrent(const vector<FDouble> & states);
    
  protected:        
    double gBar;  // uS/nF
    double gBarQ10;
    double E;     // mV
    
    double celsiusBase; // C (base temperature of model)
  
  private:
    double gBar_t; // uS/nF, instantaneous gBar
};



////////////////////////// OhmicBase Member Functions //////////////////////////

// defineStateValues() - general case
template <Power_t mPower, Power_t hPower, typename... SoluteTypes>
void
OhmicBase<mPower,hPower,SoluteTypes...>::defineStateValues(NeuronModel & model)
{
  // Define the state values of the Injector
  mIndex = DEFINE_STATE_VALUE(m, 0.0, 1.0, 0.0, 1.0, UNITLESS);
  hIndex = DEFINE_STATE_VALUE(h, 0.0, 1.0, 0.0, 1.0, UNITLESS);
}
// defineStateValues() - no inactivation
template <Power_t mPower, typename... SoluteTypes>
void
OhmicBase<mPower,0,SoluteTypes...>::defineStateValues(NeuronModel & model)
{
  // Define the state values of the Injector
  mIndex = DEFINE_STATE_VALUE(m, 0.0, 1.0, 0.0, 1.0, UNITLESS);
}
// defineStateValues() - no activation
template <Power_t hPower, typename... SoluteTypes>
void
OhmicBase<0,hPower,SoluteTypes...>::defineStateValues(NeuronModel & model)
{
  // Define the state values of the Injector
  hIndex = DEFINE_STATE_VALUE(h, 0.0, 1.0, 0.0, 1.0, UNITLESS);
}
// defineStateValues() - no activation, no inactivation
template <typename... SoluteTypes>
void
OhmicBase<0,0,SoluteTypes...>::defineStateValues(NeuronModel & model)
{
  // Define the state values of the Injector
  // no states to define
  (void) model;
}


// initialize() - general case
template <Power_t mPower, Power_t hPower, typename... SoluteTypes>
void
OhmicBase<mPower,hPower,SoluteTypes...>::initialize(void)
{
  // Initialize the CurrentInjector after parameters are set
  
  // Skip gBar q10 to match previous program. Reinstitute later:
  //const double gTempFact = pow(gBarQ10, 0.1 * (*celsius - celsiusBase));
  //gBar_t = gBar * gTempFact;
  gBar_t = gBar;
  
  mTempFact = pow(mQ10, 0.1 * (*this->celsius - celsiusBase));
  hTempFact = pow(hQ10, 0.1 * (*this->celsius - celsiusBase));
}
// initialize() - no inactivation
template <Power_t mPower, typename... SoluteTypes>
void
OhmicBase<mPower,0,SoluteTypes...>::initialize(void)
{
  // Initialize the CurrentInjector after parameters are set
  
  // Skip gBar q10 to match previous program. Reinstitute later:
  //const double gTempFact = pow(gBarQ10, 0.1 * (*celsius - celsiusBase));
  //gBar_t = gBar * gTempFact;
  gBar_t = gBar;
  
  mTempFact = pow(mQ10, 0.1 * (*this->celsius - celsiusBase));
}
// initialize() - no activation
template <Power_t hPower, typename... SoluteTypes>
void
OhmicBase<0,hPower,SoluteTypes...>::initialize(void)
{
  // Initialize the CurrentInjector after parameters are set
  
  // Skip gBar q10 to match previous program. Reinstitute later:
  //const double gTempFact = pow(gBarQ10, 0.1 * (*celsius - celsiusBase));
  //gBar_t = gBar * gTempFact;
  gBar_t = gBar;
  
  hTempFact = pow(hQ10, 0.1 * (*this->celsius - celsiusBase));
}
// initialize() - no inactivation or activation
template <typename... SoluteTypes>
void
OhmicBase<0,0,SoluteTypes...>::initialize(void)
{
  // Initialize the CurrentInjector after parameters are set
  
  // Skip gBar q10 to match previous program. Reinstitute later:
  //const double gTempFact = pow(gBarQ10, 0.1 * (*celsius - celsiusBase));
  //gBar_t = gBar * gTempFact;
  gBar_t = gBar;
}



// computeCurrent() - general case
template <Power_t mPower, Power_t hPower, typename... SoluteTypes>
void
OhmicBase<mPower,hPower,SoluteTypes...>
    ::computeCurrent(const vector<FDouble> & states)
{
  // Compute injected current
  const FDouble & v = states[this->compartmentPtr->voltageInd];

  const FDouble & m = states[mIndex];
  const FDouble & h = states[hIndex];
  this->i = gBar_t * power<mPower>(m) * power<hPower>(h) * (E - v);
}
// computeCurrent() - no inactivation
template <Power_t mPower, typename... SoluteTypes>
void
OhmicBase<mPower,0,SoluteTypes...>
    ::computeCurrent(const vector<FDouble> & states)
{
  // Compute injected current
  const FDouble & v = states[this->compartmentPtr->voltageInd];

  const FDouble & m = states[mIndex];
  this->i = gBar_t * power<mPower>(m) * (E - v);
}
// computeCurrent() - no activation
template <Power_t hPower, typename... SoluteTypes>
void
OhmicBase<0,hPower,SoluteTypes...>
    ::computeCurrent(const vector<FDouble> & states)
{
  // Compute injected current
  const FDouble & v = states[this->compartmentPtr->voltageInd];

  const FDouble & h = states[hIndex];  
  this->i = gBar_t * power<hPower>(h) * (E - v);
}
// computeCurrent() - no activation or inactivation
template <typename... SoluteTypes>
void
OhmicBase<0,0,SoluteTypes...>
    ::computeCurrent(const vector<FDouble> & states)
{
  // Compute injected current
  const FDouble & v = states[this->compartmentPtr->voltageInd];

  this->i = gBar_t * (E - v);
}

#endif
