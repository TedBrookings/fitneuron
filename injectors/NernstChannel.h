#ifndef NernstChannel_h
#define NernstChannel_h


#include "../include/CurrentInjector.h"


using namespace std;


typedef size_t Power_t;


// NernstChannel is a Hodgkin-Huxley-style channel with activation gate m raised
// to mPower and inactivation gate h raised to hPower (setting mPower or hPower
// to zero results in a channel with no activation/inactivation respectively).
// Unlike HodgkinHuxleyChannel, the changes in m and h are specified by
// steady-state value and time constant (e.g. mInf and mTau) rather than
// transition rates. It depends on at least one Solute.

// TO USE THIS CLASS:
// 1. Declare a derived class from the appropriate OhmicChannel class
// 2. Define the public function:
//    virtual void defineParameters(ParameterDescriptionList & parametersList,
//                                  const std::list<Trace> & traces,
//                                  NeuronModel & model);
// 3. If there are no SoluteType dependencies, declare the protected function:
//      virtual void dynamics(const FDouble & v, const double & t);
//    Otherwise declare the protected function:
//      virtual void dynamics(const FDouble & v, const double & t,
//                    HomogeneousTuple<numSolutes, const FDouble&> & cIntTuple,
//                    HomogeneousTuple<numSolutes, const double&> & cExtTuple);
//   (the tuples deliver internal and external Solute concentrations.

template <Power_t mPower, Power_t hPower, typename Sol, typename... SoluteTypes>
class NernstChannel
    : public OnePointInjector<Sol,SoluteTypes...>    
{
  public:
    // Return a FDouble* that points to the current injected into a given
    // Solute (return NULL if none)
    virtual FDouble const* getSoluteCurrentPtr(const Solute & solute) const {
      if(solute.name == (string)Sol::name &&
         solute.compartmentPtr == this->compartmentPtr) {
        return &this->i;
      }
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
                                    vector<FDouble> & dStates) const;
    // Compute injected current
    virtual void computeCurrent(const vector<FDouble> & states);
    
  protected:
    //         MANDATORY OVERLOAD
    // model activation and inactivation dynamics
    //   (assign mInf, mTau, hInf, and hTau)
    using IntTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const FDouble&>::Type;
    using ExtTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const double&>::Type;
    
    virtual void dynamics(const FDouble & v, const double & t,
                          IntTuple & cIntTuple, ExtTuple & cExtTuple) const = 0;

    size_t mIndex;
    mutable FDouble mInf;  // used to pass dynamics info back from derived class
    mutable FDouble mTau;  // used to pass dynamics info back from derived class
    double mQ10;
    double mTempFact;

    size_t hIndex;
    mutable FDouble hInf;  // used to pass dynamics info back from derived class
    mutable FDouble hTau;  // used to pass dynamics info back from derived class
    double hQ10;
    double hTempFact;
    
    double gBar;  // uS/nF
    double gBarQ10;
    
    double celsiusBase; // C (base temperature of model)
  
  private:
    double RTOverZF;    // mV
    double gBar_t; // uS/nF, instantaneous gBar
};



template <Power_t mPower, typename Sol, typename... SoluteTypes>
class NernstChannel<mPower, 0, Sol,SoluteTypes...>
    : public OnePointInjector<Sol,SoluteTypes...>
{
  public:
    // Return a FDouble* that points to the current injected into a given
    // Solute (return NULL if none)
    virtual FDouble const* getSoluteCurrentPtr(const Solute & solute) const {
      if(solute.name == (string)Sol::name &&
         solute.compartmentPtr == this->compartmentPtr) {
        return &this->i;
      }
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
                                    vector<FDouble> & dStates) const;
    // Compute injected current
    virtual void computeCurrent(const vector<FDouble> & states);
    
  protected:
    // model activation and inactivation dynamics
    //   (assign mInf, mTau, hInf, and hTau)
    using IntTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const FDouble&>::Type;
    using ExtTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const double&>::Type;
    
    virtual void dynamics(const FDouble & v, const double & t,
                          IntTuple & cIntTuple, ExtTuple & cExtTuple) const = 0;

    size_t mIndex;
    mutable FDouble mInf;  // used to pass dynamics info back from derived class
    mutable FDouble mTau;  // used to pass dynamics info back from derived class
    double mQ10;
    double mTempFact;

    double gBar;  // uS/nF
    double gBarQ10;
    
    double celsiusBase; // C (base temperature of model)
  
  private:
    double RTOverZF;    // mV
    double gBar_t; // uS/nF, instantaneous gBar
};


template <Power_t hPower, typename Sol, typename... SoluteTypes>
class NernstChannel<0, hPower, Sol,SoluteTypes...>
    : public OnePointInjector<Sol,SoluteTypes...>
{
  public:
    // Return a FDouble* that points to the current injected into a given
    // Solute (return NULL if none)
    virtual FDouble const* getSoluteCurrentPtr(const Solute & solute) const {
      if(solute.name == (string)Sol::name &&
         solute.compartmentPtr == this->compartmentPtr) {
        return &this->i;
      }
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
                                    vector<FDouble> & dStates) const;
    // Compute injected current
    virtual void computeCurrent(const vector<FDouble> & states);
    
  protected:
    //         MANDATORY OVERLOAD
    // model activation and inactivation dynamics
    //   (assign mInf, mTau, hInf, and hTau)
    using IntTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const FDouble&>::Type;
    using ExtTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const double&>::Type;
    
    virtual void dynamics(const FDouble & v, const double & t,
                          IntTuple & cIntTuple, ExtTuple & cExtTuple) const = 0;

    size_t hIndex;
    mutable FDouble hInf;  // used to pass dynamics info back from derived class
    mutable FDouble hTau;  // used to pass dynamics info back from derived class
    double hQ10;
    double hTempFact;

    double gBar;  // uS/nF
    double gBarQ10;
    
    double celsiusBase; // C (base temperature of model)
  
  private:
    double RTOverZF;    // mV
    double gBar_t; // uS/nF, instantaneous gBar
};


template <typename Sol, typename... SoluteTypes>
class NernstChannel<0, 0, Sol,SoluteTypes...>
    : public OnePointInjector<Sol,SoluteTypes...>
{
  public:
    // Return a FDouble* that points to the current injected into a given
    // Solute (return NULL if none)
    virtual FDouble const* getSoluteCurrentPtr(const Solute & solute) const {
      if(solute.name == (string)Sol::name &&
         solute.compartmentPtr == this->compartmentPtr) {
        return &this->i;
      }
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
                                    vector<FDouble> & dStates) const;
    // Compute injected current
    virtual void computeCurrent(const vector<FDouble> & states);
    
  protected:
    //         MANDATORY OVERLOAD
    // model activation and inactivation dynamics
    //   (assign mInf, mTau, hInf, and hTau)
    using IntTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const FDouble&>::Type;
    using ExtTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const double&>::Type;
        
    double gBar;  // uS/nF
    double gBarQ10;
    
    double celsiusBase; // C (base temperature of model)
  
  private:
    double RTOverZF;    // mV
    double gBar_t; // uS/nF, instantaneous gBar
};

//////////////////////// NernstChannel Member Functions ////////////////////////

// defineStateValues() - general case
template <Power_t mPower, Power_t hPower, typename Sol, typename... SoluteTypes>
void
NernstChannel<mPower,hPower,Sol,SoluteTypes...>
    ::defineStateValues(NeuronModel & model)
{
  // Define the state values of the Injector
  mIndex = DEFINE_STATE_VALUE(m, 0.0, 1.0, 0.0, 1.0, UNITLESS);
  hIndex = DEFINE_STATE_VALUE(h, 0.0, 1.0, 0.0, 1.0, UNITLESS);
}
template <Power_t mPower, typename Sol, typename... SoluteTypes>
void
NernstChannel<mPower,0,Sol,SoluteTypes...>
    ::defineStateValues(NeuronModel & model)
{
  // Define the state values of the Injector
  mIndex = DEFINE_STATE_VALUE(m, 0.0, 1.0, 0.0, 1.0, UNITLESS);
}
template <Power_t hPower, typename Sol, typename... SoluteTypes>
void
NernstChannel<0,hPower,Sol,SoluteTypes...>
    ::defineStateValues(NeuronModel & model)
{
  // Define the state values of the Injector
  hIndex = DEFINE_STATE_VALUE(h, 0.0, 1.0, 0.0, 1.0, UNITLESS);
}
template <typename Sol, typename... SoluteTypes>
void
NernstChannel<0,0,Sol,SoluteTypes...>
    ::defineStateValues(NeuronModel & model)
{
  // Define the state values of the Injector
}


// initialize() - general case
template <Power_t mPower, Power_t hPower, typename Sol, typename... SoluteTypes>
void
NernstChannel<mPower,hPower,Sol,SoluteTypes...>
    ::initialize(void)
{
  // Initialize the CurrentInjector after parameters are set
  
  // Skip gBar q10 to match previous program. Reinstitute later:
  //const double gTempFact = pow(gBarQ10, 0.1 * (*celsius - celsiusBase));
  //gBar_t = gBar * gTempFact;
  gBar_t = gBar;
  
  mTempFact = pow(mQ10, 0.1 * (*this->celsius - celsiusBase));
  hTempFact = pow(hQ10, 0.1 * (*this->celsius - celsiusBase));
  
  RTOverZF = 1000.0 * K_Boltzman *
    (*this->celsius + CelsiusToKelvin) / (Q_Electron * (double)Sol::Z);
}
// initialize() - no inactivation
template <Power_t mPower, typename Sol, typename... SoluteTypes>
void
NernstChannel<mPower,0,Sol,SoluteTypes...>
    ::initialize(void)
{
  // Initialize the CurrentInjector after parameters are set
  
  // Skip gBar q10 to match previous program. Reinstitute later:
  //const double gTempFact = pow(gBarQ10, 0.1 * (*celsius - celsiusBase));
  //gBar_t = gBar * gTempFact;
  gBar_t = gBar;
  
  mTempFact = pow(mQ10, 0.1 * (*this->celsius - celsiusBase));
  
  RTOverZF = 1000.0 * K_Boltzman *
    (*this->celsius + CelsiusToKelvin) / (Q_Electron * (double)Sol::Z);
}
// initialize() - no activation
template <Power_t hPower, typename Sol, typename... SoluteTypes>
void
NernstChannel<0,hPower,Sol,SoluteTypes...>
    ::initialize(void)
{
  // Initialize the CurrentInjector after parameters are set
  
  // Skip gBar q10 to match previous program. Reinstitute later:
  //const double gTempFact = pow(gBarQ10, 0.1 * (*celsius - celsiusBase));
  //gBar_t = gBar * gTempFact;
  gBar_t = gBar;
  
  hTempFact = pow(hQ10, 0.1 * (*this->celsius - celsiusBase));
  
  RTOverZF = 1000.0 * K_Boltzman *
    (*this->celsius + CelsiusToKelvin) / (Q_Electron * (double)Sol::Z);
}
// initialize() - no inactivation or activation
template <typename Sol, typename... SoluteTypes>
void
NernstChannel<0,0,Sol,SoluteTypes...>
    ::initialize(void)
{
  // Initialize the CurrentInjector after parameters are set
  
  // Skip gBar q10 to match previous program. Reinstitute later:
  //const double gTempFact = pow(gBarQ10, 0.1 * (*celsius - celsiusBase));
  //gBar_t = gBar * gTempFact;
  gBar_t = gBar;
  
  RTOverZF = 1000.0 * K_Boltzman *
    (*this->celsius + CelsiusToKelvin) / (Q_Electron * (double)Sol::Z);
}



// computeDerivatives() - general case
template <Power_t mPower, Power_t hPower, typename Sol, typename... SoluteTypes>
void
NernstChannel<mPower,hPower,Sol,SoluteTypes...>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  auto cIntTup = OnePointInjector<Sol,SoluteTypes...>::getCInt(states);
  auto cExtTup = OnePointInjector<Sol,SoluteTypes...>::getCExt();
  
  dynamics(v, time, cIntTup, cExtTup);

  dStates[mIndex] = mTempFact * (mInf - states[mIndex]) / mTau;
  dStates[hIndex] = hTempFact * (hInf - states[hIndex]) / hTau;
}
// computeDerivatives() - no inactivation
template <Power_t mPower, typename Sol, typename... SoluteTypes>
void
NernstChannel<mPower,0,Sol,SoluteTypes...>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  auto cIntTup = OnePointInjector<Sol,SoluteTypes...>::getCInt(states);
  auto cExtTup = OnePointInjector<Sol,SoluteTypes...>::getCExt();
  
  dynamics(v, time, cIntTup, cExtTup);

  dStates[mIndex] = mTempFact * (mInf - states[mIndex]) / mTau;
}
// computeDerivatives() - no activation
template <Power_t hPower, typename Sol, typename... SoluteTypes>
void
NernstChannel<0,hPower,Sol,SoluteTypes...>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  auto cIntTup = OnePointInjector<Sol,SoluteTypes...>::getCInt(states);
  auto cExtTup = OnePointInjector<Sol,SoluteTypes...>::getCExt();
  
  dynamics(v, time, cIntTup, cExtTup);

  dStates[hIndex] = hTempFact * (hInf - states[hIndex]) / hTau;
}
// computeDerivatives() - no activation or inactivation
template <typename Sol, typename... SoluteTypes>
void
NernstChannel<0,0,Sol,SoluteTypes...>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // no internal states!
}


// computeCurrent() - general case
template <Power_t mPower, Power_t hPower, typename Sol, typename... SoluteTypes>
void
NernstChannel<mPower,hPower,Sol,SoluteTypes...>
    ::computeCurrent(const vector<FDouble> & states)
{
  // Compute injected current
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  const FDouble & cInt = states[this->solutePtr->concentrationInd];
  const double & cExt = this->solutePtr->cExt;
  const auto E = RTOverZF * log(cExt / cInt);

  const FDouble & m = states[mIndex];
  const FDouble & h = states[hIndex];  
  this->i = gBar_t * power<mPower>(m) * power<hPower>(h) * (E - v);
}
// computeCurrent() - no inactivation
template <Power_t mPower, typename Sol, typename... SoluteTypes>
void
NernstChannel<mPower,0,Sol,SoluteTypes...>
    ::computeCurrent(const vector<FDouble> & states)
{
  // Compute injected current
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  const FDouble & cInt = states[this->solutePtr->concentrationInd];
  const double & cExt = this->solutePtr->cExt;
  const auto E = RTOverZF * log(cExt / cInt);

  const FDouble & m = states[mIndex];
  this->i = gBar_t * power<mPower>(m) * (E - v);
  /*
  std::cout << this->i.getValue() << " = " << gBar_t << " * "
            << power<mPower>(m).getValue() << " * (" << E.getValue()
            << " - " << v.getValue() << ")\n";
  std::cout.flush();
  */
}
// computeCurrent() - no activation
template <Power_t hPower, typename Sol, typename... SoluteTypes>
void
NernstChannel<0,hPower,Sol,SoluteTypes...>
    ::computeCurrent(const vector<FDouble> & states)
{
  // Compute injected current
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  const FDouble & cInt = states[this->solutePtr->concentrationInd];
  const double & cExt = this->solutePtr->cExt;
  const auto E = RTOverZF * log(cExt / cInt);

  const FDouble & h = states[hIndex];  
  this->i = gBar_t * power<hPower>(h) * (E - v);
}
// computeCurrent() - no activation or inactivation
template <typename Sol, typename... SoluteTypes>
void
NernstChannel<0,0,Sol,SoluteTypes...>
    ::computeCurrent(const vector<FDouble> & states)
{
  // Compute injected current
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  const FDouble & cInt = states[this->solutePtr->concentrationInd];
  const double & cExt = this->solutePtr->cExt;
  const auto E = RTOverZF * log(cExt / cInt);

  this->i = gBar_t * (E - v);
}

#endif
