#ifndef OhmicChannel_h
#define OhmicChannel_h


#include "OhmicBase.h"


// OhmicChannel is a Hodgkin-Huxley-style channel with activation gate m raised
// to mPower and inactivation gate h raised to hPower (setting mPower or hPower
// to zero results in a channel with no activation/inactivation respectively).
// Unlike HodgkinHuxleyChannel, the changes in m and h are specified by
// steady-state value and time constant (e.g. mInf and mTau) rather than
// transition rates. It may depend on one or more Solutes (or it may not).

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
//               HomogeneousTuple<numSolutes, const FDouble&> & cIntTuple,
//               HomogeneousTuple<numSolutes, const double&> & cExtTuple) const;
//   (the tuples deliver internal and external Solute concentrations.

template <Power_t mPower, Power_t hPower, typename... SoluteTypes>
class OhmicChannel;


// activation and inactivation, no SoluteTypes
template <Power_t mPower, Power_t hPower>
class OhmicChannel<mPower,hPower> : public OhmicBase<mPower,hPower>    
{
  public:
    // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const std::list<Trace> & traces,
                                  NeuronModel & model) = 0;
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const vector<FDouble> & states,
                                    const double & time,
                                    vector<FDouble> & dStates) const;
  protected:
    // model activation and inactivation dynamics
    //   (assign mInf, mTau, hInf, and hTau)
    virtual void dynamics(const FDouble & v, const double & t) const = 0;
    
    mutable FDouble mInf;  // used to pass dynamics info back from derived class
    mutable FDouble mTau;  // used to pass dynamics info back from derived class

    mutable FDouble hInf;  // used to pass dynamics info back from derived class
    mutable FDouble hTau;  // used to pass dynamics info back from derived class
};

// activation and inactivation, one or more SoluteTypes
template <Power_t mPower, Power_t hPower, typename Sol, typename... SoluteTypes>
class OhmicChannel<mPower,hPower,Sol,SoluteTypes...>
  : public OhmicBase<mPower,hPower,Sol,SoluteTypes...>    
{
  public:
    // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const std::list<Trace> & traces,
                                  NeuronModel & model) = 0;
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const vector<FDouble> & states,
                                    const double & time,
                                    vector<FDouble> & dStates) const;
  protected:
    // model activation and inactivation dynamics
    //   (assign mInf, mTau, hInf, and hTau)
    using IntTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const FDouble&>::Type;
    using ExtTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const double&>::Type;
    virtual void dynamics(const FDouble & v, const double & t,
                          IntTuple & cIntTuple, ExtTuple & cExtTuple) const = 0;

    mutable FDouble mInf;  // used to pass dynamics info back from derived class
    mutable FDouble mTau;  // used to pass dynamics info back from derived class

    mutable FDouble hInf;  // used to pass dynamics info back from derived class
    mutable FDouble hTau;  // used to pass dynamics info back from derived class
};

// activation, no inactivation, no SoluteTypes
template <Power_t mPower>
class OhmicChannel<mPower,0> : public OhmicBase<mPower,0>    
{
  public:
    // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const std::list<Trace> & traces,
                                  NeuronModel & model) = 0;
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const vector<FDouble> & states,
                                    const double & time,
                                    vector<FDouble> & dStates) const;
  protected:
    // model activation and inactivation dynamics
    //   (assign mInf, mTau, hInf, and hTau)
    virtual void dynamics(const FDouble & v, const double & t) const = 0;
    
    mutable FDouble mInf;  // used to pass dynamics info back from derived class
    mutable FDouble mTau;  // used to pass dynamics info back from derived class
};

// activation, no inactivation, one or more SoluteTypes
template <Power_t mPower, typename Sol, typename... SoluteTypes>
class OhmicChannel<mPower,0,Sol,SoluteTypes...>
  : public OhmicBase<mPower,0,Sol,SoluteTypes...>    
{
  public:
    // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const std::list<Trace> & traces,
                                  NeuronModel & model) = 0;
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const vector<FDouble> & states,
                                    const double & time,
                                    vector<FDouble> & dStates) const;
  protected:
    // model activation and inactivation dynamics
    //   (assign mInf, mTau, hInf, and hTau)
    using IntTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const FDouble&>::Type;
    using ExtTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const double&>::Type;
    virtual void dynamics(const FDouble & v, const double & t,
                          IntTuple & cIntTuple, ExtTuple & cExtTuple) const = 0;

    mutable FDouble mInf;  // used to pass dynamics info back from derived class
    mutable FDouble mTau;  // used to pass dynamics info back from derived class
};

// inactivation, no activation, no SoluteTypes
template <Power_t hPower>
class OhmicChannel<0,hPower> : public OhmicBase<0,hPower>    
{
  public:
    // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const std::list<Trace> & traces,
                                  NeuronModel & model) = 0;
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const vector<FDouble> & states,
                                    const double & time,
                                    vector<FDouble> & dStates) const;
  protected:
    // model activation and inactivation dynamics
    //   (assign mInf, mTau, hInf, and hTau)
    virtual void dynamics(const FDouble & v, const double & t) const = 0;
    
    mutable FDouble hInf;  // used to pass dynamics info back from derived class
    mutable FDouble hTau;  // used to pass dynamics info back from derived class
};

// inactivation, no activation, one or more SoluteTypes
template <Power_t hPower, typename Sol, typename... SoluteTypes>
class OhmicChannel<0,hPower,Sol,SoluteTypes...>
  : public OhmicBase<0,hPower,Sol,SoluteTypes...>    
{
  public:
    // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const std::list<Trace> & traces,
                                  NeuronModel & model) = 0;
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const vector<FDouble> & states,
                                    const double & time,
                                    vector<FDouble> & dStates) const;
  protected:
    // model activation and inactivation dynamics
    //   (assign mInf, mTau, hInf, and hTau)
    using IntTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const FDouble&>::Type;
    using ExtTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const double&>::Type;
    virtual void dynamics(const FDouble & v, const double & t,
                          IntTuple & cIntTuple, ExtTuple & cExtTuple) const = 0;

    mutable FDouble hInf;  // used to pass dynamics info back from derived class
    mutable FDouble hTau;  // used to pass dynamics info back from derived class
};

// no activation or inactivation, no SoluteTypes
template <>
class OhmicChannel<0,0> : public OhmicBase<0,0>    
{
  public:
    // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const std::list<Trace> & traces,
                                  NeuronModel & model) = 0;
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const vector<FDouble> & states,
                                    const double & time,
                                    vector<FDouble> & dStates) const;
};
// no activation or inactivation, one or more SoluteTypes
template <typename Sol, typename... SoluteTypes>
class OhmicChannel<0,0,Sol,SoluteTypes...>
  : public OhmicBase<0,0,Sol,SoluteTypes...>    
{
  public:
    // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const std::list<Trace> & traces,
                                  NeuronModel & model) = 0;
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const vector<FDouble> & states,
                                    const double & time,
                                    vector<FDouble> & dStates) const;
};


///////////////////////// OhmicChannel Member Functions ////////////////////////

// computeDerivatives() - general case, no SoluteTypes
template <Power_t mPower, Power_t hPower>
void
OhmicChannel<mPower,hPower>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  
  dynamics(v, time);
  
  dStates[this->mIndex] =
    this->mTempFact * (mInf - states[this->mIndex]) / mTau;
  dStates[this->hIndex] =
    this->hTempFact * (hInf - states[this->hIndex]) / hTau;
}
// computeDerivatives() - general case, one or more SoluteTypes
template <Power_t mPower, Power_t hPower, typename Sol, typename... SoluteTypes>
void
OhmicChannel<mPower,hPower,Sol,SoluteTypes...>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  auto cIntTup = OnePointInjector<Sol,SoluteTypes...>::getCInt(states);
  auto cExtTup = OnePointInjector<Sol,SoluteTypes...>::getCExt();
  
  dynamics(v, time, cIntTup, cExtTup);

  dStates[this->mIndex] =
    this->mTempFact * (mInf - states[this->mIndex]) / mTau;
  dStates[this->hIndex] =
    this->hTempFact * (hInf - states[this->hIndex]) / hTau;
}
// computeDerivatives() - no inactivation, no SoluteTypes
template <Power_t mPower>
void
OhmicChannel<mPower,0>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];

  dynamics(v, time);

  dStates[this->mIndex] =
    this->mTempFact * (mInf - states[this->mIndex]) / mTau;
}
// computeDerivatives() - no inactivation, one or more SoluteTypes
template <Power_t mPower, typename Sol, typename... SoluteTypes>
void
OhmicChannel<mPower,0,Sol,SoluteTypes...>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  auto cIntTup = OnePointInjector<Sol,SoluteTypes...>::getCInt(states);
  auto cExtTup = OnePointInjector<Sol,SoluteTypes...>::getCExt();
  
  dynamics(v, time, cIntTup, cExtTup);

  dStates[this->mIndex] =
    this->mTempFact * (mInf - states[this->mIndex]) / mTau;
}
// computeDerivatives() - no activation, no SoluteTypes
template <Power_t hPower>
void
OhmicChannel<0,hPower>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  dynamics(v, time);

  dStates[this->hIndex] =
    this->hTempFact * (hInf - states[this->hIndex]) / hTau;
}
// computeDerivatives() - no activation, one or more SoluteTypes
template <Power_t hPower, typename Sol, typename... SoluteTypes>
void
OhmicChannel<0,hPower,Sol,SoluteTypes...>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  auto cIntTup = OnePointInjector<Sol,SoluteTypes...>::getCInt(states);
  auto cExtTup = OnePointInjector<Sol,SoluteTypes...>::getCExt();
  
  dynamics(v, time, cIntTup, cExtTup);

  dStates[this->hIndex] =
    this->hTempFact * (hInf - states[this->hIndex]) / hTau;
}
// computeDerivatives() - no activation or inactivation, no SoluteTypes
void
OhmicChannel<0,0>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // no internal states!
  (void)states;
  (void)time;
  (void)dStates;
}
// computeDerivatives() - no activation or inactivation, one or more SoluteTypes
template <typename Sol, typename... SoluteTypes>
void
OhmicChannel<0,0,Sol,SoluteTypes...>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // no internal states!
  (void)states;
  (void)time;
  (void)dStates;
}
#endif
