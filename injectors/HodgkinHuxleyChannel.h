#ifndef HodgkinHuxleyChannel_h
#define HodgkinHuxleyChannel_h


#include "OhmicBase.h"


// HodgkinHuxleyChannel is a Hodgkin-Huxley-style channel with activation gate m
// raised to mPower and inactivation gate h raised to hPower (setting mPower or
// hPower to zero results in a channel with no activation/inactivation
// respectively). The changes in m and h are specified by transition rates (e.g.
// mAlpha and mBeta). It may depend on one or more Solutes (or it may not).

// TO USE THIS CLASS:
// 1. Declare a derived class from the appropriate HodgkinHuxleyChannel class
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
class HodgkinHuxleyChannel;


// activation and inactivation, no SoluteTypes
template <Power_t mPower, Power_t hPower>
class HodgkinHuxleyChannel<mPower,hPower> : public OhmicBase<mPower,hPower>    
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
    //   (assign mAlpha, mBeta, hAlpha, and hBeta)
    virtual void dynamics(const FDouble & v, const double & t) const = 0;
    
    mutable FDouble mAlpha;  // used to pass dynamics info back from derived class
    mutable FDouble mBeta;  // used to pass dynamics info back from derived class

    mutable FDouble hAlpha;  // used to pass dynamics info back from derived class
    mutable FDouble hBeta;  // used to pass dynamics info back from derived class
};

// activation and inactivation, one or more SoluteTypes
template <Power_t mPower, Power_t hPower, typename Sol, typename... SoluteTypes>
class HodgkinHuxleyChannel<mPower,hPower,Sol,SoluteTypes...>
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
    //   (assign mAlpha, mBeta, hAlpha, and hBeta)
    using IntTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const FDouble&>::Type;
    using ExtTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const double&>::Type;
    virtual void dynamics(const FDouble & v, const double & t,
                          IntTuple & cIntTuple, ExtTuple & cExtTuple) const = 0;

    mutable FDouble mAlpha;  // used to pass dynamics info back from derived class
    mutable FDouble mBeta;  // used to pass dynamics info back from derived class

    mutable FDouble hAlpha;  // used to pass dynamics info back from derived class
    mutable FDouble hBeta;  // used to pass dynamics info back from derived class
};

// activation, no inactivation, no SoluteTypes
template <Power_t mPower>
class HodgkinHuxleyChannel<mPower,0> : public OhmicBase<mPower,0>    
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
    //   (assign mAlpha, mBeta, hAlpha, and hBeta)
    virtual void dynamics(const FDouble & v, const double & t) const = 0;
    
    mutable FDouble mAlpha;  // used to pass dynamics info back from derived class
    mutable FDouble mBeta;  // used to pass dynamics info back from derived class
};

// activation, no inactivation, one or more SoluteTypes
template <Power_t mPower, typename Sol, typename... SoluteTypes>
class HodgkinHuxleyChannel<mPower,0,Sol,SoluteTypes...>
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
    //   (assign mAlpha, mBeta, hAlpha, and hBeta)
    using IntTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const FDouble&>::Type;
    using ExtTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const double&>::Type;
    virtual void dynamics(const FDouble & v, const double & t,
                          IntTuple & cIntTuple, ExtTuple & cExtTuple) const = 0;

    mutable FDouble mAlpha;  // used to pass dynamics info back from derived class
    mutable FDouble mBeta;  // used to pass dynamics info back from derived class
};

// inactivation, no activation, no SoluteTypes
template <Power_t hPower>
class HodgkinHuxleyChannel<0,hPower> : public OhmicBase<0,hPower>    
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
    //   (assign mAlpha, mBeta, hAlpha, and hBeta)
    virtual void dynamics(const FDouble & v, const double & t) const = 0;
    
    mutable FDouble hAlpha;  // used to pass dynamics info back from derived class
    mutable FDouble hBeta;  // used to pass dynamics info back from derived class
};

// inactivation, no activation, one or more SoluteTypes
template <Power_t hPower, typename Sol, typename... SoluteTypes>
class HodgkinHuxleyChannel<0,hPower,Sol,SoluteTypes...>
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
    //   (assign mAlpha, mBeta, hAlpha, and hBeta)
    using IntTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const FDouble&>::Type;
    using ExtTuple = typename HomogeneousTuple<1+sizeof...(SoluteTypes),
                                      const double&>::Type;
    virtual void dynamics(const FDouble & v, const double & t,
                          IntTuple & cIntTuple, ExtTuple & cExtTuple) const = 0;

    mutable FDouble hAlpha;  // used to pass dynamics info back from derived class
    mutable FDouble hBeta;  // used to pass dynamics info back from derived class
};

// no activation or inactivation, no SoluteTypes
template <>
class HodgkinHuxleyChannel<0,0>
  : public OhmicBase<0,0>    
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
class HodgkinHuxleyChannel<0,0,Sol,SoluteTypes...>
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


///////////////////////// HodgkinHuxleyChannel Member Functions ////////////////////////


// computeDerivatives() - general case, no SoluteTypes
template <Power_t mPower, Power_t hPower>
void
HodgkinHuxleyChannel<mPower,hPower>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  
  dynamics(v, time);

  dStates[this->mIndex] =
    this->mTempFact * (mAlpha - (mAlpha + mBeta) * states[this->mIndex]);
  dStates[this->hIndex] =
    this->hTempFact * (hAlpha - (hAlpha + hBeta) * states[this->hIndex]);
}
// computeDerivatives() - general case, one or more SoluteTypes
template <Power_t mPower, Power_t hPower, typename Sol, typename... SoluteTypes>
void
HodgkinHuxleyChannel<mPower,hPower,Sol,SoluteTypes...>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  auto cIntTup = OnePointInjector<Sol,SoluteTypes...>::getCInt(states);
  auto cExtTup = OnePointInjector<Sol,SoluteTypes...>::getCExt();
  
  dynamics(v, time, cIntTup, cExtTup);

  dStates[this->mIndex] =
    this->mTempFact * (mAlpha - (mAlpha + mBeta) * states[this->mIndex]);
  dStates[this->hIndex] =
    this->hTempFact * (hAlpha - (hAlpha + hBeta) * states[this->hIndex]);
}
// computeDerivatives() - no inactivation, no SoluteTypes
template <Power_t mPower>
void
HodgkinHuxleyChannel<mPower,0>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  
  dynamics(v, time);

  dStates[this->mIndex] =
    this->mTempFact * (mAlpha - (mAlpha + mBeta) * states[this->mIndex]);
}
// computeDerivatives() - no inactivation, one or more SoluteTypes
template <Power_t mPower, typename Sol, typename... SoluteTypes>
void
HodgkinHuxleyChannel<mPower,0,Sol,SoluteTypes...>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  auto cIntTup = OnePointInjector<Sol,SoluteTypes...>::getCInt(states);
  auto cExtTup = OnePointInjector<Sol,SoluteTypes...>::getCExt();
  
  dynamics(v, time, cIntTup, cExtTup);

  dStates[this->mIndex] =
    this->mTempFact * (mAlpha - (mAlpha + mBeta) * states[this->mIndex]);
}
// computeDerivatives() - no activation, no SoluteTypes
template <Power_t hPower>
void
HodgkinHuxleyChannel<0,hPower>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  dynamics(v, time);

  dStates[this->hIndex] =
    this->hTempFact * (hAlpha - (hAlpha + hBeta) * states[this->hIndex]);
}
// computeDerivatives() - no activation, one or more SoluteTypes
template <Power_t hPower, typename Sol, typename... SoluteTypes>
void
HodgkinHuxleyChannel<0,hPower,Sol,SoluteTypes...>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  const FDouble & v = states[this->compartmentPtr->voltageInd];
  auto cIntTup = OnePointInjector<Sol,SoluteTypes...>::getCInt(states);
  auto cExtTup = OnePointInjector<Sol,SoluteTypes...>::getCExt();
  
  dynamics(v, time, cIntTup, cExtTup);

  dStates[this->hIndex] =
    this->hTempFact * (hAlpha - (hAlpha + hBeta) * states[this->hIndex]);
}
// computeDerivatives() - no activation or inactivation, no SoluteTypes
void
HodgkinHuxleyChannel<0,0>
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
HodgkinHuxleyChannel<0,0,Sol,SoluteTypes...>
    ::computeDerivatives(const vector<FDouble> & states, const double & time,
                         vector<FDouble> & dStates) const
{
  // no internal states!
  (void)states;
  (void)time;
  (void)dStates;
}

#endif
