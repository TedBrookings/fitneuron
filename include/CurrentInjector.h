#ifndef CurrentInjector_h
#define CurrentInjector_h

#include <tuple>

#include "ParameterDescription.h"
#include "dynamic_class_loader.h"
#include "NeuronModel.h"
#include "NeuronException.h"



// define a macro for defining CurrentInjector state values
#define DEFINE_STATE_VALUE(parameter, minAllowed, maxAllowed, minRange, \
                           maxRange, units) \
this->defineStateValue( \
  ParameterDescription (#parameter + this->getNameSuffix(), units, minAllowed, \
                        maxAllowed, minRange, maxRange), \
  model )

// define a macro for defining CurrentInjector model parameters
#define DEFINE_PARAMETER(parameter, defaultValue, units) \
defineParameter(parameter, \
                ParameterDescription (#parameter, units, defaultValue), \
                parametersList, traces, model)

// define a macro for declaring a SoluteType that will be used by a
// OnePointInjector (e.g. ion channel)
#define DECLARE_SOLUTE(soluteName, soluteZ) \
char soluteName##Name [] = #soluteName; \
using soluteName = SoluteType<soluteName##Name, soluteZ, SoluteTypeLocation::OnePoint>

// define a macro for declaring a SoluteType that will be used on the "pre"
// side by a TwoPointInjector (e.g. synapse)
#define DECLARE_SOLUTE_PRE(soluteName, soluteZ) \
char soluteName##Name [] = #soluteName; \
using soluteName = SoluteType<soluteName##Name, soluteZ, SoluteTypeLocation::Pre>

// define a macro for declaring a SoluteType that will be used on the "post"
// side by a TwoPointInjector (e.g. synapse)
#define DECLARE_SOLUTE_POST(soluteName, soluteZ) \
char soluteName##Name [] = #soluteName; \
using soluteName = SoluteType<soluteName##Name, soluteZ, SoluteTypeLocation::Post>

// define a macro for unitless quantities
#define UNITLESS "1"



enum SoluteTypeLocation { OnePoint, Pre, Post};
template <char const* solName, int solZ, SoluteTypeLocation solLocation>
struct SoluteType
{
    constexpr static char const*const name = solName;
    constexpr static const int Z = solZ;
    constexpr static SoluteTypeLocation location = solLocation;
};
// declare a tuple with a N elements of a given type. e.g.
// HomogeneousTuple<13, double>::Type doubleTup; //declare a 13-tuple of doubles
template<int N, typename T, typename... TArgs>
struct HomogeneousTuple : HomogeneousTuple<N-1, T, T, TArgs...> {};
template<typename T, typename... TArgs>
struct HomogeneousTuple<0, T, TArgs...>{ typedef std::tuple<TArgs...> Type; };



class CurrentInjector
{
  public:
    virtual ~CurrentInjector() {}
    // Return a FDouble* that points to the current injected into a given
    // Compartment (return NULL if none)
    virtual FDouble const* getCompartmentCurrentPtr(
      const NeuronCompartment & compartment) const = 0;
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
    virtual void defineStateValues(NeuronModel & model) = 0;
    // Initialize the CurrentInjector after parameters are set
    virtual void initialize(void) {}
    // Overrideable function where CurrentInjector can request the creation
    //  of solute objects. It should call useSolute() as needed
    virtual void getSolutes(std::vector<Solute> & solutes) = 0;
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const std::vector<FDouble> & states,
                                    const double & time,
                                    std::vector<FDouble> & dStates) const = 0;
    // Compute injected current
    virtual void computeCurrent(const std::vector<FDouble> & states) = 0;
    
    
    // Injector name
    std::string name;
    // temperature
    double* celsius;
  
  protected:
    // Return a pointer to requested Solute (specified by name and compartment)
    // Create the Solute and add to solutes first if it doesn't already exist
    Solute & useSolute(const std::string & soluteName,
                             NeuronCompartment const* compartmentPtr,
                             std::vector<Solute> & solutes) {
      for(Solute & solute : solutes) {
        if(solute.name == soluteName &&
           solute.compartmentPtr == compartmentPtr) {
          return solute;
        }
      }
      Solute solute;
      solute.name = soluteName;
      solute.compartmentPtr = compartmentPtr;
      // mark the concentration has not having a state set
      solute.cIntIsState = false;
      solutes.push_back(std::move(solute));
      return solutes.back();
    }
    
    virtual size_t defineStateValue(const ParameterDescription & paramDesc,
                                    NeuronModel & model) const = 0;
    
    // Define a parameter and match it to existing parameters; or set it to
    // a default value or the intial value of a trace
    virtual void defineParameter(double & value,
                                 ParameterDescription paramInfo,
                                 ParameterDescriptionList & parametersList,
                                 const std::list<Trace> & traces,
                                 NeuronModel & model) = 0;
    
    // return a string to append to state variable names
    virtual std::string getNameSuffix(void) const = 0;
};


// create a factory for dynamically loading CurrentInjectors
extern Factory<CurrentInjector> injectorFactory;



//////////////////////////// OnePointInjector class ///////////////////////////
// OnePointInjector connects to the model at one NeuronCompartment. It is the
// natural class to use for ion channels


// OnePointInjector may depend on a variable number of Solutes
template <typename... SoluteTypes>
class OnePointInjector;


// This is the base class, which derives from CurrentInjector and depends on
// no Solutes
template <>
class OnePointInjector<> : public CurrentInjector
{
  public:
    virtual ~OnePointInjector() {}
    // Return a FDouble* that points to the current injected into a given
    // Compartment (return NULL if none)
    virtual FDouble const* getCompartmentCurrentPtr(
      const NeuronCompartment & compartment) const {
      return &compartment == compartmentPtr ? &i : NULL;
    }

    // Return a FDouble* that points to the current injected into a given
    // Solute (return NULL if none)
    virtual FDouble const* getSoluteCurrentPtr(const Solute & solute) const =0;
    // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const std::list<Trace> & traces,
                                  NeuronModel & model) = 0;
    // Define the state values of the Injector
    virtual void defineStateValues(NeuronModel & model) = 0;
    // Initialize the CurrentInjector after parameters are set
    virtual void initialize(void) {}
    // Overrideable function where CurrentInjector can connect to needed
    // Solute objects (creating them if needed). It should call useSolute()
    // for each needed Solute
    virtual void getSolutes(std::vector<Solute> & solutes) { (void)solutes; }
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const std::vector<FDouble> & states,
                                    const double & time,
                                    std::vector<FDouble> & dStates) const = 0;
    // Compute injected current
    virtual void computeCurrent(const std::vector<FDouble> & states) = 0;
    
    
    FDouble i;  // nA/nF
    NeuronCompartment const* compartmentPtr;
  
  protected:
    virtual size_t defineStateValue(const ParameterDescription & paramDesc,
                                    NeuronModel & model) const {
      return model.defineStateValue(paramDesc, compartmentPtr);
    }
  
    // Define a parameter and match it to existing parameters; or set it to
    // a default value or the intial value of a trace
    virtual void defineParameter(double & value,
                                 ParameterDescription paramInfo,
                                 ParameterDescriptionList & parametersList,
                                 const std::list<Trace> & traces,
                                 NeuronModel & model) {
      std::vector<std::string> nameSpecifiers;
      nameSpecifiers.reserve(1 + compartmentPtr->suffixes.size());
      
      for(const std::string & suffix : compartmentPtr->suffixes) {
        nameSpecifiers.push_back(name + "_" + suffix);
      }
      nameSpecifiers.push_back(name);
      
      model.defineParameter(value, paramInfo, nameSpecifiers, parametersList,
                            traces);
    }

    // return a string to append to state variable names
    virtual std::string getNameSuffix(void) const {
      return "_" + this->name;
    }
    
    // No Solutes, so return empty tuple if anyone asks for internal or
    // external concentration
    std::tuple<> getCInt(const std::vector<FDouble> & states) const {
      (void) states;
      return std::tuple<>();
    }
    std::tuple<> getCExt(void) const {
      return std::tuple<>();
    }
};


// This is a recursively-derived class that depends on a variable number of
// Solutes
template <typename Sol, typename... SoluteTypes>
class OnePointInjector<Sol, SoluteTypes...> :
    public OnePointInjector<SoluteTypes...>
{
  public:
    virtual ~OnePointInjector() {}
  protected:
    // define type of Parent class for recursion:
    using ParentType = OnePointInjector<SoluteTypes...>;
    // keep record of number of Solutes
    constexpr static size_t numSolutes = 1 + sizeof...(SoluteTypes);
    
    // Overrideable function where CurrentInjector can connect to needed
    // Solute objects (creating them if needed). It should call useSolute()
    // for each needed Solute
   virtual void getSolutes(std::vector<Solute> & solutes) {
      // call useSolute to get needed Solute
      solutePtr =
        &this->useSolute((std::string)Sol::name, this->compartmentPtr, solutes);
      // recursively call getSolutes on other needed Solute objects
      this->OnePointInjector<SoluteTypes...>::getSolutes(solutes);
    }
    
    // Recursively get internal concentration of each Solute, and return it
    // as a tuple of const FDouble references
    typename HomogeneousTuple<numSolutes, const FDouble&>::Type
    getCInt(const std::vector<FDouble> & states) const {
      const FDouble & cInt = solutePtr->cIntIsState ?
        states[solutePtr->concentrationInd] : solutePtr->cInt;
        
      return std::tuple_cat( std::forward_as_tuple(cInt),
                             ParentType::getCInt(states) );
    }
    // Recursively get external concentration of each Solute, and return it
    // as a tuple of const double references
    typename HomogeneousTuple<numSolutes, const double&>::Type
    getCExt(void) const {
      return std::tuple_cat( std::forward_as_tuple(solutePtr->cExt),
                             ParentType::getCExt() );
    }
    
    // pointer to store connection to Solute object
    Solute* solutePtr;
};



//////////////////////////// TwoPointInjector class ///////////////////////////
// TwoPointInjector connects to the model at two NeuronCompartments. It is the
// natural class to use for synapses, gap junctions, or compartment connections


// TwoPointInjector may depend on a variable number of Solutes
template <typename... SoluteTypes>
class TwoPointInjector;


// This is the base class, which derives from CurrentInjector and depends on
// no Solutes
template <>
class TwoPointInjector<> : public CurrentInjector
{
  public:
    virtual ~TwoPointInjector() {}
    // Return a FDouble* that points to the current injected into a given
    // Compartment (return NULL if none)
    virtual FDouble const* getCompartmentCurrentPtr(
                                 const NeuronCompartment & compartment) const {
      if(&compartment == preCompPtr) {
        return &iPre;
      }
      else if(&compartment == postCompPtr) {
        return &iPost;
      }
      else {
        return NULL;
      }
    }
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
    virtual void defineStateValues(NeuronModel & model) = 0;
    // Initialize the CurrentInjector after parameters are set
    virtual void initialize(void) {}
    // Overrideable function where CurrentInjector can request the creation
    //  of solute objects. It should call useSolute() as needed
    virtual void getSolutes(std::vector<Solute> & solutes) { (void) solutes; }
    // Compute derivatives of any injector internal states
    virtual void computeDerivatives(const std::vector<FDouble> & states,
                                    const double & time,
                                    std::vector<FDouble> & dStates) const = 0;
    // Compute injected current
    virtual void computeCurrent(const std::vector<FDouble> & states) = 0;
    
    
    FDouble iPre;    // nA/nF
    FDouble iPost;   // nA/nF
    
    NeuronCompartment const* preCompPtr;
    NeuronCompartment const* postCompPtr;
  
  protected:
    // Define a state value
    virtual size_t defineStateValue(const ParameterDescription & paramDesc,
                                    NeuronModel & model) const {
      return model.defineStateValue(paramDesc, NULL);
    }
    // Define a parameter and match it to existing parameters; or set it to
    // a default value or the intial value of a trace
    virtual void defineParameter(double & value,
                                 ParameterDescription paramInfo,
                                 ParameterDescriptionList & parametersList,
                                 const std::list<Trace> & traces,
                                 NeuronModel & model) {
      // Define a parameter and match it to existing parameters; or set it to
      // a default value or the intial value of a trace
      std::vector<std::string> nameSpecifiers = {
        name + "_" + this->preCompPtr->name + "_" + this->postCompPtr->name,
        name
      };
      
      model.defineParameter(value, paramInfo, nameSpecifiers, parametersList,
                            traces);
    }

    // return a string to append to state variable names
    virtual std::string getNameSuffix(void) const {
      return "_" + this->name + "_" + preCompPtr->name + "_"
             + postCompPtr->name;
    }

    // No Solutes, so return empty tuple if anyone asks for internal or
    // external concentration
    std::tuple<> getCInt(const std::vector<FDouble> & states) const {
      (void) states;
      return std::tuple<>();
    }
    std::tuple<> getCExt(void) const {
      return std::tuple<>();
    }
};


// This is a recursively-derived class that depends on a variable number of
// Solutes
template <typename Sol, typename... SoluteTypes>
class TwoPointInjector<Sol, SoluteTypes...> :
    public TwoPointInjector<SoluteTypes...>
{
  public:
    virtual ~TwoPointInjector() {}
  protected:
    // define type of Parent class for recursion:
    using ParentType = TwoPointInjector<SoluteTypes...>;
    // keep record of number of Solutes
    constexpr static size_t numSolutes = 1 + sizeof...(SoluteTypes);

    // Overrideable function where CurrentInjector can connect to needed
    // Solute objects (creating them if needed). It should call useSolute()
    // for each needed Solute
    virtual void getSolutes(std::vector<Solute> & solutes) {
      // call useSolute to get needed Solute
      if(Sol::location == SoluteTypeLocation::Pre) {
        solutePtr =
          &this->useSolute((std::string)Sol::name, this->preCompPtr, solutes);
      }
      else if (Sol::location == SoluteTypeLocation::Post) {
        solutePtr =
          &this->useSolute((std::string)Sol::name, this->postCompPtr, solutes);
      }
      else {
        throw CurrentInjectorException(
          "SoluteType for TwoPointInjector must be declared as pre or post");
      }
      
      // recursively call getSolutes on other needed Solute objects
      ParentType::getSolutes(solutes);
    }
    
    // Recursively get internal concentration of each Solute, and return it
    // as a tuple of const FDouble references
    typename HomogeneousTuple<numSolutes, const FDouble&>::Type
    getCInt(const std::vector<FDouble> & states) const {
      FDouble & cInt = solutePtr->cIntIsState ?
        states[solutePtr->concentrationInd] : solutePtr->cInt;
        
      return std::tuple_cat( std::forward_as_tuple(cInt),
                             ParentType::getCInt(states) );
    }
    // Recursively get external concentration of each Solute, and return it
    // as a tuple of const double references
    typename HomogeneousTuple<numSolutes, const double&>::Type
    getCExt(void) const {
      return std::tuple_cat( std::forward_as_tuple(solutePtr->cExt),
                             ParentType::getCExt() );
    }

    // pointer to store connection to Solute object
    Solute* solutePtr;
};



////////////////////// CurrentInjector Helper Functions ///////////////////////



template <FAD_TYPENAME(T)>
FDouble
linoid(const T & x)
{
  // handle special cases
  if(x < 1.0e-8 && x > -1.0e-8)
      // expand in Taylor series to avoid divide by *small* problems
      return 1.0 + 0.5 * x;

  // region where linoid formula is numerically well-behaved
  return x / (1.0 - exp(-x));
}
template <FLOAT_TYPENAME(T)>
T
linoid(const T & x)
{
  // handle special cases
  if(x < 1.0e-8 && x > -1.0e-8)
      // expand in Taylor series to avoid divide by *small* problems
      return 1.0 + 0.5 * x;

  // region where linoid formula is numerically well-behaved
  return x / (1.0 - exp(-x));
}


template <typename T>
auto
logistic(const T & x)
  -> decltype(1.0 / (1.0 + exp(-x)))
{
  return 1.0 / (1.0 + exp(-x));
}


template <typename S, typename T>
auto
scaleLogistic(const S & coef, const T & x)
 -> decltype(coef / (1.0 + exp(-x)))
{
  return coef / (1.0 + exp(-x));
}


#endif
