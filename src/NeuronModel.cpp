#include "../include/NeuronModel.h"
#include "../include/io_functions.h"
#include "../include/constants.h"

#include <sstream>
#include <algorithm>

using namespace NeuronGeometry;
using namespace std;



// create a factory for dynamically loading CurrentInjectors
Factory<CurrentInjector> injectorFactory;



// define a macro for defining NeuronModel parameters
#define DEFINE_MODEL_PARAMETER(parameter, defaultValue, units) \
defineParameter(parameter, \
                ParameterDescription (#parameter, units, defaultValue), \
                {}, parametersList, traces)



// define a macro for defining NeuronCompartment model parameters
#define DEFINE_COMPARTMENT_PARAMETER(parameter, defaultValue, units) \
defineParameter(compartment.parameter, \
                ParameterDescription (#parameter, units, defaultValue), \
                nameSpecifiers, parametersList, traces)




//////////////////////////////  Helper functions  /////////////////////////////



template <size_t numSum>
struct SumStruct
{
  static auto
  sum(FDouble const*const* vec)
    -> decltype(SumStruct<numSum-1>::sum(vec) + *vec[numSum-1])
  {
    return SumStruct<numSum-1>::sum(vec) + *vec[numSum-1];
  }
};
template <>
struct SumStruct<1>
{
  static const FDouble &
  sum(FDouble const*const* vec)
  {
    return **vec;
  }
};

template <size_t numSum>
inline auto
sum(FDouble const*const* vec)
 -> decltype(SumStruct<numSum>::sum(vec))
{
  return SumStruct<numSum>::sum(vec);
}



inline void
getVoltageDeriv(FDouble & dVdt,
                FDouble const*const* processCurrents, size_t numCurrents,
                const double & injectedCurrent, const double & capacitance)
{
  switch(numCurrents) {
    case 0:
      dVdt = injectedCurrent / capacitance;
      break;
    case 1:
      dVdt = *processCurrents[0] + (injectedCurrent / capacitance);
      break;
    case 2:
      dVdt = sum<2>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 3:
      dVdt = sum<3>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 4:
      dVdt = sum<4>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 5:
      dVdt = sum<5>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 6:
      dVdt = sum<6>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 7:
      dVdt = sum<7>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 8:
      dVdt = sum<8>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 9:
      dVdt = sum<9>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 10:
      dVdt = sum<10>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 11:
      dVdt = sum<11>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 12:
      dVdt = sum<12>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 13:
      dVdt = sum<13>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 14:
      dVdt = sum<14>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 15:
      dVdt = sum<15>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 16:
      dVdt = sum<16>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 17:
      dVdt = sum<17>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 18:
      dVdt = sum<18>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 19:
      dVdt = sum<19>(processCurrents) + (injectedCurrent / capacitance);
      break;
    case 20:
      dVdt = sum<20>(processCurrents) + (injectedCurrent / capacitance);
      break;
    default:
      getVoltageDeriv(dVdt, processCurrents + 20, numCurrents - 20,
                      injectedCurrent, capacitance);
      dVdt = dVdt + sum<20>(processCurrents);
  };
}


inline void
getVoltageDeriv(FDouble & dVdt, const vector<FDouble const*>& processCurrents,
                const double & injectedCurrent, const double & capacitance)
{
  getVoltageDeriv(dVdt, &processCurrents.front(), processCurrents.size(),
                  injectedCurrent, capacitance);
}


inline void
getSoluteDeriv(FDouble & dCdt, const FDouble & concentration,
               FDouble const*const* processCurrents, size_t numCurrents,
               const double & rate, const double & cEq,
               const double & coefficient)
{
  switch(numCurrents) {
    case 0:
      dCdt = rate * (cEq - concentration);
      break;
    case 1:
      dCdt = rate * (cEq - concentration + coefficient * *processCurrents[0]);
      break;
    case 2:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<2>(processCurrents) );
      break;
    case 3:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<3>(processCurrents) );
      break;
    case 4:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<4>(processCurrents) );
      break;
    case 5:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<5>(processCurrents) );
      break;
    case 6:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<6>(processCurrents) );
      break;
    case 7:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<7>(processCurrents) );
      break;
    case 8:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<8>(processCurrents) );
      break;
    case 9:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<9>(processCurrents) );
      break;
    case 10:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<10>(processCurrents) );
      break;
    case 11:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<11>(processCurrents) );
      break;
    case 12:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<12>(processCurrents) );
      break;
    case 13:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<13>(processCurrents) );
      break;
    case 14:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<14>(processCurrents) );
      break;
    case 15:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<15>(processCurrents) );
      break;
    case 16:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<16>(processCurrents) );
      break;
    case 17:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<17>(processCurrents) );
      break;
    case 18:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<18>(processCurrents) );
      break;
    case 19:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<19>(processCurrents) );
      break;
    case 20:
      dCdt = rate * ( cEq - concentration
                      + coefficient * sum<20>(processCurrents) );
      break;
    default:
      getSoluteDeriv(dCdt, concentration,
                     processCurrents + 20, numCurrents - 20,
                     rate, cEq, coefficient);
      dCdt = dCdt + rate * coefficient * sum<20>(processCurrents);
  };
}


inline void
getSoluteDeriv(FDouble & dCdt, const FDouble & concentration,
               const vector<FDouble const*> & processCurrents,
               const double & rate, const double & cEq,
               const double & coefficient)
{
  getSoluteDeriv(dCdt, concentration,
                 &processCurrents.front(), processCurrents.size(),
                 rate, cEq, coefficient);
}



inline void
connectTwoPointInjector(TwoPointInjector<>* injectorPtr,
                        NeuronCompartment & preComp,
                        NeuronCompartment & postComp)
{
  // connect two compartments together via a TwoPointInjector
  
  // specify compartments in CurrentInjector
  injectorPtr->preCompPtr = &preComp;
  injectorPtr->postCompPtr = &postComp;
  
  // Specify compartment currents (if any) in NeuronCompartments
  auto currentPtr = injectorPtr->getCompartmentCurrentPtr(preComp);
  if(currentPtr != NULL) {
    preComp.i_currentInjectors.push_back(std::move(currentPtr));
  }
  currentPtr = injectorPtr->getCompartmentCurrentPtr(postComp);
  if(currentPtr != NULL) {
    postComp.i_currentInjectors.push_back(std::move(currentPtr));
  }
}



///////////////////// NeuronModel public member functions /////////////////////

NeuronModel::NeuronModel(void) : Minimizable() , DiffEqFunctor(),
  startupFile ("startup.txt"),
  initialized(false), needsSetParameters(false)
{
  accuracyGoal = 1.0;
  minDeltaT = 1.0e-6;
}



NeuronModel::~NeuronModel()
{
  for(CurrentInjector* & injectorPtr : currentInjectors) {
    delete injectorPtr;
  }
}



void
NeuronModel::initialize(void)
{
  // prepare the NeuronModel for simulation - call after assigning startupFile
  
  if(initialized)
    return;
  
  // get the startup information from disk
  fitInfo.readStartupFile(startupFile);
  if(verbosity > 0) {
    cout << "Start information obtained.\n"; cout.flush();
  }

  // build the model
  buildModel();
  
  initialized = true;
}



void
NeuronModel::setParameters(const vector<double> & parameters)
{
  // set the parameters of the model
  
  // use the list of parameter descriptions to set all internal variables to
  //  their appropriate values, as specified by parameters
  fitInfo.parameterDescriptionList.set(parameters);
  
  // prepare the model for simulation
  initializeSimulate();
  
  // the parameters no longer need to be set
  needsSetParameters = false;
}



double
NeuronModel::getError(const vector<double> & parameters)
{
  // calculate the error associated with a set of parameters
  
  // set the model parameters to the supplied parameters
  setParameters(parameters);
  
  try {
    // simulate the model
    simulate();
    
    return getError();
  }
  catch(DifferentialEquationException const & e) {
    cout << e.what() << '\n';
    cout << "parameters:";
    for(const double & param : parameters) {
      cout << ' ' << param;
    }
    cout << '\n'; cout.flush();
    return Inf;
  }
}



double
NeuronModel::getError(void) const
{
  // return the mean squared error for the neuron model
  
  double error = 0;
  for(const Trace & trace: traces) {
    error += trace.getRMSError();
  }
  
  if(verbosity > 0) {
    cout << "Total trace error = " << error << '\n';
    cout.flush();
  }
  return error;
}



void
NeuronModel::outputRecordedTraces(const string & fileName)
{
  fitInfo.electrophysData = traces;
  fitInfo.outputRecordedTraces(fileName);
}



void
NeuronModel::outputFinalStates(const string & fileName) const
{
  const vector<double> finalStates = getFinalStates();
  const vector<string> stateNames = getStateNames();
  ofstream fileOut (fileName.c_str());
  for(size_t n = 0; n < stateNames.size(); ++n) {
    fileOut << stateNames[n] << "  " << finalStates[n] << '\n';
  }
  fileOut.close();
}

size_t
NeuronModel::defineStateValue(ParameterDescription paramDesc,
                              NeuronCompartment const* compartmentPtr)
{
  // define a state value of the NeuronModel and return its index
  
  // record the base name of the state
  const string baseName = paramDesc.name;
  // make the parameter description specific by appending compartment name
  // (otherwise won't be able to add two states with the same name!!!)
  paramDesc.name += "_" + compartmentPtr->name;
  
  // Need a way to ensure a Trace targeting this state could succeed if it
  // targets an alias.
  // loop through all the Traces:
  for(Trace & trace : fitInfo.electrophysData) {
    // loop through all the possible aliases of this state:
    for(const string & suffix : compartmentPtr->suffixes) {
      const string alias = baseName + "_" + suffix;
      // if a Trace is targeting this alias, change the Trace's targetName to
      // the true name
      if(trace.targetName == alias) {
        trace.targetDescriptiveName = trace.targetName;
        trace.targetName = paramDesc.name;
        break;
      }
    }
  }
  
  // add a description of the allowed values for this state value
  parameterDescriptionList.add(paramDesc);
  // ensure that the description sets isStateParameter = true
  parameterDescriptionList.back().isStateParameter = true;
  // store compartment information (to aid in matching state names to
  // parameter information
  stateCompartments.push_back(compartmentPtr);
  // add the state value to the differential equation and return its index
  return addStateValue(paramDesc.name, paramDesc.units);
}



void
NeuronModel::defineParameter(double & value,
                             ParameterDescription paramInfo,
                             const vector<string> & nameSpecifiers,
                             ParameterDescriptionList & parametersList,
                             const list<Trace> & traces)
{
  // Either match a model value to an existing parameter, or fit it/assign
  // default value, depending upon the value's information and model traces
  
  // keep track of whether a fit trace has matched this parameter
  bool fit = false;
  
  // check if there is a trace / parameter that matches a specified version of
  // the name from paramInfo
  for(const string & specifier : nameSpecifiers) {
    // get the specific name
    const string name = paramInfo.name + "_" + specifier;
    // try to match a trace to name
    for(const Trace & trace : traces) {
      if(trace.targetName == name) {
        DiffEqFunctor::addTraceTarget(value, name, paramInfo.units);
        if(trace.traceAction == Trace::Clamp) {
          // this state is clamped so set its initial value to the start of the
          // trace and return
          value = trace.trace.front();
          return;
        }
        else if(trace.traceAction == Trace::Fit) {
          // this state is fit so *provisionally* set its initial value to the
          // start of the trace, but allow this to be overridden by parameters
          value = trace.trace.front();
          fit = true;
          break;
        }
      }
    }
    
    // try to match a parameter to name
    if(parametersList.match(name, value, paramInfo.isStateParameter)) {
      // found a matching specific parameter
      return;
    }
  }
  
  // check if there is a trace that matches a general version of the name from
  // paramInfo
  for(const Trace & trace : traces) {
    if(trace.targetName == paramInfo.name) {
      DiffEqFunctor::addTraceTarget(value, paramInfo.name, paramInfo.units);
      if(trace.traceAction == Trace::Clamp) {
        // this state is clamped so set its initial value to the start of the
        // trace and return
        value = trace.trace.front();
        return;
      }
      else if(trace.traceAction == Trace::Fit) {
        // this state is fit so *provisionally* set its initial value to the
        // start of the trace, but allow this to be overridden by parameters
        value = trace.trace.front();
        fit = true;
        break;
      }
    }
  }

  
  // check if there is a parameter that matches a general version of the name
  // from paramInfo
  if(parametersList.match(paramInfo.name, value, paramInfo.isStateParameter)) {
    // found a matching specific parameter
    return;
  }
  
  if(!fit) {
    // this value isn't being fit and there is no pre-existing parameter
    if(paramInfo.isConstant()) {
      // set it to the default value
      value = paramInfo.defaultValue();
    }
    else {
      // mark the value as uninitialized
      value = NaN;
      // set the paramInfo to use the most specific name
      if(!nameSpecifiers.empty()) {
        paramInfo.name = paramInfo.name + "_" + nameSpecifiers.front();
      }
      // add a new parameter that must be fit
      parametersList.add(paramInfo);
      //  match the value of this variable to the new parameter description
      if(!parametersList.match(paramInfo.name, value,
                               paramInfo.isStateParameter)) {
        // (this should never happen!)
        throw NeuronModelException("Could not match state parameter: " +
                                   paramInfo.name);
      }
    }
  }
}



//////////////////// NeuronModel protected member functions ///////////////////



void
NeuronModel::buildModel(void)
{
  // use the information in fitInfo to build the neuron model
  if(verbosity > 0) {
    cout << "Building model:\n"; cout.flush();
  }
  // set time information
  DiffEqFunctor::tInitial = 0.0;
  DiffEqFunctor::tFinal = fitInfo.tFinal;
  if(verbosity > 0) {
    cout << "\tbuilding model geometry..."; cout.flush();
  }
  
  // load the CurrentInjector shared-object files
  for(const string & injectorDir : fitInfo.injectorDirs) {
    injectorFactory.loadLibrariesInDirectory(injectorDir);
  }
  
  // build the compartments and join them into the correct geometry
  buildModelGeometry();
  if(verbosity > 0) {
    cout << " done.\n\tadding current injectors..."; cout.flush();
  }
  
  // add channels to the compartments as indicated by channel map
  addCurrentInjectors();
  if(verbosity > 0) {
    cout << " done.\n\tDefining parameters..."; cout.flush();
  }
  
  // define the parameters of the model so that they can be set
  defineParameters();
  if(verbosity > 0) {
    cout << " done.\n\tInitializing differential equation"; cout.flush();
  }
  
  initializeDiffEq(fitInfo.electrophysData);
  if(verbosity > 0) {
    cout << " done.\n"; cout.flush();
  }
}



void
NeuronModel::getDerivatives(const vector<FDouble> & states,
                            const double & time,
                            vector<FDouble> & dStates) const
{
  // compute dStates
  
  // compute current from currentInjectors
  for(CurrentInjector* injectorPtr : currentInjectors) {
    injectorPtr->computeCurrent(states);
  }

  // compute derivatives of any injector internal states
  for(CurrentInjector* injectorPtr : currentInjectors) {
    injectorPtr->computeDerivatives(states, time, dStates);
  }
  
  // compute dV/dt for each compartment
  for(const NeuronCompartment & compartment : compartments) {
    FDouble & dVdt = dStates[compartment.voltageInd];
    getVoltageDeriv(dVdt, compartment.i_currentInjectors,
                    compartment.i, compartment.capacitance);
    //cout << "t=" << time << ": dV_" << compartment.name << "/dt = ";
    //disp(dVdt);
  }
}



///////////////////// NeuronModel private member functions ////////////////////



void
NeuronModel::buildModelGeometry(void)
{
  // build the compartments' basic geometry and connectivity
  
  NeuronGeometry::Geometry & geometry = fitInfo.geometry;
  
  // create the correct number of compartments
  compartments.resize(geometry.numCompartments);
  
  const double minVAllowed = -150.0;  // mV
  const double maxVAllowed = 120.0;
  const double minVRange = -80.0;
  const double maxVRange = 50.0;
  const ParameterDescription vDesc ("v", "mV", minVAllowed, maxVAllowed,
                                    minVRange, maxVRange);
  
  // loop through all the segments and
  //    -add NeuronCompartments to the model,
  //    -set their basic geometric properties
  auto compItr = compartments.begin(); // iterator to current compartment
  for(const Segment & segment : geometry.segments) {
    // Divide each segment up into the correct number of compartments
    for(size_t n = 0; n < segment.numCompartments; ++n, ++compItr) {
      // assign compartment suffixes
      compItr->suffixes = segment.compartmentNames[n];
      compItr->suffixes.insert(compItr->suffixes.end(),
                               segment.tags.begin(), segment.tags.end());
      // assign compartment name (just the first alias)
      compItr->name = compItr->suffixes.front();
      // assign the basic compartment geometry
      compItr->length = segment.length /
                             (double)segment.numCompartments;
      compItr->crossSectionArea = segment.crossSectionArea;
      compItr->surfaceArea = segment.surfaceArea /
                             (double)segment.numCompartments;
      compItr->volume = segment.volume /
                             (double)segment.numCompartments;

      // make the compartments voltage a state value (allows traces to target
      // it for clamping, recording, fitting)
      compItr->voltageInd = defineStateValue(vDesc, &(*compItr));
    }
  }
  
  // loop through all the nodes and add CompartmentConnections as specified
  // by the geometry
  for(const Node & node : geometry.nodes) {
    // loop over every pair of compartments at this node, and add connections
    // between them
    for(auto compItr1 = node.compartments.cbegin();
        compItr1 != node.compartments.cend();
        ++compItr1) {
      for(auto compItr2 = compItr1 + 1;
          compItr2 != node.compartments.cend();
          ++compItr2) {
        // add CompartmentConnection connecting comp1 and comp2
        addTwoPointCurrentInjector(compartments[*compItr1],
                                   compartments[*compItr2],
                                   "CompartmentConnection");
      }
    }
  }
}



void
NeuronModel::addCurrentInjectors(void)
{
  // add CurrentInjectors to the compartments as indicated by NeuronFitInfo
  
  // add the OnePointCurrentInjectors
  for(const auto & injectorInfo : fitInfo.getOnePointInjectorInfo()) {
    const size_t index =
      getCompartmentIndexByName(injectorInfo.compartmentName);
    addOnePointCurrentInjector(compartments[index], injectorInfo.name);
  }
  
  // add the TwoPointCurrentInjectors
  for(const auto & injectorInfo : fitInfo.getTwoPointInjectorInfo()) {
    const size_t preIndex =
      getCompartmentIndexByName(injectorInfo.preCompName);
    const size_t postIndex =
      getCompartmentIndexByName(injectorInfo.postCompName);
    addTwoPointCurrentInjector(compartments[preIndex],
                               compartments[postIndex],
                               injectorInfo.name);
  }
  
  // connect the Solutes to the CurrentInjectors
  for(Solute & solute : solutes) {
    for(CurrentInjector const* injectorPtr : currentInjectors) {
      FDouble const* iPtr = injectorPtr->getSoluteCurrentPtr(solute);
      if(iPtr != NULL) {
        solute.i_currentInjectors.push_back(std::move(iPtr));
      }
    }
  }
}



inline size_t
NeuronModel::getCompartmentIndexByName(const string & compartmentName) const
{
  // return index to compartment with specified name
  size_t index = 0;
  for(const NeuronCompartment & compartment : compartments) {
    if(compartment.name == compartmentName) {
      return index;
    }
    ++index;
  }
  throw NeuronModelException("Invalid compartment name: " + compartmentName);
  return (size_t)-1;
}



void
NeuronModel::addOnePointCurrentInjector(NeuronCompartment & compartment,
                                        const string & injectorName)
{
  // Add specified OnePointInjector to the NeuronModel
  OnePointInjector<>* injectorPtr;
  try {
    injectorPtr = (OnePointInjector<>*)injectorFactory.createNew(injectorName);
  }
  catch (InvalidDynamicClassException & e) {
    cout << "Available injectors:\n";
    for(const string & injector : injectorFactory.getDynamicNames()) {
      cout << '\t' << injector << '\n';
    }
    throw NeuronModelException("Unknown 1-point current injectory type: "
                               + injectorName);
  }

  // add OnePointInjector to model
  currentInjectors.push_back(injectorPtr);
  
  // assign any needed variables here:
  injectorPtr->name = injectorName;
  injectorPtr->celsius = &celsius;
  
  // connect OnePointInjector to NeuronCompartment
  //   -set the OnePointInjector's compartment
  injectorPtr->compartmentPtr = &compartment;
  //   -update the NeuronCompartment's i_currentInjectors
  auto currentPtr = injectorPtr->getCompartmentCurrentPtr(compartment);
  if(currentPtr != NULL) {
    compartment.i_currentInjectors.push_back(std::move(currentPtr));
  }

  // Get references to any needed Solutes, creating them if they don't exist
  injectorPtr->getSolutes(solutes);

  // define the state values for the OnePointInjector
  injectorPtr->defineStateValues(*this);

  // allow Traces to target the channel's injected current for recording
  DiffEqFunctor::addTraceTarget(injectorPtr->i.value,
                             "i_" + injectorPtr->name + "_" + compartment.name,
                             "nA/nF");
}



void
NeuronModel::addTwoPointCurrentInjector(NeuronCompartment & preComp,
                                        NeuronCompartment & postComp,
                                        const string & injectorName)
{
  // add specified TwoPointCurrentInjector to the NeuronModel
  
  // create the new TwoPointCurrentInjector
  TwoPointInjector<>* injectorPtr;
   try {
    injectorPtr = (TwoPointInjector<>*)injectorFactory.createNew(injectorName);
  }
  catch (InvalidDynamicClassException & e) {
    cout << "Available injectors:\n";
    for(const string & injector : injectorFactory.getDynamicNames()) {
      cout << '\t' << injector << '\n';
    }
    throw NeuronModelException("Unknown 2-point current injector type: "
                               + injectorName);
  }
  
  // assign any needed variables here:
  injectorPtr->name = injectorName;
  injectorPtr->celsius = &celsius;

  // add TwoPointInjector to model
  currentInjectors.push_back(injectorPtr);
  
  // connect compartments together via CompartmentConnection
  connectTwoPointInjector(injectorPtr, preComp, postComp);

  // Get references to any needed Solutes, creating them if they don't exist
  injectorPtr->getSolutes(solutes);
  
  // define the state values for the TwoPointInjector
  injectorPtr->defineStateValues(*this);
  
  // allow Traces to target the channel's injected current for recording
  DiffEqFunctor::addTraceTarget(injectorPtr->iPre.value,
                                "iPre_" + injectorPtr->name + "_"
                                  + preComp.name + "_" + postComp.name, 
                                "nA/nF");
  DiffEqFunctor::addTraceTarget(injectorPtr->iPost.value,
                                "iPost_" + injectorPtr->name + "_"
                                  + preComp.name + "_" + postComp.name, 
                                "nA/nF");  
}



void
NeuronModel::defineParameters(void)
{
  // define the parameters of the model so that they can be set
  
  ParameterDescriptionList & parametersList = fitInfo.parameterDescriptionList;
  const list<Trace> & traces = fitInfo.electrophysData;
  
  // first define the initial states as state parameters
  //std::vector<string> & stateNames = DiffEqFunctor::getStateNames();
  const size_t numStates = parameterDescriptionList.size();
  initialStates.resize(numStates);
  for(size_t n = 0; n < numStates; ++n) {
    const vector<string> & suffixes = stateCompartments[n] == NULL ?
      vector<string>() : stateCompartments[n]->suffixes;
    
    // Get the ParameterDescription
    ParameterDescription paramInfo = parameterDescriptionList[n];
    if(!suffixes.empty()) {
      // Remove the NeuronCompartment name from the end
      size_t suffixInd = paramInfo.name.find_last_of('_');
      paramInfo.name = paramInfo.name.substr(0, suffixInd);
    }
    // Define a parameter
    defineParameter(initialStates[n], paramInfo, suffixes,
                    parametersList, traces);
  }
  
  // define model-wide parameters
  DEFINE_MODEL_PARAMETER(celsius, 10.0, "C");
  
  // define NeuronCompartment parameters
  for(NeuronCompartment & compartment : compartments) {
    const vector<string> & nameSpecifiers = compartment.suffixes;
    
    DEFINE_COMPARTMENT_PARAMETER(specificCapacitance, 10.0, "nF/mm^2");
    // axial resistivity from Rabbah et al 2005
    DEFINE_COMPARTMENT_PARAMETER(axialResistivity, 0.6, "MOhm um");
    DEFINE_COMPARTMENT_PARAMETER(i, 0.0, "nA");
  }
  
  // define CurrentInjector parameters
  for(CurrentInjector* & currentInjectorPtr : currentInjectors) {
    currentInjectorPtr->defineParameters(parametersList, traces, *this);
  }
  
  // define Solute cExt parameters
  const double minCAllowed = 1.0e-9; // mM
  const double maxCAllowed = 10000.0;// mM
  const double minCRange = 0.01;     // mM
  const double maxCRange = 100.0;    // mM
  for(Solute & solute : solutes) {
    const vector<string> & nameSpecifiers = solute.compartmentPtr->suffixes;
    
    ParameterDescription paramInfo ("cExt_" + solute.name, "mM", minCAllowed,
                                    maxCAllowed, minCRange, maxCRange);
    defineParameter(solute.cExt, paramInfo, nameSpecifiers,
                    parametersList, traces);
    if(!solute.cIntIsState) {
      // cInt is not managed by an Accumulator, so make it a parameter
      
      // need to resize cInt
      //solute.cInt.setNumVariables(numStates);
      
      ParameterDescription cIntInfo ("cInt_" + solute.name, "mM", minCAllowed,
                                     maxCAllowed, minCRange, maxCRange);
      defineParameter(solute.cInt.value, std::move(cIntInfo), nameSpecifiers,
                      parametersList, traces);
    }
  }
  
  // check if the model has parameters that must be set
  needsSetParameters = (parametersList.numFitParameters() > 0);
  
  // copy parametersList to the Minimizable base (this is needed for the
  // minimization routine)
  parameterDescriptionList = parametersList;
  
  if(!needsSetParameters) {
    // prepare the model for simulation
    initializeSimulate();
  }
}




void
NeuronModel::initializeSimulate(void)
{
  // prepare the NeuronModel to call simulate()
  
  // initialize each NeuronCompartment
  for(NeuronCompartment & compartment : compartments) {
    compartment.capacitance =
      compartment.specificCapacitance * compartment.surfaceArea;
    compartment.axialResistance = compartment.axialResistivity
      * compartment.length / compartment.crossSectionArea;
  }
  
  // initialize each CurrentInjector
  for(CurrentInjector* & injectorPtr : currentInjectors) {
    injectorPtr->initialize();
  }
}
