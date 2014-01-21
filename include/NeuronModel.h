#ifndef _NEURON_MODEL_H_
#define _NEURON_MODEL_H_



#include "ParameterDescription.h"
#include "NeuronFitInfo.h"
#include "Minimizable.h"
#include "dynamic_class_loader.h"
#include "differential_equations.h"
#include "NeuronException.h"

#include <vector>
#include <list>
#include <string>



// to make NeuronModel work:
//   create Minimizable::initialize();
//     -actually, maybe change to Minimizable::initializeMinimizable()
//     -call DiffEqFunctor::addTraceTarget()
//     -call DiffEqFunctor::initializeDiffEq()
//   create DiffEqFunctor::getDerivatives()

// forward declaration of CurrentInjector class
class CurrentInjector;



struct NeuronCompartment
{
  // name of this compartment
  std::string name;
  // set of alternate names and tags that apply to this compartment
  std::vector<std::string> suffixes;
  // index to compartment's voltage in states vector
  size_t voltageInd;
  // vector of pointers to currents injected into this compartment by
  // CurrentInjector objects
  std::vector<FDouble const*> i_currentInjectors;  // pointer to nA / nF
  
  // current injected directly into this compartment
  double i;                   // nA
  
  // electrical properties:
  double capacitance;         // nF
  double axialResistance;     // MOhm  (from one endpoint to other)
  
  // raw geometric properties:
  double length;              // um
  double crossSectionArea;    // um^2
  double surfaceArea;         // mm^2
  double volume;              // mm^3
  
  // properties bridging electrical and geometric:
  double specificCapacitance; // nF/mm^2
  double axialResistivity;    // MOhm * um
};


struct Solute
{
  // Solute's name
  std::string name;
  // compartment this Solute lives in
  NeuronCompartment const* compartmentPtr;
  // external concentration
  double cExt;
  // internal concentration (when no Accumulator is adjusting concentration)
  FDouble cInt;
  // index to solute's internal concentration in states vector
  size_t concentrationInd;
  // boolean indicating if concentration is updated by differential equation
  bool cIntIsState;
  // vector of pointers to currents injected into this compartment by
  // CurrentInjector objects
  std::vector<FDouble const*> i_currentInjectors;  // pointer to nA / nF
};


// improvements:
//   NeuronCompartment class -> NeuronCompartment 
//   create outputFinalStates() functor (or similar name)
//   maybe rename IonChannel and create something like
//     OnePointProcess for channels TwoPointProcess for synapses, etc

// remove dot(), display(), orthogonalize() ?

class NeuronModel : public Minimizable, public DiffEqFunctor
{
  public:
    NeuronModel(void);
    virtual ~NeuronModel();
    
    // prepare the NeuronModel for simulation
    //   - call after assigning startupFile
    virtual void initialize(void) override;
    // get the number of parameters
    inline size_t getNumParameters(void) const
      {return fitInfo.parameterDescriptionList.size(); }
    // get the number of parameters that must be fit (i.e. are not constant)
    inline size_t getNumFitParameters(void) const
      {return fitInfo.parameterDescriptionList.numFitParameters(); }
    // get the number of regular (non-state) parameters that must be fit
    inline size_t getNumFitRegularParameters(void) const
      {return fitInfo.parameterDescriptionList.numFitRegularParameters(); }
    // get the number of state parameters that must be fit
    inline size_t getNumFitStateParameters(void) const
      {return fitInfo.parameterDescriptionList.numFitStateParameters(); }
    // set the parameters of the model
    void setParameters(const std::vector<double> & parameters);
    // calculate the error associated with a set of parameters
    double getError(const std::vector<double> & parameters);
    // return the mean squared error for the neuron model
    double getError(void) const;
    // output recorded traces to stream
    void outputRecordedTraces(const std::string & fileName);
    // output recorded traces to stream
    void outputFinalStates(const std::string & fileName) const;
    // define a state value of the NeuronModel and return its index
    size_t defineStateValue(ParameterDescription paramDesc,
                            NeuronCompartment const* compartmentPtr=NULL);
    // Either match a model value to an existing parameter, or fit it/assign
    // default value, depending upon the value's information and model traces
    void defineParameter(double & value, ParameterDescription paramInfo,
                         const std::vector<std::string> & nameSpecifiers,
                         ParameterDescriptionList & parametersList,
                         const std::list<Trace> & traces);

    // simulate model
    void simulate(void) {
      if(needsSetParameters) {
        throw NeuronModelException(
          "The model has parameters that must be set before simulating.");
      }
      
      // solve the differential equations
      DiffEqFunctor::solve(initialStates);
    }
    
    std::string startupFile;     // file specifying model information
      
  protected:
    // use the information in fitInfo to build the neuron model
    void buildModel(void);
    
    // compute d/dt of each state
    virtual void getDerivatives(const std::vector<FDouble> & states,
                                const double & time,
                                std::vector<FDouble> & dStates) const override;
    
    // inherited from Minimizable:    
    // ParameterDescriptionList parameterDescriptionList;
    NeuronFitInfo fitInfo;
    
    std::vector<NeuronCompartment> compartments;
    std::vector<CurrentInjector*> currentInjectors;
    std::vector<Solute> solutes;
    
    double celsius;        // temperature in celsius
    bool initialized;      // has Minimizable::initialize() been called?
  
  private:
    // build the compartments and join them into the correct geometry
    void buildModelGeometry(void);
    // add CurrentInjectors to the compartments as indicated by NeuronFitInfo
    void addCurrentInjectors(void);
    // return index to compartment with specified name
    size_t getCompartmentIndexByName(const std::string & compartmentName)const;
    // Add specified OnePointInjector to the NeuronModel
    void addOnePointCurrentInjector(NeuronCompartment & compartment,
                                    const std::string & injectorName);
    // add specified TwoPointCurrentInjector to the NeuronModel
    void addTwoPointCurrentInjector(NeuronCompartment & preComp,
                                    NeuronCompartment & postComp,
                                    const std::string & injectorName);
    // define the parameters of the model so that they can be set
    void defineParameters(void);
    // prepare the NeuronModel to call simulate()
    void initializeSimulate(void);
    
    bool needsSetParameters;
    std::vector<double> initialStates;
    // store a list of compartments that have alias info for matching
    // state names to parameters
    std::vector<NeuronCompartment const*> stateCompartments;
};



// create a factory for dynamically loading CurrentInjectors
extern Factory<CurrentInjector> injectorFactory;



#include "CurrentInjector.h"


 
#endif
