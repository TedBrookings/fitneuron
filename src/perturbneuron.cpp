#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include "../include/MpiParameterPopulation.h"
#include "../include/NeuronModelMpi.h"
#include "../include/io_functions.h"



using namespace std;


class PerturbNeuronException : public std::runtime_error
{
  public:
    explicit PerturbNeuronException() :
      std::runtime_error("PerturbNeuronException") {}
    explicit PerturbNeuronException(const std::string & str) :
      std::runtime_error(str) {}
 };



class PerturbMaster: public MpiParameterPopulation
{
  public:
    // this should find the original model, create the population, etc.
    virtual void initialize(void);
    // get the basic features of the perturbation from the command line args
    void parseInputArgs(int numArgs, char** args);
    // create the population of perturbed models
    void createPerturbedPopulation(void);
    // output the perturbation results so that they can be plotted/analyzed
    void outputPerturbationResults(void) const;
    
  private:
    void outputBaseParams(ostream & out) const;
  
    NeuronModelMpi model;

    size_t numPerturb;
    string startupFile;
    string baseParamsFile;
    size_t perturbPopSize;

    // names of parameters to perturb
    vector<string> perturbArgs;
    // indices (in parameterDescription) of parameters to perturb
    vector<size_t> perturbInds;
    // number of steps to take for each perturbing direction
    vector<size_t> perturbNumSteps;
    // maximum change factor along each perturbing direction
    vector<double> perturbChangeFactor;
    // base, unperturbed params
    ParameterSet baseParams;
};


class PerturbSlave : public MpiWorker
{
  protected:
    // initialize the SlaveProcess (just initialize the model)
    virtual void initialize(void);
    // evaluate parameters (get error from the model)
    virtual double evaluateParameters(const std::vector<double> & evalParams)
      { return model.getError(evalParams); }
  
  private:
    NeuronModelMpi model;
};



int main(int numArgs, char** args)
{
  MPI::Init(numArgs, args);
  const int rank = MPI::COMM_WORLD.Get_rank();

  // get the input arguments
  if(numArgs < 4) {
    MPI::Finalize();
    throw PerturbNeuronException("Incorrect number of input arguments");
  }
  
  if(rank == 0) {
    // master process
    
    // display some information about the MPI environment
    displayMpiInfo();
    
    // create a process to manage the distribution and evaluation of parameters
    PerturbMaster master;
    
    // parse input arguments to create the job
    master.parseInputArgs(numArgs, args);
    
    // initialize the population
    master.initialize();
    
    // evaluate all the perturbed models
    master.evaluatePopulation();
    
    // output results
    master.outputPerturbationResults();
  }
  else {
    // slave process
    
    // create a process to recieve parameters and evaluate them
    PerturbSlave slave;
    
    // evaluate parameter sets until shutdown
    slave.evaluatePopulation();
  }

  MPI::Finalize();
  return 0;
}



//////////////////////////////// PerturbMaster ////////////////////////////////



void
PerturbMaster::initialize(void)
{
  // prepare to do the batch evaluation
  
  if(initialized)
    // already initialized, return
    return;
  
  // don't save resume information
  resumeFileName = "/dev/null";
  // don't log to disk
  logToDisk = false;
  
  // set model startup file and initialize model
  model.startupFile = startupFile;
  model.verbosity = (int)logToTerminal;
  model.initialize();
  
  // set the ParameterPopulation's parameter descriptions to the model's
  setParameterDescriptions(model.getParameterDescriptions());
  
  // get the base model that will be perturbed
  baseParams.read(baseParamsFile, parameterDescriptions.size());
  
  // create the population of perturbed models
  createPerturbedPopulation();
  
  // done!
  initialized = true;
}



void
PerturbMaster::parseInputArgs(int numArgs, char** args)
{
  // get the basic features of the perturbation from the command line args
  
  perturbPopSize = 1;
  numPerturb = (size_t)numArgs - 3;
  startupFile = args[1];
  baseParamsFile = args[2];
  resultsFileName = "perturb_";
  
  perturbArgs.resize(numPerturb);
  perturbNumSteps.resize(numPerturb);
  perturbChangeFactor.resize(numPerturb);
  for(size_t n = 0; n < numPerturb; ++n) {
    perturbArgs[n] = args[n + 3];
    perturbNumSteps[n] = 20;
    perturbChangeFactor[n] = 5.0;
    
    perturbPopSize *= (1 + 2 * perturbNumSteps[n]);
    
    if(n > 0)
      resultsFileName += "_vs_";
    resultsFileName += perturbArgs[n];
  }
  
  resultsFileName += ".txt";
}



void
PerturbMaster::createPerturbedPopulation(void)
{
  // create the population of perturbed models
  if(model.verbosity > 0) {
    cout << "Creating perturbed population (" << perturbPopSize << ')';
    cout.flush();
  }

  // make indices to the parameters that will be perturbed
  perturbInds.resize(numPerturb, parameterDescriptions.size());
  for(size_t n = 0; n < numPerturb; n++) {
    for(size_t ind = 0; ind < parameterDescriptions.size(); ++ind)
      if(parameterDescriptions[ind].name == perturbArgs[n]) {
        perturbInds[n] = ind;
        break;
      }
    if(perturbInds[n] == parameterDescriptions.size())
      throw
        PerturbNeuronException("Unknown parameter \"" + perturbArgs[n] + "\"");
  }
  
  // create a matrix to keep track of the values of all the perturbed
  // parameters
  vector<vector<double> > perturbValues (numPerturb);
  // fill the matrix
  
  vector<vector<double> >::iterator rowItr = perturbValues.begin();
  for(size_t n = 0; n < numPerturb; ++n, ++rowItr) {
    // loop over perturb args
    
    // resize the row to the number of perturb steps
    const size_t numSteps = perturbNumSteps[n];
    rowItr->resize(1 + 2 * numSteps);
    
    // calculate how much each step changes the value
    const double delta =
      pow(perturbChangeFactor[n], 1.0 / (double)perturbNumSteps[n]);
      
    // fill the row
    //  -start in the middle, which is unchanged
    const double baseVal = baseParams.parameters[perturbInds[n]];
    double plusVal = baseVal;
    double minusVal = baseVal;
    (*rowItr)[numSteps] = plusVal;
    for(size_t step = 1; step < numSteps; ++step) {
      // moving away from the center, progressively change the value
      plusVal *= delta;
      minusVal /= delta;
      (*rowItr)[numSteps + step] = plusVal;
      (*rowItr)[numSteps - step] = minusVal;
    }
    (*rowItr)[0] = baseVal / perturbChangeFactor[n];
    (*rowItr)[2 * numSteps] = baseVal * perturbChangeFactor[n];
  }
  
  // create a list of indices to keep track of what values to use for each
  // ParameterSet in the population
  vector<size_t> valueInds (numPerturb, 0);
  
  // allocate space for the perturbed population
  setPopulationSize(perturbPopSize);
  population.resize(perturbPopSize);
  
  // fill the population with perturbed ParameterSets
  for(auto & pSet : population) {
    // loop through parameter population, creating each ParameterSet
    
    // set the value to NaN so pSet will be evaluated
    pSet.value = NaN;
    
    // initialize pSet's parameters to baseParams
    pSet.parameters = baseParams.parameters;
    
    // perturb pSet's parameters
    for(size_t n = 0; n < numPerturb; ++n) {
      pSet.parameters[perturbInds[n]] = perturbValues[n][valueInds[n]];
    }
    
    // increment the perturbation
    size_t m = 0;
    ++(valueInds[m]);
    while(valueInds[m] >= 1 + 2 * perturbNumSteps[m]) {
      valueInds[m] = 0;
      ++m;
      if(m == numPerturb)
        break;
      ++(valueInds[m]);
    }
  }

  if(model.verbosity > 0) {
    cout << "  done.\n";
    cout.flush();
  }
}



void
PerturbMaster::outputPerturbationResults(void) const
{
  // output the perturbation results so that they can be plotted/analyzed
  
  // open the output stream
  ofstream fOut (resultsFileName.c_str());
  
  // write out the parameters that were perturbed
  fOut << perturbArgs[0];
  for(size_t n = 1; n < numPerturb; ++n)
    fOut << ' ' << perturbArgs[n];
  fOut  << '\n';
  
  // write out the names and values of all the parameters of the model
  fOut << "numParameters: " << baseParams.parameters.size() << '\n';
  outputBaseParams(fOut);
  
  // write out the whole population
  fOut << "populationSize: " << population.size() << '\n';
  for(const auto & pSet : population)
    fOut << pSet;
  fOut.close();
}



void
PerturbMaster::outputBaseParams(ostream & out) const
{
  // neatly output the base parameters
  
  // calculate a few widths to ensure all the significant info in the double
  //  is output
  const streamsize doublePrecision = numeric_limits<double>::digits10 + 1;
  const streamsize expWidth = (streamsize)ceil(
    log10((double)numeric_limits<double>::max_exponent10));
  const streamsize doubleWidth = doublePrecision + expWidth + 4;
  
  // set output decimal precision
  out.precision(doublePrecision);
  
  // get length of longest parameter description name
  ParameterDescriptionList::const_iterator descItr =
    parameterDescriptions.begin();
  streamsize nameWidth = 5;
  for(; descItr != parameterDescriptions.end(); ++descItr)
    if(descItr->name.size() > (size_t)nameWidth)
      nameWidth = descItr->name.size();
  nameWidth += 1;

  // write out the function value
  out.width(nameWidth);
  out << left << "value";
  out.width(doubleWidth);
  out << right << baseParams.value << '\n';
  
  // write out the parameters
  vector<double>::const_iterator baseItr = baseParams.parameters.begin();
  descItr = parameterDescriptions.begin();
  for(; descItr != parameterDescriptions.end(); ++descItr, ++baseItr) {
    out.width(nameWidth);
    out << left << descItr->name;
    out.width(doubleWidth);
    out << right << *baseItr << '\n';
  }
}



///////////////////////////////// PerturbSlave ////////////////////////////////



void
PerturbSlave::initialize(void)
{
  if(verbosity > 0) {
    cout << "Initializing " << getRank() << '\n';
    cout.flush();
  }
  
  model.accuracyGoal = 0.1;
  model.verbosity = verbosity;
  model.initialize();
  
  if(verbosity > 0) {
    cout << "\tdone.\n";
    cout.flush();
  }
}
