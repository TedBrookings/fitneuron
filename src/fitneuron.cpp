#include <iostream>
#include <unistd.h> // for sleep

#include "../include/minimizeMpi.h"
#include "../include/ParameterDescription.h"
#include "../include/NeuronModelMpi.h"



using namespace std;



inline double
square(const double & x)
{
  return x * x;
}


int main(int numArgs, char** args)
{
  // init MPI environment
  MPI::Init(numArgs, args);
  
  int rank = MPI::COMM_WORLD.Get_rank();
  const int masterProc = 0;
  
  // create a mpi-capable neuron model
  NeuronModelMpi model;
  
  
  // Parse input arguments
  //   create a map to hold arguments for the simulation
  ArgMap argMap;
  
  //   set default values
  argMap["startup"] = "startup.txt";
  argMap["accuracy"] = "";
  argMap["tol"] = "1.0e-3";
  argMap["fTol"] = "1.0e-4";
  argMap["populationSize"] = "";
  argMap["logToDisk"] = "1";
  argMap["logToTerminal"] = "1";
  
  //   parse command-line overrides to arguments
  OverrideMap override = parseArgs(numArgs, args, argMap);
  // make sure that ~ expands to user home directory
  userExpand(argMap["startup"]);
  
  //   set values
  model.startupFile = argMap["startup"];
  if(override["accuracy"])
    model.accuracyGoal = stringToDouble(argMap["accuracy"]);
  
  const double tol = stringToDouble(argMap["tol"]);
  const double fTol = stringToDouble(argMap["fTol"]);
  const bool logToDisk = (stringToInt(argMap["logToDisk"]) == 1);
  const bool logToTerminal = (stringToInt(argMap["logToTerminal"]) == 1);

  if(rank == masterProc && logToTerminal) {
	  // display some info about the mpi environment
    displayMpiInfo();
  }

  // initialize the model (master proc from startup, slaves via mpi)
  model.verbosity = (int)stringToInt(argMap["logToTerminal"]);
  model.initialize();
  
  // chose the population size based upon the number of regular fit parameters
  //  (the number of fit parameters that aren't state parameters) or the number
  //  of fit state parameters, whichever is biggest
  //const size_t numFitParameters = model.getNumFitParameters();
  const size_t numFitRegularParameters = model.getNumFitRegularParameters();
  const size_t numFitStateParameters = model.getNumFitStateParameters();
  
  // empirical formula for best population size given a number of "hard"
  // parameters (regular) and "easy" parameters (state)
  const size_t populationSize = argMap["populationSize"].empty() ?
    (size_t)ceil( 
      pow(4.75 * (double)numFitRegularParameters +
          4.05 * (double)numFitStateParameters, 1.13) )
    :
    stringToSize(argMap["populationSize"]);
  
  if(rank == masterProc && model.verbosity > 0) {
    cout << "Fitting " << numFitRegularParameters << " regular parameters and "
         << numFitStateParameters << " state parameters.\n";
    cout << "Population size = " << populationSize << '\n';
    cout.flush();
  }
  
  
  // chose the error function
  double (NeuronModelMpi::* minFunc)(const vector<double> &) =
    &NeuronModelMpi::getError;
  
  // minimize the error function
  MinimizeResults results;
  try {
    results = minimizeMpi(model, minFunc, populationSize, tol, fTol,
                          logToDisk, logToTerminal);
  }
  catch(exception & e) {
    MPI::Finalize();
    return 0;
  }
  
  
  // output results
  if(rank == masterProc) {
    cout << results.numEvaluations << " evaluations in "
         << results.numGenerations << " generations.\n";
  
    cout << "Final value: " << results.bestParameters.value << '\n';
      
    cout << "Final best parameters:\n";
    for(size_t n = 0; n < results.bestParameters.parameters.size(); ++n)
      cout << ' ' << results.bestParameters.parameters[n];
    cout << '\n';
    cout.flush();
  }
  if(results.needsMpiAbort) {
    sleep(1);
    MPI::COMM_WORLD.Abort(-1);
  }
  else
    MPI::Finalize();
  return 0;
}
