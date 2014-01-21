#include "../include/NeuronModel.h"
#include "../include/io_functions.h"
#include "../include/Timer.h"
#include "../include/ParameterSet.h"

#include <iostream>
#include <fstream>
#include <string>
#include <map>


using namespace std;


void printHelp(ArgMap & argMap);


int main(int numArguments, char** arguments)
{
  // create a map to hold arguments for the simulation
  ArgMap argMap;
  
  // set default values
  argMap["startup"] = "startup.txt";
  argMap["outfile"] = "simulated_traces.txt";
  argMap["statesfile"] = "final_states.txt";
  argMap["accuracy"] = "";
  argMap["params"] = "";
  argMap["verbosity"] = "1";
  
  // parse command-line overrides to arguments
  OverrideMap override = parseArgs(numArguments, arguments, argMap);
  if(override["help"]) {
    printHelp(argMap);
    return 0;
  }
  // make sure that ~ expands to the home directory
  userExpand("startup");
  userExpand("outfile");
  userExpand("statesfile");
  userExpand("params");
  
  // create a NeuronModel
  NeuronModel model;
  
  // set the start up file
  model.startupFile = argMap["startup"];

  // set the verbosity
  model.verbosity = (int)stringToInt( argMap["verbosity"] );
  
  // set the accuracy goal
  if(override["accuracy"])
    model.accuracyGoal = stringToDouble( argMap["accuracy"] );
  
  // initialize
  model.initialize();
  
  if(override["params"]) {
    // if base parameters were supplied, load them in and set them

    // get the base parameters
    ParameterSet pSet;
    pSet.read(argMap["params"], model.getNumParameters());
    
    // set the model parameters
    model.setParameters(pSet.parameters);
  }
  
  Timer simTimer;
  
  // simulate
  model.simulate();
  
  if(model.verbosity > 0) {
    cout << "Finished simulation in " << simTimer << '\n';
    cout.flush();
    
    // get error from any error traces
    model.getError();
  }

  // output results
  model.outputRecordedTraces(argMap["outfile"]);
  cout << "Traces saved.\n"; cout.flush();
  model.outputFinalStates(argMap["statesfile"]);
  cout << "Final states saved.\n"; cout.flush();
  return 0;
}



void
printHelp(ArgMap & argMap)
{
  cout << "Usage:\n"
       << "simulate_neuron.cpp -help\n"
       << "   print this help\n"
       << "simulate_neuron.cpp [-arg1 val1 [-arg2 val2] ... ]\n"
       << "   simulate a neuron overriding default values specified by "
          " arguments.\n"
       << "Arguments:\n"
       << "  -startup  (" << argMap["startup"] << ") startup file name\n"
       << "  -outfile  (" << argMap["outfile"] << ") output file name\n"
       << "  -statesfile  (" << argMap["statesfile"]
       << ") output final states file name\n"
       << "  -accuracy (" << argMap["accuracy"]
       << ") override default simulation accuracy\n"
       << "  -params    (" << argMap["params"]
       << ") override startup parameters with those from params file\n"
       << "  -verbosity (" << argMap["verbosity"]
       << ") integer, 0 for silent, increase for more\n";
  cout.flush();
}
