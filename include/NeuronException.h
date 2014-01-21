#ifndef _NEURON_EXCEPTION_H_
#define _NEURON_EXCEPTION_H_



#include <stdexcept>
#include <set>



/////////// base class for all exceptions unique to this project //////////////

class NeuronException : public std::runtime_error
{
  public:
    NeuronException() : std::runtime_error("NeuronException") {}
    NeuronException(const std::string & str) : std::runtime_error(str) {}
};



//////// exceptions specific to various classes (ClassnameException) //////////



class NeuronModelException : public NeuronException
{
  public:
    NeuronModelException() : NeuronException("ModelException") {}
    NeuronModelException(const std::string & str) : NeuronException(str) {}
};


class ConvergeException : public NeuronModelException
{
  public:
    ConvergeException() :
      NeuronModelException("simulate() failed to converge") {}
    ConvergeException(const std::string & str) : NeuronModelException(str) {}
};


class SoluteException : public NeuronException
{
  public:
    SoluteException() : NeuronException("SoluteException") {}
    SoluteException(const std::string & str) : NeuronException(str) {}
};



class CurrentInjectorException : public NeuronException
{
  public:
    CurrentInjectorException() : NeuronException("CurrentInjectorException") {}
    CurrentInjectorException(const std::string & str) : NeuronException(str) {}
};



class NeuronGeometryException : public NeuronException
{
  public:
    NeuronGeometryException() : NeuronException("NeuronGeometryException") {}
    NeuronGeometryException(const std::string & str) : NeuronException(str) {}
 };
 

 
 class FitInfoException : public NeuronException
{
  public:
    FitInfoException() : NeuronException("FitInfoException") {}
    FitInfoException(const std::string & str) : NeuronException(str) {}
};



////////////// functions for frequently repeated exception code ///////////////



inline void
throwMultipleMatchException(const std::string & paramName,
                            const std::string & compartmentName,
                            std::set<std::string> & matchNames)
{
  std::string exceptStr;
  exceptStr = "Multiple fitable parameters matched \"" + paramName
            + "\" in compartment " + compartmentName + ".\n"
            + "Matching parameters:\n";
  for(const std::string & matchName : matchNames) {
    exceptStr = exceptStr + "  " + matchName + '\n';
  }
  
  throw NeuronModelException(exceptStr);
};



#endif
