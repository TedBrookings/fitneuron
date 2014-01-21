#ifndef Accumulator_h
#define Accumulator_h


#include "../include/CurrentInjector.h"



using namespace std;



template <typename Sol>
class Accumulator : public OnePointInjector<Sol>
{
  public:
    // Return a FDouble* that points to the current injected into a given
    // Compartment (return NULL if none)
    virtual FDouble const* getCompartmentCurrentPtr(
                                 const NeuronCompartment & compartment) const {
      (void)compartment;
      return NULL;
    }
    // Return a FDouble* that points to the current injected into a given
    // Solute (return NULL if none)
    virtual FDouble const* getSoluteCurrentPtr(const Solute & solute) const {
      (void) solute;
      return NULL;
    }
    // Define the parameters of the Injector
    virtual void defineParameters(ParameterDescriptionList & parametersList,
                                  const list<Trace> & traces,
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
    virtual void computeCurrent(const vector<FDouble> & states) {
      (void) states;
    }
    
  protected:
    // return a string to append to state variable names
    virtual string getNameSuffix(void) const {
      return "_" + (string)Sol::name;
    }
    
    double tau;           // ms Time constant to approach steady-state at
                          //    base temperature
    double cEquilibrium;  // mM Equilibrium internal concentration
    double domainDepth;   // um Effective depth of microdomains
    double q10;
    double celsiusBase;   // C (base temperature of model)
    
  
  private:
    double rate;
    double coef_I;
};



//////////////////////// Accumulator Member Functions ///////////////////////

template <typename Sol>
void
Accumulator<Sol>::defineStateValues(NeuronModel & model)
{
  // Define the state values of the Injector
  const double minCAllowed = 1.0e-9; // mM
  const double maxCAllowed = 10000.0;// mM
  const double minCRange = 0.01;     // mM
  const double maxCRange = 100.0;    // mM

  this->solutePtr->concentrationInd =
    DEFINE_STATE_VALUE(cInt, minCAllowed, maxCAllowed, minCRange, maxCRange,
                       "mM");
  this->solutePtr->cIntIsState = true;
}

template <typename Sol>
void
Accumulator<Sol>::initialize(void)
{
  // Initialize the CurrentInjector after parameters are set
  
  const double chargePerMole = N_Avogadro * Q_Electron * (double)Sol::Z;
  rate = pow(q10, 0.1 * (*this->celsius - celsiusBase)) / tau;
  coef_I = this->compartmentPtr->specificCapacitance
           / (rate * chargePerMole * domainDepth);
  // note: previously the modification due to temperature/q10 was calculated
  // incorrectly. Incorrect old version:
  /*
  rate = 1.0 / tau;
  coef_I = pow(q10, 0.1 * (*this->celsius - celsiusBase))
           * this->compartmentPtr->specificCapacitance
           / (rate * chargePerMole * domainDepth);
  */
}


// declare getSoluteDeriv()
inline void
getSoluteDeriv(FDouble & dCdt, const FDouble & concentration,
               const vector<FDouble const*> & processCurrents,
               const double & rate, const double & cEq,
               const double & coefficient);


template <typename Sol>
void
Accumulator<Sol>::computeDerivatives(const vector<FDouble> & states,
                                     const double & time,
                                     vector<FDouble> & dStates) const
{
  // Compute derivatives of any injector internal states
  (void)time;
  
  FDouble & dCdt = dStates[this->solutePtr->concentrationInd];
  const FDouble & concentration = states[this->solutePtr->concentrationInd];
  
  getSoluteDeriv(dCdt, concentration,
                 this->solutePtr->i_currentInjectors,
                 rate, cEquilibrium, coef_I);
}



////////////////////////////// Accumulator helper functions ///////////////////

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

#endif
