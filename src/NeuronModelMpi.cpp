#include "../include/NeuronModelMpi.h"



using namespace std;



//////////////////  NeuronModelMpi public member functions  ///////////////////



void
NeuronModelMpi::initialize(void)
{
  if(initialized)
    return;
  
  // get the MPI rank
  mpiRank = MPI::COMM_WORLD.Get_rank();

  // make sure everyone has the startup information
  if(mpiRank == 0) {
    // master process, get the startup information from disk
    fitInfo.readStartupFile(startupFile);
    if(verbosity > 0) {
      cout << "Start information obtained.\n"; cout.flush();
    }
    // dynamic cast NeuronFitInfo as NeuronFitInfoMpi
    NeuronFitInfoMpi* fitInfoPtr =
      reinterpret_cast<NeuronFitInfoMpi*>(&fitInfo);
    // use the expanded interface to send fit info via MPI Bcast
    fitInfoPtr->send(MPI::COMM_WORLD);
    if(verbosity > 0) {
      cout << "Start information sent.\n"; cout.flush();
    }
  }
  else {
    // dynamic cast NeuronFitInfo as NeuronFitInfoMpi
    NeuronFitInfoMpi* fitInfoPtr = reinterpret_cast<NeuronFitInfoMpi*>(&fitInfo);
    // use the expanded interface to receive fit info via MPI Bcast
    fitInfoPtr->receive(MPI::COMM_WORLD);
    if(verbosity > 1) {
      cout << "  " << mpiRank << ": Start information received.\n";
      cout.flush();
    }
    verbosity = 0;
  }
  
  // build the model
  buildModel();
  
  initialized = true;
}



/////////////////////  TraceMpi public member functions  //////////////////////



void
TraceMpi::send(MPI::Intracomm & comm) const
{
  // via MPI Bcast, send the information necessary to reconstitute a trace
  
  // send the target strings
  BcastArrMPIWithLength(targetName, 0, comm);
  BcastArrMPIWithLength(units, 0, comm);

  // send trace action
  BcastMPI((int)traceAction, 0, comm);
  
  // send deltaT
  BcastMPI(deltaT, 0, comm);

  // send trace length
  BcastMPI(trace.size(), 0, comm);
  // for Record traces, that's enough, other traces need additional info
  if(traceAction != Record) {
    // send the trace itself
    BcastArrMPI(trace, 0, comm);
    
    // send units
    BcastArrMPIWithLength(units, 0, comm);
    
    if(traceAction == Fit) {
      // for fit traces, send fitTau, waitTime
      BcastMPI(fitTau, 0, comm);
      BcastMPI(waitTime, 0, comm);
    }
  }
}



void
TraceMpi::receive(MPI::Intracomm & comm)
{
  // recieve the target strings
  BcastArrMPIWithLength(targetName, 0, comm);
  BcastArrMPIWithLength(units, 0, comm);
  
  // receive the trace action
  int intTraceAction;
  BcastMPI(intTraceAction, 0, comm);
  traceAction = (TraceAction)intTraceAction;
  
  // receive deltaT
  BcastMPI(deltaT, 0, comm);
  
  // receive trace length and resize
  size_t len;
  BcastMPI(len, 0, comm);
  trace.resize(len);
  // for Record traces, that's enough, other traces need additional info
  if(traceAction != Record) {
    // receive the trace itself
    BcastArrMPI(trace, 0, comm);

    // receive units
    BcastArrMPIWithLength(units, 0, comm);

    if(traceAction == Fit) {
      // for fit traces resize errorTrace
      errorTrace.resize(len);
      // and recieve fitTau, waitTime
      BcastMPI(fitTau, 0, comm);
      BcastMPI(waitTime, 0, comm);
    }
  }
}



/////////////////  NeuronFitInfoMpi public member functions  //////////////////



void
NeuronFitInfoMpi::send(MPI::Intracomm & comm) const
{
  // send all the neuron fit information via MPI Bcast
  //  (note: some of it must be rebuilt by the recipient)
  // send the startup file name
  BcastArrMPIWithLength(startupFileName, 0, comm);
  // send the startup file text
  BcastArrMPIWithLength(startupFileText, 0, comm);

  // send the geometry
  //  -dynamic cast as NeuronGeometryMpi*
  NeuronGeometryMpi const* geoPtr =
    reinterpret_cast<NeuronGeometryMpi const*>(&geometry);
  //  -use the expanded interface to send geometry via MPI bcast
  geoPtr->send(comm);
  
  // send tFinal
  BcastMPI(tFinal, 0, comm);

  // send the electrophys data traces
  // first send the number of traces
  BcastMPI(electrophysData.size(), 0, comm);
  for(const Trace & trace : electrophysData) {
    // loop through the traces
    //   -dynamic cast as ElectrophysDataMpi*
    //    NOTE: can do this because TraceMpi differs *ONLY* by adding functions
    TraceMpi const* tracePtr = reinterpret_cast<TraceMpi const*>(&trace);
    //   -use the expanded interface to send trace via MPI bcast
    tracePtr->send(comm);
  }
}



void
NeuronFitInfoMpi::receive(MPI::Intracomm & comm)
{
  // recieve all the neuron fit information via MPI Bcast
  //  (note: some of it must be rebuilt by the recipient)
  
  // receive the startup file name
  BcastArrMPIWithLength(startupFileName, 0, comm);
  
  // recieve the startup file text
  BcastArrMPIWithLength(startupFileText, 0, comm);
  
  // receive (and build) the geometry
  //  -dynamic cast as NeuronGeometryMpi*
  NeuronGeometryMpi* geoPtr = reinterpret_cast<NeuronGeometryMpi*>(&geometry);
  //  -use the expanded interface to receive geometry via MPI bcast
  geoPtr->receive(comm);
  
  // recieve tFinal
  BcastMPI(tFinal, 0, comm);
  
  // receive the electrophys data traces
  size_t len;
  BcastMPI(len, 0, comm);
  electrophysData.resize(len);
  
  for(Trace & trace : electrophysData) {
    // loop through the traces
    //   -dynamic cast as TraceMpi*
    //    NOTE: can do this because TraceMpi differs *ONLY* by adding functions
    TraceMpi* tracePtr = reinterpret_cast<TraceMpi*>(&trace);
    //   -use the expanded interface to receive trace via MPI bcast
    tracePtr->receive(comm);
    
    // shorten global tFinal if tF < tFinal
    const double tF = tracePtr->deltaT * (double)(tracePtr->trace.size() - 1);
    if(tF < tFinal)
      tFinal = tF;
  }

  dataIsLoaded = true;  // we've just received it via MPI Bcast

  // build the fit information on the recipient side
  buildFitInfo();
}



/////////////////  NeuronGeometryMpi public member functions  /////////////////



void
NeuronGeometryMpi::send(MPI::Intracomm & comm) const
{
  // send the geometry file name
  BcastArrMPIWithLength(geometryFileName, 0, comm);
  
  // send the geometry file content
  BcastArrMPIWithLength(geometryFileContent, 0, comm);
}



void
NeuronGeometryMpi::receive(MPI::Intracomm & comm)
{
  // recieve the geometry file name
  BcastArrMPIWithLength(geometryFileName, 0, comm);

  // receive the geometry file content
  BcastArrMPIWithLength(geometryFileContent, 0, comm);
  
  // build geometry member data from the file content we just received
  buildGeometry();
}

