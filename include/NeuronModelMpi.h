#ifndef _NEURON_MODEL_MPI_H_
#define _NEURON_MODEL_MPI_H_



#include "NeuronModel.h"
#include "mpi_interface.h"



class NeuronModelMpi : public NeuronModel
{
  public:
    // reimplement the initialize routine because startup data must be sent
    //   from the master process to slave procs via MPI
    virtual void initialize(void);
  protected:
    int mpiRank;     // rank in MPI environment
};



class TraceMpi : public Trace
{
  public:
    virtual ~TraceMpi() {}
    // send the ElectrophysTrace information via MPI Bcast
    void send(MPI::Intracomm & comm) const;
    // receive the ElectrophysTrace information via MPI Bcast
    void receive(MPI::Intracomm & comm);
};



class NeuronFitInfoMpi : public NeuronFitInfo
{
  public:
    virtual ~NeuronFitInfoMpi() {}
    // send all the fit information via MPI Bcast
    void send(MPI::Intracomm & comm) const;
    // receive all the fit information via MPI Bcast
    void receive(MPI::Intracomm & comm);
};



class NeuronGeometryMpi : public NeuronGeometry::Geometry
{
  public:
    virtual ~NeuronGeometryMpi() {}
    // send all the information in the Geometry via MPI Bcast
    void send(MPI::Intracomm & comm) const;
    // receive all the information in the Geometry via MPI Bcast
    void receive(MPI::Intracomm & comm);
};

#endif
