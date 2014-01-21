include Makefile.common
TARGET		= bin/simulate_neuron.bin bin/fitneuron.bin bin/perturbneuron.bin
OBJS		= bin/ParameterDescription.o bin/SignalHandler.o bin/simulate_neuron.o bin/io_functions.o bin/Trace.o bin/Timer.o bin/ParameterSet.o bin/NeuronFitInfo.o bin/NeuronModel.o bin/differential_equations.o bin/ParameterPopulation.o bin/NeuronModelMpi.o bin/perturbneuron.o bin/NeuronGeometry.o bin/probability.o bin/fitneuron.o bin/MpiParameterPopulation.o
CPPFLAGS	= $(DEFAULT_CPPFLAGS)
LDFLAGS		= $(DEFAULT_LDFLAGS)
PROJECTDIR = .

all:
	@mkdir -p bin
	@echo make $(TARGET)
	@$(MAKE) --silent $(OBJS)
	@$(MAKE) --silent $(TARGET)
	@echo Done.
	@make --silent -C injectors all

# object files:
bin/ParameterDescription.o: src/ParameterDescription.cpp include/ParameterDescription.h include/io_functions.h
	@echo "   $(CCC) bin/ParameterDescription.o"
	@$(CCC) -o bin/ParameterDescription.o -c src/ParameterDescription.cpp $(CPPFLAGS) 
bin/SignalHandler.o: include/SignalHandler.h src/SignalHandler.cpp
	@echo "   $(CCC) bin/SignalHandler.o"
	@$(CCC) -o bin/SignalHandler.o -c src/SignalHandler.cpp $(CPPFLAGS) 
bin/simulate_neuron.o: include/NeuronGeometry.h include/Minimizable.h include/NeuronFitInfo.h include/NeuronException.h include/ParameterDescription.h include/Timer.h include/io_functions.h include/constants.h include/automatic_derivatives/Fad_operators.h include/CurrentInjector.h include/automatic_derivatives/Fad_functions.h src/simulate_neuron.cpp include/automatic_derivatives/Fad_construct.h include/Trace.h include/NeuronModel.h include/dynamic_class_loader.h include/ParameterSet.h include/differential_equations.h include/automatic_derivatives/Fad.h
	@echo "   $(CCC) bin/simulate_neuron.o"
	@$(CCC) -o bin/simulate_neuron.o -c src/simulate_neuron.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/io_functions.o: include/io_functions.h include/constants.h src/io_functions.cpp
	@echo "   $(CCC) bin/io_functions.o"
	@$(CCC) -o bin/io_functions.o -c src/io_functions.cpp $(CPPFLAGS) 
bin/Trace.o: include/io_functions.h include/Trace.h src/Trace.cpp
	@echo "   $(CCC) bin/Trace.o"
	@$(CCC) -o bin/Trace.o -c src/Trace.cpp $(CPPFLAGS) 
bin/Timer.o: include/io_functions.h src/Timer.cpp include/Timer.h
	@echo "   $(CCC) bin/Timer.o"
	@$(CCC) -o bin/Timer.o -c src/Timer.cpp $(CPPFLAGS) 
bin/ParameterSet.o: include/ParameterSet.h include/io_functions.h src/ParameterSet.cpp include/constants.h
	@echo "   $(CCC) bin/ParameterSet.o"
	@$(CCC) -o bin/ParameterSet.o -c src/ParameterSet.cpp $(CPPFLAGS) 
bin/NeuronFitInfo.o: include/NeuronGeometry.h include/NeuronFitInfo.h include/ParameterDescription.h include/NeuronException.h include/io_functions.h include/constants.h include/Trace.h src/NeuronFitInfo.cpp
	@echo "   $(CCC) bin/NeuronFitInfo.o"
	@$(CCC) -o bin/NeuronFitInfo.o -c src/NeuronFitInfo.cpp $(CPPFLAGS) 
bin/NeuronModel.o: include/NeuronGeometry.h include/Minimizable.h include/NeuronFitInfo.h include/NeuronException.h include/ParameterDescription.h include/io_functions.h include/constants.h include/automatic_derivatives/Fad_operators.h include/CurrentInjector.h include/automatic_derivatives/Fad_functions.h include/automatic_derivatives/Fad_construct.h include/Trace.h include/NeuronModel.h include/dynamic_class_loader.h include/differential_equations.h src/NeuronModel.cpp include/automatic_derivatives/Fad.h
	@echo "   $(CCC) bin/NeuronModel.o"
	@$(CCC) -o bin/NeuronModel.o -c src/NeuronModel.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/differential_equations.o: include/automatic_derivatives/Fad_construct.h include/Trace.h include/automatic_derivatives/Fad.h include/io_functions.h src/differential_equations.cpp include/differential_equations.h include/constants.h include/automatic_derivatives/Fad_operators.h include/automatic_derivatives/Fad_functions.h
	@echo "   $(CCC) bin/differential_equations.o"
	@$(CCC) -o bin/differential_equations.o -c src/differential_equations.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion
bin/ParameterPopulation.o: include/SignalHandler.h src/ParameterPopulation.cpp include/ParameterPopulation.h include/ParameterDescription.h include/Timer.h include/io_functions.h include/constants.h include/probability.h include/ParameterSet.h
	@echo "   $(CCC) bin/ParameterPopulation.o"
	@$(CCC) -o bin/ParameterPopulation.o -c src/ParameterPopulation.cpp $(CPPFLAGS) 
bin/NeuronModelMpi.o: include/NeuronGeometry.h include/Minimizable.h include/NeuronFitInfo.h src/NeuronModelMpi.cpp include/ParameterDescription.h include/io_functions.h include/constants.h include/automatic_derivatives/Fad_operators.h include/CurrentInjector.h include/automatic_derivatives/Fad_functions.h include/automatic_derivatives/Fad_construct.h include/Trace.h include/NeuronModelMpi.h include/NeuronModel.h include/dynamic_class_loader.h include/mpi_interface.h include/differential_equations.h include/NeuronException.h include/automatic_derivatives/Fad.h
	@echo "   $(MPICCC) bin/NeuronModelMpi.o"
	@$(MPICCC) -o bin/NeuronModelMpi.o -c src/NeuronModelMpi.cpp $(CPPFLAGS) -Wno-shadow -fexceptions -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/perturbneuron.o: include/NeuronGeometry.h include/Minimizable.h include/ParameterPopulation.h include/ParameterDescription.h include/Timer.h include/io_functions.h include/MpiParameterPopulation.h include/probability.h include/mpi_interface.h include/automatic_derivatives/Fad_construct.h include/Trace.h include/NeuronModel.h include/dynamic_class_loader.h include/ParameterSet.h include/differential_equations.h src/perturbneuron.cpp include/SignalHandler.h include/NeuronFitInfo.h include/NeuronException.h include/constants.h include/CurrentInjector.h include/automatic_derivatives/Fad_functions.h include/NeuronModelMpi.h include/automatic_derivatives/Fad.h include/automatic_derivatives/Fad_operators.h
	@echo "   $(MPICCC) bin/perturbneuron.o"
	@$(MPICCC) -o bin/perturbneuron.o -c src/perturbneuron.cpp $(CPPFLAGS) -Wno-shadow -fexceptions -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/NeuronGeometry.o: include/NeuronGeometry.h include/io_functions.h src/NeuronGeometry.cpp include/constants.h include/NeuronException.h
	@echo "   $(CCC) bin/NeuronGeometry.o"
	@$(CCC) -o bin/NeuronGeometry.o -c src/NeuronGeometry.cpp $(CPPFLAGS) 
bin/probability.o: src/probability.cpp include/probability.h
	@echo "   $(CCC) bin/probability.o"
	@$(CCC) -o bin/probability.o -c src/probability.cpp $(CPPFLAGS) 
bin/fitneuron.o: include/NeuronGeometry.h include/Minimizable.h include/ParameterPopulation.h include/automatic_derivatives/Fad.h include/ParameterDescription.h include/Timer.h include/io_functions.h include/MpiParameterPopulation.h include/probability.h include/ParameterSet.h include/automatic_derivatives/Fad_construct.h include/Trace.h include/NeuronModel.h include/dynamic_class_loader.h include/mpi_interface.h include/differential_equations.h include/SignalHandler.h include/NeuronFitInfo.h include/NeuronException.h include/constants.h include/CurrentInjector.h include/automatic_derivatives/Fad_functions.h include/NeuronModelMpi.h include/minimizeMpi.h include/automatic_derivatives/Fad_operators.h src/fitneuron.cpp include/minimizeMpi_template_src.h
	@echo "   $(MPICCC) bin/fitneuron.o"
	@$(MPICCC) -o bin/fitneuron.o -c src/fitneuron.cpp $(CPPFLAGS) -Wno-shadow -fexceptions -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/MpiParameterPopulation.o: include/SignalHandler.h include/ParameterPopulation.h include/ParameterDescription.h include/Timer.h include/io_functions.h include/constants.h include/MpiParameterPopulation.h include/probability.h include/ParameterSet.h src/MpiParameterPopulation.cpp include/mpi_interface.h
	@echo "   $(MPICCC) bin/MpiParameterPopulation.o"
	@$(MPICCC) -o bin/MpiParameterPopulation.o -c src/MpiParameterPopulation.cpp $(CPPFLAGS) -fexceptions

# target files (executable or shared object):
bin/simulate_neuron.bin: bin/ParameterDescription.o bin/simulate_neuron.o bin/Timer.o bin/Trace.o bin/ParameterSet.o bin/NeuronFitInfo.o bin/NeuronModel.o bin/differential_equations.o bin/io_functions.o bin/NeuronGeometry.o
	@echo "   $(CCC) bin/simulate_neuron.bin"
	@$(CCC) -o bin/simulate_neuron.bin bin/ParameterDescription.o bin/simulate_neuron.o bin/Timer.o bin/Trace.o bin/ParameterSet.o bin/NeuronFitInfo.o bin/NeuronModel.o bin/differential_equations.o bin/io_functions.o bin/NeuronGeometry.o $(LDFLAGS) -rdynamic -fPIC -ldl
bin/fitneuron.bin: bin/ParameterDescription.o bin/SignalHandler.o bin/NeuronFitInfo.o bin/io_functions.o bin/Trace.o bin/Timer.o bin/ParameterSet.o bin/NeuronModel.o bin/differential_equations.o bin/ParameterPopulation.o bin/NeuronModelMpi.o bin/NeuronGeometry.o bin/probability.o bin/fitneuron.o bin/MpiParameterPopulation.o
	@echo "   $(MPICCC) bin/fitneuron.bin"
	@$(MPICCC) -o bin/fitneuron.bin bin/ParameterDescription.o bin/SignalHandler.o bin/NeuronFitInfo.o bin/io_functions.o bin/Trace.o bin/Timer.o bin/ParameterSet.o bin/NeuronModel.o bin/differential_equations.o bin/ParameterPopulation.o bin/NeuronModelMpi.o bin/NeuronGeometry.o bin/probability.o bin/fitneuron.o bin/MpiParameterPopulation.o $(LDFLAGS) -lmpi -lmpi_cxx -rdynamic -fPIC -ldl -lgsl -lgslcblas
bin/perturbneuron.bin: bin/SignalHandler.o bin/ParameterDescription.o bin/NeuronFitInfo.o bin/io_functions.o bin/Trace.o bin/Timer.o bin/ParameterSet.o bin/NeuronModel.o bin/differential_equations.o bin/ParameterPopulation.o bin/NeuronModelMpi.o bin/perturbneuron.o bin/NeuronGeometry.o bin/probability.o bin/MpiParameterPopulation.o
	@echo "   $(MPICCC) bin/perturbneuron.bin"
	@$(MPICCC) -o bin/perturbneuron.bin bin/SignalHandler.o bin/ParameterDescription.o bin/NeuronFitInfo.o bin/io_functions.o bin/Trace.o bin/Timer.o bin/ParameterSet.o bin/NeuronModel.o bin/differential_equations.o bin/ParameterPopulation.o bin/NeuronModelMpi.o bin/perturbneuron.o bin/NeuronGeometry.o bin/probability.o bin/MpiParameterPopulation.o $(LDFLAGS) -lmpi -lmpi_cxx -rdynamic -fPIC -ldl -lgsl -lgslcblas

clean:
	@make --silent -C injectors clean
	@rm -f $(OBJS)
	@rm -f $(TARGET)
	@if [ -d bin ]; then rm -r bin ; fi
