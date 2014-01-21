include ../Makefile.common
TARGET		= bin/NaCS.so bin/CaT.so bin/Na.so bin/CompartmentConnection.so bin/KCa.so bin/CAN.so bin/MI.so bin/Nap.so bin/ACS.so bin/As.so bin/KdCS.so bin/Kd.so bin/CaS.so bin/Af.so bin/Leak.so bin/CaAccumulate.so bin/h.so
OBJS		= bin/CompartmentConnection.o bin/CaAccumulate.o bin/CaS.o bin/h.o bin/KCa.o bin/KdCS.o bin/ACS.o bin/Kd.o bin/Nap.o bin/MI.o bin/NaCS.o bin/CaT.o bin/As.o bin/CAN.o bin/Na.o bin/Leak.o bin/Af.o
CPPFLAGS	= $(DEFAULT_CPPFLAGS)
LDFLAGS		= $(DEFAULT_LDFLAGS)
PROJECTDIR = ..

all:
	@mkdir -p bin
	@echo make $(TARGET)
	@$(MAKE) --silent $(OBJS)
	@$(MAKE) --silent $(TARGET)
	@echo Done.

# object files:
bin/CompartmentConnection.o: ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h ../include/differential_equations.h ../include/NeuronFitInfo.h ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h ../include/Minimizable.h CompartmentConnection.cpp ../include/ParameterDescription.h ../include/Trace.h ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h ../include/automatic_derivatives/Fad_construct.h
	@echo "   $(CCC) bin/CompartmentConnection.o"
	@$(CCC) -o bin/CompartmentConnection.o -c CompartmentConnection.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/CaAccumulate.o: ../include/NeuronFitInfo.h ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h Accumulator.h ../include/differential_equations.h ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h CaAccumulate.cpp ../include/Minimizable.h ../include/ParameterDescription.h ../include/Trace.h ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h ../include/automatic_derivatives/Fad_construct.h
	@echo "   $(CCC) bin/CaAccumulate.o"
	@$(CCC) -o bin/CaAccumulate.o -c CaAccumulate.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/CaS.o: ../include/NeuronModel.h ../include/NeuronGeometry.h ../include/automatic_derivatives/Fad.h ../include/NeuronFitInfo.h ../include/differential_equations.h ../include/automatic_derivatives/Fad_operators.h ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h NernstChannel.h ../include/Minimizable.h ../include/ParameterDescription.h ../include/Trace.h CaS.cpp ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h ../include/automatic_derivatives/Fad_construct.h
	@echo "   $(CCC) bin/CaS.o"
	@$(CCC) -o bin/CaS.o -c CaS.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/h.o: ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h ../include/NeuronFitInfo.h ../include/differential_equations.h ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h h.cpp ../include/Minimizable.h ../include/ParameterDescription.h ../include/Trace.h ../include/automatic_derivatives/Fad_construct.h OhmicBase.h ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h OhmicChannel.h
	@echo "   $(CCC) bin/h.o"
	@$(CCC) -o bin/h.o -c h.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/KCa.o: ../include/NeuronFitInfo.h ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h KCa.cpp ../include/differential_equations.h ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h ../include/Minimizable.h ../include/ParameterDescription.h ../include/Trace.h ../include/automatic_derivatives/Fad_construct.h OhmicBase.h ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h OhmicChannel.h
	@echo "   $(CCC) bin/KCa.o"
	@$(CCC) -o bin/KCa.o -c KCa.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/KdCS.o: ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h ../include/NeuronFitInfo.h ../include/differential_equations.h ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h ../include/NeuronException.h ../include/Minimizable.h ../include/ParameterDescription.h ../include/Trace.h ../include/constants.h OhmicBase.h ../include/NeuronGeometry.h KdCS.cpp ../include/dynamic_class_loader.h ../include/CurrentInjector.h HodgkinHuxleyChannel.h ../include/automatic_derivatives/Fad_construct.h
	@echo "   $(CCC) bin/KdCS.o"
	@$(CCC) -o bin/KdCS.o -c KdCS.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/ACS.o: ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h ../include/NeuronFitInfo.h ../include/differential_equations.h ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h ../include/Minimizable.h ACS.cpp ../include/ParameterDescription.h ../include/Trace.h ../include/automatic_derivatives/Fad_construct.h OhmicBase.h ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h OhmicChannel.h
	@echo "   $(CCC) bin/ACS.o"
	@$(CCC) -o bin/ACS.o -c ACS.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/Kd.o: ../include/NeuronFitInfo.h ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h ../include/differential_equations.h Kd.cpp ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h ../include/Minimizable.h ../include/ParameterDescription.h ../include/Trace.h ../include/automatic_derivatives/Fad_construct.h OhmicBase.h ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h OhmicChannel.h
	@echo "   $(CCC) bin/Kd.o"
	@$(CCC) -o bin/Kd.o -c Kd.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/Nap.o: ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h ../include/NeuronFitInfo.h ../include/differential_equations.h ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h ../include/Minimizable.h ../include/ParameterDescription.h Nap.cpp ../include/Trace.h ../include/automatic_derivatives/Fad_construct.h OhmicBase.h ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h OhmicChannel.h
	@echo "   $(CCC) bin/Nap.o"
	@$(CCC) -o bin/Nap.o -c Nap.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/MI.o: MI.cpp ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h ../include/NeuronFitInfo.h ../include/differential_equations.h ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h ../include/Minimizable.h ../include/ParameterDescription.h ../include/Trace.h ../include/automatic_derivatives/Fad_construct.h OhmicBase.h ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h OhmicChannel.h
	@echo "   $(CCC) bin/MI.o"
	@$(CCC) -o bin/MI.o -c MI.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/NaCS.o: ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h ../include/NeuronFitInfo.h ../include/differential_equations.h ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h ../include/Minimizable.h OhmicBase.h ../include/ParameterDescription.h ../include/Trace.h ../include/constants.h NaCS.cpp ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h HodgkinHuxleyChannel.h ../include/automatic_derivatives/Fad_construct.h
	@echo "   $(CCC) bin/NaCS.o"
	@$(CCC) -o bin/NaCS.o -c NaCS.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/CaT.o: ../include/automatic_derivatives/Fad_functions.h CaT.cpp ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h ../include/NeuronFitInfo.h ../include/differential_equations.h ../include/io_functions.h NernstChannel.h ../include/NeuronModel.h ../include/Minimizable.h ../include/ParameterDescription.h ../include/Trace.h ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h ../include/automatic_derivatives/Fad_construct.h
	@echo "   $(CCC) bin/CaT.o"
	@$(CCC) -o bin/CaT.o -c CaT.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/As.o: ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h ../include/NeuronFitInfo.h ../include/differential_equations.h As.cpp ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h ../include/Minimizable.h ../include/ParameterDescription.h ../include/Trace.h ../include/automatic_derivatives/Fad_construct.h OhmicBase.h ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h OhmicChannel.h
	@echo "   $(CCC) bin/As.o"
	@$(CCC) -o bin/As.o -c As.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/CAN.o: ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h ../include/NeuronFitInfo.h ../include/differential_equations.h ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h ../include/Minimizable.h OhmicBase.h ../include/ParameterDescription.h ../include/Trace.h ../include/automatic_derivatives/Fad_construct.h CAN.cpp ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h OhmicChannel.h
	@echo "   $(CCC) bin/CAN.o"
	@$(CCC) -o bin/CAN.o -c CAN.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/Na.o: ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h ../include/NeuronFitInfo.h ../include/differential_equations.h ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h ../include/Minimizable.h OhmicBase.h ../include/ParameterDescription.h ../include/Trace.h ../include/automatic_derivatives/Fad_construct.h Na.cpp ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h OhmicChannel.h
	@echo "   $(CCC) bin/Na.o"
	@$(CCC) -o bin/Na.o -c Na.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/Leak.o: ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h ../include/NeuronFitInfo.h ../include/differential_equations.h ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h ../include/Minimizable.h Leak.cpp ../include/ParameterDescription.h ../include/Trace.h ../include/automatic_derivatives/Fad_construct.h OhmicBase.h ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h OhmicChannel.h
	@echo "   $(CCC) bin/Leak.o"
	@$(CCC) -o bin/Leak.o -c Leak.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC
bin/Af.o: ../include/NeuronFitInfo.h ../include/NeuronModel.h ../include/automatic_derivatives/Fad_operators.h ../include/automatic_derivatives/Fad.h ../include/differential_equations.h Af.cpp ../include/io_functions.h ../include/automatic_derivatives/Fad_functions.h ../include/Minimizable.h ../include/ParameterDescription.h ../include/Trace.h ../include/automatic_derivatives/Fad_construct.h OhmicBase.h ../include/NeuronGeometry.h ../include/NeuronException.h ../include/dynamic_class_loader.h ../include/CurrentInjector.h ../include/constants.h OhmicChannel.h
	@echo "   $(CCC) bin/Af.o"
	@$(CCC) -o bin/Af.o -c Af.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion -fPIC

# target files (executable or shared object):
bin/NaCS.so: bin/NaCS.o
	@echo "   $(CCC) bin/NaCS.so"
	@$(CCC) -o bin/NaCS.so bin/NaCS.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/CaT.so: bin/CaT.o
	@echo "   $(CCC) bin/CaT.so"
	@$(CCC) -o bin/CaT.so bin/CaT.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/Na.so: bin/Na.o
	@echo "   $(CCC) bin/Na.so"
	@$(CCC) -o bin/Na.so bin/Na.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/CompartmentConnection.so: bin/CompartmentConnection.o
	@echo "   $(CCC) bin/CompartmentConnection.so"
	@$(CCC) -o bin/CompartmentConnection.so bin/CompartmentConnection.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/KCa.so: bin/KCa.o
	@echo "   $(CCC) bin/KCa.so"
	@$(CCC) -o bin/KCa.so bin/KCa.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/CAN.so: bin/CAN.o
	@echo "   $(CCC) bin/CAN.so"
	@$(CCC) -o bin/CAN.so bin/CAN.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/MI.so: bin/MI.o
	@echo "   $(CCC) bin/MI.so"
	@$(CCC) -o bin/MI.so bin/MI.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/Nap.so: bin/Nap.o
	@echo "   $(CCC) bin/Nap.so"
	@$(CCC) -o bin/Nap.so bin/Nap.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/ACS.so: bin/ACS.o
	@echo "   $(CCC) bin/ACS.so"
	@$(CCC) -o bin/ACS.so bin/ACS.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/As.so: bin/As.o
	@echo "   $(CCC) bin/As.so"
	@$(CCC) -o bin/As.so bin/As.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/KdCS.so: bin/KdCS.o
	@echo "   $(CCC) bin/KdCS.so"
	@$(CCC) -o bin/KdCS.so bin/KdCS.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/Kd.so: bin/Kd.o
	@echo "   $(CCC) bin/Kd.so"
	@$(CCC) -o bin/Kd.so bin/Kd.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/CaS.so: bin/CaS.o
	@echo "   $(CCC) bin/CaS.so"
	@$(CCC) -o bin/CaS.so bin/CaS.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/Af.so: bin/Af.o
	@echo "   $(CCC) bin/Af.so"
	@$(CCC) -o bin/Af.so bin/Af.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/Leak.so: bin/Leak.o
	@echo "   $(CCC) bin/Leak.so"
	@$(CCC) -o bin/Leak.so bin/Leak.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/CaAccumulate.so: bin/CaAccumulate.o
	@echo "   $(CCC) bin/CaAccumulate.so"
	@$(CCC) -o bin/CaAccumulate.so bin/CaAccumulate.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl
bin/h.so: bin/h.o
	@echo "   $(CCC) bin/h.so"
	@$(CCC) -o bin/h.so bin/h.o -shared $(LDFLAGS) -rdynamic -fPIC -ldl

clean:
	@rm -f $(OBJS)
	@rm -f $(TARGET)
	@if [ -d bin ]; then rm -r bin ; fi
