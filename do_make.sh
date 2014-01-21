#!/bin/bash
#./python/CreateMakefile.py fitneuron.bin simulate_neuron.bin perturbneuron.bin -shared ions ions_AT ions_mismatch
./python/CreateMakefile.py simulate_neuron.bin fitneuron.bin perturbneuron.bin -shared injectors
make
