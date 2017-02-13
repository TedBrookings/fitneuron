#!/bin/bash

# can pass accuracy as an argument
if [ -n "$1" ]; then
  ../../fitneuron/bin/simulate_neuron.bin -startup startGenModelTrace.txt -accuracy "$1"
else
  ../../fitneuron/bin/simulate_neuron.bin -startup startGenModelTrace.txt -accuracy 0.01
fi
  ../../scripts/neuron_plot_trace.py ModelTrace.txt simulated_traces.txt
