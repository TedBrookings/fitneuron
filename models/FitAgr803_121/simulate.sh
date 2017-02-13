#!/bin/bash

# can pass accuracy as an argument
if [ -n "$1" ]; then
  acc="$1"
else
  acc="0.1"
fi
../../fitneuron/bin/simulate_neuron.bin -startup startSimResults.txt -params results.txt -accuracy "${acc}"
../../scripts/neuron_plot_trace.py PinkNoise.txt simulated_traces.txt
