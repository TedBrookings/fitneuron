#!/bin/bash

# must pass number of processes as an argument
if [ -z "$1" ]; then
  echo "Must pass number of processes"
  exit -1
fi

# can optionally pass the population size
if [ -n "$2" ] ; then
  pSize="$2"
else
  pSize="214"
fi

mpirun -np "$1" ../../fitneuron/bin/fitneuron.bin -startup startFit.txt -accuracy 1.0 -populationSize "$pSize"
