#!/usr/bin/python

import os, sys

program = '../fitneuron/bin/perturbneuron.bin'
plotProgram = '../scripts/neuron_plot_perturb.py'
startfile = 'startPerturb.txt'
modelParams = 'ModelParams.txt'


def visualizeErrorLandscape(perturbArgs):
  perturbFile = 'perturb_' + '_vs_'.join(perturbArgs) + '.txt'

  # first ensure that the perturb file exists
  if not os.path.isfile(perturbFile):
    # run perturbneuron.bin to create perturb file
    varyParams = ' '.join(perturbArgs)
    sysCmd = 'mpirun -np 8 %s %s %s %s' % \
             (program, startfile, modelParams, varyParams)
    
    os.system(sysCmd)
  
  # perturb file should exist now, if there was no error
  if not os.path.isfile(perturbFile):
    raise RuntimeError("Could not create perturb file")

  # plot perturb
  sysCmd = '%s %s' % (plotProgram, perturbFile)
  os.system(sysCmd)


def _parseArguments():
  if len(sys.argv) != 3:
    sys.tracebacklimit=0
    raise ValueError('Must pass the name of two parameters from the model')
  
  param1 = sys.argv[1]
  param2 = sys.argv[2]
  with open(modelParams, 'r') as fIn:
    params = [line.split()[0] for line in fIn]
  if param1 not in params:
    sys.tracebacklimit=0
    raise ValueError('%s is not a valid parameter of the model' % param1)
  if param2 not in params:
    sys.tracebacklimit=0
    raise ValueError('%s is not a valid parameter of the model' % param2)
  
  return sys.argv[1:]


if __name__ == '__main__':
  perturbArgs = _parseArguments()
  visualizeErrorLandscape(perturbArgs)
