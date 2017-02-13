#!/usr/bin/python
# -*- coding: utf-8 -*-



_usageStr=\
"""usage: make_control_demo_figure.py
         draws control demo
"""



import sys, os
import scipy
from scipy import linalg
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.patches as patches
from plotXY import *



# Make output fonts compatible with Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42

# period of target data
_target_period = 1.0
# sample time for target (correct) data
_target_dt = _target_period / 16.0
# final simulation time for close views
_t_final_short = _target_period * 1.0 + _target_dt
# final simulation time for long views
_t_final_long = 5 * _target_period
# timescale for control
_control_tau = _target_dt * 2.0
# sample time for fit (candidate) model
_test_dt = _target_dt / 3.0
# number of mini steps to take for each recorded step
_iterate_N = 100
# number of candidate frequencies to test
_num_frequencies = 101


###############################################################################
def simulateTarget(dt, omega, tFinal=_t_final_long):
  # create the list of times
  targetN = int(1 + round(tFinal / dt))
  target_t = scipy.array([dt * n for n in range(targetN)])
  # create the list of values
  target_f = scipy.ndarray((targetN,))
  
  # set the initial state
  target_state = scipy.array([[1.0], [0.0]])
  # record initial value
  target_f[0] = target_state[0][0]
  
  # while simulation is not complete, iterate
  cosWT = scipy.cos(omega * dt)
  sinWT = scipy.sin(omega * dt)
  iterMat = scipy.array([[cosWT, sinWT], [-sinWT, cosWT]])
  for n in range(1, targetN):
    target_state = iterMat.dot(target_state)
    target_f[n] = target_state[0][0]
  
  return (target_t, target_f)


###############################################################################
def simulateFit(dt, testDt, omega, target_t, target_f, tFinal=_t_final_long):
  # create the list of times
  fitN = int(1 + round(tFinal / testDt))
  targetN = int(1 + round(tFinal / dt))
  nPoints = fitN + targetN - 1# -1 because the first point will not be adjusted
  # create the list of values
  controlled_t = scipy.ndarray((targetN,))
  controlled_f = scipy.ndarray((targetN,))
  discarded_t = scipy.ndarray((fitN-1,))
  discarded_f = scipy.ndarray((fitN-1,))
  fit_t = scipy.ndarray((nPoints,))
  fit_f = scipy.ndarray((nPoints,))
  
  # set the initial state
  t = 0.0
  fitInd = 0
  targetInd = 1
  discardInd = 0
  fit_state = scipy.array([[1.0], [0.0]])
  # record initial value
  fit_t[fitInd] = t
  fit_f[fitInd] = fit_state[0][0]
  controlled_t[0] = t
  controlled_f[0] = fit_state[0][0]
  
  # while simulation is not complete, iterate
  cosWT = scipy.cos(omega * testDt)
  sinWT = scipy.sin(omega * testDt)
  iterMat = scipy.array([[cosWT, sinWT], [-sinWT, cosWT]])
  controlFact = 1 - scipy.exp(-dt / _control_tau)
  while targetInd < targetN:
    fit_state = iterMat.dot(fit_state)
    fitInd += 1
    t += testDt
    fit_t[fitInd] = t
    fit_f[fitInd] = fit_state[0][0]
    
    if t + 1.0e-10 >= target_t[targetInd]:
      discarded_t[discardInd] = t
      discarded_f[discardInd] = fit_state[0][0]

      fitError = fit_state[0][0] - target_f[targetInd]
      fit_state[0][0] -= controlFact * fitError
      controlled_t[targetInd] = t
      controlled_f[targetInd] = fit_state[0][0]
      
      targetInd += 1
      fitInd += 1
      discardInd += 1
      
      fit_t[fitInd] = t
      fit_f[fitInd] = fit_state[0][0]
    else:
      discarded_t[discardInd] = t
      discarded_f[discardInd] = fit_state[0][0]
      discardInd += 1
  
  return (controlled_t, controlled_f, discarded_t, discarded_f, fit_t, fit_f)



###############################################################################
def plotProcess(target_t, target_f, controlled_t, controlled_f, \
                discarded_t, discarded_f, xMin, xMax, \
                titleStr="Process of Control Adjustments", xLabel="Time", \
                yLabel="y", markerSize=10.0, alphaLevel=0.5):
  
  # Make a plot that demonstrates the control adjustment process

  # Calculate a few things first
  # 0. Slightly broaden the plot window
  xMin -= 0.5 * _test_dt
  xMax += 0.5 * _test_dt
  
  # 1. Calculate the indices corresponding to target/controlled points that
  #    are within the plot window
  tInd = [n for n in range(len(target_t)) if target_t[n] >= xMin and \
                                             target_t[n] <= xMax]
  # 2. Calculate the indices corresponding to discarded (working) points
  dInd = [n for n in range(len(discarded_t)) if discarded_t[n] >= xMin and \
                                                discarded_t[n] <= xMax]
  # 3. Compute the size of adjustments during this window. This is necessary
  #    to choose an appropriate size for arrows on the plot
  deltaYVec = []
  yBase = []
  for n in tInd:
    t = controlled_t[n]
    yControl = controlled_f[n]
    discardInd = \
      min(((abs(discarded_t[m]-t), m) for m in range(len(discarded_t))),
          key = lambda z:z[0])[1]
    yDiscarded = discarded_f[discardInd]
    yBase.append(yDiscarded)
    deltaYVec.append(yControl - yDiscarded)
  
  arrowHeadLen = 0.2 * abs(max(deltaYVec, key=lambda x: abs(x)))
  arrowHeadWidth = 0.5 * arrowHeadLen
  
  # Make the figures and plot
  fig = pyplot.figure()
  
  plotXY(target_t[tInd], target_f[tInd], 'v', color='r', \
         markersize=markerSize, markeredgecolor='none', xLabel=xLabel, \
         yLabel=yLabel, title=titleStr, legendLabel="target", figure=fig)
  plotXY(controlled_t[tInd], controlled_f[tInd], 'o', color='g', \
         alpha=alphaLevel, markersize=markerSize, markeredgecolor='none', \
         xLabel=xLabel, yLabel=yLabel, title=titleStr, \
         legendLabel="controlled", figure=fig)
  plotXY(discarded_t[dInd], discarded_f[dInd], 'o', color='k', \
         alpha=alphaLevel, markersize=markerSize, markeredgecolor='none', \
         xLabel=xLabel, yLabel=yLabel, title=titleStr, \
         legendLabel="working values", figure=fig)

  # draw arrows on the plot annotating the action of the control
  
  for (n, yDiscarded, deltaY) in zip(tInd, yBase, deltaYVec):
    # shrink beginning and end of arrow
    yDiscarded += 0.05 * deltaY;
    deltaY *= 0.9
    # plot arrow
    pyplot.arrow(controlled_t[n], yDiscarded, 0, deltaY, \
                 length_includes_head=True, \
                 head_width=arrowHeadWidth, head_length=arrowHeadLen,
                 facecolor='k', edgecolor='none')

  pyplot.xlim([xMin, xMax])
  yMin = min([min(target_f[tInd]), min(discarded_f[dInd])])
  yMax = max([max(target_f[tInd]), max(discarded_f[dInd])])
  yBoarder = 0.05 * (yMax - yMin)
  yMin -= yBoarder
  yMax += yBoarder
  pyplot.ylim([yMin, yMax])
  
  pyplot.legend(loc=0)
  return (yMin, yMax)


###############################################################################
def plotSimErrors(target_t, target_f, fit_raw_t, fit_raw_f, \
                  controlled_t, controlled_f, xMin, xMax, \
                  titleStr="Error With and Without Control", xLabel="Time", \
                  yLabel="y", markerSize=10.0, lineWidth=6.0, alphaLevel=0.5):
  # Calculate the indices corresponding to target/controlled points that are
  # within the plot window
  xMin -= 0.5 * _test_dt
  xMax += 0.5 * _test_dt
  tInd = [n for n in range(len(target_t)) if target_t[n] >= xMin and \
                                             target_t[n] <= xMax]
  
  fig = pyplot.figure()
  
  plotXY(target_t[tInd], target_f[tInd], 'v', color='r', \
         markersize=markerSize, markeredgecolor='none', xLabel=xLabel, \
         yLabel=yLabel, title=titleStr, \
         legendLabel="target", figure=fig)
  plotXY(fit_raw_t[tInd], fit_raw_f[tInd], 'd', color='b', alpha=alphaLevel, \
         markersize=markerSize, markeredgecolor='none', xLabel=xLabel, \
         yLabel=yLabel, title=titleStr, legendLabel="uncontrolled", \
         figure=fig)
  plotXY(controlled_t[tInd], controlled_f[tInd], 'o', color='g', \
         alpha=alphaLevel, markersize=markerSize, markeredgecolor='none', \
         xLabel=xLabel, yLabel=yLabel, title=titleStr, \
         legendLabel="controlled", figure=fig)
  
  eps = 0.1 * _test_dt
  ax = pyplot.gca()
  for n in tInd:
    t = target_t[n]
    h = fit_raw_f[n] - target_f[n]
    ax.add_patch(\
      patches.Rectangle((t - eps, target_f[n]), -eps, h, facecolor='b', \
                        alpha=0.5*alphaLevel, ec='none')
                )
    h = controlled_f[n] - target_f[n]
    ax.add_patch(\
      patches.Rectangle((t + eps, target_f[n]), eps, h, facecolor='g', \
                        alpha=0.5*alphaLevel, ec='none')
                )
  
  pyplot.legend(loc=0)
  pyplot.xlim([xMin, xMax])
  yMin = min([min(target_f[tInd]), min(fit_raw_f[tInd])])
  yMax = max([max(target_f[tInd]), max(fit_raw_f[tInd])])
  yBoarder = 0.05 * (yMax - yMin)
  yMin -= yBoarder
  yMax += yBoarder
  pyplot.ylim([yMin, yMax])


###############################################################################
def plotSims(target_t, target_f, fit_raw_t, fit_raw_f, \
             controlled_t, controlled_f, xMin, xMax, yMin, yMax, \
             titleStr="Control Promotes Alignment", xLabel="Time", \
             yLabel="y", markerSize=10.0, alphaLevel=0.5):
  
  fig = pyplot.figure()
  
  plotXY(target_t, target_f, 'v', color='r', markersize=markerSize, \
         markeredgecolor='none', xLabel=xLabel, yLabel=yLabel, title=titleStr,\
         legendLabel="target", figure=fig)
  plotXY(fit_raw_t, fit_raw_f, 'd', color='b', alpha=alphaLevel, \
         markersize=markerSize, markeredgecolor='none', xLabel=xLabel, \
         yLabel=yLabel, title=titleStr, legendLabel="uncontrolled", figure=fig)
  plotXY(controlled_t, controlled_f, 'o', color='g', alpha=alphaLevel, \
         markersize=markerSize, markeredgecolor='none', xLabel=xLabel, \
         yLabel=yLabel, title=titleStr, legendLabel="controlled", \
         figure=fig)
  xMin -= 0.5 * _test_dt
  xMax += 0.5 * _test_dt
  pyplot.plot([xMin, xMax, xMax, xMin, xMin], [yMin, yMin, yMax, yMax, yMin], \
              'k--', linewidth=4.0)

  pyplot.legend(loc=0)
  pyplot.xlim([-0.5 * _test_dt, 0.5 * _test_dt + target_t[-1]])
  pyplot.ylim([-1.05, 1.05])



###############################################################################
def make_demo_figure():
  omegaTarget = 2 * scipy.pi / _target_period
  numCycle = _t_final_long / _target_period
  omegaFit = omegaTarget * (1.0 + 0.5 / numCycle)
  
  # make a figure to demonstrate the control adjustment process:
  tFinal = _t_final_short  # short time
  dt = _target_dt
  testDt = _test_dt
  
  # get synthetic target data
  (target_t, target_f) = simulateTarget(dt, omegaTarget, tFinal)
  # get fit data with no control
  (fit_raw_t, fit_raw_f) = simulateTarget(dt, omegaFit, tFinal)
  # get fit data with control
  (controlled_t, controlled_f, discarded_t, discarded_f, fit_t, fit_f) = \
    simulateFit(dt, testDt, omegaFit, target_t, target_f, tFinal)

  # plot target and fit
  xMax = _t_final_short
  xMin = xMax - 3.0 * _target_dt
  yMin, yMax = plotProcess(target_t, target_f, controlled_t, controlled_f, \
                           discarded_t, discarded_f, xMin, xMax)


  # visualize errors with and without control
  plotSimErrors(target_t, target_f, fit_raw_t, fit_raw_f, \
                controlled_t, controlled_f, xMin, xMax)


  # make a figure to demonstrate the results of the control adjustment
  tFinal = _t_final_long  # long time
  dt = 0.25 * _target_dt # short dt to connect the curves better
  testDt = dt
  
  # get synthetic target data
  (target_t, target_f) = simulateTarget(dt, omegaTarget, tFinal)
  # get fit data with no control
  (fit_raw_t, fit_raw_f) = simulateTarget(dt, omegaFit, tFinal)
  # get fit data with control
  (controlled_t, controlled_f, discarded_t, discarded_f, fit_t, fit_f) = \
    simulateFit(dt, testDt, omegaFit, target_t, target_f, tFinal)

  # plot target and fit
  plotSims(target_t, target_f, fit_raw_t, fit_raw_f, \
           controlled_t, controlled_f, xMin, xMax, yMin, yMax)



###############################################################################
def plotErrors(frequencies, errUncontrol, errControl, \
               titleStr="Error vs Frequency", \
               xLabel="Test frequency", yLabel="RMS Error", \
               markerSize=10.0, alphaLevel=0.5):
  fig = pyplot.figure()
  
  plotXY(frequencies, errUncontrol, 'd', color='b', markersize=markerSize, \
         markeredgecolor='none', xLabel=xLabel, yLabel=yLabel, title=titleStr,\
         legendLabel="uncontrolled", figure=fig, xScale='log', \
         alpha=alphaLevel)
  plotXY(frequencies, errControl, 'o', color='g', markersize=markerSize, \
         markeredgecolor='none', xLabel=xLabel, yLabel=yLabel, title=titleStr,\
         legendLabel="controlled", figure=fig, xScale='log', alpha=alphaLevel)
  
  pyplot.legend(loc=0)
  #dF = scipy.log(frequencies[-1] / frequencies[0]) / len(freqencies);
  dF = (frequencies[-1] / frequencies[0])**(1.0/len(frequencies))
  
  pyplot.xlim([frequencies[0] / dF, frequencies[-1] * dF])
  pyplot.ylim([-0.05, 0.05 + max(errUncontrol)])
  
  a = pyplot.gca()
  a.set_xticks([0.5, 1.0, 2.0])
  a.set_xticklabels(["0.5", "1.0", "2.0"])


###############################################################################
def rmsError(f1, f2):
  if len(f1) != len(f2):
    raise RuntimeError("Trace lengths don't match")
  mse = sum((yn-xn)**2 for (yn,xn) in zip(f1, f2))
  return (mse / len(f1))**0.5


###############################################################################
def make_error_figure():
  # make a figure showing the effect of control adjustment on RMS error of
  # oscillators with the incorrect frequency

  # get synthetic target data
  tFinal = _t_final_long
  omegaTarget = 2 * scipy.pi / _target_period
  dt = _target_dt
  
  (target_t, target_f) = simulateTarget(_target_dt, omegaTarget, tFinal=tFinal)
  targetN = int(1 + round(_t_final_short / _target_dt))
  
  # make test omegas
  frequencies = \
    [2**x / _target_period for x in scipy.linspace(-1, 1, _num_frequencies)]
  omegas = [2 * scipy.pi * f for f in frequencies]
  
  numOmegas = len(omegas)
  errUncontrol = scipy.ndarray((numOmegas,))
  errControl = scipy.ndarray((numOmegas,))
  
  ind = 0
  for omega in omegas:
    (fit_raw_t, fit_raw_f) = simulateTarget(_target_dt, omega, tFinal=tFinal)
    errUncontrol[ind] = rmsError(fit_raw_f, target_f)
    
    # get fit data with control
    (controlled_t, controlled_f, discarded_t, discarded_f, fit_t, fit_f) = \
      simulateFit(dt, dt, omega, target_t, target_f, tFinal=tFinal)
    errControl[ind] = rmsError(controlled_f, target_f)
    ind += 1
  
  plotErrors(frequencies, errUncontrol, errControl)


###############################################################################
def _parseArguments():
  arguments = sys.argv
  if len(arguments) >  1:
    print(_usageStr)
    sys.tracebacklimit = 1
    raise TypeError('Incorrect number of arguments.')

  return None


  
###############################################################################
if __name__ == "__main__":
  _parseArguments()
  
  make_demo_figure()
  make_error_figure()
  
  # wait until figures are closed
  pyplot.show()

  sys.exit(0)
