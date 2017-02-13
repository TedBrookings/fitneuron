#!/usr/bin/python


_usageStr = \
  """usage: neuron_plot_perturb perturbFileName
           visualize perturbation landscape from file of perturbed neuron
           model"""



import sys, os, numpy

_useMatplotlib = True
_useMayavi = not _useMatplotlib
_useFakeData = False



###############################################################################
if _useMatplotlib:
  # import necessary modules
  import mpl_toolkits.mplot3d.axes3d as axes3d
  import matplotlib
  from matplotlib import pyplot
  from matplotlib import cm


  # update font properties for illustrator compatibility
  matplotlib.rcParams['pdf.fonttype'] = 42
  #matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
  #matplotlib.rc('font', **{'sans-serif' : 'Arial', 'family' : 'sans-serif'})


  
  
  def plotTrace(trace, labelSize=30, tickSize=24, loglog=False):
    """plot a trace using matplotlib"""
    mSize = 10.0
  
    fig = pyplot.figure()
    if loglog:
      pyplot.loglog(trace['x'], trace['value'], 'bo', ms = mSize)
    else:
      pyplot.semilogx(trace['x'], trace['value'], 'bo', ms = mSize)
    pyplot.title(trace['title'], fontsize=labelSize)
    pyplot.xlabel(trace['xName'], fontsize=labelSize)
    pyplot.ylabel('rms error (mV)', fontsize=labelSize)
    pyplot.xticks(fontsize=tickSize)
    pyplot.yticks(fontsize=tickSize)


  def plotSurface(surface, labelSize=30, tickSize=24):
    """plot a surface using matplotlib"""
    
    # create a figure
    fig = pyplot.figure()
    # create an axis
    axes = axes3d.Axes3D(fig)

    x = numpy.log10(surface['x'])
    y = numpy.log10(surface['y'])
    value = -numpy.log10(surface['value'])
    # create a surface plot
    # had cmap=cm.jet_r
    # hsv is better
    surfPlot = axes.plot_surface(x, y, value, \
                                 rstride=1, cstride=1, cmap=cm.gist_ncar, \
                                 linewidth=0, antialiased=False, edgecolors='none')
    axes.set_title(surface['title'], fontsize=labelSize)
    axes.set_xlabel(surface['xName'], fontsize=labelSize)
    axes.set_ylabel(surface['yName'], fontsize=labelSize)
    axes.set_zlabel('rms error (mV)', fontsize=labelSize)
    axes.view_init(45, 255)
    
    def _unLog(_t):
      """
      Convert _t back to 10^t, keeping 2 digits of precision
      """
      if t >= 1:
        if t > 2:
          if t >= 5:
            return '%.2g' % 10**t
          else:
            return '%d' % (100 * round(10**(t-2)))
        else:
          return '%d' % round(10**t)
      else:
        if t < 0:
          return '%.2f' % 10**t
        else:
          return '%.2f' % 10**t
    
    axes.set_xticklabels([_unLog(t) for t in axes.get_xticks()], \
                         fontsize=tickSize)
    axes.set_yticklabels([_unLog(t) for t in axes.get_yticks()], \
                         fontsize=tickSize)
    axes.set_zticklabels([_unLog(-t) for t in reversed(axes.get_zticks())], \
                         fontsize=tickSize)
    
    # add a colorbar
    #fig.colorbar(surfPlot, shrink = 0.5, aspect = 5)


  
  def showFigures():
    """after all figures are created, show them and wait"""
    pyplot.show()



if _useMayavi:
  # import necessary modules
  from enthought.mayavi import mlab
  
  def plotTrace(trace):
    # make a blank figure
    mlab.figure()
  
  
  def plotSurface(surface):
    """plot a surface using mayavi"""
    
    # create a figure
    mlab.figure()
    # create a surface plot
    surfPlot = mlab.mesh(surface['x'], surface['y'], surface['value'])
  
  
  def showFigures():
    """after all figures are created, show them and wait"""
    mlab.show()



###############################################################################
def _makeFakeData():
  """make a couple surfaces to test plotting routines"""
  
  # choose the x and y points that will be evaluated
  x = numpy.arange(-numpy.pi, numpy.pi, 0.05 * numpy.pi)
  y = numpy.arange(-numpy.pi, numpy.pi, 0.05 * numpy.pi)
  
  # form x and y into a grid
  (x, y) = numpy.meshgrid(x, y)
  # compute the radius of each x,y pair
  r = numpy.sqrt(x**2 + y**2)
  
  # make the first surface
  z1 = numpy.sin(r**2) / (1 + r**2)
  # make the second surface
  z2 = r
  
  # put into a list
  surfaces = [\
    {'x' : x, 'y' : y, 'value' : z1, \
     'title' : 'ripple', 'xName' : 'x', 'yName' : 'y'},
    {'x' : x, 'y' : y, 'value' : z2, \
     'title' : 'cone', 'xName' : 'x', 'yName' : 'y'}]
  
  return surfaces



###############################################################################
def _getNextLine(fIn):
  """get the next line from the file, comments removed, split into words"""
  
  splitLine = []
  while not splitLine:
    # get the next line from the file
    line = fIn.next()
    
    # remove comments and endline
    line = line.split('#', 1)[0].strip('\n')
  
    # split line into words separated by white space      
    splitLine = line.split(None)
  return splitLine



###############################################################################
def loadPerturbInfo(perturbFile):
  with open(perturbFile, 'r') as fIn:
    # get the parameters that are perturbed
    perturbParams = _getNextLine(fIn)
    # make a list of indices to perturbed params, set to NaN for now
    perturbInds = [float('NaN') for n in range(len(perturbParams))]
    
    # get the number of total parameters
    numParameters = int(_getNextLine(fIn)[1])
    
    # get the value of the unperturbed parameters
    splitLine = _getNextLine(fIn)
    unperturbedValue = float(splitLine[1])
    
    # read in the parameter descriptions (starting with the value)
    parameterNames = []
    parameterVals = []
    for n in range(numParameters):
      # get the description
      splitLine = _getNextLine(fIn)
      name = splitLine[0]
      parameterNames.append(name)
      parameterVals.append(float(splitLine[1]))
      
      if name in perturbParams:
        # this is one of the parameters that is perturbed
        perturbInds[perturbParams.index(name)] = n
    
    # get the number of perturbed parameter sets
    numParamSets = int(_getNextLine(fIn)[1])
    
    # get the perturbed parameter sets
    values = []
    perturbedParamSets = []
    for n in range(numParamSets):
      splitLine = _getNextLine(fIn)
      values.append(float(splitLine[0]))
      perturbedParamSets.append([float(x) for x in splitLine[1:]])
  
  perturbInfo = { \
    'numPerturbed' : len(perturbParams), \
    'perturbParams' : perturbParams, \
    'perturbInds' : perturbInds, \
    'unperturbedValue' : unperturbedValue, \
    'numParameters' : numParameters, \
    'parameterNames' : parameterNames, \
    'parameterVals' : parameterVals, \
    'numParamSets' : numParamSets, \
    'values' : values, \
    'perturbedParamSets' : perturbedParamSets \
  }
  
  return perturbInfo



###############################################################################
def getPerturbTraces(perturbInfo):
  """load in the data from a perturb file and turn it into surfaces"""
  
  perturbTraces = []
  perturbInds = perturbInfo['perturbInds']
  for n in range(perturbInfo['numPerturbed']):
    # loop over each perturb param and get its trace
    
    # get the parameter index of the parameter that's changing in this trace
    ind_n = perturbInds[n]
    
    # make a list of what parameters aren't changing in this trace
    constInds = [ind for ind in perturbInds if ind != ind_n]
        
    # create a dictionary obj to hold trace info
    paramName = perturbInfo['perturbParams'][n]
    trace = {'x' : [], 'value' : [], \
             'title' : 'Perturbing ' + paramName, 'xName' : paramName}

    # loop over population of perturbed parameter sets, adding appropriate
    # parameter sets to trace
    for m in range(perturbInfo['numParamSets']):
      paramSet = perturbInfo['perturbedParamSets'][m]
      if wantedSet(paramSet, constInds, perturbInfo):
        trace['x'].append(paramSet[ind_n])
        trace['value'].append(perturbInfo['values'][m])
    
    perturbTraces.append(trace)
  return perturbTraces



###############################################################################
def getPerturbSurfaces(perturbInfo):
  """load in the data from a perturb file and turn it into surfaces"""
  
  # get all the pairs of perturbed parameters
  pairs = []
  for n in range(perturbInfo['numPerturbed']):
    for m in range(n + 1, perturbInfo['numPerturbed']):
      pairs.append((n, m))
  
  # loop over all possible pairs
  perturbInds = perturbInfo['perturbInds']
  perturbSurfaces = []
  for pair in pairs:
    # create an error landscape for this pair of perturbed parameters
  
    # get the parameter index of the parameters that are changing
    ind_x = perturbInds[pair[0]]
    ind_y = perturbInds[pair[1]]
    
    # make a list of what parameters aren't changing in this trace
    constInds = [ind for ind in perturbInds if ind not in [ind_x, ind_y]]
        
    # create a dictionary obj to hold trace info
    xName = perturbInfo['perturbParams'][pair[0]]
    yName = perturbInfo['perturbParams'][pair[1]]
    surface = {'x' : [], 'y' : [], 'value' : [], \
             'title' : 'Perturbing ' + xName + ' vs ' + yName, \
             'xName' : xName, 'yName' : yName}

    # loop over population of perturbed parameter sets, adding appropriate
    # parameter sets to trace
    for m in range(perturbInfo['numParamSets']):
      paramSet = perturbInfo['perturbedParamSets'][m]
      if wantedSet(paramSet, constInds, perturbInfo):
        surface['x'].append(paramSet[ind_x])
        surface['y'].append(paramSet[ind_y])
        surface['value'].append(perturbInfo['values'][m])
    
    # reorder/reform the x/y/values arrays to conform to the format necessary
    # for plotting (produced by numpy.meshgrid)
    reorderSurface(surface)
    
    perturbSurfaces.append(surface)
  
  return perturbSurfaces



###############################################################################
def wantedSet(paramSet, constInds, perturbInfo):
  """return True if all the constInds are unchanged from their base values,
     False otherwise"""
  
  for n in range(len(constInds)):
    ind = constInds[n]
    if paramSet[ind] != perturbInfo['parameterVals'][ind]:
      return False
  
  return True



###############################################################################
def reorderSurface(surface):
  """reorder/reform the x/y/values arrays to conform to the format necessary
     for plotting (produced by numpy.meshgrid)"""
  
  oldX = surface['x'][:]
  uniqueX = list(set(oldX))
  uniqueX.sort()
  numX = len(uniqueX)
  
  oldY = surface['y'][:]
  uniqueY = list(set(oldY))
  uniqueY.sort()
  numY = len(uniqueY)
  
  oldZ = surface['value'][:]
  
  shape = (numY, numX)
  meshX = numpy.ndarray(shape)
  meshY = numpy.ndarray(shape)
  meshZ = numpy.ndarray(shape)
  meshZ[:] = 1.0e10
  
  for (x_n, y_n, z_n) in zip(oldX, oldY, oldZ):
    xInd = uniqueX.index(x_n)
    yInd = uniqueY.index(y_n)
    if meshZ[yInd, xInd] != 1.0e10:
      print('Duplicate of %g, %g' % (y_n, x_n))
      print('indices: %d, %d' % (yInd, xInd))
    
    meshX[yInd, xInd] = x_n
    meshY[yInd, xInd] = y_n
    meshZ[yInd, xInd] = z_n
  
  for yInd in range(numY):
    for xInd in range(numX):
      if meshZ[yInd, xInd] == 1.0e10:
        print('Missing %g, %g' % (uniqueY[yInd], uniqueX[xInd]))
        print('Indices: %d, %d' % (yInd, xInd))
      
  
  surface['x'] = meshX
  surface['y'] = meshY
  surface['value'] = meshZ



###############################################################################
def _parseArguments():
  arguments = sys.argv
  
  if len(arguments) != 2:
    print(_usageStr)
    if len(arguments) > 1:
      raise TypeError('Incorrect number of arguments.')
    sys.exit(0)
  
  perturbFile = arguments[1]
  return perturbFile



###############################################################################
if __name__ == "__main__":
  # load in data (or make fake data if only testing
  if _useFakeData:
    surfaces = _makeFakeData()
  else:
    perturbFile = _parseArguments()
    perturbInfo = loadPerturbInfo(perturbFile)
    perturbTraces = getPerturbTraces(perturbInfo)
    perturbSurfaces = getPerturbSurfaces(perturbInfo)
  
  # plot each perturb trace
  #for trace in perturbTraces:
  #  plotTrace(trace)

  # plot each perturb surface
  for surface in perturbSurfaces:
    plotSurface(surface)
  
  # show the figures
  showFigures()
  
  sys.exit(0)
