#!/usr/bin/python


import matplotlib
from matplotlib import pyplot
import scipy



_usageStr = """
Not meant to be called from command line, but imported and used
"""



# update font properties for illustrator compatibility
matplotlib.rcParams['pdf.fonttype'] = 42



###############################################################################
def plotXY(x, y, markerStyle, color=None, xLabel='x', yLabel='y', title=None,
           labelSize=30, tickSize=24, titleSize=None, figure=None,
           legendLabel=None, xScale=None, yScale=None, **kwargs):
           
  """
  Plot typical x vs y plot with requested options.
  Return figure object
  """
  if figure is None:
    fig = pyplot.figure()
    axes = fig.gca()
  else:
    fig = figure
    if fig.axes:
      axes = fig.axes[0]
      pyplot.sca(axes)
    else:
      #axes = fig.gca()
      axes = fig.add_subplot(111)
  
  
  if legendLabel is None:
    legendLabel = yLabel
  
  if color is None:
    axes.plot(x, y, markerStyle, label=legendLabel, **kwargs)
  else:
    axes.plot(x, y, markerStyle, color=color, label=legendLabel, **kwargs)

  if xScale is not None:
    axes.set_xscale(xScale)
  if yScale is not None:
    axes.set_yscale(yScale)
  
  override = {
   'fontsize'            : labelSize,
   'verticalalignment'   : 'center',
   'horizontalalignment' : 'center',
   'rotation'            : 'vertical'}
  if not axes.get_ylabel():
    override['horizontalalignment'] = 'right'
    pyplot.ylabel(yLabel, override)
    override['horizontalalignment'] = 'center'
  override['rotation'] = 'horizontal'
  if not axes.get_xlabel():
    pyplot.xlabel(xLabel, override)
  
  if not axes.get_title():
    if title is None:
      title = '%s vs %s' % (xLabel, yLabel)
    if titleSize is None:
      titleSize = labelSize

    override['fontsize'] = titleSize
    override['verticalalignment'] = 'bottom'
    pyplot.title(title, override)
  
  pyplot.setp(axes.get_xticklabels(), rotation='horizontal', fontsize=tickSize)
  pyplot.setp(axes.get_yticklabels(), rotation='horizontal', fontsize=tickSize)
  
  return fig



if __name__ == "__main__":
  print(_usageStr)
  # maybe put in a demo?
  
