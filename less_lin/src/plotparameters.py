params = {
 	  #'text.usetex': True,
 	  #'figure.figsize': [16, 10], # instead of 4.5, 4.5
 	  'figure.figsize': [8, 5], # instead of 4.5, 4.5
	  #'text.latex.preamble': ['\usepackage{gensymb}'],
	  'font.family' : 'serif',
	  'font.sans-serif'    : 'Helvetica',#, Avant Garde, Computer Modern Sans serif
          'font.size' : 15,
          #dont work with latex
	  #'font.weight' : 'normal',
	  'font.weight' : 500,
	  'lines.linewidth':2,
	 #note that font.size controls default text sizes.  To configure
	# special text sizes tick labels, axes, labels, title, etc, see the rc
	# settings for axes and ticks. Special text sizes can be defined
	# relative to font.size, using the following values: xx-small, x-small,
	# small, medium, large, x-large, xx-large, larger, or smaller
 	  #'legend.fontsize': 'large',
          'legend.fontsize': 18,
          'legend.handlelength': 2,
 	  #'axes.labelsize' : 'medium',
 	  'axes.labelweight' :'normal',
 	  'axes.labelweight' :'500',
 	  'xtick.labelsize': 'large',
 	  'ytick.labelsize': 'large',
          #'xtick.major.size': 4,      # major tick size in points
          #'xtick.major.size': 4,      # major tick size in points
          #'xtick.minor.size': 4,   # minor tick size in points
          #'xtick.major.pad': 4,      # distance to major tick label in points
          #'xtick.minor.pad': 4,      # distance to the minor tick label in points
          #'xtick.color': 'k',      # color of the tick labels
          'xtick.direction': 'out',     # direction: in or out
	  #'mathtext.bf'  : 'serif:bold'
        }

#def set_size(w,h, ax=None):
#
#""" w, h: width, height in inches """
#
#if not ax: ax=plt.gca()
#    l = ax.figure.subplotpars.left
#    r = ax.figure.subplotpars.right
#    t = ax.figure.subplotpars.top
#    b = ax.figure.subplotpars.bottom
#    figw = float(w)/(r-l)
#    figh = float(h)/(t-b)
#    ax.figure.set_size_inches(figw, figh)
#
#    fig, ax=plt.subplots()
#
#    ax.plot([1,3,2])
#
#    set_size(8,5)
