
def figure_format(fig_width='3.4',fig_height='2.1'):
    """
    Parameter
    ---------
    fig_width : float
        figure width in inches
        e.g. 3.4
        
    fig_height : float
        figure height in inches
        e.g. 3.4
   
    Returns
    -------
    A tuple with:
    - Float corresponding to figure width and height
    - A dictionary containing several parameters, such as fontsize
      , etc. If fig_width is not set, the default value is 3.4 inch, 
      corresponding to the width of a column in a double colomn paper.
    """
    golden_ratio = 1.618      # Aesthetic ratio
    fig_size     = [fig_width,fig_height]
    fontsize     = 15
    linewidth    = 1.5
    params = {'backend': 'ps',
              'axes.labelsize': fontsize,
              'axes.titlesize': fontsize,
              'font.size': fontsize,
              'legend.frameon': False,
              'legend.fontsize': fontsize,
              'legend.loc': 'best',
              'lines.linewidth': linewidth,
              'xtick.labelsize': fontsize,
              'ytick.labelsize': fontsize,
              'xtick.direction': 'in',
              'ytick.direction': 'in',
              'xtick.top': True,
              'ytick.right': True,
              'xtick.major.size': linewidth*4,
              'xtick.major.top': True,
              'xtick.major.width': linewidth,
              'xtick.minor.size': linewidth*2,
              'xtick.minor.top': True,
              'xtick.minor.width': linewidth,
              'ytick.major.size': linewidth*4,
              'ytick.major.width': linewidth,
              'ytick.minor.size': linewidth*2,
              'ytick.minor.width': linewidth,
              'figure.figsize': fig_size,
		  'pgf.texsystem': 'pdflatex',
        	  'font.family': 'serif',
        	  'text.usetex': True,
        	  'pgf.rcfonts': False,
        	  'ps.usedistiller': 'xpdf'}
    return(fig_width,fig_height,params)

