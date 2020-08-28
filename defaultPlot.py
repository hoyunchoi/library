import matplotlib
import matplotlib.pyplot as plt

new_rc_params = {
    # Default parameters for plot
    #* text
    'text.usetex' : True,
    'font.size' : 50,

    #* figure
    'figure.figsize' : (10,10),
    'axes.linewidth' : 3,

    #* line
    'lines.linewidth' : 4,
    'lines.markersize' : 15,

    #* legend
    'legend.frameon' : False,
    'legend.fontsize' : 40,
    'legend.handlelength' : 1.0,
    'legend.handletextpad' : 0.5,
    'legend.labelspacing' : 0.2,

    #* label
    'axes.labelsize' : 50,

    #* tick
    'xtick.direction' : 'in',
    'ytick.direction' : 'in',
    'xtick.major.width' : 3,
    'ytick.major.width' : 3,
    'xtick.major.size' : 10,
    'ytick.major.size' : 10,
    'xtick.minor.width' : 0.5,
    'ytick.minor.width' : 0.5,

    #* tick label
    'xtick.labelsize' : 36,
    'ytick.labelsize' : 36,
    'xtick.major.pad' : 15,
    'ytick.major.pad' : 15,

    #* save
    'savefig.dpi' : 500,
    'savefig.transparent' : True,
    # 'savefig.facecolor' : "#ffffff",
    'savefig.bbox' : 'tight',
    'savefig.format' : 'pdf'
}


matplotlib.rcParams.update(new_rc_params)

if __name__=="__main__":
    print("This is a module draw.py")