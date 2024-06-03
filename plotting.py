from matplotlib import pyplot as plt
from matplotlib import dates as mdates
import math
import numpy as np
import sys
import json

def readJSON(filename):
    f = open(filename, 'r')
    return json.load(f)

def plotJSON(filename, fileout = 'plot.png'):

    # Reading JSON file
    data = readJSON(filename)
    plot_name = data['Plot name']
    plot_x_type = data['X type']
    plot_x_min = data['X min']
    plot_x_max = data['X max']
    plot_x_label = data['X label']
    plot_y_label = data['Y label']

    # Plotting
    if plot_x_type == 'date':
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d-%Y'))
        plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.xlim((plot_x_min, plot_x_max))
    plt.plot(data['Data'])
    plt.title(plot_name)
    plt.xlabel(plot_x_label)
    plt.ylabel(plot_y_label)
    plt.savefig(fileout)


if __name__ == '__main__':

    # Taking input
    if len(sys.argv) == 1:
        assert False, 'Filename unknown'
    filename = sys.argv[1]
    fileout = sys.argv[2] if len(sys.argv) > 2 else 'plot.png'

    # Plotting
    plotJSON(filename, fileout)

