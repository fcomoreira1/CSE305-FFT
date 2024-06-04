from matplotlib import pyplot as plt
from matplotlib import dates as mdates
from datetime import datetime
import math
import pandas as pd
import numpy as np
import sys
import json

def readJSON(filename):
    f = open(filename, 'r')
    return json.load(f)

def plotJSON(filename, fileout = 'plot.png'):

    # Reading JSON file
    with open(filename, 'r') as file:
        data = json.load(file)
        plot_name = data['Plot name']
        plot_x_type = data['X type']
        plot_x_min = data['X min']
        plot_x_max = data['X max']
        plot_x_label = data['X label']
        plot_y_label = data['Y label']

    # Plotting
    if plot_x_type == 'date':
        x = pd.date_range(start = plot_x_min, end=plot_x_max, periods = int(input_size))
        plt.xticks(rotation=30)
        plt.gca().xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[4,7,10]))
    else:
        x = np.linspace(plot_x_min, plot_x_max, data['Size'])
    plt.plot(data['Data'])
    plt.title(plot_name)
    plt.xlabel(plot_x_label)
    plt.ylabel(plot_y_label)
    plt.savefig(fileout)

def plotFile(filename, fileout = 'plot.png'):

    # Reading file
    with open(filename, 'r') as file:
        data = []
        line = file.readline()
        while (line not in {"","\n"}):
            data.append(float(line))
            line = file.readline()
        input_size, plot_name, plot_x_type, plot_x_min, plot_x_max, plot_x_label, plot_y_label = file.readline().split(',')

    # Plotting
    if plot_x_type == 'date':
        x = pd.date_range(start = plot_x_min, end=plot_x_max, periods = int(input_size))
        plt.xticks(rotation=30)
        plt.gca().xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[4,7,10]))
    else:
        x = np.linspace(plot_x_min, plot_x_max, int(input_size))
    plt.plot(x, data)
    plt.title(plot_name)
    plt.xlabel(plot_x_label)
    plt.ylabel(plot_y_label)
    plt.savefig(fileout)


if __name__ == '__main__':

    # # Taking input
    # if len(sys.argv) == 1:
    #     assert False, 'Filename unknown'
    # filename = sys.argv[1]
    # fileout = sys.argv[2] if len(sys.argv) > 2 else 'plot.png'

    # # Plotting
    # plotFile(filename, fileout)

    # inp = sys.argv
    # if len(inp) < 10:
    #     assert False, 'Not enough args'

    f = open(sys.argv[-1], 'r')
    line = f.read()[:-1]
    
    inp = line.split(' ')

    data = []
    for i in range(1, len(inp)-8):
        data.append(float(inp[i]))

    input_size, plot_name, plot_x_type, plot_x_min, plot_x_max, plot_x_label, plot_y_label = inp[-8:-1]

    # Plotting
    if plot_x_type == 'date':
        x = pd.date_range(start = plot_x_min, end=plot_x_max, periods = int(input_size)-1)
        plt.xticks(rotation=30)
        plt.gca().xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[4,7,10]))
    else:
        x = np.linspace(plot_x_min, plot_x_max, int(input_size))
    plt.plot(x, data)
    plt.title(plot_name)
    plt.xlabel(plot_x_label)
    plt.ylabel(plot_y_label)
    plt.savefig(inp[-1])

