# Synapse Modelers Workshop
A GUI to model synaptic recordings

This tool is a simple GUI that helps simulate synaptic electrophysiology recordings efficiently. Currently, it only supports Tsodyks Markram short-term plasticity model. Both voltage-clamp and current-clamp traces are supported. We also allow the generation of pseudo-traces knowing measures of synaptic activity. Optimization of multiple traces is allowed at once. Modeling and optimization for most traces take a few seconds using parallel processing, analog mathematical solutions and just-in-time compilation of some of the functions.

![image](https://user-images.githubusercontent.com/18602635/59128584-74aa4300-8938-11e9-9b16-85b3f4b221f6.png)

## Installation
The program is tested on standard python3 conda distribution. If you use conda all the dependencies will be taken care of. Then, just copy Main.py and my_scatter_matrix.py files in a folder on your computer, and create two subfolders named "csvs" and "jsons" in it and run the Main.py.

## Dependencies:
psutil, sys, os, multiprocessing, queue, tkinter, matplotlib, pandas, math, scipy, re, time, numba, and numpy

## How to use
Synaptic traces should be digitized in a specific way. Each synaptic event in a trace should have strictly three points: (1) start, (2) peak, and (3) half-decay points. Digitized traces should be saved as CSV files and saved in the csvs folder you made in the installation step. The first column should be the time in milliseconds and the second column synaptic potential in millivolts or current in picoampere. We encourage using free and open source [Engauge digitizer](https://github.com/markummitchell/engauge-digitizer/releases) for digitization purpose.

CSV files should be imported to the program. Then users need to set parameters like synaptic reversal potential, postsynaptic membrane potential, input resistance, and capacitance. After pressing optimize button parameters of Tsodyks Markram are found by optimization the model parameters. You can save the results as JSON. Saving of multiple optimization results is allowed for bootstrapping purpose.
