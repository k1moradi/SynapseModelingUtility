# Synapse Modelers Workshop
A GUI to model synaptic recordings

This tool is a simple GUI that helps simulate synaptic electrophysiology recordings efficiently. Currently, it only supports Tsodyks Markram short-term plasticity model. Both voltage-calmp and current-clamp traces are supported. We also allow generation of psudo-traces knowning measures of synaptic activity. Optimization of multiple traces are allowed at once. Using parallel processing, analog mathematical solutions and just-in-time compilation of some of the functions optimization for most traces takes a few seconds.

![image](https://user-images.githubusercontent.com/18602635/59127417-96ee9180-8935-11e9-908a-3f4a15db42e7.png)

## Installation
The program is tested on standard python3 conda distribution. If you use conda you will not need to install anything, just copy Main.py and my_scatter_matrix.py files in a folder on your conputer, and create two subfolders named "csvs" and "jsons" in it and run the Main.py.

## Dependencies:
psutil, sys, os, multiprocessing, queue, tkinter, matplotlib, pandas, math, scipy, re, time, numba, and numpy

