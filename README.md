# Synapse Modelers Workshop
A GUI to model synaptic recordings

This GUI helps simulate synaptic electrophysiology recordings efficiently. Currently, it only supports Tsodyks and Markram's short-term plasticity model. Both voltage-clamp and current-clamp traces are supported. We also allow the generation of pseudo-traces knowing measures of synaptic activity like the rise and decay times, and paired-pulse ratios. Parallel optimization of multiple traces is allowed. Modeling and optimization for most traces take a few seconds using parallel processing, analog mathematical solutions and just-in-time compilation of some of the functions.

![image](https://user-images.githubusercontent.com/18602635/62956427-ed0e0400-bdc0-11e9-99f7-3a21b997141f.png)

## Installation
The program is tested on standard python 3 [anaconda distribution](https://www.anaconda.com/distribution/). If you use anaconda all the dependencies will be taken care of. Then, just copy Main.py and my_scatter_matrix.py files in a folder on your computer, and create two subfolders named "csvs" and "jsons" in it and run the Main.py.

## Module Dependencies:
psutil, matplotlib, pandas, math, scipy, numba, numpy, sys, os multiprocessing, queue, tkinter, re, and time.

|module|tested version|
|---|---|
|Python|3.7.3|
|psutil|5.6.2|
|matplotlib|3.1.0|
|pandas|0.24.2|
|scipy|1.2.1|
|numba|0.43.1|
|numpy|1.16.4|
|tk|8.6.8|
|PyQt|5|


## How to use
Synaptic traces should be digitized in a specific way. Each synaptic event in a trace should have strictly three points: (1) initiation, (2) peak, and (3) decay points. Digitized traces should be saved as CSV files and put in the csvs folder you made earlier. The first column of CSV file should be the time in milliseconds and the second column synaptic potential in millivolts or current in picoamperes. We encourage using free and open source [Engauge digitizer](https://github.com/markummitchell/engauge-digitizer/releases) for digitization purpose.

An example of a digitized trace in Engauge digitizer:
![image](https://user-images.githubusercontent.com/18602635/59129236-3ca3ff80-893a-11e9-858d-bb6e74625ea6.png)

CSV files should be imported to the program (Start -> Open CSV files). Then users need to set parameters like synaptic reversal potential (Erev), postsynaptic membrane potential (Vm), input resistance (Rin), and capacitance (Cm). You can use entries with blue background to enter the lower bount of the search space and those with pink backround to set the upper bound. After pressing optimize button parameters of Tsodyks Markram are found by optimization techniques. You can save the results in panda compatile JSON format. Saving of multiple optimization results is allowed for bootstrapping purpose.
