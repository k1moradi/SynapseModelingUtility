# Synapse Modeling Utility
A simulator to reconstruct and model synaptic recordings.

This program efficiently simulates synaptic electrophysiology recordings. Currently, we only incorporated [Tsodyks, Pawelzik and Markram's short-term plasticity model](https://pubmed.ncbi.nlm.nih.gov/9573407/). Both voltage-clamp and current-clamp traces are supported. Synaptic current is injected to a resistor-capcitor circuit, a simple model of biomembranes, in order to simulate current-clmap condition. This software also allows the reconstruction of traces knowing synaptic measures like the rise and decay times, and paired-pulse ratios. Parallel optimization of multiple traces is allowed. Modeling and optimization for most traces take a few seconds using an efficient genetic algorithm, analog mathematical modeling, parallel processing technique, and just-in-time compilation of costly functions.

![97099398-da4abc80-165e-11eb-997a-2930a680dffa](https://user-images.githubusercontent.com/18602635/108758308-f0875d00-7518-11eb-97f0-6cdbb34c52c0.png)

## Installation
The program is tested on standard python 3 [anaconda distribution](https://www.anaconda.com/distribution/). If you use anaconda all the dependencies will be taken care of. Then, just copy Main.py and Data.rar files in a folder on your computer. Extract the data.rar file, which should add two subfolders named "csvs" and "jsons". Run the Main.py via command line with the following command.

```{cmd}
> python Main.py
```

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
Synaptic traces should be digitized in a specific way. Each synaptic event in a trace should have strictly three points: (1) initiation, (2) peak, and (3) decay points. Digitized traces should be saved as CSV files and put in the csvs folder you made earlier. The first column of CSV file should be time in milliseconds and the second column synaptic potential in millivolts or synaptic current in picoamperes. We encourage using free, open source, and multiplatform [Engauge digitizer](https://github.com/markummitchell/engauge-digitizer/releases) for digitization purpose.

An example of a digitized trace in Engauge Digitizer:
![image](https://user-images.githubusercontent.com/18602635/59129236-3ca3ff80-893a-11e9-858d-bb6e74625ea6.png)

You can also reconstruct a trace, knowing the synaptic measures (synaptimetrics). Go to "start" and select "Reconstruct a trace".

![Reconstructor](https://user-images.githubusercontent.com/18602635/110882512-0da98300-82b0-11eb-90c7-5c3a3cad1eda.png)

You need to enter the values of amplitude, rise time, and decay time and choose the right synaptometric type for each. If you add ISI and paired-pulse ratio (PPR) values a trace containing multiple synaptic events will be created. On pressing "Save CSV" a synaptic trace is generated. You should save the CSV file in the "csvs" folder, then you can open it. If you enter more than one ISI value separated by commas, then you need to add more than one PPR value, separated by commas as well. In this case, more than one trace is generated, with one-minute inactivity separating them. You may also read a CSV file, which estimates synaptometrics 0%-100% rise time, decay time constant, amplitude, ISI, and PPRs of the first trace. Upon re-reconstruction, some of the slow membrane fluctuations affecting the signals would be eliminated.

CSV files should be imported to the program (Start -> Open CSV files). Then users need to set parameters like synaptic reversal potential (Erev), postsynaptic membrane potential (Vm), input resistance (Rin), and capacitance (Cm). You can use entries with blue background to enter the lower bounds of the search space and those with pink backround to set the upper bounds. After pressing optimize button synaptic parameters are found by optimization techniques (differential evolution algorithm). You can save the results in pandas compatile JSON format. Saving multiple optimization results is allowed for bootstrapping purpose, which can be done automatically by checking the box next to the optimize button. After optimization, you can press the summarize button to get an average of saved resutls. You may also press the correct button to correct the signals based on optimization results.

### Note:
* To get biological conductance (g), the g_syn value found by optimization should be multiplied by U value, i.e., g = g_syn * U. For instance, g is 716.684 * 0.036 = 25.8 (nS) in the above figure.
* Shortcuts: 
Open: CTRL + O
Close: CTRL + w
