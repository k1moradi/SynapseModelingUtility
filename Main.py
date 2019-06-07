# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 16:59:35 2019

@author: Keivan Moradi
"""
from psutil import cpu_count
from sys import platform, setrecursionlimit
from os import listdir, path, getcwd, environ, remove, system
from multiprocessing import freeze_support, Process, Queue
from queue import Empty
from tkinter import Tk, Menu, filedialog, Canvas, Scrollbar, Toplevel, Scale, messagebox, TclError
from tkinter import Frame, Button, Checkbutton, Label, Entry, OptionMenu, StringVar, BooleanVar, IntVar
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure, rcParams
import matplotlib.pyplot as plt
from pandas import DataFrame, read_csv, read_json
from my_scatter_matrix import scatter_matrix
from math import exp, fabs, isnan, log
from scipy.integrate import odeint
from scipy.optimize import differential_evolution
from re import match, findall, split
from time import time, sleep
from numba import jit
import numpy as np
environ['NUMBA_ENABLE_AVX'] = '1'
rcParams['font.family'] = 'monospace'  # set the font family used in
setrecursionlimit(6000)
if platform.startswith("win"):
    from ctypes import windll, c_int, byref
    # Query DPI Awareness (Windows 10 and 8)
    awareness = c_int()
    windll.shcore.GetProcessDpiAwareness(0, byref(awareness))
    # Set DPI Awareness  (Windows 10 and 8)
    windll.shcore.SetProcessDpiAwareness(1)
    # the argument is the awareness level, which can be 0, 1 or 2:
    # for 1-to-1 pixel control it seems to need it to be non-zero (I'm using level 2)
    # Set DPI Awareness  (Windows 7 and Vista)
    windll.user32.SetProcessDPIAware()
    # behaviour on later OSes is undefined, although when I run it on my Windows 10 machine,
    # it seems to work with effects identical to SetProcessDpiAwareness(1)


class Experiment:
    rowCount = -1
    CSVs_FOLDER = 'csvs'
    JSONs_FOLDER = 'jsons'
    KEYS = ['Erev', 'Vm', 'Rin', 'Cm', 'g_syn', 'tau_d', 'tau_r', 'tau_f', 'U']
    UNITS = ['mV', 'mV', 'M'+u'\u2126', 'PF', 'nS', 'ms', 'ms', 'ms', '']
    OPTIMIZATION = ['g_syn', 'tau_d', 'tau_r', 'tau_f', 'U', 'Rin', 'Cm']  # keep this order
    MAY_NEED_OPTIMIZATION = ['Cm', 'Rin']
    WORKING_DIRECTORY = '.'
    DECIMAL_POINTS = 3
    default_mode = 'voltage-clamp'
    summary_window = None

    @staticmethod
    def is_float(p):
        if match(r'^\s*-?(?:\d*\.\d+|\d+\.?\d*)(?:e-?\d*)?\s*$', p) or p == '' or p == '-':
            return True
        else:
            return False

    @staticmethod
    def set_initial_parameters_values(file_name):
        parameters = dict()
        parameters['Mode'] = Experiment.default_mode
        parameters['Erev'] = 0.0
        parameters['Vm'] = -80.0
        parameters['Data'] = []

        parameters['Rin'] = 100.0
        parameters['Cm'] = 100.0
        parameters['g_syn'] = 3.0
        parameters['tau_d'] = 5.0
        parameters['tau_f'] = 1000.0
        parameters['tau_r'] = 100.0
        parameters['U'] = 0.5

        parameters['Rin_min'] = 10.0
        parameters['Cm_min'] = 10.0
        parameters['g_syn_min'] = 0.01
        parameters['tau_d_min'] = 0.3
        parameters['tau_f_min'] = 1.0
        parameters['tau_r_min'] = 1.0
        parameters['U_min'] = 0.001

        parameters['Rin_max'] = 400.0
        parameters['Cm_max'] = 400.0
        parameters['g_syn_max'] = 30.0
        parameters['tau_d_max'] = 30.0
        parameters['tau_f_max'] = 30000.0
        parameters['tau_r_max'] = 30000.0
        parameters['U_max'] = 1 - 10**-Experiment.DECIMAL_POINTS

        json_file = Experiment.WORKING_DIRECTORY + '/' + Experiment.JSONs_FOLDER + '/' + file_name + '.json'
        if path.isfile(json_file):
            loaded_parameters = Experiment.get_json_file_as_dict_list(json_file)[0]
            for key in parameters:
                if key in loaded_parameters:
                    parameters[key] = loaded_parameters[key]

        return parameters

    @staticmethod
    def get_json_file_as_dict_list(json_file):
        return Experiment.get_json_file_as_data_frame(json_file).to_dict(orient='records')

    @staticmethod
    def get_json_file_as_data_frame(json_file):
        with open(json_file, 'r') as file_pointer:
            return read_json(file_pointer)

    @staticmethod
    def get_csv_file_as_data_frame(file_name):
        csv_file = Experiment.WORKING_DIRECTORY + '/' + Experiment.CSVs_FOLDER + '/' + file_name + '.csv'
        with open(csv_file, 'r') as file_pointer:
            data_frame = read_csv(file_pointer)
            data_frame = data_frame.rename(columns={
                data_frame.columns[0]: "time",
                data_frame.columns[1]: "signal"})
        return data_frame

    @staticmethod
    def correct_digitization(data_frame, mode, membrane_potential):
        signal_correction = data_frame.signal[0]
        signal_correction -= 0.0 if mode == "voltage-clamp" else membrane_potential
        return data_frame.assign(
            time=data_frame.time - data_frame.time[0],
            signal=data_frame.signal - signal_correction)

    @staticmethod
    def initiation_points(mode, time_list, signal_list):
        # return time_list[1:len(time_list) - 1:3], signal_list[1:len(signal_list) - 1:3]
        if mode == "current-clamp":
            return time_list[0:len(time_list) - 1:3], signal_list[0:len(signal_list) - 1:3]
        elif mode == "voltage-clamp":
            return time_list[1:len(time_list) - 1:3], signal_list[1:len(signal_list) - 1:3]
        else:
            print("Error in update_initiation_points function: experimentType is not current-clamp or voltage-clamp")

    @staticmethod
    def signal_label(mode):
        return 'Signal (' + ('pA' if mode == 'voltage-clamp' else 'mV') + ')'

    @staticmethod
    def delete_file_(file_name, file_extension, message):
        if file_extension == '.json':
            file_ = Experiment.WORKING_DIRECTORY + '/' + Experiment.JSONs_FOLDER + '/' + file_name + file_extension
        elif file_extension == '.csv':
            file_ = Experiment.WORKING_DIRECTORY + '/' + Experiment.CSVs_FOLDER + '/' + file_name + file_extension
        if path.isfile(file_):
            user_wants_to_delete = messagebox.askyesno(file_extension + 'File:', message)
            if user_wants_to_delete:
                remove(file_)

    @staticmethod
    def plot_saved_results_(file_name, window, amplitude, ppr, isi):
        row, width = 0, 70
        font =('Courier New', 18)
        ref_id, location = findall(r'^(\d+)-?(.*)$', file_name)[0]
        json_file = Experiment.WORKING_DIRECTORY + '/' + Experiment.JSONs_FOLDER + '/' + file_name + ".json"
        if location:
            referencing_stuff = '(n=1)@%s{%s:digitized:modeled already}' % (ref_id, location)
        else:
            referencing_stuff = '(n=1)@%s{digitized:modeled already}' % ref_id
        if Experiment.summary_window:
            Experiment.summary_window.destroy()
        Experiment.summary_window = window
        window.title('Summary')
        window.geometry('+%d+%d' % (55, 60))
        window.columnconfigure(1, weight=1)
        msp = StringVar(window)  # measure of spread
        msp.set(global_msp)
        mct = StringVar(window)  # measure of central tendency
        mct.set(global_mct)
        Label(window, text='Stats:', font=font).grid(row=0, column=0, sticky="NEWS")
        frame = Frame(window)
        OptionMenu(frame, mct, *{'Mean', 'Median', 'Mode'}).grid(row=0, column=0, sticky='W')
        OptionMenu(frame, msp, *{'SEM', 'SD', 'IQR'}).grid(row=0, column=1, sticky='W')
        frame.grid(row=0, column=1, sticky='W')
        entries = dict()
        keys_colors = [('File', 'DarkOrchid1'), ('error', 'red'),
                       ('g_syn', 'blue'), ('tau_d', 'blue'), ('tau_r', 'blue'), ('tau_f', 'blue'), ('U', 'blue'),
                       ('Peak', 'green')]
        if ppr:
            keys_colors += [('2/1PPR', 'green'), ('ISI', 'green'), ('ST-P', 'green')]
        for row, key_color in enumerate(keys_colors, start=1):
            key, color = key_color
            Label(window, text=key+':', font=font, anchor="e").grid(row=row, column=0, sticky="NEWS")
            entries[key] = EntryWithPlaceholder(master=window, color=color, state='readonly', width=width, font=font)
            entries[key].grid(row=row, column=1, sticky="WE")
        entries['error'].config(state='disabled', disabledforeground='red')
        entries['File'].put_placeholder(file_name)
        entries['Peak'].put_placeholder("% 9.2f%s" % (amplitude, referencing_stuff))
        if ppr:
            entries['2/1PPR'].put_placeholder("% 9.2f%s" % (ppr, referencing_stuff))
            entries['ISI'].put_placeholder("% 6.0f%s" % (isi, referencing_stuff))
            entries['ST-P'].put_placeholder('Facilitating (F), Pseudolinear (F▶D) or No ST-P' if ppr > 1
                                            else 'Depressing (D), Either D or F▶D or No ST-P')

        def set_print_values(*args):
            global global_msp, global_mct, bootstrap_max_iteration
            global_msp, global_mct = msp.get(), mct.get()
            if platform.startswith("win"):
                system('cls')
            elif platform.startswith("linux"):
                system('clear')
            console_print = '\n\33[36m  File\33[0m:\33[31m%s\33[0m\n' % file_name
            try:
                n = bootstrap_max_iteration.get()+1
                df = Experiment.get_json_file_as_data_frame(json_file)
                if len(df) > n:
                    df = df.sort_values(by=['error']).head(n)
                keys = ['g_syn', 'tau_d', 'tau_r', 'tau_f', 'U', 'error']
                df = df[keys]
                df_description = df.describe()
                annotation = '\n File: %s (n=%d)\n' % (file_name, len(df))
                for key in keys:
                    count, mean, std, v_min, p25, p50, p75, v_max = df_description[key]
                    if global_msp == 'SEM':
                        std = std / count**0.5
                    elif global_msp == 'IQR':
                        std = p75 - p25
                    if global_mct == 'Median':
                        mean = df[key].median()
                    elif global_mct == 'Mode':
                        mean = df[key].mode()[0]
                    console_print += "\33[36m%6s\33[0m:% 11.4f\u00B1%.4f [%.4f to %.4f](n=%d)\n" % (
                        key, mean, std, v_min, v_max, count)
                    annotation += "%6s:% 10.3f \u00B1% 11.3f [% 10.3f to% 10.3f]\n" % (key, mean, std, v_min, v_max)
                    entries[key].put_placeholder(
                        "% 11.3f\u00B1%.3f [%.3f to %.3f](n=%d)" % (mean, std, v_min, v_max, count))
                if show_matrix_plot.get():
                    scatter_matrix(df, alpha=0.5, figsize=(12, 9), diagonal='kde', grid=True)
                    plt.get_current_fig_manager().window.showMaximized()
                    plt.annotate(annotation, xy=(0.365, 0.67), xycoords='figure fraction',
                                 size=24 if window.master.winfo_screenwidth() > 1280 else 16)
                    plt.show(block=False)
            except FileNotFoundError:
                console_print += '  Optimization results available yet to be summarized!\n'
                pass
            console_print += '\33[36m  Peak\33[0m:% 11.2f%s\n' % (amplitude, referencing_stuff)
            if ppr:
                console_print += '\33[36m2/1PPR\33[0m:% 11.2f%s\n' % (ppr, referencing_stuff)
                console_print += '\33[36m   ISI\33[0m:% 11.0f%s\n' % (isi, referencing_stuff)
                console_print += '\33[36m  ST-P\33[0m:% 11s' % (
                    'Facilitating (F) or Pseudolinear (F▶D)' if ppr > 1 else 'Depressing (D)')
            print(console_print)
        set_print_values.entries = entries
        set_print_values.file_name = file_name
        set_print_values.json_file = json_file
        set_print_values()
        msp.trace('w', set_print_values)
        mct.trace('w', set_print_values)

    def __init__(self, parent, file_name):
        # each experiment is one row, count the number of experiments created to put them on the next row
        Experiment.rowCount += 1
        # make a Tcl wrapper around is_float function so that tk can use it as a command
        Experiment.IS_FLOAT = parent.register(Experiment.is_float)

        self.FILE_NAME = file_name

        self.parameters = Experiment.set_initial_parameters_values(file_name)

        data = Experiment.get_csv_file_as_data_frame(file_name)
        data = Experiment.correct_digitization(data, self.parameters['Mode'], self.parameters['Vm'])

        self.experimentFrame = Frame(parent)
        self.experimentFrame.grid(row=Experiment.rowCount, column=0, sticky='NESW')

        self.topButtonFrame = Frame(self.experimentFrame)
        width = 16
        Button(self.topButtonFrame, text="Run", command=self.run_model, width=width).grid(
            row=0, column=0)
        self.optimizeFrame = Frame(self.topButtonFrame)
        self.optimizeButton = Button(self.optimizeFrame, text="Optimize", command=self.optimize, width=width-4)
        self.optimizeButton.grid(row=0, column=0)
        self.bootstrap_mode = BooleanVar()
        Checkbutton(self.optimizeFrame, variable=self.bootstrap_mode).grid(row=0, column=1)
        self.optimizeFrame.grid(
            row=0, column=1)
        self.buttonCorrect = Button(
            self.topButtonFrame, text="Correct", command=self.correct_data, width=width, state='disabled')
        self.buttonCorrect.grid(
            row=0, column=2)
        Button(self.topButtonFrame, text="Reload", command=self.reload_data, width=width).grid(
            row=1, column=0)
        self.buttonSave = Button(self.topButtonFrame, text="Save", command=self.save, width=width, state='disabled')
        self.buttonSave.grid(
            row=1, column=1)
        Button(self.topButtonFrame, text="Summarize", command=self.plot_saved_results, width=width).grid(
            row=1, column=2)
        self.topButtonFrame.grid(row=0, column=0)

        self.lowerBoxFrame = Frame(self.experimentFrame)
        Label(self.lowerBoxFrame, text="Mode:").grid(
            row=0, column=0)
        self.experimentModeChoice = StringVar(self.lowerBoxFrame, self.parameters['Mode'])
        self.experimentModeChoice.trace('w', self.mode_change)
        OptionMenu(self.lowerBoxFrame, self.experimentModeChoice, *{'voltage-clamp', 'current-clamp'}).grid(
            row=0, column=1, columnspan=2, sticky='WENS')
        Button(self.lowerBoxFrame, text="delete file", command=self.delete_file).grid(row=0, column=3, sticky="NEWS")
        self.entries = dict()
        self.checkButton = dict()
        self.parameterMayNeedOptimization = dict()
        row = 1
        col_width = int(width * 3 / 4)
        for row, key in enumerate(self.KEYS, start=1):
            unit = ' ('+Experiment.UNITS[row-1]+')' if Experiment.UNITS[row-1] != '' else ''
            Label(self.lowerBoxFrame, text=key+unit, height=1, justify='left').grid(row=row, column=0, sticky="W")
            if key in Experiment.OPTIMIZATION:
                for column, (suffix, color) in enumerate(
                        [('', 'white'), ('_min', 'aliceblue'), ('_max', 'lavenderblush')], start=1):
                    self.entries[key+suffix] = EntryWithPlaceholder(
                        master=self.lowerBoxFrame, bd=3, width=col_width, bg=color,
                        validate='all', validatecommand=(Experiment.IS_FLOAT, '%P'))
                    self.entries[key+suffix].set(round(self.parameters[key+suffix], Experiment.DECIMAL_POINTS))
                    self.entries[key+suffix].grid(row=row, column=column, sticky='WENS')
            else:
                self.entries[key] = EntryWithPlaceholder(
                    master=self.lowerBoxFrame, bd=3, width=col_width*3,
                    validate='all', validatecommand=(Experiment.IS_FLOAT, '%P'))
                self.entries[key].set(round(self.parameters[key], Experiment.DECIMAL_POINTS))
                self.entries[key].grid(row=row, column=1, columnspan=3, sticky='WENS')
            if key in Experiment.MAY_NEED_OPTIMIZATION:
                self.parameterMayNeedOptimization[key] = BooleanVar()
                self.parameterMayNeedOptimization[key].set(False)
                self.checkButton[key] = Checkbutton(self.lowerBoxFrame, variable=self.parameterMayNeedOptimization[key])
                self.checkButton[key].grid(row=row, column=4)
        self.entries['Vm'].bind("<FocusOut>",
                                lambda event: self.reload_data()if self.parameters['Mode'] == 'current-clamp' else None)
        self.entries['Vm'].bind("<Return>",
                                lambda event: self.reload_data()if self.parameters['Mode'] == 'current-clamp' else None)
        self.error, self.bootstrap_counter = StringVar(), IntVar()
        self.bootstrap_counter.set(0)
        Label(self.lowerBoxFrame, textvariable=self.bootstrap_counter, height=1, justify='right').grid(
            row=row + 1, column=0)
        Label(self.lowerBoxFrame, textvariable=self.error, height=1, justify='right').grid(
            row=row+1, column=1, columnspan=3)
        self.lowerBoxFrame.grid(row=1, column=0, pady=(0, 0))
        self.toggle_entries()

        # setup the initial plots
        self.fig = Figure(figsize=(6, 4), dpi=100 if self.experimentFrame.winfo_screenwidth() > 1280 else 72)
        self.fig.subplots_adjust(left=0.07, right=0.98, top=0.93, bottom=0.12)
        self.subplot = self.fig.add_subplot(111)  # n_rows, n_cols, index
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.experimentFrame)
        self.canvas.get_tk_widget().grid(row=0, column=1, rowspan=len(self.KEYS), sticky='NSEW')
        self.experimentFrame.columnconfigure(1, weight=1)
        self.experimentFrame.rowconfigure(1, weight=1)
        [t_init, y_init] = Experiment.initiation_points(self.parameters['Mode'], data.time, data.signal)
        self.plotData, = self.subplot.plot(data.time, data.signal, '--bo', color='blue', label='Data')
        self.plotModel, = self.subplot.plot([], [], '-X', color='red', label='Model')
        self.plotCorrectedSignal, = self.subplot.plot([], [], '-.bo', color='green', label='Corrected Data')
        self.plotInitTimes, = self.subplot.plot(t_init, y_init, 'P', color='magenta', label='Init')
        self.subplot.set_title(self.FILE_NAME, fontsize=14, x=0.8, y=1)
        self.subplot.set_autoscale_on(True)
        self.subplot.set_ylabel(Experiment.signal_label(self.parameters['Mode']), fontsize=12)
        self.subplot.set_xlabel("time (ms)", fontsize=12)
        self.subplot.grid(True)
        self.subplot.spines['right'].set_position('zero')
        self.subplot.spines['top'].set_position('zero')
        self.subplot.legend(loc='center right', ncol=1, borderaxespad=0.)
        # add toolbar
        self.toolbarFrame = Frame(self.experimentFrame)
        self.toolbarFrame.place(
            relx=0.0, rely=0.0, anchor="se", width=500, y=33,
            x=width * 3 * 19 + 3 if self.experimentFrame.winfo_screenwidth() > 1280 else width * 3 * 18 + 3)
        NavigationToolbar2Tk(self.canvas, self.toolbarFrame).configure(bg='white')

        # needed for multi-threading
        self.queue = Queue()

    def plot(self, plot_obj, x, y):
        plot_obj.set_xdata(x)
        plot_obj.set_ydata(y)
        self.subplot.relim()  # recompute the axis size
        self.subplot.autoscale_view()  # update axis view size
        self.canvas.draw()

    def get_ty_data(self):
        return self.plotData.get_xydata().tolist()

    def get_t_data(self):
        return list(self.plotData.get_xdata())

    def get_init_times(self):
        return list(self.plotInitTimes.get_xdata())

    def reload_data(self):
        self.update_params()
        data = Experiment.get_csv_file_as_data_frame(self.FILE_NAME)
        data = Experiment.correct_digitization(data, self.parameters['Mode'], self.parameters['Vm'])
        [t, y] = Experiment.initiation_points(self.parameters['Mode'], data.time, data.signal)
        # since the data stored in the plots are accessible we do not need to store is somewhere else
        self.plot(self.plotCorrectedSignal, [], [])
        self.plot(self.plotInitTimes, t, y)
        self.plot(self.plotData, data.time, data.signal)

    def toggle_entries(self):
        if self.parameters['Mode'] == "voltage-clamp":
            Experiment.default_mode = "voltage-clamp"
            for key in Experiment.MAY_NEED_OPTIMIZATION:
                self.entries[key].config(state='readonly')
                self.entries[key+'_min'].config(state='readonly')
                self.entries[key+'_max'].config(state='readonly')
                self.checkButton[key].config(state='disabled')
        else:
            Experiment.default_mode = "current-clamp"
            for key in Experiment.MAY_NEED_OPTIMIZATION:
                self.entries[key].config(state='normal')
                self.entries[key + '_min'].config(state='normal')
                self.entries[key + '_max'].config(state='normal')
                self.checkButton[key].config(state='normal')

    def mode_change(self, *args):
        new_mode = self.experimentModeChoice.get()
        if new_mode != self.parameters['Mode']:
            self.parameters['Mode'] = new_mode
            Experiment.delete_file_(
                self.FILE_NAME, '.json',
                'Changing mode makes saved files irrelevant.\n'
                'Do you want to delete the file?')
            self.reload_data()
            self.toggle_entries()
            self.plot(self.plotModel, [], [])
            self.subplot.set_ylabel(Experiment.signal_label(self.parameters['Mode']))

    def delete_file(self):
        Experiment.delete_file_(self.FILE_NAME, '.json', 'Do you want to delete the file?')

    def update_params(self):
        vm_before = self.parameters['Vm']
        for key in Experiment.KEYS:
            self.parameters[key] = float(self.entries[key].get())
            if key in Experiment.OPTIMIZATION:
                self.parameters[key+'_min'] = float(self.entries[key+'_min'].get())
                self.parameters[key+'_max'] = float(self.entries[key+'_max'].get())
        if round(vm_before, 2) != round(self.parameters['Vm'], 2) and self.parameters['Mode'] == 'current-clamp':
            self.reload_data()

    def run_model(self):
        self.update_params()
        model_input = [self.parameters[key] for key in self.OPTIMIZATION[0:-2]]
        bounds = tuple((self.parameters[key+'_min'], self.parameters[key+'_max']) for key in self.OPTIMIZATION[0:-2])
        if self.parameters['Mode'] == "voltage-clamp":
            model = ExperimentVoltageClamp(
                self.get_ty_data(), self.get_init_times(), self.parameters['Erev'], self.parameters['Vm'], bounds)
        elif self.parameters['Mode'] == "current-clamp":
            input_tuple = (
                self.get_ty_data(), self.get_init_times(), self.parameters['Erev'],
                self.parameters['Vm'], self.parameters['Rin'], self.parameters['Cm'])
            model = ExperimentCurrentClamp(*(input_tuple + (bounds,)))
            if self.parameterMayNeedOptimization['Cm'].get() and self.parameterMayNeedOptimization['Rin'].get():
                model_input += [self.parameters['Rin'], self.parameters['Cm']]
                bounds += tuple((self.parameters[key+'_min'], self.parameters[key+'_max']) for key in ['Rin', 'Cm'])
                model = ExperimentCurrentClampRinCm(*(input_tuple + (bounds,)))
            elif self.parameterMayNeedOptimization['Cm'].get():
                model_input += [self.parameters['Cm']]
                bounds += tuple((self.parameters[key + '_min'], self.parameters[key + '_max']) for key in ['Cm'])
                model = ExperimentCurrentClampCapacitance(*(input_tuple + (bounds,)))
            elif self.parameterMayNeedOptimization['Rin'].get():
                model_input += [self.parameters['Rin']]
                bounds += tuple((self.parameters[key + '_min'], self.parameters[key + '_max']) for key in ['Rin'])
                model = ExperimentCurrentClampResistance(*(input_tuple + (bounds,)))
        else:
            raise Exception("optimize function: experimentType should be either voltage-clamp or current-clamp")
        model.simulate(model_input)
        self.plot(self.plotModel, self.get_t_data(), model.simulatedSignal)
        return model

    def process_queue(self):
        try:  # check the queue for the optimization results then show result
            self.show_results(*self.queue.get(block=False))
        except Empty:
            self.lowerBoxFrame.after(100, self.process_queue)

    def optimize(self):
        model = self.run_model()
        global thread_numbers, population_size, print_results_when_optimization_ends
        MultiProcessOptimization(
            self.queue, model, self.parameters['Mode'], thread_numbers.get(), population_size.get(),
            self.bootstrap_counter.get(), self.bootstrap_mode.get(), self.FILE_NAME,
            print_results_when_optimization_ends.get()
        ).start()
        self.lowerBoxFrame.after(100, self.process_queue)
        self.optimizeButton.configure(state='disabled')
        for frame in self.experimentFrame.winfo_children():
            if frame.winfo_class() == 'Toplevel':
                continue
            for child in frame.winfo_children():
                if child.winfo_class() not in ['Frame', 'Label'] and child.cget('text') != 'Summarize':
                    child.configure(state='disabled') if child.winfo_class() != 'Entry' \
                        else child.configure(state='readonly')

    def show_results(self, results, corrected_signal, optimization_time):
        global ring_bell_when_optimization_ends
        # set the entries
        self.error.set('Optimized in %.1f seconds. Error is %.5f' % (optimization_time, results.fun))
        keys = Experiment.KEYS[4:]
        # oder maters here since in results.x[-2] is Rin and results.x[-1] is Cm
        if self.parameterMayNeedOptimization['Rin'].get():
            keys.append('Rin')
        if self.parameterMayNeedOptimization['Cm'].get():
            keys.append('Cm')
        for key, result in zip(keys, results.x):
            self.parameters[key] = result
            self.entries[key].set(round(result, Experiment.DECIMAL_POINTS))

        if self.bootstrap_mode.get() and self.bootstrap_counter.get() < bootstrap_max_iteration.get():
            self.bootstrap_counter.set(self.bootstrap_counter.get()+1)
            self.save()
            self.optimize()
        else:
            # enable GUI widgets
            self.optimizeButton.configure(state='normal')
            for frame in self.experimentFrame.winfo_children():  # enable button if optimization had results
                if frame.winfo_class() == 'Toplevel':
                    continue
                for child in frame.winfo_children():
                    if child.winfo_class() != 'Frame':
                        child.configure(state='normal')
            self.toggle_entries()
            fake_high_res_t_data = sorted([float(x) for x in range(int(self.get_t_data()[-1]))] + self.get_t_data())
            fake_high_res_data = [[t, 0] for t in fake_high_res_t_data]
            model = self.run_model()  # model should run after setting the entries
            model.dataList = fake_high_res_data
            model.simulate(results.x)
            self.plot(self.plotModel, [row[0] for row in fake_high_res_data], model.simulatedSignal)
            self.bootstrap_counter.set(0)
            if ring_bell_when_optimization_ends.get():
                self.experimentFrame.bell()
            if corrected_signal:
                self.plot(self.plotCorrectedSignal, self.get_t_data(), corrected_signal)

    def correct_data(self):
        times = self.plotCorrectedSignal.get_xdata()
        signals = self.plotCorrectedSignal.get_ydata()
        self.plot(self.plotData, times, signals)
        init_times = self.get_init_times()
        init_time = init_times.pop(0)
        init_signals = []
        for time_, signal in zip(times, signals):
            if time_ == init_time:
                init_signals.append(signal)
                if init_times:
                    init_time = init_times.pop(0)
        self.plot(self.plotInitTimes, self.get_init_times(), init_signals)
        self.plot(self.plotCorrectedSignal, [], [])
        self.buttonCorrect.config(state='disabled')

    def save(self):
        model_info = dict()
        model_info['Mode'] = [self.experimentModeChoice.get()]
        for key in Experiment.KEYS:
            model_info[key] = [float(self.entries[key].get())]
            if key in Experiment.OPTIMIZATION:
                model_info[key+'_min'] = [float(self.entries[key+'_min'].get())]
                model_info[key+'_max'] = [float(self.entries[key+'_max'].get())]
        model_info['error'] = [float(findall(r'\d+(?:[.]\d+)?', self.error.get())[1])]
        # model_info['Data'] = [self.plotCorrectedSignal.get_xydata().tolist()]
        current_data = DataFrame.from_dict(model_info)
        json_file = Experiment.WORKING_DIRECTORY + '/jsons/' + self.FILE_NAME + ".json"
        if path.isfile(json_file):
            current_data = current_data.append(
                DataFrame.from_dict(Experiment.get_json_file_as_dict_list(json_file)),
                sort=False)
        current_data.to_json(path_or_buf=json_file, orient='records')
        self.buttonSave.config(state='disabled')

    def plot_saved_results(self):
        data = self.get_ty_data()
        (t0, y0), (t1, y1) = data[0], data[1]
        amplitude, ppr,  isi = y1-y0, None, None
        if len(data) > 4 and amplitude != 0.0:
            (t3, y3), (t4, y4) = data[3], data[4]
            ppr, isi = abs((y4 - y3) / amplitude), t4 - t1
        Experiment.plot_saved_results_(self.FILE_NAME, Toplevel(self.experimentFrame), amplitude, ppr, isi)


class MultiProcessOptimization(Process):
    def __init__(self, queue_, model, mode, thread_numbers_, population_size_, counter, bootstrap_mode, file_name,
                 print_results_when_optimization_ends_):
        Process.__init__(self)
        self.queue = queue_
        self.model = model
        self.mode = mode
        self.thread_numbers = thread_numbers_
        self.population_size = population_size_
        self.bootstrap_counter = counter
        self.bootstrap_mode = bootstrap_mode
        self.FILE_NAME = file_name
        self.print_results_when_optimization_ends = print_results_when_optimization_ends_

    def run(self):
        t = time()
        results = differential_evolution(
            self.model.simulate,
            bounds=self.model.BOUNDS,  # differential_evolution params
            strategy='best2bin',
            polish=True,
            updating='deferred' if (self.mode == 'current-clamp' and self.thread_numbers > 1) else 'immediate',
            workers=self.thread_numbers if self.mode == 'current-clamp' else 1,
            popsize=self.population_size,
            maxiter=10000000,
            # mutation=(0.1, 1),
            tol=10**-Experiment.DECIMAL_POINTS
        )
        if self.bootstrap_mode:
            corrected_signal = None
        else:
            self.model.simulate(results.x)
            corrected_signal = self.model.correct_data()
        if self.print_results_when_optimization_ends:
            print(
                '\33[34m error \33[0m % 9.5f' % results.fun,
                '\33[32m result\33[0m'+str(['% 8.2f' % elem for elem in results.x]),
                self.FILE_NAME,
                ('number %d' % self.bootstrap_counter) if self.bootstrap_mode else ''
            )
        self.queue.put([results, corrected_signal, round(time() - t, 2)])


class ExperimentVoltageClamp:
    def __init__(self, data_list, init_times, synaptic_reversal_potential, holding_potential, bounds):
        self.BOUNDS = bounds
        self.dataList = data_list  # format [[t in ms, i in pA], ]
        self.initTimes = init_times  # [synTimes in ms]
        self.num_events = len(init_times) + 1
        self.drivingForce = float(holding_potential) - float(synaptic_reversal_potential)  # in mV
        self.simulatedSignal = [0.0]
        amp = fabs(self.dataList[0][1] - self.dataList[1][1])
        self.normalization_value = (amp if amp != 0.0 else 1.0)*(self.num_events if self.num_events != 0 else 1.0)

    def simulate(self, input_vec):
        g0, tau_d, tau_r, tau_f, u = input_vec
        x0, y0, u0, delta_t = 1.0, 0.0, 0.0, 0.0
        g, x0, y0, u0 = synaptic_event(delta_t, g0, tau_d, tau_r, tau_f, u, u0, x0, y0)
        init_times = self.initTimes.copy()
        init_time = init_times[0]
        del init_times[0]
        signal = []
        error = 0.0
        error_weight = 1.0
        for data_time, data_signal in self.dataList:
            delta_t = data_time - init_time
            if delta_t < 0:
                signal.append(0)
                continue
            if init_times and data_time >= init_times[0]:
                g, x0, y0, u0 = synaptic_event(delta_t, g0, tau_d, tau_r, tau_f, u, u0, x0, y0)
                init_time, delta_t, error_weight = init_times[0], 0.0, 1.0
                del init_times[0]
            signal.append(synaptic_current(g, delta_t, tau_d, self.drivingForce))
            if not isnan(data_signal):
                error = soft_l1_loss(error, signal[-1], data_signal, error_weight)
            error_weight /= self.num_events
            self.simulatedSignal = signal
        error /= self.normalization_value
        if self.num_events == 2:
            error += (1 - u/self.BOUNDS[4][1]) * error
        return error

    def correct_data(self):
        # make an array of difference between simulated and recorded initiation times
        # self.initTimes are at the signal pick because of the limitations of current synaptic model
        # I assume one event before pick (i.e. init) is the real initiation time signal
        real_init_times = []
        init_times = self.initTimes.copy()
        init_time = init_times[0]
        del init_times[0]
        for idx, dataPoint in enumerate(self.dataList):
            if dataPoint[0] == init_time:
                real_init_times.append(self.dataList[idx-1][0])
                if init_times:
                    init_time = init_times[0]
                    del init_times[0]
        # the following code is almost the same as the one in the current clamp class
        init_times = real_init_times.copy()
        init_times.append(self.dataList[-1][0])
        init_time = init_times[0]
        del init_times[0]
        simulated_vs_recorded_diff_at_initiation_points = []
        for t, signal, simulated_signal in zip(*zip(*self.dataList), self.simulatedSignal):
            if t == init_time:
                simulated_vs_recorded_diff_at_initiation_points.append(simulated_signal - signal)
                if init_times:
                    init_time = init_times[0]
                    del init_times[0]

        # start correcting data
        corrected_signal = []
        corrected_data = []
        init_times = real_init_times.copy()
        init_times.append(self.dataList[-1][0])
        init_time = init_times[0]
        del init_times[0]
        delta_init_t = init_times[0] - init_time
        delta_init_signal = simulated_vs_recorded_diff_at_initiation_points[0]
        del simulated_vs_recorded_diff_at_initiation_points[0]
        delta_init_signal_next = simulated_vs_recorded_diff_at_initiation_points[0]
        for t, signal in self.dataList:
            delta_t = t - init_time
            if t < init_time:
                corrected_signal.append(signal)
                corrected_data.append([t, corrected_signal[-1]])
                continue
            corrected_signal.append(
                signal +
                (delta_init_signal * (delta_init_t - delta_t) + delta_init_signal_next * delta_t) / delta_init_t
            )
            corrected_data.append([t, corrected_signal[-1]])
            if len(init_times) > 1 and t >= init_times[0]:
                init_time = init_times[0]
                del init_times[0]
                delta_init_t = init_times[0] - init_time
                delta_init_signal = simulated_vs_recorded_diff_at_initiation_points[0]
                del simulated_vs_recorded_diff_at_initiation_points[0]
                delta_init_signal_next = simulated_vs_recorded_diff_at_initiation_points[0]

        if len(corrected_signal) == len(self.dataList):
            self.dataList = corrected_data
        else:
            raise Exception
        return corrected_signal


class ExperimentCurrentClamp:
    def __init__(self, data_list, init_times, synaptic_reversal_potential, membrane_potential,
                 input_resistance, membrane_capacitance, bounds):
        self.BOUNDS = bounds
        self.dataList = data_list  # format [[t in ms, membrane potential in pA], ]
        self.initTimes = init_times  # [synTimes in ms]
        self.num_events = len(init_times) + 1
        self.synapticReversalPotential = float(synaptic_reversal_potential)  # in mV
        self.V0 = self.leakReversalPotential = float(membrane_potential)  # in mV
        # input_resistance in MegaOhm -> leak in nanoSiemens -> normalize by membrane_capacitance
        self.leakConductanceBar = 1/float(input_resistance) * 1e3/float(membrane_capacitance)
        self.membraneCapacitance = float(membrane_capacitance)  # in pF
        self.simulatedSignal = []
        self.g0, self.tau_d, self.tau_d, self.tau_f, self.tau_r, self.U = None, None, None, None, None, None
        amp = fabs(self.dataList[0][1] - self.dataList[1][1])
        self.normalization_value = (amp if amp != 0.0 else 1.0) * (self.num_events if self.num_events != 0 else 1.0)

    def input_parser(self, input_vec):
        # This function separated for overloading purpose
        self.g0 = input_vec[0] / self.membraneCapacitance
        self.tau_d = input_vec[1]
        self.tau_r = input_vec[2]
        self.tau_f = input_vec[3]
        self.U = input_vec[4]

    def interevent_signal(self, signal_init, delta_ts, g):
        signal = odeint(cell, signal_init, delta_ts,
                        args=(g, self.tau_d, self.leakConductanceBar, self.leakReversalPotential,
                              self.synapticReversalPotential)
                        ).ravel().tolist()
        del signal[0]
        return signal

    def simulate(self, input_vec):
        self.input_parser(input_vec)
        error, x0, y0, u0, delta_t = 0.0, 1.0, 0.0, 0.0, 0.0
        g, x0, y0, u0 = synaptic_event(delta_t, self.g0, self.tau_d, self.tau_r, self.tau_f, self.U, u0, x0, y0)
        init_times, data = self.initTimes.copy(), self.dataList.copy()
        init_time, signal, delta_ts, epsilon = init_times[0], [self.V0], [], (data[1][0] - data[0][0]) / 10  # ms
        del init_times[0]
        # add two time points before and after the pick of first signal to check the overshooting of simulatedSignal
        data.insert(2, [data[1][0] + epsilon, None])
        data.insert(1, [data[1][0] - epsilon, None])
        for dataPoint in data:
            delta_t = dataPoint[0] - init_time
            if delta_t < 0:
                signal.append(self.V0)
                continue
            delta_ts.append(delta_t)
            if init_times and dataPoint[0] >= init_times[0]:
                signal += self.interevent_signal(signal[-1], delta_ts, g)
                g, x0, y0, u0 = synaptic_event(delta_t, self.g0, self.tau_d, self.tau_r, self.tau_f, self.U, u0, x0, y0)
                init_time, delta_ts = init_times[0], [0.0]
                del init_times[0]
        signal += self.interevent_signal(signal[-1], delta_ts, g)
        # calculate error
        init_times = self.initTimes.copy()
        init_time, error_weight = init_times[0], 1.0
        del init_times[0]
        if len(signal) == len(data):
            for data_time, data_signal, model_signal in zip(*zip(*data), signal):
                if data_signal is None or isnan(data_signal):
                    continue
                if init_times and (init_time - 0.05) <= data_time < (init_time + 0.05):
                    init_time, error_weight = init_times[0], 1.0
                    del init_times[0]
                else:
                    error = soft_l1_loss(error, data_signal, model_signal, error_weight)
                    error_weight /= self.num_events
        else:
            raise Exception
        # penalize the data pick not match signal pick
        error += 0.0 if signal[1] <= signal[2] >= signal[3] else 1.0
        del signal[3], signal[1]
        error /= self.normalization_value
        if self.num_events == 2:
            error += (1 - self.U/self.BOUNDS[4][1]) * error
        self.simulatedSignal = signal
        return error

    def correct_data(self):
        # make an array of difference between simulated and recorded initiation times
        init_times = self.initTimes.copy()
        init_times.append(self.dataList[-1][0])
        init_time = init_times[0]
        del init_times[0]
        simulated_vs_recorded_diff_at_initiation_points = []
        for t, signal, simulated_signal in zip(*zip(*self.dataList), self.simulatedSignal):
            if t == init_time:
                simulated_vs_recorded_diff_at_initiation_points.append(simulated_signal - signal)
                if init_times:
                    init_time = init_times[0]
                    del init_times[0]

        # start correcting data
        corrected_signal = []
        corrected_data = []
        init_times = self.initTimes.copy()
        init_times.append(self.dataList[-1][0])
        init_time = init_times[0]
        del init_times[0]
        delta_init_t = init_times[0] - init_time
        delta_init_signal = simulated_vs_recorded_diff_at_initiation_points[0]
        del simulated_vs_recorded_diff_at_initiation_points[0]
        delta_init_signal_next = simulated_vs_recorded_diff_at_initiation_points[0]
        for t, signal in self.dataList:
            delta_t = t - init_time

            if t < init_time:
                corrected_signal.append(signal)
                corrected_data.append([t, corrected_signal[-1]])
                continue

            corrected_signal.append(
                signal +
                (delta_init_signal * (delta_init_t - delta_t) + delta_init_signal_next * delta_t) / delta_init_t
            )
            corrected_data.append([t, corrected_signal[-1]])

            if len(init_times) > 1 and t >= init_times[0]:
                init_time = init_times[0]
                del init_times[0]
                delta_init_t = init_times[0] - init_time
                delta_init_signal = simulated_vs_recorded_diff_at_initiation_points[0]
                del simulated_vs_recorded_diff_at_initiation_points[0]
                delta_init_signal_next = simulated_vs_recorded_diff_at_initiation_points[0]

        if len(corrected_signal) == len(self.dataList):
            self.dataList = corrected_data
        else:
            raise Exception
        return corrected_signal


class ExperimentCurrentClampCapacitance(ExperimentCurrentClamp):
    def __init__(self, data_list, init_times, synaptic_reversal_potential, membrane_potential,
                 input_resistance, membrane_capacitance, bounds):
        ExperimentCurrentClamp.__init__(self, data_list, init_times, synaptic_reversal_potential, membrane_potential,
                                        input_resistance, membrane_capacitance, bounds)
        self.BOUNDS = bounds

    def input_parser(self, input_vec):
        self.g0 = input_vec[0] / input_vec[5]
        self.tau_d = input_vec[1]
        self.tau_r = input_vec[2]
        self.tau_f = input_vec[3]
        self.U = input_vec[4]
        self.leakConductanceBar = self.leakConductanceBar * self.membraneCapacitance / input_vec[5]
        self.membraneCapacitance = input_vec[5]


class ExperimentCurrentClampResistance(ExperimentCurrentClamp):
    def __init__(self, data_list, init_times, synaptic_reversal_potential, membrane_potential,
                 input_resistance, membrane_capacitance, bounds):
        ExperimentCurrentClamp.__init__(self, data_list, init_times, synaptic_reversal_potential, membrane_potential,
                                        input_resistance, membrane_capacitance, bounds)
        self.BOUNDS = bounds
        self.inputResistance = float(input_resistance)

    def input_parser(self, input_vec):
        self.g0 = input_vec[0] / self.membraneCapacitance
        self.tau_d = input_vec[1]
        self.tau_r = input_vec[2]
        self.tau_f = input_vec[3]
        self.U = input_vec[4]
        self.leakConductanceBar = self.leakConductanceBar * self.inputResistance / input_vec[5]
        self.inputResistance = input_vec[5]


class ExperimentCurrentClampRinCm(ExperimentCurrentClamp):
    def __init__(self, data_list, init_times, synaptic_reversal_potential, membrane_potential,
                 input_resistance, membrane_capacitance, bounds):
        ExperimentCurrentClamp.__init__(self, data_list, init_times, synaptic_reversal_potential, membrane_potential,
                                        input_resistance, membrane_capacitance, bounds)
        self.BOUNDS = bounds
        self.inputResistance = float(input_resistance)

    def input_parser(self, input_vec):
        self.g0 = input_vec[0] / input_vec[6]
        self.tau_d = input_vec[1]
        self.tau_r = input_vec[2]
        self.tau_f = input_vec[3]
        self.U = input_vec[4]
        self.leakConductanceBar = \
            self.leakConductanceBar * self.inputResistance / input_vec[5] * self.membraneCapacitance / input_vec[6]
        self.inputResistance = input_vec[5]
        self.membraneCapacitance = input_vec[6]


class ScrollableFrame(Tk):
    # this class adds the "mainFrame" which is scrollable
    def __init__(self, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)  # make a window
        self.mainCanvas = Canvas(self, borderwidth=0)
        self.vertical_scroll_bar = Scrollbar(self, orient='vertical', command=self.mainCanvas.yview)
        self.mainCanvas.configure(yscrollcommand=self.vertical_scroll_bar.set)
        self.mainFrame = Frame(self.mainCanvas)
        self.mainWindow = self.mainCanvas.create_window(
            (0, 0), window=self.mainFrame, anchor='nw', tags="self.mainFrame")
        self.mainFrame.bind("<Configure>", self._on_main_frame_configure)
        self.mainCanvas.bind("<Configure>", self._on_main_canvas_configure)
        self.mainCanvas.bind_all("<MouseWheel>", self._on_mousewheel)

    def display_scrollable_frame(self):
        self.mainCanvas.grid(row=0, column=0, sticky="NSWE")
        self.vertical_scroll_bar.grid(row=0, column=1, sticky='NS')
        self.mainFrame.columnconfigure(0, weight=1)

    def _on_main_frame_configure(self, event):
        # Reset the scroll region to encompass the inner frame
        self.mainCanvas.configure(yscrollcommand=self.vertical_scroll_bar.set, scrollregion=self.mainCanvas.bbox("all"))

    def _on_main_canvas_configure(self, event):
        # Resize the inner frame to match the canvas
        min_height = self.mainFrame.winfo_reqheight()
        if self.winfo_height() >= min_height:
            new_height = self.winfo_height()
            # Hide the scrollbar when not needed
            self.vertical_scroll_bar.grid_remove()
        else:
            new_height = min_height
            # Show the scrollbar when needed
            self.vertical_scroll_bar.grid()

        self.mainCanvas.itemconfig(self.mainWindow, height=new_height, width=self.winfo_width())

    def _on_mousewheel(self, event):
        self.mainCanvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

    def update_scrollbar(self):
        self.mainCanvas.update()
        self._on_main_canvas_configure(0)
        self._on_main_frame_configure(0)
        self.mainCanvas.yview_moveto(1)


class EntryWithPlaceholder(Entry):
    def __init__(self, master=None, placeholder="PLACEHOLDER", color='blue', *args, **kwargs):
        # super(EntryWithPlaceholder, self).__init__(master)
        Entry.__init__(self, master, *args, **kwargs)
        self.placeholder = placeholder
        self.placeholder_color = color
        self.default_fg_color = self['fg']
        self.bind("<FocusIn>", self.foc_in)
        self.bind("<FocusOut>", self.foc_out)
        self.bind("<ButtonRelease-3>", self.copy_to_clipboard)
        self.put_placeholder()

    def put_placeholder(self, *placeholder):
        if placeholder:
            self.placeholder = placeholder[0]
        self.set(self.placeholder)
        self['fg'] = self.placeholder_color

    def foc_in(self, *args):
        if self['fg'] == self.placeholder_color:
            self.delete('0', 'end')
            self['fg'] = self.default_fg_color

    def foc_out(self, *args):
        if not self.get():
            self.put_placeholder()

    def copy_to_clipboard(self, event=None):
        try:
            self.clipboard_clear()
            self['fg'] = self.default_fg_color
            sleep(0.1)
            self.clipboard_append(self.get_().strip())
        except TclError:
            pass

    def set(self, txt):
        state = self['state']
        self['state'] = 'normal'
        self.delete('0', 'end')
        self.insert(0, txt)
        self['fg'] = self.default_fg_color
        self['state'] = state

    def get_(self):
        if self['fg'] == self.placeholder_color:
            return ''
        else:
            return self.get()


class Main(ScrollableFrame):
    CSVs_FOLDER = 'csvs'

    @staticmethod
    def chunks(l, n):
        """Yield successive n-sized chunks from l."""
        for i in range(0, len(l), n):
            yield l[i:i + n]

    def __init__(self, *args, **kwargs):
        Experiment.CSVs_FOLDER = Main.CSVs_FOLDER
        # make a scrollable Frame
        ScrollableFrame.__init__(self, *args, **kwargs)
        # Window size and position
        self.title('Synapse Modelers Workshop')
        self.state('zoomed')  # maximize the window
        pad = 200  # the pad size of restored window with respect to full screen
        self.geometry("{0}x{1}+{2}+{3}".format(
            self.winfo_screenwidth()-pad,   # restored window width
            self.winfo_screenheight()-pad,  # restored window height
            int(pad/2),                     # restored window position on the center
            int(pad/2))                     # restored window position on the center
        )
        # make content match the frame size
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        # add a menu bar and file menu
        menu = Menu(self)
        start_menu = Menu(menu, tearoff=0)
        start_menu.add_command(label="Open CSV files", command=self.open)
        start_menu.add_command(label="Create a pseudo-trace", command=self.pseudo_trace)
        menu.add_cascade(label="Start", menu=start_menu)
        menu.add_command(label="Options", command=self.options)
        menu.add_command(label="Calculator", command=self.calculator)
        menu.add_command(label="Close-the-last", command=self.close_last)
        self.bind('<Control-w>', self.close_last)

        self.config(menu=menu)

        self.working_directory = getcwd()
        while not path.isdir(self.working_directory + '/' + Main.CSVs_FOLDER):
            self.working_directory = filedialog.askdirectory(
                initialdir=".",
                title="Choose Working Directory: should have a \"" + Main.CSVs_FOLDER + "\" folder")

        # setup the experiments
        Experiment.WORKING_DIRECTORY = self.working_directory
        Experiment.CSVs_FOLDER = Main.CSVs_FOLDER
        self.experiments = []
        # csv_file_names_list = self._get_all_csv_file_names()
        # for file_name in csv_file_names_list:
        #     self.append_experiment(file_name)

        # Show the frame and scrollbar
        self.display_scrollable_frame()

    def _get_all_csv_file_names(self):
        data = []
        csv_files = listdir(self.working_directory + '/' + Main.CSVs_FOLDER)
        for csv_file in csv_files:
            file_name, file_ext = path.splitext(csv_file)
            if file_ext.lower() != '.csv':
                continue
            data.append(file_name)
        return data

    def append_experiment(self, file_name):
        try:
            self.experiments.append(Experiment(self.mainFrame, file_name))
        except TypeError:
            ok = messagebox.askokcancel(
                title='Error',
                message='The %s.csv file is probably containing multiple traces\n' % file_name +
                        'Do you want to try to split it into multiple files?')
            if ok:
                try:
                    path_ = Experiment.WORKING_DIRECTORY + '\\' + Experiment.CSVs_FOLDER + '\\'
                    for sections in filter(None, split(r'\n{2,}', open(path_ + file_name + '.csv').read())):
                        col_name = findall(r'^.+?,(.+?)\n', sections)[0]
                        new_file_name = file_name + '-' + col_name
                        open(path_ + new_file_name + '.csv', 'w').write(sections)
                        self.experiments.append(Experiment(self.mainFrame, new_file_name))
                    self.update_scrollbar()
                    Experiment.delete_file_(
                        file_name, '.csv',
                        'File splitting was successful.\n'
                        'Do you want to delete %s.csv file?' % file_name)
                except Exception as e:
                    print(e)
                    messagebox.showerror(
                        title='Error',
                        message='Please check the file:\n%s.csv' % file_name
                    )

    def open(self):
        path_file_names = filedialog.askopenfilenames(
            initialdir=Main.CSVs_FOLDER,
            title="Choose CSV files",
            filetypes=[('csv', '*.csv')]
        )
        open_file_names = [e.FILE_NAME for e in self.experiments]
        for chunk in Main.chunks(path_file_names, 50):
            for path_file_name in chunk:
                file_name, extension = path.splitext(path.basename(path_file_name))
                if file_name in open_file_names:
                    messagebox.showwarning('warning', 'File %s was already open' % file_name)
                else:
                    self.append_experiment(file_name)
                self.update_scrollbar()

    def pseudo_trace(self):
        window = Toplevel(self)
        window.title('Pseudo Trace Generator')
        window.columnconfigure(0, weight=1)
        types, entries, row = dict(), dict(), 0
        row_values = [
            {'type': 'Rise', 'title': 'time', 'unit': 'ms',
             'options': {'10-90%', '0-100%', '20-80%', 'time constant', 'unspecified'}},
            {'type': 'Decay', 'title': 'time', 'unit': 'ms',
             'options': {'100-0%', '100-37%', '100-50%', '100-63%', '90-37%', 'unspecified'}},
            {'type': 'Signal', 'title': 'potency', 'unit': 'mV or pA', 'options': None},
            {'type': 'ISI', 'title': 'time', 'unit': 'ms', 'options': None},
            {'type': 'PPR', 'title': '2/1 amplitude', 'unit': 'ratio', 'options': None}
        ]
        for row, value in enumerate(row_values):
            entries[value['type']] = EntryWithPlaceholder(
                master=window,
                placeholder='%s %s value (%s)' % (value['type'], value['title'], value['unit']),
                font=('Courier New', 16),
                width=35
            )
            entries[value['type']].grid(row=row, column=0, sticky='NEWS')
            if value['options']:
                types[value['type']] = StringVar(window)
                types[value['type']].set(next(iter(value['options'])))
                OptionMenu(window, types[value['type']], *value['options']).grid(row=row, column=1, sticky='NEWS')

        def make_save_pseudo_trace():
            inter_stimulus_intervals, paired_pulse_ratios = [None], [None]
            all_dfs = DataFrame({'time': [], 'signal': []})
            if entries['ISI'].get_() and entries['PPR'].get_():
                inter_stimulus_intervals = map(float, split(r'\s*[,; ]+\s*', entries['ISI'].get_()))
                paired_pulse_ratios = map(float, split(r'\s*[,; ]+\s*', entries['PPR'].get_()))
            for idx, (isi, ppr) in enumerate(zip(inter_stimulus_intervals, paired_pulse_ratios), start=0):
                df = DataFrame({'time': [0], 'signal': [0]})
                t_rise, a_rise, t_decay, a_decay = None, None, None, None
                rise_conversion_factor = {
                    '10-90%': 0.8, '0-100%': 1.0, '20-80%': 0.6, 'time constant': 0.632, 'unspecified': 0.6}
                t_rise = float(entries['Rise'].get_()) / rise_conversion_factor.get(types['Rise'].get())
                a_rise = float(entries['Signal'].get_())
                df.loc[1] = [t_rise, a_rise]
                decay_conversion_factor = {
                    '100-0%': 0.0, '100-37%': 0.37, '100-50%': 0.5, '100-63%': 0.63, '90-37%': 0.37, 'unspecified': 0.0}
                t_decay = t_rise + float(entries['Decay'].get_()) * (
                    1 if types['Decay'].get() != '90-37%' else 0.63 / 0.57)
                a_decay = a_rise * decay_conversion_factor.get(types['Decay'].get())
                df.loc[2] = [t_decay, a_decay]
                if entries['ISI'].get_() and entries['PPR'].get_():
                    tau = (t_decay - t_rise) / (log(fabs(a_rise)) - log(fabs(a_decay)))
                    df.loc[2] = [(isi + t_rise)/2, a_rise * exp(-(isi - t_rise) / 2 / tau)]
                    df.loc[3] = [isi, a_rise * exp(-(isi-t_rise)/tau)]
                    df.loc[4] = [isi + t_rise, df.loc[3].signal + a_rise * ppr]
                    df.loc[5] = [isi + t_decay, df.loc[4].signal * exp(-(t_decay - t_rise) / tau)]
                df.time = df.time + 60000 * idx
                print(df)
                all_dfs = all_dfs.append(df)

            f = filedialog.asksaveasfile(defaultextension=".csv")
            if f is None:  # asksaveasfile return `None` if dialog closed with "cancel".
                return
            all_dfs.to_csv(f, index=False)
            # window.destroy()
        make_save_pseudo_trace.window = window
        make_save_pseudo_trace.entries = entries
        make_save_pseudo_trace.types = types

        Button(window, command=make_save_pseudo_trace, text='Save CSV').grid(row=row, column=1, sticky='NEWS')

    def options(self):
        window = Toplevel(self)
        window.grab_set()  # make the main window unclickable until closing the settings window
        Label(window, text='Current-clamp process numbers:').grid(
            row=0, column=0, sticky="W")
        Scale(window, from_=1, to=8, variable=thread_numbers, orient='horizontal').grid(
            row=0, column=1, sticky="W")
        Label(window, text='Optimizer population size:').grid(
            row=1, column=0, sticky="W")
        Scale(window, from_=15, to=150, variable=population_size, orient='horizontal').grid(
            row=1, column=1, sticky="W")
        Label(window, text='Bootstrap max iteration').grid(
            row=2, column=0, sticky="W")
        Scale(window, from_=2, to=100, variable=bootstrap_max_iteration, orient='horizontal').grid(
            row=2, column=1, sticky="W")
        Label(window, text='Enable bell sound:').grid(
            row=3, column=0, sticky="W")
        Checkbutton(window, variable=ring_bell_when_optimization_ends).grid(
            row=3, column=1, sticky="W")
        Label(window, text='Print each optimization results:').grid(
            row=4, column=0, sticky="W")
        Checkbutton(window, variable=print_results_when_optimization_ends).grid(
            row=4, column=1, sticky="W")
        Label(window, text='Show matrix plot in summary:').grid(
            row=5, column=0, sticky="W")
        Checkbutton(window, variable=show_matrix_plot).grid(
            row=5, column=1, sticky="W")

    def calculator(self):

        def stats(list_, type_='SEM', comment=''):
            list_ = list(
                map(float, split(
                    r"(?:(?:\s*\u00B1\s*\d+(?:\.\d+)?(?:e\d+)?\s*)?(?:\s*\[\d+(?:\.\d+)?(?:e\d+)? to \d+(?:\.\d+)?(?:e"
                    r"\d+)?\]\s*)?(?:\s*\(n=\d+\)\s*)?(?:\s*{.*?}\s*)?(?:\s*@\d+\s*(?:&\d+)*)?(?:\s*{.*?}\s*)?)?\s*[;,]"
                    r"\s*", list_)))
            a, n = np.array(list_), len(list_)
            (iqr25, iqr75), mst = np.percentile(a, [75, 25]), a.std()
            if type_ == 'IQR':
                mst = iqr75 - iqr25
            elif type_ == 'SEM':
                mst = a.std() / n ** 0.5
            comment += ':' if comment else ''
            stat_ = '%.2f\u00B1%.2f [%.2f to %.2f] (n=%d) {%scalculated=(%s)/%d}' % (
                np.median(a) if type_ == 'IQR' else a.mean(), mst, a.min(), a.max(), n, comment,
                "+".join(map(str, a)), n)
            return stat_

        def u_calc(failure_, mse):
            inputs = tuple(map(
                float,
                filter(None, findall(
                    r'^\s*(\d+(?:\.\d+)?(?:e-?\d+)?)(?:\s*\u00B1\s*(\d+(?:\.\d+)?(?:e-?\d+)?).*?\(n=(\d+)\))?\s*$',
                    failure_)[0])))
            if len(inputs) == 3:
                mean, mse_v, n = inputs
            else:
                mean, mse_v, n = inputs[0], 0.0, 1.0
            if mse == 'SEM':
                mse_v = mse_v / 100
            elif mse == 'SD':
                mse_v = mse_v / n ** 0.5 / 100
            elif mse == 'IQR':
                mse_v = mse_v / n ** 0.5 / 1.35 / 100
            mean = 1-mean/100
            conf_80 = 1.282 * mse_v
            conf_95 = 1.960 * mse_v
            return '80%% = %.2f [%.2f to %.2f], 95%% = %.2f [%.2f to %.2f]' % (
                mean, mean-conf_80, mean+conf_80, mean, mean-conf_95, mean+conf_95)

        window = Toplevel(self.mainFrame)
        msp_frame = Frame(window)
        Label(msp_frame, text='Measure of Spread (MS):').grid(row=0, column=0, sticky='W')
        msp = StringVar(window)
        msp.set(global_msp)
        msp.trace('w', lambda *args: globals().update(global_msp=msp.get()))
        OptionMenu(msp_frame, msp, *{'SEM', 'SD', 'IQR'}).grid(row=0, column=1, sticky='W')
        Label(msp_frame, text='Stats Inputs:', fg='green').grid(row=1, column=0, sticky='W')
        msp_frame.grid(row=0, column=0, sticky='W')

        width = 100
        half_width = int(width / 2)
        csv_ = EntryWithPlaceholder(window, placeholder='Comma separated values', width=width)
        csv_.grid(row=1, column=0, sticky='NEWS')
        csv_.grid_columnconfigure(0, weight=1)
        note = EntryWithPlaceholder(window, placeholder='Endnote', width=width)
        note.grid(row=2, column=0, sticky='NEWS')
        note.grid_columnconfigure(0, weight=1)
        stat_frame = Frame(window)
        Button(stat_frame, text="Stats=", command=lambda: stat.set(stats(csv_.get(), msp.get(), note.get_()))).grid(
            row=0, column=0, sticky='NEWS')
        stat = EntryWithPlaceholder(stat_frame, state='readonly', placeholder='')
        stat.grid(row=0, column=1, sticky='NEWS')
        stat_frame.grid(row=3, column=0, sticky='NEWS')
        stat_frame.grid_columnconfigure(1, weight=1)

        cm_frame_upper = Frame(window)
        Label(cm_frame_upper, text='Membrane Properties Inputs:', fg='green').grid(row=0, column=0, sticky='W')
        tau_ = EntryWithPlaceholder(cm_frame_upper, placeholder='Membrane time constant (ms)', width=half_width)
        tau_.grid(row=1, column=0, sticky='NEWS')
        r_in = EntryWithPlaceholder(cm_frame_upper, placeholder='Membrane input resistance (M\u2126)', width=half_width)
        r_in.grid(row=1, column=1, sticky='NEWS')
        cm_frame_upper.grid(row=4, column=0)
        cm_frame_lower = Frame(window)
        cm__ = EntryWithPlaceholder(cm_frame_lower, state='readonly', placeholder='')
        cm__.grid(row=0, column=1, sticky='NEWS')
        Button(cm_frame_lower, text="Membrane capacitance (PF)=", width=22,
               command=lambda: cm__.set(
                   round(float(tau_.get())/float(r_in.get())*1000, Experiment.DECIMAL_POINTS))).grid(row=0, column=0)
        cm_frame_lower.grid(row=5, column=0, sticky='NEWS')
        cm_frame_lower.grid_columnconfigure(1, weight=1)

        u_frame = Frame(window)
        Label(u_frame, text='U Calculator Inputs:', fg='green').grid(row=0, column=0, sticky='W')
        failure = EntryWithPlaceholder(master=u_frame, placeholder='% Failure MCT\u00B1MS(n=)', width=22)
        failure.grid(row=1, column=0, sticky='NEWS')
        u_value = EntryWithPlaceholder(master=u_frame, placeholder='', state='readonly')
        u_value.grid(row=1, column=2, sticky='NEWS')
        Button(u_frame, text="U=", command=lambda: u_value.set(u_calc(failure.get(), msp.get()))).grid(row=1, column=1)
        u_frame.grid(row=6, column=0, sticky='NEWS')
        u_frame.grid_columnconfigure(2, weight=1)

    def close_last(self, event=None):
        if messagebox.askokcancel('Close', 'Do you want to close the last experiment'):
            self.experiments[-1].experimentFrame.grid_forget()
            self.experiments[-1].experimentFrame.destroy()
            del self.experiments[-1]
            Experiment.rowCount -= 1
            self.update_scrollbar()


@jit(nopython=True, fastmath=True, cache=True)
def synaptic_event(delta_t, g0, tau_d, tau_r, tau_f, u, u0, x0, y0):
    tau1r = tau_d / ((tau_d - tau_r) if tau_d != tau_r else 1e-13)
    y_ = y0 * exp(-delta_t / tau_d)
    x_ = 1 + (x0 - 1 + tau1r * y0) * exp(-delta_t / tau_r) - tau1r * y_
    u_ = u0 * exp(-delta_t / tau_f)
    u0 = u_ + u * (1 - u_)
    y0 = y_ + u0 * x_
    x0 = x_ - u0 * x_
    g = g0 * y0
    return g, x0, y0, u0


@jit(nopython=True, fastmath=True, cache=True)
def synaptic_current(g, delta_t, tau_d, e_syn):
    return g * exp(-delta_t / tau_d) * e_syn


@jit(nopython=True, fastmath=True, cache=True)
def cell(v, t, g, tau_d, g_leak, e_leak, e_syn):
    return g_leak * (e_leak - v) + g * exp(-t / tau_d) * (e_syn - v)


@jit(nopython=True, fastmath=True, cache=True)
def soft_l1_loss(error, data_signal, model_signal, error_weight):  # smoothed least square function
    return error + (2.0*((1.0+(data_signal - model_signal)**2.0)**0.5)-1.0)*error_weight


# Required for multi-threading on windows
if __name__ == '__main__':
    freeze_support()
    app = Main()
    thread_numbers = IntVar()
    thread_numbers.set((cpu_count(logical=False) - 1) if cpu_count(logical=False) == cpu_count(logical=True) else
                       cpu_count(logical=False) + cpu_count(logical=False)/2)
    population_size = IntVar()
    population_size.set(15)
    bootstrap_max_iteration = IntVar()
    bootstrap_max_iteration.set(29)
    ring_bell_when_optimization_ends = BooleanVar()
    ring_bell_when_optimization_ends.set(True)
    print_results_when_optimization_ends = BooleanVar()
    print_results_when_optimization_ends.set(False)
    show_matrix_plot = BooleanVar()
    show_matrix_plot.set(False)
    global_msp = 'SEM'
    global_mct = 'Mean'
    app.mainloop()