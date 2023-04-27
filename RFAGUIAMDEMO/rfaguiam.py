#!/usr/bin/python3
# rfaguiam.py
# D. Gibbon
# Rhythm Formant Analysis
# Amplitude Demodulation
# Alpha version (partially tested)
# 2023-03-20 V01
"""
Description

The "rfaguiam.py" application illustrates the acoustic analyses which are used by Rhythm Formant Analysis, the methodology which accompanies Rhythm Formant Theory. The application is designed for teaching purposes, both for illustration and for first steps in computational phonetics with Python. Functional programming is used. The GUI uses a TkInter loop.
The software is distributed as FOSS (Free and Open Source Software) and may be freely modified on condition that this file is cited as source.
For additional software see
	https://github.com/dafyddg/RFA
	
	https://github.com/dafyddg/RFAGUI
Source:
D. Gibbon, “The rhythms of rhythm”,
	Online: JIPA First View, pp. 1–33, 2021.
	Print and online: JIPA Volume 53 , Issue 1 , April 2023 , pp. 233 - 265
	DOI: https://doi.org/10.1017/S0025100321000086
© The Author(s), 2021. Published by Cambridge University Press on behalf of the International Phonetic Association
@article{gibbon_2023, title={The rhythms of rhythm}, volume={53}, DOI={10.1017/S0025100321000086}, number={1}, journal={Journal of the International Phonetic Association}, publisher={Cambridge University Press}, author={Gibbon, Dafydd}, year={2023}, pages={233–265}}

Implementation:
The RFA implementation is design
"""

import re, sys, datetime
import numpy as np
from scipy.io import wavfile
from scipy.signal import butter, lfilter, medfilt, find_peaks
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk

####################################################################
# Internally fixed values

argvlen = len(sys.argv)

figwidth = 10
figheight = 6
rows = 4
columns = 1
fontsize = 10
linewid = 1
envskip = 400

specheatmaptype = "YlOrRd"
heatmapdotsize = 60
waterfalldotsize = 20

speccolor = "g"
peakcolor = "r"

application = re.sub(".*[/]","",sys.argv[0])

datetimeInstance = datetime.datetime.today()
dateInstance = datetimeInstance.date()
date = re.sub("-","",str(dateInstance))

####################################################################
# create Tkinter window

root = tk.Tk()
root.title("WAV File Signal Display [%s]"%application)

####################################################################
# My rounding function

def roundint(num):
    x = round(num)
    y = int(num)
    z = y + 0.5
    if num >= z:
        return y + 1
    else:
        return y

####################################################################
# Input WAV file (mono)

def read_wav():

    global argvlen, file_path, filebase, fs, signal, siglen, sigdur, nykvist, resolution, hzperpt, signalmessage

    if argvlen > 1:
        file_path = sys.argv[1]
    else:
        file_path = tk.filedialog.askopenfilename(filetypes=[("Audio files","*.wav")])

    filebase = re.sub(".*/","",file_path)

    filebox.delete(0, tk.END)
    filebox.insert(0, filebase)

    # read in WAV file using scipy
    try:
        fs, signal = wavfile.read(file_path)
    except:
        print("WAV filename error."); exit()

    if len(signal.shape) > 1:
        print("Mono WAV file required."); exit()

    siglen = len(signal)

    nykvist = int(fs/2)
    sigdur = siglen/fs
    resolution = int(round(sigdur))
    hzperpt = 1/resolution
    statusboxrefresh("Check/set parameters, then click RFA. Then wait...")

    if argvlen > 1:
        argvlen = 0
        run_program()

################################################################

def process_input():

    global signal, envelope, signalmessage, envsmooth, specsmooth, envsmo

###################################################
# Envelope smoothing

    envsmo = int(envsmooth.get())
    if envsmo < 0: envsmo = 0-envsmo
    if envsmo % 2 == 0: envsmo += 1

# AM demodulation
    envelope = np.abs(signal)

# Signal scaling to -1...0...1
    envmin = np.min(envelope); envmax = np.max(envelope)
    signal = (signal - envmin) / (envmax - envmin)

# Moving median filter    
    envelope = medfilt(envelope, envsmo)

# Envelope scaling to 0...1
    envmin = np.min(envelope); envmax = np.max(envelope)
    envelope = (envelope - envmin) / (envmax - envmin)


################################################################

def createspecmatrix():
    global lfspecmatrix, lffrequencies, sigwin, sampskip, maxfreq, minfreq

    minfreq = float(minfrequency.get())
    maxfreq = float(maxfrequency.get())

    sigwin = int(sigwindowing.get())
    sigwindow = sigwin * fs

    sampskip = float(sampskipping.get()) * fs
    sampskip = int(sampskip)

    # Extend the signal with a buffer the length of the window.
    siglist = list(signal) + [0]*sigwindow
    sig = np.array(siglist)

    siglen = len(sig)            # get new length; new duration not needed?
    nykvist = roundint(fs/2)
    lfspecmatrix = []
    for i in range(0,len(sig)-sigwindow, sampskip):
        siglet = sig[i:i+sigwindow]
#        signal = signal**2
        envlet = abs(siglet)
        siglet = siglet / np.max(envlet)    # scale -1...0...1

        spectrum = abs(np.fft.rfft(envlet))**2    # FFT power spectrum
        spectrum[0] = spectrum[1]                 # reduce DC
        speclen = len(spectrum)

        frequencies = np.linspace(0,nykvist,speclen)

        lflen = roundint(maxfreq*speclen/nykvist)    # LF spectrum
        lfspectrum = spectrum[:lflen]
        lffrequencies = frequencies[:lflen]

        lfspecmin = np.min(lfspectrum)
        lfspecmax = np.max(lfspectrum)
        lfspectrum = (lfspectrum-lfspecmin)/(lfspecmax-lfspecmin) # Scale 0...1
        lfspecmatrix += [ lfspectrum ]

    lfspecmatrix = np.asarray(lfspecmatrix)

    return


####################################################################

def createlfspectrum():

    global spectrum, lfspectrum, speclen, lflen, peaks, peakheight, peakmags, specsmo
    print("... createlfspectrum")
    minfreq = float(minfrequency.get())
    maxfreq = float(maxfrequency.get())

    specsmo = int(float(specsmooth.get()))	# Correct to odd values
    if specsmo < 0: specsmo = 0-specsmo
    if specsmo % 2 == 0: specsmo += 1
    specsmooth.delete(0,tk.END)
    specsmooth.insert(0, specsmo)

    peakheight = float(peakheightval.get())

    spectrum = abs(np.fft.rfft(envelope))**2    # FFT power spectrum
    spectrum[0] = spectrum[1]    # remove DC
#    spectrum = spectrum**2                # use instead of abs!
    speclen = len(spectrum)

    lflen = int(round(maxfreq*speclen/nykvist))    # LF spectrum
    lfspectrum = spectrum[:lflen]

    print("... filtering LF spectrum")
    lfspectrum = medfilt(lfspectrum, specsmo)
    print("... done filtering LF spectrum")

    lfspecmin = np.min(lfspectrum)
    lfspecmax = np.max(lfspectrum)
    lfspectrum = (lfspectrum-lfspecmin)/(lfspecmax-lfspecmin)
    peaks, peakheights = find_peaks(lfspectrum, height=peakheight)
    peaks = np.array(peaks)/resolution
    peakmags = peakheights["peak_heights"]

####################################################################
####################################################################
# Figure

def waveformplot():

    axwave.set_title("RFA [file: %s]"%filebase, fontsize=fontsize)

    print("... waveformplot")
    x = np.linspace(0,sigdur, siglen)
    axwave.plot(x, signal, linewidth=linewid, color="lightgrey")

    # Envelope thinning for display purposes
    envskip = int(20 * sigdur)
    xx = x[::envskip]
    env = envelope[::envskip]
    env = (env - np.min(env))/(np.max(env)-np.min(env))
    axwave.plot(xx, env, color="r")

    axwavexticks = axwave.get_xticks()[1:-1]    # "grid" didn't work
    for i in axwavexticks:
        axwave.axvline(x=i, linewidth=linewid)

    axwavexlabels = [0] + [ "%.1fs"%x for x in axwavexticks ]
    axwave.set_xticklabels(axwavexlabels)

    axwaveyticks = axwave.get_yticks()
    for i in axwaveyticks:
        axwave.axhline(y=i, linewidth=linewid)

    axwave.set_ylim(-1.1, 1.1)
    axwave.set_ylabel('Waveform\noscillogram &\nAM envelope\n$Amplitude$')
    print("... waveformplot done")

#################################################################
# Heatmap spectrogram

def heatmapspectrogramplot():

    global sigwindow, sampskipper

    print("... heatmap spectrogram plot")

    sigwin = int(sigwindowing.get())
    sigwindow = sigwin * fs

    sampskipper = float(sampskipping.get())
    sampskip = int(sampskipper*fs)

    lfspecmatrix01 = lfspecmatrix[:-int(sigwindow/sampskip)] # calculate this

# Heatmap spectrogram
    if True:
    	# x-axis as signal time range, number of spectra in spectrogram
    	x = np.linspace(0,sigdur,len(lfspecmatrix01))
    	# Colormap is derived from magnitudes at each frequency/spectrum
    	for t, magvals in zip(x, lfspecmatrix01):
    		xx = [t] * len(lffrequencies)
    		axheat.scatter(
	    		xx,lffrequencies, c=magvals, cmap=specheatmaptype,
	    		marker="s", s=heatmapdotsize)

    axheatxticks = axheat.get_xticks()[1:-1]    # "grid" didn't work
    for i in axheatxticks:
        axheat.axvline(x=i, linewidth=linewid)

    axheatxlabels = [0] + [ "%.1fs"%x for x in axheatxticks ]
    axheat.set_xticklabels(axheatxlabels)

    axheat.set_yticks(np.arange(0, maxfreq+0.5, maxfreq/5))
    axheat.set_ylabel("Heatmap\nspectrogram\n\n$Frequency$ (Hz)")

    lfspecmatrix01 = []
    
    axheat.grid()

    print("... heatmap spectrogram plot done")

#################################################################
# Waterfall spectrogram

def waterfallspectrogramplot():

    global waterfallskip

    print("... waterfall spectrogram plot")

    waterfallskip = int(waterfallskipping.get())
    waterfallmultiplier = float(waterfallmultiplying.get())

	# Create shorter matrix by skipping waterfallskip intervals
    lfspecmatrix02 = lfspecmatrix[::waterfallskip]
    x = np.linspace(0,maxfreq, lfspecmatrix02.shape[1])
    linewid02 = float(linewid)/2

    specpairs = []
    for i, spec in enumerate(lfspecmatrix02):
        displacement = spec-i*waterfallmultiplier
        axwater.plot(x, displacement, linewidth=linewid02)

        specmax = np.max(displacement)
        freq = list(displacement).index(specmax)
        axwater.scatter(x[freq], specmax, s=waterfalldotsize, color="r")

    xx = np.linspace(0, maxfreq, lflen)
    axwater.set_xticks(np.arange(0, maxfreq+0.2, float("%.2f"%(maxfreq/10))))

    axwaterxticks = axwater.get_xticks()
    for t in axwaterxticks:
        axwater.axvline(t)

    axwaterxlabels = [ "%.2fHz"%x for x in axwaterxticks ]
    axwater.set_xticklabels(axwaterxlabels)

    axwater.set_yticks([])
    axwater.set_ylabel("Waterfall\nspectrogram\n\n<--$Time$--\n\n%.2fs"%sigdur)

    print("... waterfall spectrogram plot done")

####################################################################
# Spectrum plot

def lfspectrumplot():

    xx = np.linspace(0, maxfreq, lflen)
    axspect.set_xticks(np.arange(0, maxfreq+0.2, float("%.2f"%(maxfreq/5))))

    axspectxticks = axspect.get_xticks()
    for t in axspectxticks:
        axspect.axvline(t)

    axspectxlabels = [ "%.2fHz"%x for x in axspectxticks ]

    axspectyticks = axspect.get_yticks()
    for i in axspectyticks:
        axspect.axhline(y=i, linewidth=1, color="g")

    speccolor = str(speccoloring.get())
    axspect.plot(xx, lfspectrum, linewidth=linewid, color=speccolor)

    axspect.set_xticklabels(axspectxlabels)

    for p, m in zip(peaks,peakmags):
        if p >= 0.5:
            axspect.scatter(p+0.015, m, color=peakcolor)    # peak correction
            axspect.text(p-0.05, m-0.15, "%.2fHz"%p)
    
    axspect.set_ylabel('Global\nspectrum\n$Energy$\n(rescaled)')
    axspect.axhline(y=0.7, linewidth=1, linestyle=":", color="blue")

####################################################################

def creategraphics():

    global axwave, axf0, axheat, axwater, axspect

    plt.clf()

    plotcount = 0

    plotcount += 1
    axwave = fig.add_subplot(rows, columns, plotcount)
    waveformplot()
    """
    plotcount += 1
    axf0 = fig.add_subplot(rows, columns, plotcount)
    f0plot()
    """
    plotcount += 1
    axheat = fig.add_subplot(rows, columns, plotcount)
    heatmapspectrogramplot()

    plotcount += 1
    axwater = fig.add_subplot(rows, columns, plotcount)
    waterfallspectrogramplot()

    plotcount += 1
    axspect = fig.add_subplot(rows, columns, plotcount)
    lfspectrumplot()

# Save graphics as file
    plt.tight_layout(pad=1, w_pad=0, h_pad=0.5)

    figfullname = "FIGURES/%s-%s-%d-%d-%d-%d-%d-%d-%d-%.1f"%(filebase, date, envsmo, minfreq, maxfreq, sigwin, sampskipper, waterfallskip, specsmo, peakheight)

    paramnames = """The filenames contain a value list for the following parameters:
    filebase, date,
    envsmo (envelope smoothing),
    minfreq, maxfreq (LF spectrum frequency bounds),
    sigwin, sampskipper (LF spectrogram parameters),
    waterfallskip, specsmo (spectrum smoothing),
    peakheight (threshold for peak marking)
    This file is written with the figures after each RFA call.
    D. Gibbon, rfaguiam.py"""

    with open('FIGURES/paramnames.txt', 'w') as f:
        f.write(paramnames)

    plt.savefig(figfullname+".png")
    plt.savefig(figfullname+".pdf")

    canvas.draw()

####################################################################

def create_figure():

    global fig, canvas

#    plt.clf()

    fig = plt.figure(figsize=(figwidth, figheight))

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack(side=tk.BOTTOM)


####################################################################
# Main program

def run_program():

    global f0array, f0framerate, f0frameduration

    signalmessage = "Run RFA: takes time, signal duration %.3fs, sampling %dHz)."%(sigdur, fs)
    statusboxrefresh(signalmessage)
    print("... process input")
    process_input()
    print("... process input done")

    signalmessage = "... analysing spectrum (FFT) ..."
    statusboxrefresh(signalmessage)
    print("... createlfspectrum")
    createlfspectrum()
    print("... createlfspectrum done")

    signalmessage = "... creating spectrogram matrix ..."
    statusboxrefresh(signalmessage)
    print("... create spec matrix")
    createspecmatrix()
    print("... create spec matrix done")

    signalmessage = "... creating graphics ..."
    statusboxrefresh(signalmessage)
    print("... creategraphics")
    creategraphics()
    print("... creategraphics done")

    statusboxrefresh("Waveform - LF spectrograms (heatmap, waterfall) - LF spectrum")

####################################################################
####################################################################
# Widget helpers

def temp_text(e):
    infobox.delete(0,"end")

def onHelp():
    helpmessage = """RFAGUI
D. Gibbon, alpha-01, 2023-03-20

RFAGUI (Rhythm Formant Analysis Graphical User Interface) is a visualiser for exploratory research on acoustic rhythm patterns (cf. Gibbon 2021). Python3 with NumPy and MatPlotLib is required.

"File": load data from file explorer, show filename; "RFA": analyse; "Status": show status; "Defaults": reset parameters; "Help": show help text; "Quit": exit. Values shown (except and status) are modifiable.

Spectrogram values refer to sample count.

Gibbon, Dafydd. 2021. The rhythms of rhythm. Journal of the International Phonetic Association, First View, pp. 1-33.
DOI={10.1017/S0025100321000086}"""
    tk.messagebox.showinfo("RFA help message",  helpmessage)

def onResults():
    resultmessage = """Results (located in RESULTS directory)
1. AM: AM envelope, LF Spectrum, LF Spectrum Peak vector, LF Spectrogram mag traj
2. Output file formats: filebase, date, envsmo, minfreq, maxfreq, sigwin, sampskipper, waterfallskip, specsmo, peakheight"""
    tk.messagebox.showinfo("RFA help message",  resultmessage)

####################################################################
# create input fields for figure variables

def create_widgets():

    global filebox, envsmooth, minfrequency, maxfrequency, sigwindowing, sampskipping, waterfallskipping, waterfallmultiplying, specsmooth, speccoloring, peakheightval, infobox

    print("... create_widgets")

#    frame01 = tk.Frame(root, highlightbackground="blue", highlightthickness=1)
    frame01 = tk.Frame(root)
    frame01.pack(side=tk.TOP)

#    frame02 = tk.LabelFrame(root, text="AM", relief=tk.FLAT, highlightbackground="blue", highlightthickness=1)
    frame02 = tk.LabelFrame(root)
    frame02.pack(side=tk.LEFT)

#############################################################
# Parameter frames

# AM

    frame03 = tk.LabelFrame(frame02, text="   Waveform  ", relief=tk.FLAT, width=60, highlightbackground="blue", highlightthickness=1)
    frame03.pack(padx=1, pady=5)

    frame04 = tk.LabelFrame(frame02, text="    Heatmap   ", relief=tk.FLAT, width=60, highlightbackground="blue", highlightthickness=1)
    frame04.pack(padx=1, pady=5)

    frame05 = tk.LabelFrame(frame02, text="   Waterfall   ", relief=tk.FLAT, width=60, highlightbackground="blue", highlightthickness=1)
    frame05.pack(padx=1, pady=5)

    frame06 = tk.LabelFrame(frame02, text="LF spectrum", relief=tk.FLAT, width=60, highlightbackground="blue", highlightthickness=1)
    frame06.pack(padx=1, pady=5)

##############################################################
# AM

    envsmooth_label = tk.Label(frame03, text="LP medfilt:")
    envsmooth_label.pack()
    envsmooth = tk.Entry(frame03, width=4)
    envsmooth.pack()


    minfrequency_label = tk.Label(frame03, text="Min freq:")
    minfrequency_label.pack()
    minfrequency = tk.Entry(frame03, width=4)
    minfrequency.pack()

    maxfrequency_label = tk.Label(frame03, text="Max freq:")
    maxfrequency_label.pack()
    maxfrequency = tk.Entry(frame03, width=4)
    maxfrequency.pack()

    sigwindowing_label = tk.Label(frame04, text="Window (s):")
    sigwindowing_label.pack()
    sigwindowing = tk.Entry(frame04, width=4)
    sigwindowing.pack()

    sampskipping_label = tk.Label(frame04, text="skip (s):")
    sampskipping_label.pack()
    sampskipping = tk.Entry(frame04, width=4)
    sampskipping.pack()

    waterfallskipping_label = tk.Label(frame05, text="skip (s):")
    waterfallskipping_label.pack()
    waterfallskipping = tk.Entry(frame05, width=4)
    waterfallskipping.pack()

    waterfallmultiplying_label = tk.Label(frame05, text="Multiplier:")
    waterfallmultiplying_label.pack()
    waterfallmultiplying = tk.Entry(frame05, width=4)
    waterfallmultiplying.pack()

    specsmooth_label = tk.Label(frame06, text="LP medfilt:")
    specsmooth_label.pack()
    specsmooth = tk.Entry(frame06, width=4)
    specsmooth.pack()

    speccoloring_label = tk.Label(frame06, text="Colour:")
    speccoloring_label.pack()
    speccoloring = tk.Entry(frame06, width=4)
    speccoloring.pack()

    peakheightval_label = tk.Label(frame06, text="Peak min:")
    peakheightval_label.pack()
    peakheightval = tk.Entry(frame06, width=4)
    peakheightval.pack()

####################################################################
####################################################################
# Buttons etc.

    dir_button = tk.Button(frame01, text="File", command=read_wav)
    dir_button.pack(side=tk.LEFT)

# Filebox
#    fileboxlabel = tk.Label(frame01, text="Status:")
#    fileboxlabel.pack(side=tk.LEFT, padx=10, pady=10)
    filebox = tk.Entry(frame01, bg="lightgrey", width=25, borderwidth=2)
    filebox.pack(side=tk.LEFT, pady=35)
    filebox.bind("<FocusIn>", temp_text)

# Run
    run_button = tk.Button(frame01, text="RFA", command=run_program)
    run_button.pack(side=tk.LEFT)

# Infobox
    infoboxlabel = tk.Label(frame01, text="Status:")
    infoboxlabel.pack(side=tk.LEFT, padx=10, pady=10)
    infobox = tk.Entry(frame01, bg="lightgrey", width=60, borderwidth=2)
    infobox.pack(side=tk.LEFT, pady=35)

# Add a button to return to defaults
    default_button = tk.Button(frame01, text="Defaults", command=param_default)
    default_button.pack(side=tk.LEFT)

    help_button = tk.Button(frame01, text="Help", command=onHelp)
    help_button.pack(side=tk.LEFT)

    results_button = tk.Button(frame01, text="Results", command=onResults)
    results_button.pack(side=tk.LEFT)

    quit_button = tk.Button(frame01, text="Quit", command=root.destroy)
    quit_button.pack(side=tk.LEFT)

####################################################################
# Set params

def param_default():

    param_delete()

    # AM

    envsmooth.insert(tk.END, 1)
    minfrequency.insert(tk.END, 0) # Max LF spectrum frequency
    maxfrequency.insert(tk.END, 3) # Max LF spectrum frequency
    sigwindowing.insert(tk.END, 3)
    sampskipping.insert(tk.END, 0.05)
    waterfallskipping.insert(tk.END, 6)
    waterfallmultiplying.insert(tk.END, 0.05)

    specsmooth.insert(tk.END, 9) # Default medfilt value
    speccoloring.insert(tk.END, speccolor)
    peakheightval.insert(tk.END, 1.1)

####################################################################

def param_delete():

    # AM
    envsmooth.delete(0, tk.END)
    minfrequency.delete(0, tk.END)
    maxfrequency.delete(0, tk.END)
    sigwindowing.delete(0, tk.END)
    sampskipping.delete(0, tk.END)
    waterfallskipping.delete(0, tk.END)
    waterfallmultiplying.delete(0, tk.END)
#    linewidth.delete(0, tk.END)
    specsmooth.delete(0, tk.END)
    speccoloring.delete(0, tk.END)
    peakheightval.delete(0, tk.END)

####################################################################

def statusboxrefresh(signalmessage):

    print(signalmessage)

    infobox.delete(0, tk.END)
    infobox.insert(0, signalmessage)

####################################################################
# Workarounds - recheck program logic

def resetvariables():

    plt.clf()

    axwave.set_xticks([]); axwave.set_yticks([])
    axwater.set_xticks([]); axwater.set_yticks([])
    axheat.set_xticks([]); axheat.set_yticks([])
    axspect.set_xticks([]); axspect.set_yticks([])

    axwave.set_xticklabels(""); axwave.set_yticklabels("")
    axwater.set_xticklabels(""); axwater.set_yticklabels("")
    axheat.set_xticklabels(""); axheat.set_yticklabels("")
    axspect.set_xticklabels(""); axspect.set_yticklabels("")

    axwave.set_xlabel(""); axwave.set_ylabel("")
    axwater.set_xlabel(""); axwater.set_ylabel("")
    axheat.set_xlabel(""); axheat.set_ylabel("")
    axspect.set_xlabel(""); axspect.set_ylabel("")

####################################################################
####################################################################
# One-time procedures.

print("RUNTIME: CREATING WIDGETS")

create_widgets()

param_default()

infobox.delete(0, tk.END)
infobox.insert(0, "Expecting mono WAV file.")

filebox.delete(0, tk.END)
filebox.insert(0, "<-- search filename")

print("RUNTIME: CREATING FIGURE")
create_figure()

if argvlen > 1:
    read_wav()

####################################################################
# run main event loop

root.mainloop()

####################################################################

