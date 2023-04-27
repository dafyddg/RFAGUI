# RFAGUI

## Code for RFA demo, teaching and parameter experimentation, with GUI.

There is currently one RFAGUI application, rfaguiam.py, illustrating low frequency spectral analysis of the amplitude modulation of speech as a correlate of perceived speech rhythms. A companion application for Low Frequency Frequency Modulation of speech, rfaguiamfm.py, will be added later. The applications provide graphical user interfaces (GUI) for rhythm formant analysis (RFA). A rhythm formant is a higher magnitude frequency band in the low frequency spectrum (usually below 10Hz) which corresponds to a perceived speech rhythm.

The application rfaguiam.py takes a mono WAV audio file input, often several minutes in length, preferably with 16kHz or less sampling frequency for speed in procedding. Relevant processing parameters can be adjusted on-screen in upper and left-hand bars. The output is a four panel on-screen figure (also saved as PNG and PDF files):
  - Waveform and AM envelope
  - LF spectrogram (heatmap format)
  - LF spectrogram (waterfall format)
  - LF spectrum with entire recording as FFT window

The basic principles are described in detail in the following article:
@article{gibbon_2023,
title={The rhythms of rhythm},
volume={53},
DOI={10.1017/S0025100321000086},
number={1},
journal={Journal of the International Phonetic Association},
publisher={Cambridge University Press},
author={Gibbon, Dafydd},
year={2023},
pages={233–265}}

** Note: in addition to the RFAGUI project, the main RFA implementation with a more complete CLI application suite is in the RFA project repository.**
