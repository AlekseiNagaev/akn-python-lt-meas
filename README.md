# python-lt-meas
Python code for LT measurements

analysis.py includes vol2ppr that checks the input data for pulses above specified level

conv.py includes API for unit conversions

instruments.py includes API for waveform signal generator, multimeter, resistance bridge, DC power source and lock-in amplifier:
  - Base class instr with GPIB/COM/LAN connection setup 
  - Waveform signal generator class wfsg
  - Multimeter class mult
  - Resistance bridge class rbrg
  - DC power source dcpw
  - Lock-in amplifier class liao
  
main.py is the main file

measurements.py includes programs for IV and Pulse measurements

nic.py includes the program for data aquisition with NI DAQ card

nuf.py, plfit.py and plfit2.py include programs for numerical fitting of the measurement data with theoretical formulas

setup.py sets measurement options for further use 

tef.py includes a function for loading measurement data from Matlab file
