# ALPHA-g online analyzer

# To build

git clone git clone https://bitbucket.org/teamalphag/agana.git
cd agana
git clone https://bitbucket.org/teamalphag/agcfmdb.git
make
./agana.exe -h

# To run online from midas

- build agana.exe
- ./agana.exe ### should connect to MIDAS, if a run is in progress, should see events
- in MIDAS ODB:
-- in /Programs/rootana:
--- set "required" to "y"
--- set "start command" to "cd /home/agmini/online/agana; ./agana.exe -R8081 -Hlocalhost -- --wfexport"
-- in /Webserver/Prooxy
--- add TID_STRING "rootana" with value "http://localhost:8081" ### "8081" is the port number in -R on the agana command line
-- in /Alias
--- add TID_STRING "rootana" with value "/proxy/rootana"
- reload the midas status page
- start analyzer from the "programs" page "start rootana" button
- stop analyzer from the "programs" page
- view analyzer plots by following link "rootana" in the midas left-hand side menu

# To run online or offline

Online: ./agana.exe -Hlocalhost -- ...module arguments...
Offline: ./agana.exe /daq/alpha_agmini/data/run903333sub*.mid.lz4 -- ...module arguments...

Generated files:

run903333.root - output histograms

# Special diagnostic functions:

## Waveform viewer

./agana.exe ... -- --wfexport

- view the ADC waveforms under manalyzer/adc_waveforms
- view the PWB waveforms under manalyzer/pwb_waveforms

## Diagnostics of bad PWB waveform

Bad PWB waveforms are reported with prefix XXX, they can
be saved into the output ROOT file:

./agana.exe ... -- --pwb-wf-save-bad | grep ^XXX

## Calibration of PWB data suppression

The pwb_module has an exact simulator of the PWB data suppression algorithm. It is used
to calibrate the data suppression (measure the best threshold value) and
to validate the FPGA implementation of the algorithm.

To calibrate the data suppression, record some uncomressed (ch_force=y) PWB data
and run the analyzer:

./agana.exe ... -- --pwb-wf-suppress | grep ^TTT

Tag TTT will reports all mismatches between fpga and software data suppression. The suppression
threshold ch_force is taken from ODB. To specify a different threshold, use command
line switch"--pwb-wf-threshold TTT".

Calibration produces histograms in the pads/wfsuppress subdirectory:

- ADC amp - hit pulse height amplitude
- ADC amp pos, neg - amplitude of positive and negative polarity pulses
- cumul keep - number of kept data, the higher the threshold the fewer data is kept (100% -> 0%)
- cumul drop - inverse of the "keep" plot, number of data dropped, the higher the threshold, the more is dropped (0% -> 100%)
- cumul keep map - the 2D "keep" plot to see if some PWBs need different thresholds
- pwbNN_keep - per-pwb "keep" plots
- pwbNN_drop - inverse of the "keep" plots

The best data suppression threshold is selected by looking at the "cumul_keep" plot
and selecting the ch_threshold value where data is sufficiently reduced. Lowest possible
value should be used to avoid suppression of valid TPC hits that have a small pulse height.

In addition, plots are made to calibrate the old "min adc" data suppression (no longer used):

- adc_min_map - 2D plot of adc_min for each PWB
- pwbNN_adc_min - per-pwb plot of adc_min

# End
