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

view the ADC waveforms under manalyzer/adc_waveforms
view the PWB waveforms under manalyzer/pwb_waveforms

# End
